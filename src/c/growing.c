/*
 * growing.c : a stack whose top layer grows: the monitoring curve at one
 * wavelength, the spectra after each deposited layer, and a growing top layer
 * sampled over a wavelength grid. Conventions are in tmm_internal.h.
 */

#include "tmm_internal.h"

/* ── Shared tails for the growing-stack kernels ──────────────────────────────
 * Forward tail: reproduces the [B,C] → r,t → R,T block of tmm_core on an
 * already-built characteristic matrix. A is not formed here; the incoherent
 * slab combination that consumes these results computes it after the sum, and
 * a semi-infinite caller forms 1−R−T itself.
 *
 * Reverse tail: reflectance of the same stack seen from the substrate side.
 * The characteristic matrix of the reversed stack is the anti-transpose of the
 * forward one (every layer matrix is invariant under anti-transposition), so
 * the reverse pass costs no second product, and the common rescale factor
 * cancels in the amplitude ratio.
 *
 * The substrate enters as its vector [bs, cs] (see substrateVector), so its
 * admittance cs/bs is the incident one of the reverse pass. */

static inline void growingTailFwd(mat2 M, cx eta0, vec2 substrate, double logScale,
                                  double *outR, double *outT) {
    vec2 BC = cmatvec(M, substrate);
    cx B = BC.x, C = BC.y;
    cx eta0B = cmul(eta0, B);
    cx r = cdiv(csub(eta0B, C), cadd(eta0B, C));
    cx t = cmul(cdiv(cmul(cmk(2.0, 0.0), eta0), cadd(eta0B, C)), substrate.x);
    double R = cabs2(r);
    double T = transmittedFlux(substrate) / eta0.re * cabs2(t) * exp(-2.0 * logScale);
    if (T < 0.0) T = 0.0;
    *outR = R; *outT = T;
}

static inline double growingTailRev(mat2 M, cx eta0, vec2 substrate) {
    cx B = cadd(M.d, cmul(M.b, eta0));
    cx C = cadd(M.c, cmul(M.a, eta0));
    cx etaSB = cmul(substrate.y, B);
    cx bsC = cmul(substrate.x, C);
    cx r = cdiv(csub(etaSB, bsC), cadd(etaSB, bsC));
    return cabs2(r);
}

/* ── Exported: monitor curve of one growing layer ────────────────────────────
 * The signal of one layer as it grows on a completed stack, at one wavelength,
 * both polarizations: the incremental control algorithm of Tikhonravov &
 * Trubetskov, Appl. Opt. 44, 6877 (2005). The completed stack's characteristic
 * matrix is built once; each sample then costs one layer matrix, one 2×2
 * multiply and the tails, so a sweep is linear in samples rather than in
 * samples × stack depth.
 *
 * The growing layer faces the incident medium, so its matrix multiplies the
 * completed product from the LEFT.
 *
 *   base  : NB triples [n_re, n_im, d_nm], the completed stack outermost first
 *           (zero-thickness entries skipped, as everywhere)
 *   ng    : growing layer ñ at this wavelength
 *   dArr  : nD sample thicknesses of the growing layer (nm); d ≤ 0 evaluates
 *           the bare completed stack
 * Outputs (each nD):
 *   outRs/outTs, outRp/outTp : forward R and T of the coated surface
 *   outRrs/outRrp            : its reflectance from the substrate side, for
 *                              callers that model the witness as an incoherent
 *                              slab with a bare back face */

TMM_EXPORT("tmm_monitor_curve")
void tmm_monitor_curve(double lambda_nm, double theta_deg,
                       double n0_re, double n0_im, double ns_re, double ns_im,
                       const double *base, int NB,
                       double ng_re, double ng_im,
                       const double *dArr, int nD,
                       double *outRs, double *outTs,
                       double *outRp, double *outTp,
                       double *outRrs, double *outRrp) {
    cx n0 = cmk(n0_re, n0_im);
    cx ns = cmk(ns_re, ns_im);
    cx ng = cmk(ng_re, ng_im);
    incidence_t a = incidence(theta_deg);
    cx cosTheta0 = incidentCosTheta(n0, a);

    for (int pol = 0; pol < 2; pol++) {
        cx eta0 = layerAdmittance(n0, cosTheta0, pol);
        vec2 substrate = substrateVector(ns, snellCosTheta(n0, a, ns), pol);

        mat2 Mb;
        Mb.a = cmk(1.0, 0.0); Mb.b = cmk(0.0, 0.0);
        Mb.c = cmk(0.0, 0.0); Mb.d = cmk(1.0, 0.0);
        double logScaleB = 0.0;
        for (int i = 0; i < NB; i++) {
            cx n = cmk(base[3 * i + 0], base[3 * i + 1]);
            double d = base[3 * i + 2];
            if (!(d > 0.0)) continue;
            cx cosThetaJ = snellCosTheta(n0, a, n);
            Mb = matmul(Mb, layerMatrix(n, d, lambda_nm, cosThetaJ, pol));
            logScaleB += layerLogScale(n, d, lambda_nm, cosThetaJ) + rescaleMatrix(&Mb);
        }
        cx cosThetaG = snellCosTheta(n0, a, ng);

        double *oR  = pol ? outRp  : outRs;
        double *oT  = pol ? outTp  : outTs;
        double *oRr = pol ? outRrp : outRrs;

        for (int k = 0; k < nD; k++) {
            mat2 M = Mb;
            double logScale = logScaleB;
            double d = dArr[k];
            if (d > 0.0) {
                M = matmul(layerMatrix(ng, d, lambda_nm, cosThetaG, pol), Mb);
                logScale += layerLogScale(ng, d, lambda_nm, cosThetaG) + rescaleMatrix(&M);
            }
            growingTailFwd(M, eta0, substrate, logScale, &oR[k], &oT[k]);
            oRr[k] = growingTailRev(M, eta0, substrate);
        }
    }
}

/* ── Exported: spectra of a growing stack, one per deposited layer ───────────
 * The front-surface passes of a deposition run in a single call: layers arrive
 * in DEPOSITION order (first deposited, substrate-adjacent, first), each is
 * folded into the running product from the left, since the newest layer faces
 * the incident medium, and the spectrum after every fold is written out. The
 * same incremental structure as tmm_monitor_curve over a wavelength grid, so
 * all N step spectra together cost what the final one costs alone.
 *
 * Memory layout mirrors tmm_rta_spectrum:
 *   matNK : N × nLam × 2, [layer][λ][re, im], deposition order
 *   thick : N thicknesses (nm); a step of zero thickness repeats the previous
 *           spectrum, matching the d ≤ 0 skip of the reference loop
 * Outputs, each N × nLam, step-major ([step][λ]):
 *   outRs/outTs, outRp/outTp : forward R and T after each step
 *   outRrs/outRrp            : reflectance from the substrate side after each
 *                              step, for the incoherent slab combination */

/* Returns 1 on success; 0, with the outputs untouched, when the working
 * state could not be allocated. */
TMM_EXPORT("tmm_deposition_spectra")
int tmm_deposition_spectra(const double *lambdas, int nLam,
                           const double *theta_deg, int nTheta,
                           const double *n0arr, const double *nsarr,
                           const double *matNK, const double *thick, int N,
                           double *outRs, double *outTs,
                           double *outRp, double *outTp,
                           double *outRrs, double *outRrp) {
    /* Running product and admittances per (λ, pol); pol-major blocks. */
    size_t cells = (size_t)nLam * 2;
    mat2 *M = (mat2 *)malloc(sizeof(mat2) * cells);
    double *logScale = (double *)malloc(sizeof(double) * cells);
    cx *eta0v = (cx *)malloc(sizeof(cx) * cells);
    vec2 *substrateV = (vec2 *)malloc(sizeof(vec2) * cells);
    incidence_t *angle = (incidence_t *)malloc(sizeof(incidence_t) * (nLam > 0 ? nLam : 1));
    void *const blocks[] = { M, logScale, eta0v, substrateV, angle };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        free(M); free(logScale); free(eta0v); free(substrateV); free(angle);
        return 0;
    }

    for (int li = 0; li < nLam; li++) {
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        incidence_t a = incidence(theta_at(theta_deg, nTheta, li));
        angle[li] = a;
        cx cosTheta0 = incidentCosTheta(n0, a);
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * nLam + li;
            eta0v[at] = layerAdmittance(n0, cosTheta0, pol);
            substrateV[at] = substrateVector(ns, snellCosTheta(n0, a, ns), pol);
            M[at].a = cmk(1.0, 0.0); M[at].b = cmk(0.0, 0.0);
            M[at].c = cmk(0.0, 0.0); M[at].d = cmk(1.0, 0.0);
            logScale[at] = 0.0;
        }
    }

    for (int k = 0; k < N; k++) {
        double d = thick[k];
        for (int li = 0; li < nLam; li++) {
            double lam = lambdas[li];
            cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
            cx n = cmk(matNK[((size_t)k * nLam + li) * 2 + 0],
                       matNK[((size_t)k * nLam + li) * 2 + 1]);
            for (int pol = 0; pol < 2; pol++) {
                size_t at = (size_t)pol * nLam + li;
                if (d > 0.0) {
                    cx cosThetaJ = snellCosTheta(n0, angle[li], n);
                    M[at] = matmul(layerMatrix(n, d, lam, cosThetaJ, pol), M[at]);
                    logScale[at] += layerLogScale(n, d, lam, cosThetaJ) + rescaleMatrix(&M[at]);
                }
                size_t out = (size_t)k * nLam + li;
                double *oR  = pol ? outRp  : outRs;
                double *oT  = pol ? outTp  : outTs;
                double *oRr = pol ? outRrp : outRrs;
                growingTailFwd(M[at], eta0v[at], substrateV[at], logScale[at],
                               &oR[out], &oT[out]);
                oRr[out] = growingTailRev(M[at], eta0v[at], substrateV[at]);
            }
        }
    }

    free(M); free(logScale); free(eta0v); free(substrateV); free(angle);
    return 1;
}

/* ── Exported: growing top layer over a wavelength grid ──────────────────────
 * The stateful counterpart of tmm_monitor_curve, batched the other way: one
 * thickness of the growing layer per call, the WHOLE wavelength grid at once.
 * This is the shape a broadband monitor scan needs, where every scan reads a
 * full spectrum of the growing stack and the completed layers beneath do not
 * change until the layer is cut.
 *
 * Lifecycle: create() folds the completed stack's characteristic matrices for
 * every (λ, pol) once and keeps them; set_top() declares the growing layer's
 * ñ(λ); sample() then answers one thickness across the grid, costing one layer
 * matrix, one 2×2 multiply and the tails per (λ, pol). free() releases the
 * state. Handles are opaque pointers; the caller owns their lifetime.
 *
 * Layout matches tmm_deposition_spectra: matNK is NB × nLam × 2 layer-major
 * with the completed stack OUTERMOST FIRST (the fold is M = M · M_k);
 * zero-thickness entries are skipped. Outputs of sample(), each nLam: forward
 * R and T plus the substrate-side reflectance, for the incoherent slab
 * combination. */

typedef struct {
    int nLam;
    int topSet;
    incidence_t *angle; /* nLam */
    double *lam;        /* nLam */
    cx *n0;             /* nLam */
    cx *ng;             /* nLam, set_top */
    cx *cosThetaG;      /* nLam, set_top (polarization-independent) */
    cx *eta0;           /* 2·nLam, pol-major */
    vec2 *substrate;    /* 2·nLam */
    mat2 *Mb;           /* 2·nLam completed-stack product */
    double *logScaleB;  /* 2·nLam */
} growing_stack;

TMM_EXPORT("tmm_growing_stack_free")
void tmm_growing_stack_free(growing_stack *h) {
    if (!h) return;
    free(h->angle); free(h->lam); free(h->n0); free(h->ng); free(h->cosThetaG);
    free(h->eta0); free(h->substrate); free(h->Mb); free(h->logScaleB);
    free(h);
}

TMM_EXPORT("tmm_growing_stack_create")
growing_stack *tmm_growing_stack_create(const double *lambdas, int nLam,
                                        const double *theta_deg, int nTheta,
                                        const double *n0arr, const double *nsarr,
                                        const double *matNK, const double *thick,
                                        int NB) {
    growing_stack *h = (growing_stack *)malloc(sizeof(growing_stack));
    if (!h) return 0;
    size_t count = (size_t)(nLam > 0 ? nLam : 1);
    h->nLam = nLam;
    h->topSet = 0;
    h->angle = (incidence_t *)malloc(sizeof(incidence_t) * count);
    h->lam = (double *)malloc(sizeof(double) * count);
    h->n0 = (cx *)malloc(sizeof(cx) * count);
    h->ng = (cx *)malloc(sizeof(cx) * count);
    h->cosThetaG = (cx *)malloc(sizeof(cx) * count);
    h->eta0 = (cx *)malloc(sizeof(cx) * 2 * count);
    h->substrate = (vec2 *)malloc(sizeof(vec2) * 2 * count);
    h->Mb = (mat2 *)malloc(sizeof(mat2) * 2 * count);
    h->logScaleB = (double *)malloc(sizeof(double) * 2 * count);
    void *const blocks[] = { h->angle, h->lam, h->n0, h->ng, h->cosThetaG,
                             h->eta0, h->substrate, h->Mb, h->logScaleB };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        tmm_growing_stack_free(h);
        return 0;
    }

    for (int li = 0; li < nLam; li++) {
        h->lam[li] = lambdas[li];
        h->angle[li] = incidence(theta_at(theta_deg, nTheta, li));
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        h->n0[li] = n0;
        cx cosTheta0 = incidentCosTheta(n0, h->angle[li]);
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * nLam + li;
            h->eta0[at] = layerAdmittance(n0, cosTheta0, pol);
            h->substrate[at] = substrateVector(ns, snellCosTheta(n0, h->angle[li], ns), pol);
            mat2 M;
            M.a = cmk(1.0, 0.0); M.b = cmk(0.0, 0.0);
            M.c = cmk(0.0, 0.0); M.d = cmk(1.0, 0.0);
            double logScale = 0.0;
            for (int k = 0; k < NB; k++) {
                double d = thick[k];
                if (!(d > 0.0)) continue;
                cx n = cmk(matNK[((size_t)k * nLam + li) * 2 + 0],
                           matNK[((size_t)k * nLam + li) * 2 + 1]);
                cx cosThetaJ = snellCosTheta(n0, h->angle[li], n);
                M = matmul(M, layerMatrix(n, d, lambdas[li], cosThetaJ, pol));
                logScale += layerLogScale(n, d, lambdas[li], cosThetaJ) + rescaleMatrix(&M);
            }
            h->Mb[at] = M;
            h->logScaleB[at] = logScale;
        }
    }
    return h;
}

TMM_EXPORT("tmm_growing_stack_set_top")
void tmm_growing_stack_set_top(growing_stack *h, const double *ngNK) {
    if (!h) return;
    for (int li = 0; li < h->nLam; li++) {
        cx ng = cmk(ngNK[2 * li + 0], ngNK[2 * li + 1]);
        h->ng[li] = ng;
        h->cosThetaG[li] = snellCosTheta(h->n0[li], h->angle[li], ng);
    }
    h->topSet = 1;
}

/* Returns 1 on success; 0, with the outputs untouched, for a null handle or
 * a d > 0 sample before set_top declared the growing layer. The status is
 * what keeps an unwritten buffer from being read back as data. */
TMM_EXPORT("tmm_growing_stack_sample")
int tmm_growing_stack_sample(growing_stack *h, double d,
                             double *outRs, double *outTs,
                             double *outRp, double *outTp,
                             double *outRrs, double *outRrp) {
    if (!h || (d > 0.0 && !h->topSet)) return 0;
    for (int li = 0; li < h->nLam; li++) {
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * h->nLam + li;
            mat2 M = h->Mb[at];
            double logScale = h->logScaleB[at];
            if (d > 0.0) {
                M = matmul(layerMatrix(h->ng[li], d, h->lam[li], h->cosThetaG[li], pol),
                           h->Mb[at]);
                logScale += layerLogScale(h->ng[li], d, h->lam[li], h->cosThetaG[li])
                          + rescaleMatrix(&M);
            }
            double *oR  = pol ? outRp  : outRs;
            double *oT  = pol ? outTp  : outTs;
            double *oRr = pol ? outRrp : outRrs;
            growingTailFwd(M, h->eta0[at], h->substrate[at], logScale, &oR[li], &oT[li]);
            oRr[li] = growingTailRev(M, h->eta0[at], h->substrate[at]);
        }
    }
    return 1;
}
