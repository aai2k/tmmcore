/*
 * tmm_kernel.c : transfer-matrix method for multilayer thin films.
 *
 * Computes reflectance, transmittance and absorptance for a stack of absorbing,
 * dispersive layers at arbitrary angle of incidence, in s and p polarization,
 * together with exact analytic derivatives: the thickness Jacobian, the
 * thickness Hessian, and the needle-insertion P-function.
 *
 * Also computes phase, group delay, group delay dispersion and third-order
 * dispersion, by carrying the same characteristic matrix in third-order Taylor
 * arithmetic with angular frequency as the differentiation variable.
 *
 * Conventions:
 *   ñ = n + ik              k > 0 for absorbing media
 *   time factor exp(−iωt)   a wave exp(i(kz − ωt)) decays for k > 0
 *   off-diagonals of the characteristic matrix carry −i
 *   wavelengths and thicknesses in nm, angles in degrees from normal
 *
 * All arithmetic is double precision. The JavaScript implementation shipped
 * alongside this file is the behavioural reference; the two agree to float64
 * round-off (see tests/).
 *
 * Builds with any C99 compiler; no dependencies beyond libm.
 *
 * References:
 *   • Macleod, Thin-Film Optical Filters 5th ed., §2.4, Eqs. 2.111, 2.123–2.125;
 *     Eq. 2.83; §10.2, Eqs. 10.11–10.13 (tilted admittances, real invariant)
 *   • Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6)
 *   • Tikhonravov, Trubetskov & DeBell, Appl. Opt. 35, 5493 (1996)
 *   • Birge & Kärtner, Appl. Opt. 45, 1478 (2006)   [phase dispersion]
 */

#include <math.h>
#include <stdlib.h>

/* WebAssembly builds export by name; native builds need no annotation. */
#if defined(__wasm__) || defined(__EMSCRIPTEN__)
#define TMM_EXPORT(name) __attribute__((export_name(name)))
#else
#define TMM_EXPORT(name)
#endif

#define PI 3.14159265358979323846
#define MAX_IM_DELTA 50.0
#define MATRIX_RESCALE_THRESHOLD 1e100

/* ── Complex number arithmetic (cx = {re, im}) ───────────────────────────────
 * Mirrors thinFilmMath.js cadd/csub/cmul/cdiv/cabs2/cconj/csqrt/ccos/csin. */

typedef struct { double re, im; } cx;

static inline cx cmk(double re, double im) { cx z; z.re = re; z.im = im; return z; }

static inline cx cadd(cx a, cx b) { return cmk(a.re + b.re, a.im + b.im); }
static inline cx csub(cx a, cx b) { return cmk(a.re - b.re, a.im - b.im); }
static inline cx cmul(cx a, cx b) {
    return cmk(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}
static inline cx cdiv(cx a, cx b) {
    double d = b.re * b.re + b.im * b.im;
    return cmk((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
}
static inline double cabs2(cx a) { return a.re * a.re + a.im * a.im; }
static inline cx cconj(cx a) { return cmk(a.re, -a.im); }
static inline cx csqrt_(cx a) {
    /* r = sqrt(sqrt(re²+im²)); theta = atan2(im,re)/2  : JS csqrt verbatim */
    double r = sqrt(sqrt(a.re * a.re + a.im * a.im));
    double theta = atan2(a.im, a.re) / 2.0;
    return cmk(r * cos(theta), r * sin(theta));
}
static inline cx ccos_(cx a) {
    return cmk(cos(a.re) * cosh(a.im), -sin(a.re) * sinh(a.im));
}
static inline cx csin_(cx a) {
    return cmk(sin(a.re) * cosh(a.im), cos(a.re) * sinh(a.im));
}

/* ── 2×2 complex matrix [[a,b],[c,d]] and 2-vector {x,y} ──────────────────── */

typedef struct { cx a, b, c, d; } mat2;
typedef struct { cx x, y; } vec2;

static inline mat2 matmul(mat2 A, mat2 B) {
    mat2 M;
    M.a = cadd(cmul(A.a, B.a), cmul(A.b, B.c));
    M.b = cadd(cmul(A.a, B.b), cmul(A.b, B.d));
    M.c = cadd(cmul(A.c, B.a), cmul(A.d, B.c));
    M.d = cadd(cmul(A.c, B.b), cmul(A.d, B.d));
    return M;
}
/* Common scaling cancels from reflectance (Macleod Eq. 2.123, p. 45);
 * tmm_core retains the scale required by transmittance (Eq. 2.125). */
static inline double rescaleMatrix(mat2 *M) {
    double scale = 0.0;
    scale = fmax(scale, fmax(fabs(M->a.re), fabs(M->a.im)));
    scale = fmax(scale, fmax(fabs(M->b.re), fabs(M->b.im)));
    scale = fmax(scale, fmax(fabs(M->c.re), fabs(M->c.im)));
    scale = fmax(scale, fmax(fabs(M->d.re), fabs(M->d.im)));
    if (scale <= MATRIX_RESCALE_THRESHOLD) return 0.0;
    double inverse = 1.0 / scale;
    M->a.re *= inverse; M->a.im *= inverse;
    M->b.re *= inverse; M->b.im *= inverse;
    M->c.re *= inverse; M->c.im *= inverse;
    M->d.re *= inverse; M->d.im *= inverse;
    return log(scale);
}
static inline vec2 cmatvec(mat2 M, vec2 v) {
    vec2 o;
    o.x = cadd(cmul(M.a, v.x), cmul(M.b, v.y));
    o.y = cadd(cmul(M.c, v.x), cmul(M.d, v.y));
    return o;
}

/* ── Snell's law: cosθ_j from incident (n0, sinθ0) into medium nj ───────────
 * The transverse invariant is Re(n0) sinθ0 (Macleod 5th ed., §10.2, Eqs.
 * 10.11–10.13): an absorbing incident medium carries a wave whose amplitude
 * falls along the normal only, and every medium's cosθ follows from that one
 * real invariant. Mirrors snellCosTheta / incidentCosTheta in tmm.js. */

static inline cx snellCosTheta(cx n0, cx sinTheta0, cx nj) {
    cx sinThetaJ = cdiv(cmul(cmk(n0.re, 0.0), sinTheta0), nj);
    return csqrt_(csub(cmk(1.0, 0.0), cmul(sinThetaJ, sinThetaJ)));
}

/* cosθ0 of the incident medium: the plain cosine for a transparent medium,
 * otherwise from the same real invariant as the layers. */
static inline cx incidentCosTheta(cx n0, cx sinTheta0) {
    if (n0.im == 0.0) return csqrt_(csub(cmk(1.0, 0.0), cmul(sinTheta0, sinTheta0)));
    return snellCosTheta(n0, sinTheta0, n0);
}

/* ── Layer characteristic matrix (pol: 0 = s, 1 = p) ─────────────────────── */

static inline mat2 layerMatrix(cx nj, double dj_nm, double lambda_nm, cx cosTheta_j, int pol) {
    double k0 = (2.0 * PI) / lambda_nm;
    cx delta = cmul(cmul(nj, cmk(k0 * dj_nm, 0.0)), cosTheta_j);
    if (delta.im > MAX_IM_DELTA) delta.im = MAX_IM_DELTA;
    else if (delta.im < -MAX_IM_DELTA) delta.im = -MAX_IM_DELTA;
    cx cosD = ccos_(delta);
    cx sinD = csin_(delta);

    cx eta = (pol == 0) ? cmul(nj, cosTheta_j)   /* s: n cosθ */
                        : cdiv(nj, cosTheta_j);  /* p: n / cosθ */

    cx negI = cmk(0.0, -1.0);
    cx iSinD_div_eta = cmul(negI, cdiv(sinD, eta));
    cx iEta_sinD     = cmul(negI, cmul(eta, sinD));

    mat2 M;
    M.a = cosD;          M.b = iSinD_div_eta;
    M.c = iEta_sinD;     M.d = cosD;
    return M;
}

/* ── Core TMM for one wavelength / polarization ──────────────────────────────
 * Faithful port of tmm() in thinFilmMath.js. `layers` is N triples
 * [n_re, n_im, d_nm]. Writes R,T,A to out[0..2]. (Zero-thickness layers skipped,
 * exactly like the JS reference.) */

static void tmm_core(double lambda_nm, double theta_deg, int pol,
                     cx n0, cx ns, const double *layers, int N,
                     double *outR, double *outT, double *outA) {
    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);
    cx cosTheta0 = incidentCosTheta(n0, sinTheta0);

    cx eta0 = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);

    cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    cx etaS = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    /* identity */
    mat2 M;
    M.a = cmk(1.0, 0.0); M.b = cmk(0.0, 0.0);
    M.c = cmk(0.0, 0.0); M.d = cmk(1.0, 0.0);
    double logScale = 0.0;

    for (int i = 0; i < N; i++) {
        cx n = cmk(layers[3 * i + 0], layers[3 * i + 1]);
        double d = layers[3 * i + 2];
        if (d <= 0.0) continue;
        cx cosThetaJ = snellCosTheta(n0, sinTheta0, n);
        mat2 Mj = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
        M = matmul(M, Mj);
        logScale += rescaleMatrix(&M);
    }

    cx B = cadd(M.a, cmul(M.b, etaS));
    cx C = cadd(M.c, cmul(M.d, etaS));
    cx eta0B = cmul(eta0, B);
    cx r = cdiv(csub(eta0B, C), cadd(eta0B, C));
    cx t = cdiv(cmul(cmk(2.0, 0.0), eta0), cadd(eta0B, C));

    double R = cabs2(r);
    double T = etaS.re / eta0.re * cabs2(t) * exp(-2.0 * logScale);
    if (T < 0.0) T = 0.0;
    double A = 1.0 - R - T;
    if (A < 0.0) A = 0.0;

    *outR = R; *outT = T; *outA = A;
}

/* ── Exported: single (λ, θ, pol) evaluation ─────────────────────────────────
 * out[0..2] = R, T, A. pol: 0 = s, 1 = p. Drop-in for tmm() in thinFilmMath.js. */

TMM_EXPORT("tmm_one")
void tmm_one(double lambda_nm, double theta_deg, int pol,
             double n0_re, double n0_im, double ns_re, double ns_im,
             const double *layers, int N, double *out) {
    tmm_core(lambda_nm, theta_deg, pol,
             cmk(n0_re, n0_im), cmk(ns_re, ns_im), layers, N,
             &out[0], &out[1], &out[2]);
}

/* ── Exported: batched spectrum over a wavelength grid ───────────────────────
 * One WASM call evaluates the full stack across all `nLam` wavelengths for BOTH
 * polarizations, amortizing the JS↔WASM boundary (the actual speed win). This
 * is the batched analogue of looping tmmAvg() in evaluateSpectrum().
 *
 * Memory layout (all f64, caller-owned, written from JS):
 *   lambdas : nLam                    wavelengths (nm)
 *   n0arr   : nLam × 2  [re, im]      incident-medium ñ per λ
 *   nsarr   : nLam × 2  [re, im]      substrate ñ per λ
 *   matNK   : N × nLam × 2            per-layer ñ per λ, layout [layer][λ][re,im]
 *   thick   : N                       layer thicknesses (nm)
 * Outputs (each nLam):
 *   outRs,outTs,outAs  (s-pol)   outRp,outTp,outAp  (p-pol)
 * Caller forms avg = (s+p)/2 and selects the requested polarization in JS,
 * exactly as tmmAvg() does, so the JS-side semantics are unchanged. */

TMM_EXPORT("tmm_spectrum")
void tmm_spectrum(const double *lambdas, int nLam,
                  const double *n0arr, const double *nsarr,
                  const double *matNK, const double *thick, int N,
                  double theta_deg,
                  double *outRs, double *outTs, double *outAs,
                  double *outRp, double *outTp, double *outAp) {
    /* Reusable per-λ layer scratch [n_re, n_im, d] × N. */
    double *layers = (N > 0) ? (double *)malloc(sizeof(double) * 3 * N) : NULL;

    for (int li = 0; li < nLam; li++) {
        double lam = lambdas[li];
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);

        for (int k = 0; k < N; k++) {
            /* matNK[k][li][.] */
            long base = ((long)k * nLam + li) * 2;
            layers[3 * k + 0] = matNK[base + 0];
            layers[3 * k + 1] = matNK[base + 1];
            layers[3 * k + 2] = thick[k];
        }

        tmm_core(lam, theta_deg, 0, n0, ns, layers, N, &outRs[li], &outTs[li], &outAs[li]);
        tmm_core(lam, theta_deg, 1, n0, ns, layers, N, &outRp[li], &outTp[li], &outAp[li]);
    }

    if (layers) free(layers);
}

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
 * cancels in the amplitude ratio. */

static inline void growingTailFwd(mat2 M, cx eta0, cx etaS, double logScale,
                                  double *outR, double *outT) {
    cx B = cadd(M.a, cmul(M.b, etaS));
    cx C = cadd(M.c, cmul(M.d, etaS));
    cx eta0B = cmul(eta0, B);
    cx r = cdiv(csub(eta0B, C), cadd(eta0B, C));
    cx t = cdiv(cmul(cmk(2.0, 0.0), eta0), cadd(eta0B, C));
    double R = cabs2(r);
    double T = etaS.re / eta0.re * cabs2(t) * exp(-2.0 * logScale);
    if (T < 0.0) T = 0.0;
    *outR = R; *outT = T;
}

static inline double growingTailRev(mat2 M, cx eta0, cx etaS) {
    cx B = cadd(M.d, cmul(M.b, eta0));
    cx C = cadd(M.c, cmul(M.a, eta0));
    cx etaSB = cmul(etaS, B);
    cx r = cdiv(csub(etaSB, C), cadd(etaSB, C));
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
 * completed product from the LEFT, exactly as createMonitorTmmEvaluator in the
 * JS reference.
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
    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);
    cx cosTheta0 = incidentCosTheta(n0, sinTheta0);

    for (int pol = 0; pol < 2; pol++) {
        cx eta0 = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
        cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
        cx etaS = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

        mat2 Mb;
        Mb.a = cmk(1.0, 0.0); Mb.b = cmk(0.0, 0.0);
        Mb.c = cmk(0.0, 0.0); Mb.d = cmk(1.0, 0.0);
        double logScaleB = 0.0;
        for (int i = 0; i < NB; i++) {
            cx n = cmk(base[3 * i + 0], base[3 * i + 1]);
            double d = base[3 * i + 2];
            if (d <= 0.0) continue;
            cx cosThetaJ = snellCosTheta(n0, sinTheta0, n);
            Mb = matmul(Mb, layerMatrix(n, d, lambda_nm, cosThetaJ, pol));
            logScaleB += rescaleMatrix(&Mb);
        }
        cx cosThetaG = snellCosTheta(n0, sinTheta0, ng);

        double *oR  = pol ? outRp  : outRs;
        double *oT  = pol ? outTp  : outTs;
        double *oRr = pol ? outRrp : outRrs;

        for (int k = 0; k < nD; k++) {
            mat2 M = Mb;
            double logScale = logScaleB;
            double d = dArr[k];
            if (d > 0.0) {
                M = matmul(layerMatrix(ng, d, lambda_nm, cosThetaG, pol), Mb);
                logScale += rescaleMatrix(&M);
            }
            growingTailFwd(M, eta0, etaS, logScale, &oR[k], &oT[k]);
            oRr[k] = growingTailRev(M, eta0, etaS);
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
 * Memory layout mirrors tmm_spectrum:
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
                           const double *n0arr, const double *nsarr,
                           const double *matNK, const double *thick, int N,
                           double theta_deg,
                           double *outRs, double *outTs,
                           double *outRp, double *outTp,
                           double *outRrs, double *outRrp) {
    /* Running product and admittances per (λ, pol); pol-major blocks. */
    size_t cells = (size_t)nLam * 2;
    mat2 *M = (mat2 *)malloc(sizeof(mat2) * cells);
    double *logScale = (double *)malloc(sizeof(double) * cells);
    cx *eta0v = (cx *)malloc(sizeof(cx) * cells);
    cx *etaSv = (cx *)malloc(sizeof(cx) * cells);
    if (!M || !logScale || !eta0v || !etaSv) {
        free(M); free(logScale); free(eta0v); free(etaSv);
        return 0;
    }

    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);

    for (int li = 0; li < nLam; li++) {
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        cx cosTheta0 = incidentCosTheta(n0, sinTheta0);
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * nLam + li;
            eta0v[at] = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
            cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
            etaSv[at] = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);
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
                    cx cosThetaJ = snellCosTheta(n0, sinTheta0, n);
                    M[at] = matmul(layerMatrix(n, d, lam, cosThetaJ, pol), M[at]);
                    logScale[at] += rescaleMatrix(&M[at]);
                }
                size_t out = (size_t)k * nLam + li;
                double *oR  = pol ? outRp  : outRs;
                double *oT  = pol ? outTp  : outTs;
                double *oRr = pol ? outRrp : outRrs;
                growingTailFwd(M[at], eta0v[at], etaSv[at], logScale[at],
                               &oR[out], &oT[out]);
                oRr[out] = growingTailRev(M[at], eta0v[at], etaSv[at]);
            }
        }
    }

    free(M); free(logScale); free(eta0v); free(etaSv);
    return 1;
}

/* ── Exported: persistent growing-layer evaluator over a wavelength grid ─────
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
 * with the completed stack OUTERMOST FIRST (the fold is M = M · M_k, exactly
 * the JS reference createMonitorTmmEvaluator); zero-thickness entries are
 * skipped. Outputs of sample(), each nLam: forward R and T plus the
 * substrate-side reflectance, for the incoherent slab combination. */

typedef struct {
    int nLam;
    cx sinTheta0;
    int topSet;
    double *lam;        /* nLam */
    cx *n0;             /* nLam */
    cx *ng;             /* nLam, set_top */
    cx *cosThetaG;      /* nLam, set_top (polarization-independent) */
    cx *eta0;           /* 2·nLam, pol-major */
    cx *etaS;           /* 2·nLam */
    mat2 *Mb;           /* 2·nLam completed-stack product */
    double *logScaleB;  /* 2·nLam */
} growing_eval;

TMM_EXPORT("tmm_growing_eval_free")
void tmm_growing_eval_free(growing_eval *h) {
    if (!h) return;
    free(h->lam); free(h->n0); free(h->ng); free(h->cosThetaG);
    free(h->eta0); free(h->etaS); free(h->Mb); free(h->logScaleB);
    free(h);
}

TMM_EXPORT("tmm_growing_eval_create")
growing_eval *tmm_growing_eval_create(const double *lambdas, int nLam,
                                      double theta_deg,
                                      const double *n0arr, const double *nsarr,
                                      const double *matNK, const double *thick,
                                      int NB) {
    growing_eval *h = (growing_eval *)malloc(sizeof(growing_eval));
    if (!h) return 0;
    h->nLam = nLam;
    h->topSet = 0;
    h->lam = (double *)malloc(sizeof(double) * nLam);
    h->n0 = (cx *)malloc(sizeof(cx) * nLam);
    h->ng = (cx *)malloc(sizeof(cx) * nLam);
    h->cosThetaG = (cx *)malloc(sizeof(cx) * nLam);
    h->eta0 = (cx *)malloc(sizeof(cx) * 2 * nLam);
    h->etaS = (cx *)malloc(sizeof(cx) * 2 * nLam);
    h->Mb = (mat2 *)malloc(sizeof(mat2) * 2 * nLam);
    h->logScaleB = (double *)malloc(sizeof(double) * 2 * nLam);
    if (!h->lam || !h->n0 || !h->ng || !h->cosThetaG
        || !h->eta0 || !h->etaS || !h->Mb || !h->logScaleB) {
        tmm_growing_eval_free(h);
        return 0;
    }

    h->sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);

    for (int li = 0; li < nLam; li++) {
        h->lam[li] = lambdas[li];
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        h->n0[li] = n0;
        cx cosTheta0 = incidentCosTheta(n0, h->sinTheta0);
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * nLam + li;
            h->eta0[at] = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
            cx cosThetaS = snellCosTheta(n0, h->sinTheta0, ns);
            h->etaS[at] = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);
            mat2 M;
            M.a = cmk(1.0, 0.0); M.b = cmk(0.0, 0.0);
            M.c = cmk(0.0, 0.0); M.d = cmk(1.0, 0.0);
            double logScale = 0.0;
            for (int k = 0; k < NB; k++) {
                double d = thick[k];
                if (d <= 0.0) continue;
                cx n = cmk(matNK[((size_t)k * nLam + li) * 2 + 0],
                           matNK[((size_t)k * nLam + li) * 2 + 1]);
                cx cosThetaJ = snellCosTheta(n0, h->sinTheta0, n);
                M = matmul(M, layerMatrix(n, d, lambdas[li], cosThetaJ, pol));
                logScale += rescaleMatrix(&M);
            }
            h->Mb[at] = M;
            h->logScaleB[at] = logScale;
        }
    }
    return h;
}

TMM_EXPORT("tmm_growing_eval_set_top")
void tmm_growing_eval_set_top(growing_eval *h, const double *ngNK) {
    if (!h) return;
    for (int li = 0; li < h->nLam; li++) {
        cx ng = cmk(ngNK[2 * li + 0], ngNK[2 * li + 1]);
        h->ng[li] = ng;
        h->cosThetaG[li] = snellCosTheta(h->n0[li], h->sinTheta0, ng);
    }
    h->topSet = 1;
}

/* Returns 1 on success; 0, with the outputs untouched, for a null handle or
 * a d > 0 sample before set_top declared the growing layer. The status is
 * what keeps an unwritten buffer from being read back as data. */
TMM_EXPORT("tmm_growing_eval_sample")
int tmm_growing_eval_sample(growing_eval *h, double d,
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
                logScale += rescaleMatrix(&M);
            }
            double *oR  = pol ? outRp  : outRs;
            double *oT  = pol ? outTp  : outTs;
            double *oRr = pol ? outRrp : outRrs;
            growingTailFwd(M, h->eta0[at], h->etaS[at], logScale, &oR[li], &oT[li]);
            oRr[li] = growingTailRev(M, h->eta0[at], h->etaS[at]);
        }
    }
    return 1;
}

/* ── Exported: analytic thickness Jacobian for one (λ, θ, pol) ────────────────
 * Faithful port of tmmThicknessJacobian() in thinFilmMath.js. Returns the exact
 * analytic dR/dd_k, dT/dd_k, dA/dd_k for every layer at one sample : the DLS
 * refiner's per-step gradient (2·N fewer evals than central differences).
 *
 * `layers` is N triples [n_re, n_im, d]; layers used AS-IS (no d>0 filter) for
 * index parity with the caller's design, exactly like the JS reference.
 * Outputs (each length N): dRdd, dTdd, dAdd. Also writes base R,T,A to base[0..2]. */

TMM_EXPORT("tmm_jacobian")
void tmm_jacobian(double lambda_nm, double theta_deg, int pol,
                  double n0_re, double n0_im, double ns_re, double ns_im,
                  const double *layers, int N,
                  double *dRdd, double *dTdd, double *dAdd, double *base) {
    cx n0 = cmk(n0_re, n0_im);
    cx ns = cmk(ns_re, ns_im);

    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);
    cx cosTheta0 = incidentCosTheta(n0, sinTheta0);
    cx eta0 = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    cx etaS = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    /* Per-layer cosθ and characteristic matrices. */
    cx   *cosThJ = (cx *)  malloc(sizeof(cx)   * (N > 0 ? N : 1));
    mat2 *Ms     = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    for (int k = 0; k < N; k++) {
        cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        double d = layers[3 * k + 2];
        if (d < 0.0) d = 0.0;
        cosThJ[k] = snellCosTheta(n0, sinTheta0, n);
        Ms[k] = layerMatrix(n, d, lambda_nm, cosThJ[k], pol);
    }

    /* Pre[j] = M_0·…·M_{j-1};  Post[j] = M_j·…·M_{N-1}·[1, ηs]. */
    mat2 *Pre = (mat2 *)malloc(sizeof(mat2) * (N + 1));
    vec2 *Post = (vec2 *)malloc(sizeof(vec2) * (N + 1));
    Pre[0].a = cmk(1.0, 0.0); Pre[0].b = cmk(0.0, 0.0);
    Pre[0].c = cmk(0.0, 0.0); Pre[0].d = cmk(1.0, 0.0);
    for (int j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    Post[N].x = cmk(1.0, 0.0); Post[N].y = etaS;
    for (int j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    cx Bv = Post[0].x, Cv = Post[0].y;
    cx den = cadd(cmul(eta0, Bv), Cv);
    cx den2 = cmul(den, den);
    cx r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    cx t = cdiv(cmul(cmk(2.0, 0.0), eta0), den);
    double R = cabs2(r);
    double Tfac = etaS.re / eta0.re;
    double T = Tfac * cabs2(t); if (T < 0.0) T = 0.0;
    double A = 1.0 - R - T; if (A < 0.0) A = 0.0;
    base[0] = R; base[1] = T; base[2] = A;

    cx f = cdiv(cmul(cmk(2.0, 0.0), eta0), den2);
    double k0 = (2.0 * PI) / lambda_nm;

    for (int k = 0; k < N; k++) {
        cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        double d = layers[3 * k + 2]; if (d < 0.0) d = 0.0;
        cx cth = cosThJ[k];
        cx etaK = (pol == 0) ? cmul(n, cth) : cdiv(n, cth);
        cx Q = cmul(cmul(n, cmk(k0, 0.0)), cth);                 /* (2π/λ) n cosθ */
        cx delta = cmul(cmul(n, cmk(k0 * d, 0.0)), cth);
        cx cD = ccos_(delta), sD = csin_(delta);

        /* dM_k/dd_k = Q · [[ −sinδ, −i cosδ/η ], [ −i η cosδ, −sinδ ]] */
        cx negI = cmk(0.0, -1.0), neg1 = cmk(-1.0, 0.0);
        mat2 dMk;
        dMk.a = cmul(Q, cmul(neg1, sD));
        dMk.b = cmul(Q, cmul(negI, cdiv(cD, etaK)));
        dMk.c = cmul(Q, cmul(negI, cmul(etaK, cD)));
        dMk.d = cmul(Q, cmul(neg1, sD));

        vec2 dV = cmatvec(Pre[k], cmatvec(dMk, Post[k + 1]));
        cx dB = dV.x, dC = dV.y;

        /* metrics(dB,dC) : verbatim from validated tmmNeedleScan.metrics */
        cx dr = cmul(f, csub(cmul(Cv, dB), cmul(Bv, dC)));
        double dR = 2.0 * (cmul(cconj(r), dr)).re;
        cx dt = cmul(neg1, cmul(f, cadd(cmul(eta0, dB), dC)));
        double dT = Tfac * 2.0 * (cmul(cconj(t), dt)).re;
        dRdd[k] = dR; dTdd[k] = dT; dAdd[k] = -(dR + dT);
    }

    free(cosThJ); free(Ms); free(Pre); free(Post);
}

/* ── Analytic needle P-function scan ─────────────────────────────────────────
 * Faithful port of tmmNeedleScan() in thinFilmMath.js : the d→0 limit of
 * Sullivan's pre/post method (Tikhonravov's analytic P-function). Returns the
 * merit-gradient ingredients {dR,dT,dA} of inserting an infinitesimal needle of
 * each candidate index at every gap position (0..N) and, optionally, at intra-
 * layer split fractions. The {dR,dT,dA} chain rule is identical to tmm_jacobian.
 *
 * `layers` used AS-IS (no d>0 filter), index parity with the caller's stack.
 * Outputs (caller-owned f64):
 *   base  : 3                      [R,T,A]
 *   gaps  : (N+1)*nCand*3          layout [pos][cand][dR,dT,dA]
 *   intra : N*nFrac*nCand*3        layout [layer][frac][cand][dR,dT,dA] (nFrac>0)
 */

/* {dR,dT,dA} from d[B,C]/dd : verbatim from tmm_jacobian.metrics. */
static void needle_metrics(cx Bv, cx Cv, cx eta0, cx f, cx r, cx t, double Tfac,
                           cx dB, cx dC, double *o) {
    cx dr = cmul(f, csub(cmul(Cv, dB), cmul(Bv, dC)));
    double dR = 2.0 * (cmul(cconj(r), dr)).re;
    cx dt = cmul(cmk(-1.0, 0.0), cmul(f, cadd(cmul(eta0, dB), dC)));
    double dT = Tfac * 2.0 * (cmul(cconj(t), dt)).re;
    o[0] = dR; o[1] = dT; o[2] = -(dR + dT);
}

TMM_EXPORT("tmm_needle_scan")
void tmm_needle_scan(double lambda_nm, double theta_deg, int pol,
                     double n0_re, double n0_im, double ns_re, double ns_im,
                     const double *layers, int N,
                     const double *candNs, int nCand,
                     const double *fracs, int nFrac,
                     double *base, double *gaps, double *intra) {
    cx n0 = cmk(n0_re, n0_im), ns = cmk(ns_re, ns_im);
    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);
    cx cosTheta0 = incidentCosTheta(n0, sinTheta0);
    cx eta0 = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    cx etaS = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    double k0 = (2.0 * PI) / lambda_nm;

    cx   *cosThJ = (cx *)  malloc(sizeof(cx)   * (N > 0 ? N : 1));
    mat2 *Ms     = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    for (int k = 0; k < N; k++) {
        cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        double d = layers[3 * k + 2]; if (d < 0.0) d = 0.0;
        cosThJ[k] = snellCosTheta(n0, sinTheta0, n);
        Ms[k] = layerMatrix(n, d, lambda_nm, cosThJ[k], pol);
    }

    mat2 *Pre  = (mat2 *)malloc(sizeof(mat2) * (N + 1));
    vec2 *Post = (vec2 *)malloc(sizeof(vec2) * (N + 1));
    Pre[0].a = cmk(1.0, 0.0); Pre[0].b = cmk(0.0, 0.0);
    Pre[0].c = cmk(0.0, 0.0); Pre[0].d = cmk(1.0, 0.0);
    for (int j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    Post[N].x = cmk(1.0, 0.0); Post[N].y = etaS;
    for (int j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    cx Bv = Post[0].x, Cv = Post[0].y;
    cx den = cadd(cmul(eta0, Bv), Cv);
    cx den2 = cmul(den, den);
    cx r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    cx t = cdiv(cmul(cmk(2.0, 0.0), eta0), den);
    double R = cabs2(r);
    double Tfac = etaS.re / eta0.re;
    double T = Tfac * cabs2(t); if (T < 0.0) T = 0.0;
    double A = 1.0 - R - T; if (A < 0.0) A = 0.0;
    base[0] = R; base[1] = T; base[2] = A;

    cx f = cdiv(cmul(cmk(2.0, 0.0), eta0), den2);

    /* Needle derivative-matrix A per candidate: A = [[0,−iQ/ηₐ],[−iηₐQ,0]]. */
    mat2 *Ac = (mat2 *)malloc(sizeof(mat2) * (nCand > 0 ? nCand : 1));
    for (int c = 0; c < nCand; c++) {
        cx nA = cmk(candNs[2 * c + 0], candNs[2 * c + 1]);
        cx cthA = snellCosTheta(n0, sinTheta0, nA);
        cx etaA = (pol == 0) ? cmul(nA, cthA) : cdiv(nA, cthA);
        cx Q = cmul(cmul(nA, cmk(k0, 0.0)), cthA);
        cx negI = cmk(0.0, -1.0);
        mat2 Am;
        Am.a = cmk(0.0, 0.0);                 Am.b = cmul(negI, cdiv(Q, etaA));
        Am.c = cmul(negI, cmul(etaA, Q));     Am.d = cmk(0.0, 0.0);
        Ac[c] = Am;
    }

    /* Gaps: pos = 0..N, every candidate. */
    for (int pos = 0; pos <= N; pos++) {
        vec2 post = Post[pos];
        mat2 pre  = Pre[pos];
        for (int c = 0; c < nCand; c++) {
            vec2 dV = cmatvec(pre, cmatvec(Ac[c], post));
            needle_metrics(Bv, Cv, eta0, f, r, t, Tfac, dV.x, dV.y,
                           &gaps[((long)pos * nCand + c) * 3]);
        }
    }

    /* Intra-layer splits (host-split), if requested. */
    if (nFrac > 0) {
        for (int k = 0; k < N; k++) {
            cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
            double d = layers[3 * k + 2];     /* raw d (matches JS intra) */
            cx cth = cosThJ[k];
            for (int fi = 0; fi < nFrac; fi++) {
                double frac = fracs[fi];
                double dl = frac * d;        if (dl < 1e-9) dl = 1e-9;
                double dr_ = (1.0 - frac) * d; if (dr_ < 1e-9) dr_ = 1e-9;
                mat2 Mleft  = layerMatrix(n, dl,  lambda_nm, cth, pol);
                mat2 Mright = layerMatrix(n, dr_, lambda_nm, cth, pol);
                mat2 preIn  = matmul(Pre[k], Mleft);
                vec2 postIn = cmatvec(Mright, Post[k + 1]);
                for (int c = 0; c < nCand; c++) {
                    vec2 dV = cmatvec(preIn, cmatvec(Ac[c], postIn));
                    long off = (((long)k * nFrac + fi) * nCand + c) * 3;
                    needle_metrics(Bv, Cv, eta0, f, r, t, Tfac, dV.x, dV.y, &intra[off]);
                }
            }
        }
    }

    free(cosThJ); free(Ms); free(Pre); free(Post); free(Ac);
}

/* ── Analytic thickness-Hessian kernel ───────────
 * LINE-BY-LINE port of tmmThicknessHessian() in thinFilmMath.js : the EXACT
 * analytic second derivatives ∂²{R,T,A}/∂dᵢ∂dⱼ (full N×N symmetric) plus the
 * first derivatives, at one (λ,θ,pol). Used by the bounded-SQP / Newton inner
 * refiner; the JS remains the oracle (tests/wasm_hessian_equivalence.mjs).
 *
 * Outputs (caller-owned f64):
 *   base   : 3        [R,T,A]
 *   dRdd/dTdd/dAdd : N each      first derivatives (verbatim tmm_jacobian)
 *   d2Rdd/d2Tdd/d2Add : N*N each row-major, FULL symmetric (both triangles set)
 */
TMM_EXPORT("tmm_hessian")
void tmm_hessian(double lambda_nm, double theta_deg, int pol,
                 double n0_re, double n0_im, double ns_re, double ns_im,
                 const double *layers, int N,
                 double *dRdd, double *dTdd, double *dAdd,
                 double *d2Rdd, double *d2Tdd, double *d2Add, double *base) {
    cx n0 = cmk(n0_re, n0_im);
    cx ns = cmk(ns_re, ns_im);

    cx sinTheta0 = cmk(sin(theta_deg * PI / 180.0), 0.0);
    cx cosTheta0 = incidentCosTheta(n0, sinTheta0);
    cx eta0 = (pol == 0) ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    cx cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    cx etaS = (pol == 0) ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    cx   *cosThJ = (cx *)  malloc(sizeof(cx)   * (N > 0 ? N : 1));
    mat2 *Ms     = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    for (int k = 0; k < N; k++) {
        cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        double d = layers[3 * k + 2]; if (d < 0.0) d = 0.0;
        cosThJ[k] = snellCosTheta(n0, sinTheta0, n);
        Ms[k] = layerMatrix(n, d, lambda_nm, cosThJ[k], pol);
    }

    mat2 *Pre  = (mat2 *)malloc(sizeof(mat2) * (N + 1));
    vec2 *Post = (vec2 *)malloc(sizeof(vec2) * (N + 1));
    Pre[0].a = cmk(1.0, 0.0); Pre[0].b = cmk(0.0, 0.0);
    Pre[0].c = cmk(0.0, 0.0); Pre[0].d = cmk(1.0, 0.0);
    for (int j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    Post[N].x = cmk(1.0, 0.0); Post[N].y = etaS;
    for (int j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    cx Bv = Post[0].x, Cv = Post[0].y;
    cx den  = cadd(cmul(eta0, Bv), Cv);
    cx den2 = cmul(den, den);
    cx den3 = cmul(den2, den);
    cx r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    cx t = cdiv(cmul(cmk(2.0, 0.0), eta0), den);
    double R = cabs2(r);
    double Tfac = etaS.re / eta0.re;
    double T = Tfac * cabs2(t); if (T < 0.0) T = 0.0;
    double A = 1.0 - R - T; if (A < 0.0) A = 0.0;
    base[0] = R; base[1] = T; base[2] = A;
    cx f = cdiv(cmul(cmk(2.0, 0.0), eta0), den2);

    double k0 = (2.0 * PI) / lambda_nm;

    /* Per-layer first-derivative pieces: dM[k], d2M[k], v[k]=dMₖ·Post[k+1],
     * [dB,dC] and the first-derivative metrics. */
    mat2 *dM  = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    mat2 *d2M = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    vec2 *v   = (vec2 *)malloc(sizeof(vec2) * (N > 0 ? N : 1));
    cx   *dBa = (cx *)  malloc(sizeof(cx)   * (N > 0 ? N : 1));
    cx   *dCa = (cx *)  malloc(sizeof(cx)   * (N > 0 ? N : 1));
    cx negI = cmk(0.0, -1.0), posI = cmk(0.0, 1.0), neg1 = cmk(-1.0, 0.0);
    for (int k = 0; k < N; k++) {
        cx n = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        double d = layers[3 * k + 2]; if (d < 0.0) d = 0.0;
        cx cth = cosThJ[k];
        cx etaK = (pol == 0) ? cmul(n, cth) : cdiv(n, cth);
        cx Q  = cmul(cmul(n, cmk(k0, 0.0)), cth);     /* (2π/λ) n cosθ */
        cx Q2 = cmul(Q, Q);
        cx delta = cmul(cmul(n, cmk(k0 * d, 0.0)), cth);
        cx cD = ccos_(delta), sD = csin_(delta);
        /* dMₖ/ddₖ = Q·[[ −sinδ, −i cosδ/η ], [ −i η cosδ, −sinδ ]] */
        dM[k].a = cmul(Q, cmul(neg1, sD));
        dM[k].b = cmul(Q, cmul(negI, cdiv(cD, etaK)));
        dM[k].c = cmul(Q, cmul(negI, cmul(etaK, cD)));
        dM[k].d = cmul(Q, cmul(neg1, sD));
        /* d²Mₖ/ddₖ² = Q²·[[ −cosδ, i sinδ/η ], [ i η sinδ, −cosδ ]] */
        d2M[k].a = cmul(Q2, cmul(neg1, cD));
        d2M[k].b = cmul(Q2, cmul(posI, cdiv(sD, etaK)));
        d2M[k].c = cmul(Q2, cmul(posI, cmul(etaK, sD)));
        d2M[k].d = cmul(Q2, cmul(neg1, cD));
        v[k] = cmatvec(dM[k], Post[k + 1]);
        vec2 dVk = cmatvec(Pre[k], v[k]);
        dBa[k] = dVk.x; dCa[k] = dVk.y;
        cx dr = cmul(f, csub(cmul(Cv, dBa[k]), cmul(Bv, dCa[k])));
        double dR = 2.0 * (cmul(cconj(r), dr)).re;
        cx dt = cmul(neg1, cmul(f, cadd(cmul(eta0, dBa[k]), dCa[k])));
        double dT = Tfac * 2.0 * (cmul(cconj(t), dt)).re;
        dRdd[k] = dR; dTdd[k] = dT; dAdd[k] = -(dR + dT);
    }

    mat2 I; I.a = cmk(1.0, 0.0); I.b = cmk(0.0, 0.0); I.c = cmk(0.0, 0.0); I.d = cmk(1.0, 0.0);
    for (int i = 0; i < N; i++) {
        mat2 Wmat_i = matmul(Pre[i], dM[i]);   /* Pre[i]·dMᵢ (used for j>i) */
        mat2 Cmid = I;                          /* M_{i+1}···M_{j-1}, empty at j=i+1 */
        for (int j = i; j < N; j++) {
            cx d2Bv, d2Cv;
            if (j == i) {
                vec2 w = cmatvec(Pre[i], cmatvec(d2M[i], Post[i + 1]));
                d2Bv = w.x; d2Cv = w.y;
            } else {
                vec2 w = cmatvec(Wmat_i, cmatvec(Cmid, v[j]));
                d2Bv = w.x; d2Cv = w.y;
            }
            cx dBi = dBa[i], dCi = dCa[i], dBj = dBa[j], dCj = dCa[j];
            cx dr_i = cmul(f, csub(cmul(Cv, dBi), cmul(Bv, dCi)));
            cx dr_j = cmul(f, csub(cmul(Cv, dBj), cmul(Bv, dCj)));
            cx dden_i = cadd(cmul(eta0, dBi), dCi);
            cx dden_j = cadd(cmul(eta0, dBj), dCj);
            /* d²r_ij = f(dCᵢdBⱼ + C d²B − dBᵢdCⱼ − B d²C) − 2 drⱼ ddenᵢ/den */
            cx innerR = csub(
                cadd(cmul(dCi, dBj), cmul(Cv, d2Bv)),
                cadd(cmul(dBi, dCj), cmul(Bv, d2Cv)));
            cx d2r = csub(cmul(f, innerR),
                          cdiv(cmul(cmul(cmk(2.0, 0.0), dr_j), dden_i), den));
            double d2Rij = 2.0 * ((cmul(cconj(dr_i), dr_j)).re + (cmul(cconj(r), d2r)).re);
            cx dt_i = cmul(neg1, cmul(f, dden_i));
            cx dt_j = cmul(neg1, cmul(f, dden_j));
            cx d2den = cadd(cmul(eta0, d2Bv), d2Cv);
            cx d2t = cadd(
                cmul(cmul(cmk(-2.0, 0.0), eta0), cdiv(d2den, den2)),
                cmul(cmul(cmk(4.0, 0.0), eta0), cdiv(cmul(dden_i, dden_j), den3)));
            double d2Tij = Tfac * 2.0 * ((cmul(cconj(dt_i), dt_j)).re + (cmul(cconj(t), d2t)).re);
            double d2Aij = -(d2Rij + d2Tij);
            d2Rdd[(long)i * N + j] = d2Rdd[(long)j * N + i] = d2Rij;
            d2Tdd[(long)i * N + j] = d2Tdd[(long)j * N + i] = d2Tij;
            d2Add[(long)i * N + j] = d2Add[(long)j * N + i] = d2Aij;
            if (j >= i + 1) Cmid = matmul(Cmid, Ms[j]);   /* advance middle: include M_j */
        }
    }

    free(cosThJ); free(Ms); free(Pre); free(Post);
    free(dM); free(d2M); free(v); free(dBa); free(dCa);
}

/* ── Third-order Taylor jets ─────────────────────────────────────────────────
 * A jet holds [f, f', f''/2!, f'''/3!], each entry complex. Ordinary power-
 * series algebra on these differentiates a function exactly, with no finite
 * differences. Port of taylorJet.js; where the JS multiplies by a reciprocal
 * rather than dividing, so does this, since the two are not bit-identical. */

#define JET_N 4

typedef struct { cx c[JET_N]; } jet;

static inline jet jconst(double re, double im) {
    jet j;
    j.c[0] = cmk(re, im);
    j.c[1] = cmk(0.0, 0.0); j.c[2] = cmk(0.0, 0.0); j.c[3] = cmk(0.0, 0.0);
    return j;
}
static inline jet jread(const double *p) {
    jet j;
    for (int i = 0; i < JET_N; i++) j.c[i] = cmk(p[2 * i], p[2 * i + 1]);
    return j;
}
static inline jet jadd(jet a, jet b) {
    jet o; for (int i = 0; i < JET_N; i++) o.c[i] = cadd(a.c[i], b.c[i]); return o;
}
static inline jet jsub(jet a, jet b) {
    jet o; for (int i = 0; i < JET_N; i++) o.c[i] = csub(a.c[i], b.c[i]); return o;
}
static inline jet jscale(jet a, double s) {
    jet o; for (int i = 0; i < JET_N; i++) o.c[i] = cmk(a.c[i].re * s, a.c[i].im * s); return o;
}
static inline jet jmul(jet a, jet b) {
    jet o;
    for (int order = 0; order < JET_N; order++) {
        cx sum = cmk(0.0, 0.0);
        for (int i = 0; i <= order; i++) sum = cadd(sum, cmul(a.c[i], b.c[order - i]));
        o.c[order] = sum;
    }
    return o;
}
static inline jet jrecip(jet a) {
    jet o;
    o.c[0] = cdiv(cmk(1.0, 0.0), a.c[0]);
    for (int order = 1; order < JET_N; order++) {
        cx sum = cmk(0.0, 0.0);
        for (int i = 1; i <= order; i++) sum = cadd(sum, cmul(a.c[i], o.c[order - i]));
        cx q = cdiv(sum, a.c[0]);
        o.c[order] = cmk(-q.re, -q.im);
    }
    return o;
}
static inline jet jdiv(jet a, jet b) { return jmul(a, jrecip(b)); }

static inline jet jsqrt_j(jet a) {
    jet o;
    o.c[0] = csqrt_(a.c[0]);
    cx twiceRoot = cmk(o.c[0].re * 2.0, o.c[0].im * 2.0);
    for (int order = 1; order < JET_N; order++) {
        cx known = cmk(0.0, 0.0);
        for (int i = 1; i < order; i++) known = cadd(known, cmul(o.c[i], o.c[order - i]));
        o.c[order] = cdiv(csub(a.c[order], known), twiceRoot);
    }
    return o;
}

static void jsincos(jet a, jet *sine, jet *cosine) {
    double re = a.c[0].re, im = a.c[0].im;
    sine->c[0]   = cmk(sin(re) * cosh(im),  cos(re) * sinh(im));
    cosine->c[0] = cmk(cos(re) * cosh(im), -sin(re) * sinh(im));
    for (int order = 1; order < JET_N; order++) {
        cx sineSum = cmk(0.0, 0.0), cosineSum = cmk(0.0, 0.0);
        for (int i = 1; i <= order; i++) {
            cx ts = cmul(a.c[i], cosine->c[order - i]);
            cx tc = cmul(a.c[i], sine->c[order - i]);
            sineSum   = cadd(sineSum,   cmk(ts.re * i, ts.im * i));
            cosineSum = cadd(cosineSum, cmk(tc.re * i, tc.im * i));
        }
        double inv = 1.0 / (double)order;
        sine->c[order]   = cmk( sineSum.re   * inv,  sineSum.im   * inv);
        cosine->c[order] = cmk(-cosineSum.re * inv, -cosineSum.im * inv);
    }
}

/* Past the limit the layer is opaque: the derivatives are zero to machine
 * precision and dropping them keeps cosh from overflowing the whole product. */
static inline jet jclampim(jet a, double limit) {
    if (a.c[0].im > limit || a.c[0].im < -limit) {
        double held = (a.c[0].im > limit) ? limit : -limit;
        jet o;
        o.c[0] = cmk(a.c[0].re, held);
        for (int i = 1; i < JET_N; i++) o.c[i] = cmk(a.c[i].re, 0.0);
        return o;
    }
    return a;
}

/* [f, f', f'', f'''] from the stored [f, f', f''/2!, f'''/3!]. */
static inline void jderivs(jet a, cx *out) {
    out[0] = a.c[0];
    out[1] = a.c[1];
    out[2] = cmk(a.c[2].re * 2.0, a.c[2].im * 2.0);
    out[3] = cmk(a.c[3].re * 6.0, a.c[3].im * 6.0);
}

/* λ(ω) = 2πc/ω. Needs no value for c: with λ and ω given, λ' = −λ/ω. */
static inline jet jwavelength(double lambda, double omega) {
    jet o;
    o.c[0] = cmk(lambda, 0.0);
    o.c[1] = cmk(-lambda / omega, 0.0);
    o.c[2] = cmk(lambda / (omega * omega), 0.0);
    o.c[3] = cmk(-lambda / (omega * omega * omega), 0.0);
    return o;
}

/* ── Jet-valued 2×2 matrices ─────────────────────────────────────────────── */

typedef struct { jet a, b, c, d; } jmat2;

static jmat2 jmatmul(jmat2 A, jmat2 B) {
    jmat2 M;
    M.a = jadd(jmul(A.a, B.a), jmul(A.b, B.c));
    M.b = jadd(jmul(A.a, B.b), jmul(A.b, B.d));
    M.c = jadd(jmul(A.c, B.a), jmul(A.d, B.c));
    M.d = jadd(jmul(A.c, B.b), jmul(A.d, B.d));
    return M;
}
static jmat2 jidentity(void) {
    jmat2 M;
    M.a = jconst(1.0, 0.0); M.b = jconst(0.0, 0.0);
    M.c = jconst(0.0, 0.0); M.d = jconst(1.0, 0.0);
    return M;
}
static jmat2 jzero(void) {
    jmat2 M;
    M.a = jconst(0.0, 0.0); M.b = jconst(0.0, 0.0);
    M.c = jconst(0.0, 0.0); M.d = jconst(0.0, 0.0);
    return M;
}
/* The order-0 matrix controls overflow in the physical coefficient. Once
 * selected, one plain scalar rescales every jet order and cancels from r. */
static double jrescale(jmat2 *M, double threshold) {
    jet *e[4] = { &M->a, &M->b, &M->c, &M->d };
    double scale = 0.0;
    for (int i = 0; i < 4; i++) {
        scale = fmax(scale, fabs(e[i]->c[0].re));
        scale = fmax(scale, fabs(e[i]->c[0].im));
    }
    if (scale <= threshold) return 0.0;
    double inverse = 1.0 / scale;
    for (int i = 0; i < 4; i++) *e[i] = jscale(*e[i], inverse);
    return log(scale);
}
static double jmatmag(jmat2 M) {
    jet *e[4] = { &M.a, &M.b, &M.c, &M.d };
    double magnitude = 0.0;
    for (int i = 0; i < 4; i++)
        for (int o = 0; o < JET_N; o++)
            magnitude = fmax(magnitude, fmax(fabs(e[i]->c[o].re), fabs(e[i]->c[o].im)));
    return magnitude;
}

/* The real part of a jet, order by order. */
static inline jet jreal(jet a) {
    jet o;
    for (int i = 0; i < JET_N; i++) o.c[i] = cmk(a.c[i].re, 0.0);
    return o;
}
static inline int jhas_imag(jet a) {
    for (int i = 0; i < JET_N; i++) if (a.c[i].im != 0.0) return 1;
    return 0;
}

/* Same real invariant Re(n0) sinθ0 as snellCosTheta, in jet arithmetic. */
static inline jet jsnell_cos(jet n0, jet sin0, jet nj) {
    jet s = jdiv(jmul(jreal(n0), sin0), nj);
    return jsqrt_j(jsub(jconst(1.0, 0.0), jmul(s, s)));
}
/* cosθ0 of the incident medium: from the same real invariant as the layers
 * when the medium absorbs at any order, otherwise the plain cosine, taken
 * from the sine jet when the caller supplied one. Mirrors incidentCosine in
 * phase.js. */
static inline jet jincident_cos(jet n0, jet incidentSine, int fromSineJet, double theta_deg) {
    if (jhas_imag(n0)) return jsnell_cos(n0, incidentSine, n0);
    return fromSineJet
        ? jsqrt_j(jsub(jconst(1.0, 0.0), jmul(incidentSine, incidentSine)))
        : jconst(cos(theta_deg * PI / 180.0), 0.0);
}
static inline jet jadmittance(jet n, jet cosv, int pol) {
    return (pol == 0) ? jmul(n, cosv) : jdiv(n, cosv);
}

static jmat2 jlayer_matrix(jet index, double thickness, jet wavelength, jet cosine, int pol) {
    jet phase = jclampim(jscale(jdiv(jmul(index, cosine), wavelength),
                                2.0 * PI * thickness), MAX_IM_DELTA);
    jet sine, cosinePhase;
    jsincos(phase, &sine, &cosinePhase);
    jet eta = jadmittance(index, cosine, pol);
    jet minusI = jconst(0.0, -1.0);
    jmat2 M;
    M.a = cosinePhase;
    M.b = jmul(minusI, jdiv(sine, eta));
    M.c = jmul(minusI, jmul(eta, sine));
    M.d = cosinePhase;
    return M;
}

static void jlayer_matrix_dd(jet index, double thickness, jet wavelength, jet cosine, int pol,
                             jmat2 *M, jmat2 *dM) {
    jet phasePerUnit = jscale(jdiv(jmul(index, cosine), wavelength), 2.0 * PI);
    jet rawPhase = jscale(phasePerUnit, thickness);
    jet phase = jclampim(rawPhase, MAX_IM_DELTA);
    jet phaseDerivative;
    if (rawPhase.c[0].im == phase.c[0].im) {
        phaseDerivative = phasePerUnit;
    } else {
        for (int i = 0; i < JET_N; i++)
            phaseDerivative.c[i] = cmk(phasePerUnit.c[i].re, 0.0);
    }
    jet sine, cosinePhase;
    jsincos(phase, &sine, &cosinePhase);
    jet sineDerivative = jmul(cosinePhase, phaseDerivative);
    jet cosineDerivative = jscale(jmul(sine, phaseDerivative), -1.0);
    jet eta = jadmittance(index, cosine, pol);
    jet minusI = jconst(0.0, -1.0);
    M->a = cosinePhase;
    M->b = jmul(minusI, jdiv(sine, eta));
    M->c = jmul(minusI, jmul(eta, sine));
    M->d = cosinePhase;
    dM->a = cosineDerivative;
    dM->b = jmul(minusI, jdiv(sineDerivative, eta));
    dM->c = jmul(minusI, jmul(eta, sineDerivative));
    dM->d = cosineDerivative;
}

typedef struct { jet reflection, transmission, denominator; } jcoef;

static jcoef jcoef_from_matrix(jmat2 M, jet incidentEta, jet substrateEta, double logScale) {
    jet boundaryB = jadd(M.a, jmul(M.b, substrateEta));
    jet boundaryC = jadd(M.c, jmul(M.d, substrateEta));
    jet incidentB = jmul(incidentEta, boundaryB);
    jcoef o;
    o.denominator = jadd(incidentB, boundaryC);
    o.reflection = jdiv(jsub(incidentB, boundaryC), o.denominator);
    o.transmission = jdiv(jscale(incidentEta, 2.0), o.denominator);
    if (logScale != 0.0) o.transmission = jscale(o.transmission, exp(-logScale));
    return o;
}

static void jcoef_thickness(jmat2 dMat, jcoef base, jet incidentEta, jet substrateEta,
                            jet *dReflection, jet *dTransmission) {
    jet dB = jadd(dMat.a, jmul(dMat.b, substrateEta));
    jet dC = jadd(dMat.c, jmul(dMat.d, substrateEta));
    jet dIncidentB = jmul(incidentEta, dB);
    jet dDenominator = jadd(dIncidentB, dC);
    jet dNumerator = jsub(dIncidentB, dC);
    *dReflection = jdiv(jsub(dNumerator, jmul(base.reflection, dDenominator)), base.denominator);
    *dTransmission = jscale(jdiv(jmul(base.transmission, dDenominator), base.denominator), -1.0);
}

/* ── Phase quantities from a coefficient jet ─────────────────────────────────
 * Writes [phaseRad, GD, GDD, TOD, |coefficient|²]; all NaN where the
 * coefficient is exactly zero and the phase is undefined.
 *
 *   GD  = Im(r'/r)
 *   GDD = Im(r''/r − (r'/r)²)
 *   TOD = Im(r'''/r − 3 r'r''/r² + 2 (r'/r)³)                Birge & Kärtner
 *
 * GD comes out in the reciprocal of the caller's ω unit, GDD in its square and
 * TOD in its cube. */

static void jphase(jet coefficient, double *out5) {
    cx d[4];
    jderivs(coefficient, d);
    cx value = d[0];
    double magnitudeSquared = value.re * value.re + value.im * value.im;
    if (magnitudeSquared == 0.0 || !isfinite(magnitudeSquared)) {
        for (int i = 0; i < 5; i++) out5[i] = NAN;
        return;
    }
    cx inverse = cdiv(cmk(1.0, 0.0), value);
    cx firstRatio  = cmul(d[1], inverse);
    cx secondRatio = cmul(d[2], inverse);
    cx thirdRatio  = cmul(d[3], inverse);
    cx squareFirst = cmk(firstRatio.re * firstRatio.re - firstRatio.im * firstRatio.im,
                         2.0 * firstRatio.re * firstRatio.im);
    cx firstTimesSecond = cmk(
        firstRatio.re * secondRatio.re - firstRatio.im * secondRatio.im,
        firstRatio.re * secondRatio.im + firstRatio.im * secondRatio.re);
    cx cubeFirst = cmk(squareFirst.re * firstRatio.re - squareFirst.im * firstRatio.im,
                       squareFirst.re * firstRatio.im + squareFirst.im * firstRatio.re);
    out5[0] = -atan2(value.im, value.re);
    out5[1] = firstRatio.im;
    out5[2] = secondRatio.im - squareFirst.im;
    out5[3] = thirdRatio.im - 3.0 * firstTimesSecond.im + 2.0 * cubeFirst.im;
    out5[4] = magnitudeSquared;
}

/* ── Phase core: one wavelength, both coefficients ───────────────────────────
 * `sinJet` is the incident-side sine as a jet, for a stack embedded in a
 * dispersive medium at a fixed external angle; NULL uses the constant
 * sin(theta_deg). out10 = [r: phaseRad, GD, GDD, TOD, |r|²][t: same]. */

static void jphase_core(double lambda, double omega, double theta_deg, int pol,
                        jet n0, jet ns, const jet *layerN, const double *thick, int N,
                        const jet *sinJet, double *out10) {
    jet wavelength = jwavelength(lambda, omega);
    jet incidentSine = sinJet ? *sinJet : jconst(sin(theta_deg * PI / 180.0), 0.0);
    jet incidentCosine = jincident_cos(n0, incidentSine, sinJet != NULL, theta_deg);
    jet incidentEta = jadmittance(n0, incidentCosine, pol);
    jet substrateCosine = jsnell_cos(n0, incidentSine, ns);
    jet substrateEta = jadmittance(ns, substrateCosine, pol);

    jmat2 M = jidentity();
    double logScale = 0.0;
    for (int k = 0; k < N; k++) {
        if (!(thick[k] > 0.0)) continue;
        jet cosine = jsnell_cos(n0, incidentSine, layerN[k]);
        M = jmatmul(M, jlayer_matrix(layerN[k], thick[k], wavelength, cosine, pol));
        logScale += jrescale(&M, MATRIX_RESCALE_THRESHOLD);
    }
    jcoef coefficients = jcoef_from_matrix(M, incidentEta, substrateEta, logScale);
    jphase(coefficients.reflection, &out10[0]);
    jphase(coefficients.transmission, &out10[5]);
}

/* ── Exported: phase dispersion at one wavelength ─────────────────────────────
 * Mirrors tmmPhaseDispersion() in phase.js. Index jets are 8 doubles each,
 * [re,im] per order. `sinJet` may be NULL. */

TMM_EXPORT("tmm_phase_one")
void tmm_phase_one(double lambda, double omega, double theta_deg, int pol,
                   const double *n0jet, const double *nsjet,
                   const double *layerJets, const double *thick, int N,
                   const double *sinJet, double *out) {
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    for (int k = 0; k < N; k++) layerN[k] = jread(&layerJets[8 * k]);
    jet sine;
    if (sinJet) sine = jread(sinJet);
    jphase_core(lambda, omega, theta_deg, pol,
                jread(n0jet), jread(nsjet), layerN, thick, N,
                sinJet ? &sine : NULL, out);
    free(layerN);
}

/* ── Exported: batched phase dispersion over a wavelength grid ────────────────
 * One call evaluates the whole grid, amortizing the JS↔WASM boundary the same
 * way tmm_spectrum does. Polarization is an argument rather than both-at-once,
 * because this kernel is an order of magnitude dearer per sample than the plain
 * spectrum and callers at normal incidence would pay twice for nothing.
 *
 * Memory layout (all f64, caller-owned):
 *   lambdas  : nLam
 *   omegas   : nLam                   angular frequency per λ; sets the unit
 *   n0jets   : nLam × 8               incident-medium index jet per λ
 *   nsjets   : nLam × 8               substrate index jet per λ
 *   matJets  : N × nLam × 8           per-layer index jet, layout [layer][λ]
 *   thick    : N
 *   sinJets  : nLam × 8, or NULL
 *   out      : nLam × 10              [r: phaseRad,GD,GDD,TOD,|r|²][t: same] */

TMM_EXPORT("tmm_phase_spectrum")
void tmm_phase_spectrum(const double *lambdas, const double *omegas, int nLam,
                        const double *n0jets, const double *nsjets,
                        const double *matJets, const double *thick, int N,
                        double theta_deg, int pol,
                        const double *sinJets, double *out) {
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    for (int li = 0; li < nLam; li++) {
        for (int k = 0; k < N; k++) {
            long base = ((long)k * nLam + li) * 8;
            layerN[k] = jread(&matJets[base]);
        }
        jet sine;
        if (sinJets) sine = jread(&sinJets[8 * (long)li]);
        jphase_core(lambdas[li], omegas[li], theta_deg, pol,
                    jread(&n0jets[8 * (long)li]), jread(&nsjets[8 * (long)li]),
                    layerN, thick, N,
                    sinJets ? &sine : NULL, &out[10 * (long)li]);
    }
    free(layerN);
}

/* ── Exported: phase dispersion plus exact thickness derivatives ──────────────
 * Mirrors tmmPhaseThicknessJacobian() in phase.js. Frequency stays the Taylor
 * variable, so each thickness derivative is itself a third-order frequency jet.
 * Zero-thickness layers are retained so derivative indices line up with the
 * caller's design array. Negative and NaN thicknesses are skipped with a zero
 * derivative, matching the point evaluator's base result.
 *
 *   out   : 10        as tmm_phase_one
 *   deriv : 8 × N     [side][quantity][layer], side 0 = r, 1 = t,
 *                     quantity 0 = dPhaseDeg, 1 = dGD, 2 = dGDD, 3 = dTOD
 *
 * The prefix/suffix decomposition cannot carry a rescaling, so if the matrix
 * product overflows, `out` is still filled from the plain path and every entry
 * of `deriv` is set to NaN. */

TMM_EXPORT("tmm_phase_jacobian")
void tmm_phase_jacobian(double lambda, double omega, double theta_deg, int pol,
                        const double *n0jet, const double *nsjet,
                        const double *layerJets, const double *thick, int N,
                        const double *sinJet, double *out, double *deriv) {
    jet n0 = jread(n0jet), ns = jread(nsjet);
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    for (int k = 0; k < N; k++) layerN[k] = jread(&layerJets[8 * k]);
    jet sine;
    if (sinJet) sine = jread(sinJet);
    const jet *sinePtr = sinJet ? &sine : NULL;

    jet wavelength = jwavelength(lambda, omega);
    jet incidentSine = sinePtr ? *sinePtr : jconst(sin(theta_deg * PI / 180.0), 0.0);
    jet incidentCosine = jincident_cos(n0, incidentSine, sinePtr != NULL, theta_deg);
    jet incidentEta = jadmittance(n0, incidentCosine, pol);
    jet substrateCosine = jsnell_cos(n0, incidentSine, ns);
    jet substrateEta = jadmittance(ns, substrateCosine, pol);

    int M = (N > 0 ? N : 1);
    jmat2 *layerM  = (jmat2 *)malloc(sizeof(jmat2) * M);
    jmat2 *layerDM = (jmat2 *)malloc(sizeof(jmat2) * M);
    for (int k = 0; k < N; k++) {
        /* Match jphase_core's skip rule for invalid negative/NaN thicknesses,
         * while retaining the useful derivative of a zero-thickness layer. */
        if (!(thick[k] >= 0.0)) {
            layerM[k] = jidentity();
            layerDM[k] = jzero();
            continue;
        }
        jet cosine = jsnell_cos(n0, incidentSine, layerN[k]);
        jlayer_matrix_dd(layerN[k], thick[k], wavelength, cosine, pol,
                         &layerM[k], &layerDM[k]);
    }

    jmat2 *prefix = (jmat2 *)malloc(sizeof(jmat2) * (N + 1));
    jmat2 *suffix = (jmat2 *)malloc(sizeof(jmat2) * (N + 1));
    prefix[0] = jidentity();
    int overflowed = 0;
    for (int k = 0; k < N; k++) {
        prefix[k + 1] = jmatmul(prefix[k], layerM[k]);
        if (jmatmag(prefix[k + 1]) > MATRIX_RESCALE_THRESHOLD) { overflowed = 1; break; }
    }

    if (overflowed) {
        jphase_core(lambda, omega, theta_deg, pol, n0, ns, layerN, thick, N, sinePtr, out);
        for (long i = 0; i < 8L * N; i++) deriv[i] = NAN;
        free(layerN); free(layerM); free(layerDM); free(prefix); free(suffix);
        return;
    }

    suffix[N] = jidentity();
    for (int k = N - 1; k >= 0; k--) suffix[k] = jmatmul(layerM[k], suffix[k + 1]);

    jcoef coefficients = jcoef_from_matrix(prefix[N], incidentEta, substrateEta, 0.0);
    jphase(coefficients.reflection, &out[0]);
    jphase(coefficients.transmission, &out[5]);

    for (int k = 0; k < N; k++) {
        jmat2 matrixDerivative = jmatmul(jmatmul(prefix[k], layerDM[k]), suffix[k + 1]);
        jet dReflection, dTransmission;
        jcoef_thickness(matrixDerivative, coefficients, incidentEta, substrateEta,
                        &dReflection, &dTransmission);
        const jet sides[2] = { dReflection, dTransmission };
        const jet base[2] = { coefficients.reflection, coefficients.transmission };
        for (int s = 0; s < 2; s++) {
            cx dd[4];
            jderivs(jdiv(sides[s], base[s]), dd);
            deriv[((long)s * 4 + 0) * N + k] = -dd[0].im * 180.0 / PI;
            deriv[((long)s * 4 + 1) * N + k] = dd[1].im;
            deriv[((long)s * 4 + 2) * N + k] = dd[2].im;
            deriv[((long)s * 4 + 3) * N + k] = dd[3].im;
        }
    }

    free(layerN); free(layerM); free(layerDM); free(prefix); free(suffix);
}
