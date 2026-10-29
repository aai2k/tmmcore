/*
 * rta.c : reflectance, transmittance and absorptance of a stack, at one
 * wavelength and over a wavelength grid. Conventions are in tmm_internal.h.
 */

#include "tmm_internal.h"

/* ── Core TMM for one wavelength / polarization ──────────────────────────────
 * Faithful port of the R, T, A kernel in src/js/rta.js. `layers` is N triples [n_re, n_im, d_nm].
 * Writes R,T,A to out[0..2]. A layer whose thickness is not positive, NaN
 * included, is absent, exactly like the JS reference. */

static void tmm_core(double lambda_nm, double theta_deg, int pol,
                     cx n0, cx ns, const double *layers, int N,
                     double *outR, double *outT, double *outA) {
    incidence_t a = incidence(theta_deg);
    cx eta0 = layerAdmittance(n0, incidentCosTheta(n0, a), pol);
    vec2 substrate = substrateVector(ns, snellCosTheta(n0, a, ns), pol);

    /* identity */
    mat2 M;
    M.a = cmk(1.0, 0.0); M.b = cmk(0.0, 0.0);
    M.c = cmk(0.0, 0.0); M.d = cmk(1.0, 0.0);
    double logScale = 0.0;

    for (int i = 0; i < N; i++) {
        cx n = cmk(layers[3 * i + 0], layers[3 * i + 1]);
        double d = layers[3 * i + 2];
        if (!(d > 0.0)) continue;
        cx cosThetaJ = snellCosTheta(n0, a, n);
        mat2 Mj = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
        M = matmul(M, Mj);
        logScale += layerLogScale(n, d, lambda_nm, cosThetaJ) + rescaleMatrix(&M);
    }

    vec2 BC = cmatvec(M, substrate);
    cx B = BC.x, C = BC.y;
    cx eta0B = cmul(eta0, B);
    cx r = cdiv(csub(eta0B, C), cadd(eta0B, C));
    cx t = cmul(cdiv(cmul(cmk(2.0, 0.0), eta0), cadd(eta0B, C)), substrate.x);

    double R = cabs2(r);
    double T = transmittedFlux(substrate) / eta0.re * cabs2(t) * exp(-2.0 * logScale);
    if (T < 0.0) T = 0.0;
    double A = layerAbsorptance(R, T, eta0, r);
    if (A < 0.0) A = 0.0;

    *outR = R; *outT = T; *outA = A;
}

/* ── Exported: R, T, A at one wavelength ─────────────────────────────────────
 * `layers` is N triples [n_re, n_im, d_nm], incident side first.
 * out[0..2] = R, T, A. */

TMM_EXPORT("tmm_rta")
void tmm_rta(double lambda_nm, double theta_deg, int pol,
             double n0_re, double n0_im, double ns_re, double ns_im,
             const double *layers, int N, double *out) {
    tmm_core(lambda_nm, theta_deg, pol,
             cmk(n0_re, n0_im), cmk(ns_re, ns_im), layers, N,
             &out[0], &out[1], &out[2]);
}

/* ── Exported: R, T, A over a wavelength grid, s and p ───────────────────────
 * The whole grid in one call, both polarizations, so a caller crossing a
 * language boundary pays for it once per spectrum.
 *
 * Memory layout (all f64, caller-owned):
 *   lambdas   : nLam                    wavelengths (nm)
 *   theta_deg : nTheta, 1 or nLam       angle of incidence (degrees)
 *   n0arr     : nLam × 2  [re, im]      incident-medium ñ per λ
 *   nsarr     : nLam × 2  [re, im]      substrate ñ per λ
 *   matNK     : N × nLam × 2            per-layer ñ per λ, layout [layer][λ][re,im]
 *   thick     : N                       layer thicknesses (nm)
 * Outputs (each nLam):
 *   outRs,outTs,outAs  (s-pol)   outRp,outTp,outAp  (p-pol) */

TMM_EXPORT("tmm_rta_spectrum")
void tmm_rta_spectrum(const double *lambdas, int nLam,
                      const double *theta_deg, int nTheta,
                      const double *n0arr, const double *nsarr,
                      const double *matNK, const double *thick, int N,
                      double *outRs, double *outTs, double *outAs,
                      double *outRp, double *outTp, double *outAp) {
    const grid_layout g = { lambdas, nLam, theta_deg, nTheta, n0arr, nsarr, matNK, thick, N };
    /* Reusable per-λ layer scratch [n_re, n_im, d] × N. */
    double *layers = (N > 0) ? (double *)malloc(sizeof(double) * 3 * N) : NULL;
    if (N > 0 && !layers) {
        double *outs[6] = { outRs, outTs, outAs, outRp, outTp, outAp };
        for (int o = 0; o < 6; o++) fill_nan(outs[o], nLam);
        return;
    }

    for (int li = 0; li < nLam; li++) {
        double lam = lambdas[li];
        double theta = theta_at(theta_deg, nTheta, li);
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        layers_at(&g, li, layers);
        tmm_core(lam, theta, 0, n0, ns, layers, N, &outRs[li], &outTs[li], &outAs[li]);
        tmm_core(lam, theta, 1, n0, ns, layers, N, &outRp[li], &outTp[li], &outAp[li]);
    }

    if (layers) free(layers);
}
