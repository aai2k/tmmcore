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
 *   p admittance n/cosθ     so r_p = r_s at normal incidence (Macleod §9.7.1)
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

/* An exported kernel that cannot allocate its working state fills its outputs
 * with NaN, so an unwritten buffer is never read back as data. */
static void fill_nan(double *out, long count) {
    for (long i = 0; i < count; i++) out[i] = NAN;
}

/* 1 when every one of `count` allocations succeeded. */
static int all_allocated(void *const *blocks, int count) {
    for (int i = 0; i < count; i++) if (!blocks[i]) return 0;
    return 1;
}

/* ── Snell's law: cosθ_j from incident (n0, θ0) into medium nj ──────────────
 * The transverse invariant is Re(n0) sinθ0 (Macleod 5th ed., §10.2, Eqs.
 * 10.11–10.13): an absorbing incident medium carries a wave whose amplitude
 * falls along the normal only, and every medium's cosθ follows from that one
 * real invariant. Mirrors incidence / snellCosTheta / incidentCosTheta in
 * tmm.js. */

/* sinθ0 and cosθ0 of the angle of incidence. The cosine is taken directly:
 * sqrt(1 − sin²θ0) keeps none of it near grazing incidence. */
typedef struct { cx sin0, cos0; } incidence_t;

static inline incidence_t incidence(double theta_deg) {
    double rad = theta_deg * PI / 180.0;
    incidence_t a;
    a.sin0 = cmk(sin(rad), 0.0);
    a.cos0 = cmk(cos(rad), 0.0);
    return a;
}

/* cos²θ = 1 − (a/nj)² with a = n0 sinθ0, or the same number as
 * ((nj − n0)(nj + n0) + (n0 cosθ0)²)/nj², which keeps its digits at grazing
 * incidence into a medium of the incident index. The second is taken where
 * |nj² − n0²| + 2 (n0 cosθ0)² < n0², where its rounding error is the smaller;
 * never at normal incidence, so cosθ = 1 stays exact there. */
static inline cx snellCosTheta(cx n0, incidence_t a, cx nj) {
    double nr = n0.re;
    double nc = nr * a.cos0.re;
    double q = nc * nc;
    /* Positive only past 45°; |nj² − n0²| ≥ ||nj|² − n0²| rules most media
     * out before (nj − n0)(nj + n0) is formed. */
    double room = nr * nr - 2.0 * q;
    if (room > 0.0 && fabs(cabs2(nj) - nr * nr) < room) {
        cx m = cmul(csub(nj, cmk(nr, 0.0)), cadd(nj, cmk(nr, 0.0)));
        if (cabs2(m) < room * room) return csqrt_(cdiv(cadd(m, cmk(q, 0.0)), cmul(nj, nj)));
    }
    cx sinThetaJ = cdiv(cmul(cmk(n0.re, 0.0), a.sin0), nj);
    return csqrt_(csub(cmk(1.0, 0.0), cmul(sinThetaJ, sinThetaJ)));
}

/* cosθ0 of the incident medium: its own cosine for a transparent medium,
 * otherwise from the same real invariant as the layers. */
static inline cx incidentCosTheta(cx n0, incidence_t a) {
    if (n0.im != 0.0) return snellCosTheta(n0, a, n0);
    return a.cos0;
}

/* ── Layer characteristic matrix (pol: 0 = s, 1 = p) ───────────────────────
 * MAX_IM_DELTA bounds |Im δ| so that cosh and sinh stay finite. Past it the
 * held matrix is the true one divided by e^{|Im δ| − MAX_IM_DELTA}, to within
 * e^{−100} relative; that factor cancels from r and every admittance, and is
 * carried as a log scale for t, T and fields. The bound tests Im δ whatever k
 * is, so a thick lossless layer past the critical angle reaches it too.
 * Mirrors layerPhase / layerLogScale / layerAdmittance / phaseMatrix /
 * layerMatrix in tmm.js. */

/* Phase thickness δ = (2π/λ) n d cosθ (Macleod 5th ed., Eq. 9.2) with Im δ
 * held to ±MAX_IM_DELTA. *excess, when given, is |Im δ| − MAX_IM_DELTA where
 * the bound held it and 0 elsewhere. */
static inline cx layerPhase(cx nj, double dj_nm, double lambda_nm, cx cosTheta_j, double *excess) {
    double k0 = (2.0 * PI) / lambda_nm;
    cx delta = cmul(cmul(nj, cmk(k0 * dj_nm, 0.0)), cosTheta_j);
    if (excess) *excess = fabs(delta.im) > MAX_IM_DELTA ? fabs(delta.im) - MAX_IM_DELTA : 0.0;
    if (delta.im > MAX_IM_DELTA) delta.im = MAX_IM_DELTA;
    else if (delta.im < -MAX_IM_DELTA) delta.im = -MAX_IM_DELTA;
    return delta;
}

/* The log scale layerMatrix leaves out for the same arguments; added to the
 * log scale a product carries for t, like the return of rescaleMatrix. */
static inline double layerLogScale(cx nj, double dj_nm, double lambda_nm, cx cosTheta_j) {
    double excess;
    layerPhase(nj, dj_nm, lambda_nm, cosTheta_j, &excess);
    return excess;
}

/* Tilted admittance: n cosθ for s, n / cosθ for p (Macleod 5th ed., §9.2). */
static inline cx layerAdmittance(cx nj, cx cosTheta_j, int pol) {
    return (pol == 0) ? cmul(nj, cosTheta_j) : cdiv(nj, cosTheta_j);
}

/* M = [[cosδ, −i sinδ/η], [−i η sinδ, cosδ]]; matrices of one admittance
 * compose by adding phases, M(a)·M(b) = M(a + b). */
static inline mat2 phaseMatrix(cx delta, cx eta) {
    cx cosD = ccos_(delta);
    cx sinD = csin_(delta);
    cx negI = cmk(0.0, -1.0);
    cx iSinD_div_eta = cmul(negI, cdiv(sinD, eta));
    cx iEta_sinD     = cmul(negI, cmul(eta, sinD));

    mat2 M;
    M.a = cosD;          M.b = iSinD_div_eta;
    M.c = iEta_sinD;     M.d = cosD;
    return M;
}

/* ── The critical angle ───────────────────────────────────────────────────────
 * Where cosθ is exactly zero δ and the s admittance vanish together and the p
 * admittance diverges; phaseMatrix is 0/0 while its limit is I + d·G with
 * G = [[0, −iP], [−iS, 0]], P = Q/η, S = Qη, Q = dδ/dd, and (P, S) =
 * (2π/λ, 0) for s and (0, (2π/λ) n²) for p. Mirrors the section of the same
 * name in tmm.js. */

typedef struct { cx P, S; } generator_t;

static inline int atCriticalAngle(cx c) { return c.re == 0.0 && c.im == 0.0; }

static inline generator_t criticalGenerator(cx nj, double k0, int pol) {
    generator_t g;
    if (pol == 0) { g.P = cmk(k0, 0.0); g.S = cmk(0.0, 0.0); }
    else {
        cx n2 = cmul(nj, nj);
        g.P = cmk(0.0, 0.0); g.S = cmk(n2.re * k0, n2.im * k0);
    }
    return g;
}

/* [[0, −iP], [−iS, 0]] */
static inline mat2 generatorMatrix(generator_t g) {
    mat2 A;
    A.a = cmk(0.0, 0.0);                 A.b = cmul(cmk(0.0, -1.0), g.P);
    A.c = cmul(cmk(0.0, -1.0), g.S);     A.d = cmk(0.0, 0.0);
    return A;
}

/* I + d·G */
static inline mat2 criticalMatrix(generator_t g, double dj_nm) {
    mat2 M;
    M.a = cmk(1.0, 0.0);                 M.b = cmul(cmk(0.0, -dj_nm), g.P);
    M.c = cmul(cmk(0.0, -dj_nm), g.S);   M.d = cmk(1.0, 0.0);
    return M;
}

static inline mat2 layerMatrix(cx nj, double dj_nm, double lambda_nm, cx cosTheta_j, int pol) {
    if (atCriticalAngle(cosTheta_j)) {
        return criticalMatrix(criticalGenerator(nj, (2.0 * PI) / lambda_nm, pol), dj_nm);
    }
    return phaseMatrix(layerPhase(nj, dj_nm, lambda_nm, cosTheta_j, NULL),
                       layerAdmittance(nj, cosTheta_j, pol));
}

/* ── Boundaries ───────────────────────────────────────────────────────────────
 * The substrate closes the product as [B, C] = M·[1, ηs] (Macleod Eq. 2.113);
 * at its critical angle in p the vector is [0, 1], r keeps only the direction
 * of [B, C], and t = 2η0 bs/(η0B + C) and the transmitted flux Re(cs·conj(bs))
 * vanish. Mirrors substrateVector / transmittedFlux in tmm.js. */

static inline vec2 substrateVector(cx ns, cx cosThetaS, int pol) {
    vec2 v;
    if (pol == 1 && atCriticalAngle(cosThetaS)) {
        v.x = cmk(0.0, 0.0); v.y = cmk(1.0, 0.0);
    } else {
        v.x = cmk(1.0, 0.0); v.y = layerAdmittance(ns, cosThetaS, pol);
    }
    return v;
}

static inline double transmittedFlux(vec2 s) { return cmul(s.y, cconj(s.x)).re; }

/* Absorptance of the layers: the net irradiance entering the front surface
 * (Macleod Eq. 2.120) less the transmitted, per unit incident irradiance,
 *   A = 1 − R − T + 2 (Im η0 / Re η0) Im r,
 * whose last term, the mixed Poynting vector of an absorbing incident medium
 * (Macleod §2.13, Eq. 2.163), is zero for a transparent one. Mirrors
 * layerAbsorptance in tmm.js. */
static inline double layerAbsorptance(double R, double T, cx eta0, cx r) {
    return 1.0 - R - T + 2.0 * (eta0.im / eta0.re) * r.im;
}

/* ── Core TMM for one wavelength / polarization ──────────────────────────────
 * Faithful port of tmm() in tmm.js. `layers` is N triples [n_re, n_im, d_nm].
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
    if (N > 0 && !layers) {
        double *outs[6] = { outRs, outTs, outAs, outRp, outTp, outAp };
        for (int o = 0; o < 6; o++) fill_nan(outs[o], nLam);
        return;
    }

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
    vec2 *substrateV = (vec2 *)malloc(sizeof(vec2) * cells);
    if (!M || !logScale || !eta0v || !substrateV) {
        free(M); free(logScale); free(eta0v); free(substrateV);
        return 0;
    }

    incidence_t a = incidence(theta_deg);

    for (int li = 0; li < nLam; li++) {
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
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
                    cx cosThetaJ = snellCosTheta(n0, a, n);
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

    free(M); free(logScale); free(eta0v); free(substrateV);
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
    incidence_t angle;
    int topSet;
    double *lam;        /* nLam */
    cx *n0;             /* nLam */
    cx *ng;             /* nLam, set_top */
    cx *cosThetaG;      /* nLam, set_top (polarization-independent) */
    cx *eta0;           /* 2·nLam, pol-major */
    vec2 *substrate;    /* 2·nLam */
    mat2 *Mb;           /* 2·nLam completed-stack product */
    double *logScaleB;  /* 2·nLam */
} growing_eval;

TMM_EXPORT("tmm_growing_eval_free")
void tmm_growing_eval_free(growing_eval *h) {
    if (!h) return;
    free(h->lam); free(h->n0); free(h->ng); free(h->cosThetaG);
    free(h->eta0); free(h->substrate); free(h->Mb); free(h->logScaleB);
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
    h->substrate = (vec2 *)malloc(sizeof(vec2) * 2 * nLam);
    h->Mb = (mat2 *)malloc(sizeof(mat2) * 2 * nLam);
    h->logScaleB = (double *)malloc(sizeof(double) * 2 * nLam);
    void *const blocks[] = { h->lam, h->n0, h->ng, h->cosThetaG,
                             h->eta0, h->substrate, h->Mb, h->logScaleB };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        tmm_growing_eval_free(h);
        return 0;
    }

    h->angle = incidence(theta_deg);

    for (int li = 0; li < nLam; li++) {
        h->lam[li] = lambdas[li];
        cx n0 = cmk(n0arr[2 * li + 0], n0arr[2 * li + 1]);
        cx ns = cmk(nsarr[2 * li + 0], nsarr[2 * li + 1]);
        h->n0[li] = n0;
        cx cosTheta0 = incidentCosTheta(n0, h->angle);
        for (int pol = 0; pol < 2; pol++) {
            size_t at = (size_t)pol * nLam + li;
            h->eta0[at] = layerAdmittance(n0, cosTheta0, pol);
            h->substrate[at] = substrateVector(ns, snellCosTheta(n0, h->angle, ns), pol);
            mat2 M;
            M.a = cmk(1.0, 0.0); M.b = cmk(0.0, 0.0);
            M.c = cmk(0.0, 0.0); M.d = cmk(1.0, 0.0);
            double logScale = 0.0;
            for (int k = 0; k < NB; k++) {
                double d = thick[k];
                if (!(d > 0.0)) continue;
                cx n = cmk(matNK[((size_t)k * nLam + li) * 2 + 0],
                           matNK[((size_t)k * nLam + li) * 2 + 1]);
                cx cosThetaJ = snellCosTheta(n0, h->angle, n);
                M = matmul(M, layerMatrix(n, d, lambdas[li], cosThetaJ, pol));
                logScale += layerLogScale(n, d, lambdas[li], cosThetaJ) + rescaleMatrix(&M);
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
        h->cosThetaG[li] = snellCosTheta(h->n0[li], h->angle, ng);
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

/* ── Shared pieces of the derivative kernels ─────────────────────────────────
 * Port of the section of the same name in tmm.js. Every prefix
 * Pre[j] = M_0·…·M_{j−1} and suffix Post[j] = M_j·…·M_{N−1}·[1, ηs] is stored
 * as a mantissa times 2^e with the integer e kept alongside; a power of two
 * scales a double exactly, and every derivative, one prefix times one suffix
 * over den = η0 B + C, takes the exponent e_pre + e_post − e_0. */

#define BINARY_RESCALE_LIMIT 18446744073709551616.0   /* 2^64 */

static inline cx cscale(cx a, double s) { return cmk(a.re * s, a.im * s); }

/* Divides the values in place by a power of two when the largest part exceeds
 * BINARY_RESCALE_LIMIT, and returns that power's exponent. */
static int rescaleBinary(cx **values, int count) {
    double largest = 0.0;
    for (int i = 0; i < count; i++)
        largest = fmax(largest, fmax(fabs(values[i]->re), fabs(values[i]->im)));
    if (!(largest > BINARY_RESCALE_LIMIT) || isinf(largest)) return 0;
    int exponent = (int)floor(log2(largest));
    double factor = ldexp(1.0, -exponent);
    for (int i = 0; i < count; i++) {
        values[i]->re *= factor;
        values[i]->im *= factor;
    }
    return exponent;
}
static int rescaleMat(mat2 *M) {
    cx *e[4] = { &M->a, &M->b, &M->c, &M->d };
    return rescaleBinary(e, 4);
}
static int rescaleVec(vec2 *v) {
    cx *e[2] = { &v->x, &v->y };
    return rescaleBinary(e, 2);
}

/* The pair (P, S) = (Q/η, Qη) of a medium, at its critical-angle limit where
 * cosθ = 0. Mirrors generator in tmm.js. */
static inline generator_t generatorOf(cx nj, cx cosTheta, cx Q, double k0, int pol) {
    if (atCriticalAngle(cosTheta)) return criticalGenerator(nj, k0, pol);
    cx eta = layerAdmittance(nj, cosTheta, pol);
    generator_t g;
    g.P = cdiv(Q, eta);
    g.S = cmul(Q, eta);
    return g;
}

/* Everything the derivative kernels share at one (λ, θ, pol). `layers` is N
 * triples [n_re, n_im, d] used as given, with no d > 0 filter, for index parity
 * with the caller's design. A thickness that is not a number ≥ 0 leaves the
 * layer out: identity matrix, and Q = P = S = 0 so its derivative is zero. Q is
 * dδ/dd in full: past the bound the held matrix times e^{excess} moves with d
 * as the true one does, and Q with the held δ is its derivative. `logScale`
 * sums the excesses, which t carries. Mirrors derivativeLayer and
 * derivativeStack in tmm.js. */
typedef struct {
    int N;
    double k0;
    incidence_t angle;
    cx eta0;
    vec2 substrate;
    cx *n, *cosTheta, *delta, *eta, *Q;
    generator_t *G;
    int *critical;
    double *thickness, *excess;
    double logScale;
    mat2 *M, *Pre;
    vec2 *Post;
    int *preExp, *postExp;
} deriv_stack;

/* The arguments every derivative kernel shares. */
typedef struct {
    double lambda_nm, theta_deg;
    int pol;
    cx n0, ns;
    const double *layers;
    int N;
} deriv_request;

static void deriv_stack_free(deriv_stack *s) {
    free(s->n); free(s->cosTheta); free(s->delta); free(s->eta); free(s->Q);
    free(s->G); free(s->critical);
    free(s->thickness); free(s->excess); free(s->M); free(s->Pre); free(s->Post);
    free(s->preExp); free(s->postExp);
}

/* The matrix of `thickness` of layer k whose phase over it is `delta`. */
static inline mat2 partMatrix(const deriv_stack *s, int k, cx delta, double thickness) {
    return s->critical[k] ? criticalMatrix(s->G[k], thickness) : phaseMatrix(delta, s->eta[k]);
}

/* Returns 1, or 0 with nothing left allocated when the working state could not
 * be allocated. */
static int deriv_stack_build(const deriv_request *q, deriv_stack *out) {
    deriv_stack s;
    int N = q->N, pol = q->pol;
    double lambda_nm = q->lambda_nm;
    const double *layers = q->layers;
    cx n0 = q->n0;
    size_t count = (size_t)(N > 0 ? N : 1);
    s.N = N;
    s.k0 = (2.0 * PI) / lambda_nm;
    s.angle = incidence(q->theta_deg);
    s.eta0 = layerAdmittance(n0, incidentCosTheta(n0, s.angle), pol);
    s.substrate = substrateVector(q->ns, snellCosTheta(n0, s.angle, q->ns), pol);
    s.n         = (cx *)         malloc(sizeof(cx) * count);
    s.cosTheta  = (cx *)         malloc(sizeof(cx) * count);
    s.delta     = (cx *)         malloc(sizeof(cx) * count);
    s.eta       = (cx *)         malloc(sizeof(cx) * count);
    s.Q         = (cx *)         malloc(sizeof(cx) * count);
    s.G         = (generator_t *)malloc(sizeof(generator_t) * count);
    s.critical  = (int *)        malloc(sizeof(int) * count);
    s.thickness = (double *)     malloc(sizeof(double) * count);
    s.excess    = (double *)     malloc(sizeof(double) * count);
    s.M         = (mat2 *)       malloc(sizeof(mat2) * count);
    s.Pre       = (mat2 *)       malloc(sizeof(mat2) * ((size_t)N + 1));
    s.Post      = (vec2 *)       malloc(sizeof(vec2) * ((size_t)N + 1));
    s.preExp    = (int *)        malloc(sizeof(int) * ((size_t)N + 1));
    s.postExp   = (int *)        malloc(sizeof(int) * ((size_t)N + 1));
    void *const blocks[] = { s.n, s.cosTheta, s.delta, s.eta, s.Q, s.G, s.critical,
                             s.thickness, s.excess, s.M, s.Pre, s.Post, s.preExp, s.postExp };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        deriv_stack_free(&s);
        return 0;
    }

    s.logScale = 0.0;
    for (int k = 0; k < N; k++) {
        double d = layers[3 * k + 2];
        int present = d >= 0.0;
        if (!present) d = 0.0;
        s.n[k] = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        s.thickness[k] = d;
        s.cosTheta[k] = snellCosTheta(n0, s.angle, s.n[k]);
        s.delta[k] = layerPhase(s.n[k], d, lambda_nm, s.cosTheta[k], &s.excess[k]);
        s.logScale += s.excess[k];
        s.eta[k] = layerAdmittance(s.n[k], s.cosTheta[k], pol);
        cx Q = cmul(cmul(s.n[k], cmk(s.k0, 0.0)), s.cosTheta[k]);   /* (2π/λ) n cosθ */
        s.Q[k] = present ? Q : cmk(0.0, 0.0);
        if (present) {
            s.G[k] = generatorOf(s.n[k], s.cosTheta[k], s.Q[k], s.k0, pol);
        } else {
            s.G[k].P = cmk(0.0, 0.0); s.G[k].S = cmk(0.0, 0.0);
        }
        s.critical[k] = atCriticalAngle(s.cosTheta[k]);
        s.M[k] = partMatrix(&s, k, s.delta[k], d);
    }

    s.Pre[0].a = cmk(1.0, 0.0); s.Pre[0].b = cmk(0.0, 0.0);
    s.Pre[0].c = cmk(0.0, 0.0); s.Pre[0].d = cmk(1.0, 0.0);
    s.preExp[0] = 0;
    for (int j = 0; j < N; j++) {
        s.Pre[j + 1] = matmul(s.Pre[j], s.M[j]);
        s.preExp[j + 1] = s.preExp[j] + rescaleMat(&s.Pre[j + 1]);
    }
    s.Post[N] = s.substrate;
    s.postExp[N] = 0;
    for (int j = N - 1; j >= 0; j--) {
        s.Post[j] = cmatvec(s.M[j], s.Post[j + 1]);
        s.postExp[j] = s.postExp[j + 1] + rescaleVec(&s.Post[j]);
    }
    *out = s;
    return 1;
}

/* R, T, A from [B, C] = Post[0] (Macleod Eqs. 2.123–2.125), with r, the true
 * t, b = B/den and c = C/den, and g = 2 Im η0 / Re η0, the weight of Im r in
 * the absorptance. Only t carries the exponent of [B, C] and the log scale of
 * the held layers. */
typedef struct {
    double R, T, A, Tfac, g;
    cx eta0, r, t, b, c, inverseDen;
    int exponent;
} stack_response;

static stack_response stack_response_of(const deriv_stack *s) {
    stack_response o;
    cx B = s->Post[0].x, C = s->Post[0].y;
    cx den = cadd(cmul(s->eta0, B), C);
    o.inverseDen = cdiv(cmk(1.0, 0.0), den);
    o.r = cdiv(csub(cmul(s->eta0, B), C), den);
    o.t = cscale(cmul(cdiv(cmul(cmk(2.0, 0.0), s->eta0), den), s->substrate.x),
                 ldexp(1.0, -s->postExp[0]));
    if (s->logScale > 0.0) o.t = cscale(o.t, exp(-s->logScale));
    o.Tfac = transmittedFlux(s->substrate) / s->eta0.re;
    o.R = cabs2(o.r);
    o.T = o.Tfac * cabs2(o.t); if (o.T < 0.0) o.T = 0.0;
    o.A = layerAbsorptance(o.R, o.T, s->eta0, o.r); if (o.A < 0.0) o.A = 0.0;
    o.g = 2.0 * (s->eta0.im / s->eta0.re);
    o.eta0 = s->eta0;
    o.b = cmul(B, o.inverseDen);
    o.c = cmul(C, o.inverseDen);
    o.exponent = s->postExp[0];
    return o;
}

/* A derivative of [B, C], a mantissa and its exponent, divided by den. */
static inline vec2 per_den(const stack_response *s, vec2 v, int exponent) {
    double scale = ldexp(1.0, exponent - s->exponent);
    vec2 o;
    o.x = cscale(cmul(v.x, s->inverseDen), scale);
    o.y = cscale(cmul(v.y, s->inverseDen), scale);
    return o;
}

/* {dR,dT,dA} from (u, w) = d[B,C]/den: dr = 2η0 (c·u − b·w), dt = −t (η0·u + w),
 * the chain rule with f = 2η0/den² and den divided out once; dA = −(dR + dT)
 * + g Im dr, as layerAbsorptance. */
static void response_derivative(const stack_response *s, vec2 uw, double *o) {
    cx dr = cmul(cmul(cmk(2.0, 0.0), s->eta0), csub(cmul(s->c, uw.x), cmul(s->b, uw.y)));
    double dR = 2.0 * (cmul(cconj(s->r), dr)).re;
    cx dt = cmul(cmk(-1.0, 0.0), cmul(s->t, cadd(cmul(s->eta0, uw.x), uw.y)));
    double dT = s->Tfac * 2.0 * (cmul(cconj(s->t), dt)).re;
    o[0] = dR; o[1] = dT; o[2] = -(dR + dT) + s->g * dr.im;
}

/* dM/dd = [[ −Q sinδ, −iP cosδ ], [ −iS cosδ, −Q sinδ ]], P = Q/η, S = Qη */
static inline mat2 layerMatrixDerivative(cx delta, cx Q, generator_t G) {
    cx cD = ccos_(delta), sD = csin_(delta);
    cx negI = cmk(0.0, -1.0), neg1 = cmk(-1.0, 0.0);
    mat2 dM;
    dM.a = cmul(neg1, cmul(Q, sD));
    dM.b = cmul(negI, cmul(G.P, cD));
    dM.c = cmul(negI, cmul(G.S, cD));
    dM.d = cmul(neg1, cmul(Q, sD));
    return dM;
}

/* ── Exported: analytic thickness Jacobian for one (λ, θ, pol) ────────────────
 * Port of tmmThicknessJacobian() in tmm.js. Returns the exact analytic
 * dR/dd_k, dT/dd_k, dA/dd_k for every layer at one sample.
 *
 * `layers` is N triples [n_re, n_im, d] used as given (see deriv_stack).
 * Outputs (each length N): dRdd, dTdd, dAdd. Also writes base R,T,A to base[0..2]. */

TMM_EXPORT("tmm_jacobian")
void tmm_jacobian(double lambda_nm, double theta_deg, int pol,
                  double n0_re, double n0_im, double ns_re, double ns_im,
                  const double *layers, int N,
                  double *dRdd, double *dTdd, double *dAdd, double *base) {
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    deriv_stack s;
    if (!deriv_stack_build(&q, &s)) {
        fill_nan(base, 3); fill_nan(dRdd, N); fill_nan(dTdd, N); fill_nan(dAdd, N);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    for (int k = 0; k < N; k++) {
        mat2 dMk = layerMatrixDerivative(s.delta[k], s.Q[k], s.G[k]);
        vec2 dV = cmatvec(s.Pre[k], cmatvec(dMk, s.Post[k + 1]));
        double o[3];
        response_derivative(&resp, per_den(&resp, dV, s.preExp[k] + s.postExp[k + 1]), o);
        dRdd[k] = o[0]; dTdd[k] = o[1]; dAdd[k] = o[2];
    }

    deriv_stack_free(&s);
}

/* ── Analytic needle P-function scan ─────────────────────────────────────────
 * Port of tmmNeedleScan() in tmm.js : the d→0 limit of
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

TMM_EXPORT("tmm_needle_scan")
void tmm_needle_scan(double lambda_nm, double theta_deg, int pol,
                     double n0_re, double n0_im, double ns_re, double ns_im,
                     const double *layers, int N,
                     const double *candNs, int nCand,
                     const double *fracs, int nFrac,
                     double *base, double *gaps, double *intra) {
    cx n0 = cmk(n0_re, n0_im);
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    deriv_stack s;
    mat2 *Ac = (mat2 *)malloc(sizeof(mat2) * (nCand > 0 ? nCand : 1));
    if (!Ac || !deriv_stack_build(&q, &s)) {
        free(Ac);
        fill_nan(base, 3);
        fill_nan(gaps, (long)(N + 1) * nCand * 3);
        if (nFrac > 0) fill_nan(intra, (long)N * nFrac * nCand * 3);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    /* Needle derivative-matrix A = G per candidate: [[0,−iQ/ηₐ],[−iηₐQ,0]]. */
    for (int c = 0; c < nCand; c++) {
        cx nA = cmk(candNs[2 * c + 0], candNs[2 * c + 1]);
        cx cthA = snellCosTheta(n0, s.angle, nA);
        cx Q = cmul(cmul(nA, cmk(s.k0, 0.0)), cthA);
        Ac[c] = generatorMatrix(generatorOf(nA, cthA, Q, s.k0, pol));
    }

    /* Gaps: pos = 0..N, every candidate. */
    for (int pos = 0; pos <= N; pos++) {
        for (int c = 0; c < nCand; c++) {
            vec2 dV = cmatvec(s.Pre[pos], cmatvec(Ac[c], s.Post[pos]));
            response_derivative(&resp, per_den(&resp, dV, s.preExp[pos] + s.postExp[pos]),
                                &gaps[((long)pos * nCand + c) * 3]);
        }
    }

    /* Intra-layer positions. Short of the bound the host splits into halves
     * whose product is its own matrix: the front half keeps its own phase and
     * the back half takes the rest. Past it each half is held on its own, and
     * the pair takes e^{front excess + back excess − host excess}, between
     * e^{−50} and 1, the ratio of what the halves and the host leave out. */
    if (nFrac > 0) {
        for (int k = 0; k < N; k++) {
            for (int fi = 0; fi < nFrac; fi++) {
                double frontExcess, backExcess, factor = 1.0;
                cx front = layerPhase(s.n[k], fracs[fi] * s.thickness[k], lambda_nm,
                                      s.cosTheta[k], &frontExcess);
                cx back = csub(s.delta[k], front);
                if (s.excess[k] > 0.0) {
                    back = layerPhase(s.n[k], (1.0 - fracs[fi]) * s.thickness[k], lambda_nm,
                                      s.cosTheta[k], &backExcess);
                    factor = exp(frontExcess + backExcess - s.excess[k]);
                }
                mat2 preIn  = matmul(s.Pre[k], partMatrix(&s, k, front, fracs[fi] * s.thickness[k]));
                vec2 postIn = cmatvec(partMatrix(&s, k, back, (1.0 - fracs[fi]) * s.thickness[k]),
                                      s.Post[k + 1]);
                for (int c = 0; c < nCand; c++) {
                    vec2 dV = cmatvec(preIn, cmatvec(Ac[c], postIn));
                    if (factor != 1.0) { dV.x = cscale(dV.x, factor); dV.y = cscale(dV.y, factor); }
                    long off = (((long)k * nFrac + fi) * nCand + c) * 3;
                    response_derivative(&resp, per_den(&resp, dV, s.preExp[k] + s.postExp[k + 1]),
                                        &intra[off]);
                }
            }
        }
    }

    free(Ac);
    deriv_stack_free(&s);
}

/* d²M/dd² = Q²·[[ −cosδ, i sinδ/η ], [ i η sinδ, −cosδ ]], with Q²/η = QP and
 * Q²η = QS */
static inline mat2 layerMatrixSecondDerivative(cx delta, cx Q, generator_t G) {
    cx cD = ccos_(delta), sD = csin_(delta);
    cx Q2 = cmul(Q, Q);
    cx posI = cmk(0.0, 1.0), neg1 = cmk(-1.0, 0.0);
    mat2 d2M;
    d2M.a = cmul(neg1, cmul(Q2, cD));
    d2M.b = cmul(posI, cmul(cmul(Q, G.P), sD));
    d2M.c = cmul(posI, cmul(cmul(Q, G.S), sD));
    d2M.d = cmul(neg1, cmul(Q2, cD));
    return d2M;
}

/* d²R, d²T and d²A from (uᵢ,wᵢ), (uⱼ,wⱼ) and the mixed partial (U,W), each a
 * derivative of [B,C] divided by den. With gₖ = η0 uₖ + wₖ and G = η0 U + W:
 *   d²r = 2η0 (wᵢuⱼ + cU − uᵢwⱼ − bW) − 2 drⱼ gᵢ,   d²t = −t G + 2t gᵢgⱼ
 * Port of responseSecondDerivative in tmm.js. */
typedef struct { double R, T, A; } second_derivative;

static second_derivative response_second_derivative(const stack_response *s,
                                                    vec2 fi, vec2 fj, vec2 second) {
    second_derivative o;
    cx twoEta0 = cmul(cmk(2.0, 0.0), s->eta0);
    cx neg1 = cmk(-1.0, 0.0);
    cx dr_i = cmul(twoEta0, csub(cmul(s->c, fi.x), cmul(s->b, fi.y)));
    cx dr_j = cmul(twoEta0, csub(cmul(s->c, fj.x), cmul(s->b, fj.y)));
    cx g_i = cadd(cmul(s->eta0, fi.x), fi.y);
    cx g_j = cadd(cmul(s->eta0, fj.x), fj.y);
    cx innerR = csub(cadd(cmul(fi.y, fj.x), cmul(s->c, second.x)),
                     cadd(cmul(fi.x, fj.y), cmul(s->b, second.y)));
    cx d2r = csub(cmul(twoEta0, innerR), cmul(cmul(cmk(2.0, 0.0), dr_j), g_i));
    o.R = 2.0 * ((cmul(cconj(dr_i), dr_j)).re + (cmul(cconj(s->r), d2r)).re);
    cx dt_i = cmul(neg1, cmul(s->t, g_i));
    cx dt_j = cmul(neg1, cmul(s->t, g_j));
    cx G = cadd(cmul(s->eta0, second.x), second.y);
    cx d2t = cadd(cmul(neg1, cmul(s->t, G)), cmul(cmul(cmk(2.0, 0.0), s->t), cmul(g_i, g_j)));
    o.T = s->Tfac * 2.0 * ((cmul(cconj(dt_i), dt_j)).re + (cmul(cconj(s->t), d2t)).re);
    o.A = -(o.R + o.T) + s->g * d2r.im;
    return o;
}

/* ── Analytic thickness-Hessian kernel ───────────
 * Port of tmmThicknessHessian() in tmm.js : the EXACT analytic second
 * derivatives ∂²{R,T,A}/∂dᵢ∂dⱼ (full N×N symmetric) plus the first
 * derivatives, at one (λ,θ,pol). The middle product M_{i+1}···M_{j-1} carries
 * its own binary exponent, like the prefixes and suffixes.
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
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    deriv_stack s;
    /* Per-layer first-derivative pieces: dM[k], v[k] = dMₖ·Post[k+1] (exponent
     * postExp[k+1]) and d[B,C]/den. */
    mat2 *dM    = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    vec2 *v     = (vec2 *)malloc(sizeof(vec2) * (N > 0 ? N : 1));
    vec2 *first = (vec2 *)malloc(sizeof(vec2) * (N > 0 ? N : 1));
    if (!dM || !v || !first || !deriv_stack_build(&q, &s)) {
        free(dM); free(v); free(first);
        fill_nan(base, 3); fill_nan(dRdd, N); fill_nan(dTdd, N); fill_nan(dAdd, N);
        fill_nan(d2Rdd, (long)N * N); fill_nan(d2Tdd, (long)N * N); fill_nan(d2Add, (long)N * N);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    for (int k = 0; k < N; k++) {
        dM[k] = layerMatrixDerivative(s.delta[k], s.Q[k], s.G[k]);
        v[k] = cmatvec(dM[k], s.Post[k + 1]);
        first[k] = per_den(&resp, cmatvec(s.Pre[k], v[k]), s.preExp[k] + s.postExp[k + 1]);
        double o[3];
        response_derivative(&resp, first[k], o);
        dRdd[k] = o[0]; dTdd[k] = o[1]; dAdd[k] = o[2];
    }

    mat2 I; I.a = cmk(1.0, 0.0); I.b = cmk(0.0, 0.0); I.c = cmk(0.0, 0.0); I.d = cmk(1.0, 0.0);
    for (int i = 0; i < N; i++) {
        mat2 Wmat_i = matmul(s.Pre[i], dM[i]); /* Pre[i]·dMᵢ (used for j>i) */
        mat2 Cmid = I;                          /* M_{i+1}···M_{j-1}, empty at j=i+1 */
        int midExp = 0;
        for (int j = i; j < N; j++) {
            vec2 second;
            if (j == i) {
                mat2 d2M = layerMatrixSecondDerivative(s.delta[i], s.Q[i], s.G[i]);
                second = per_den(&resp, cmatvec(s.Pre[i], cmatvec(d2M, s.Post[i + 1])),
                                 s.preExp[i] + s.postExp[i + 1]);
            } else {
                second = per_den(&resp, cmatvec(Wmat_i, cmatvec(Cmid, v[j])),
                                 s.preExp[i] + midExp + s.postExp[j + 1]);
            }
            second_derivative h = response_second_derivative(&resp, first[i], first[j], second);
            double d2Rij = h.R, d2Tij = h.T, d2Aij = h.A;
            d2Rdd[(long)i * N + j] = d2Rdd[(long)j * N + i] = d2Rij;
            d2Tdd[(long)i * N + j] = d2Tdd[(long)j * N + i] = d2Tij;
            d2Add[(long)i * N + j] = d2Add[(long)j * N + i] = d2Aij;
            if (j >= i + 1) {                   /* advance middle: include M_j */
                Cmid = matmul(Cmid, s.M[j]);
                midExp += rescaleMat(&Cmid);
            }
        }
    }

    free(dM); free(v); free(first);
    deriv_stack_free(&s);
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

/* Past the limit the imaginary part is held at it at every order and the real
 * part, which carries the phase, is kept: the matrix is then the true one over
 * a real factor, which callers carry for |t|. Mirrors jetClampImaginary. */
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
static void jmat_scale(jmat2 *M, double factor) {
    M->a = jscale(M->a, factor);
    M->b = jscale(M->b, factor);
    M->c = jscale(M->c, factor);
    M->d = jscale(M->d, factor);
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
    jmat_scale(M, 1.0 / scale);
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
/* Divides every order by a power of two, which is exact, when the largest part
 * at any order exceeds BINARY_RESCALE_LIMIT, and returns that power's exponent.
 * Mirrors rescaleMatrixBinary in phase.js. */
static int jrescale_binary(jmat2 *M) {
    double magnitude = jmatmag(*M);
    if (!(magnitude > BINARY_RESCALE_LIMIT) || isinf(magnitude)) return 0;
    int exponent = (int)floor(log2(magnitude));
    jmat_scale(M, ldexp(1.0, -exponent));
    return exponent;
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

/* Same real invariant Re(n0) sinθ0 as snellCosTheta, and the same choice of
 * form made on the values, in jet arithmetic. `cos0` is NULL when the incident
 * sine arrives as a jet, and then the form 1 − (a/n)² is used. Mirrors
 * snellCosine in phase.js. */
static inline jet jsnell_cos(jet n0, jet sin0, const jet *cos0, jet nj) {
    jet realIndex = jreal(n0);
    if (cos0) {
        double nr = n0.c[0].re;
        double nc = nr * cos0->c[0].re;
        double q = nc * nc;
        double room = nr * nr - 2.0 * q;
        if (room > 0.0 && fabs(cabs2(nj.c[0]) - nr * nr) < room) {
            jet m = jmul(jsub(nj, realIndex), jadd(nj, realIndex));
            if (m.c[0].re * m.c[0].re + m.c[0].im * m.c[0].im < room * room) {
                jet normal = jmul(realIndex, *cos0);
                return jsqrt_j(jdiv(jadd(m, jmul(normal, normal)), jmul(nj, nj)));
            }
        }
    }
    jet s = jdiv(jmul(realIndex, sin0), nj);
    return jsqrt_j(jsub(jconst(1.0, 0.0), jmul(s, s)));
}
/* cosθ0 of the incident medium: from the same real invariant as the layers
 * when the medium absorbs at any order, otherwise its own cosine, or
 * sqrt(1 − sin²θ0) from the sine jet when the caller supplied one. Mirrors
 * incidentCosine in phase.js. */
static inline jet jincident_cos(jet n0, jet incidentSine, const jet *cos0) {
    if (jhas_imag(n0)) return jsnell_cos(n0, incidentSine, cos0, n0);
    return cos0 ? *cos0 : jsqrt_j(jsub(jconst(1.0, 0.0), jmul(incidentSine, incidentSine)));
}
static inline jet jadmittance(jet n, jet cosv, int pol) {
    return (pol == 0) ? jmul(n, cosv) : jdiv(n, cosv);
}

/* |Im δ| − MAX_IM_DELTA where the bound held the phase, else 0; see
 * heldExcess in phase.js. */
static inline double jheld_excess(jet rawPhase, jet phase) {
    return rawPhase.c[0].im == phase.c[0].im ? 0.0 : fabs(rawPhase.c[0].im) - MAX_IM_DELTA;
}

static jmat2 jlayer_matrix(jet index, double thickness, jet wavelength, jet cosine, int pol,
                           double *excess) {
    jet rawPhase = jscale(jdiv(jmul(index, cosine), wavelength), 2.0 * PI * thickness);
    jet phase = jclampim(rawPhase, MAX_IM_DELTA);
    *excess = jheld_excess(rawPhase, phase);
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

/* The full phase derivative with the held sines and cosines is the derivative
 * of the held matrix times its factor, as in layerMatrixWithThicknessDerivative
 * in phase.js. */
static void jlayer_matrix_dd(jet index, double thickness, jet wavelength, jet cosine, int pol,
                             jmat2 *M, jmat2 *dM, double *excess) {
    jet phasePerUnit = jscale(jdiv(jmul(index, cosine), wavelength), 2.0 * PI);
    jet rawPhase = jscale(phasePerUnit, thickness);
    jet phase = jclampim(rawPhase, MAX_IM_DELTA);
    *excess = jheld_excess(rawPhase, phase);
    jet phaseDerivative = phasePerUnit;
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

/* `transmissionScale` restores the factor a rescaled matrix gave up: it
 * cancels from r but not from t. */
static jcoef jcoef_from_matrix(jmat2 M, jet incidentEta, jet substrateEta, double transmissionScale) {
    jet boundaryB = jadd(M.a, jmul(M.b, substrateEta));
    jet boundaryC = jadd(M.c, jmul(M.d, substrateEta));
    jet incidentB = jmul(incidentEta, boundaryB);
    jcoef o;
    o.denominator = jadd(incidentB, boundaryC);
    o.reflection = jdiv(jsub(incidentB, boundaryC), o.denominator);
    o.transmission = jdiv(jscale(incidentEta, 2.0), o.denominator);
    if (transmissionScale != 1.0) o.transmission = jscale(o.transmission, transmissionScale);
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

/* The incident sine as a jet and, at a fixed angle of incidence, the constant
 * cosθ0; `cosine` is NULL when the caller gave the sine as a jet. Mirrors
 * incidentAngle in phase.js. */
typedef struct { jet sine, cos0; const jet *cosine; } jangle;

static inline void jincidence(double theta_deg, const jet *sinJet, jangle *a) {
    if (sinJet) {
        a->sine = *sinJet;
        a->cosine = NULL;
    } else {
        a->sine = jconst(sin(theta_deg * PI / 180.0), 0.0);
        a->cos0 = jconst(cos(theta_deg * PI / 180.0), 0.0);
        a->cosine = &a->cos0;
    }
}

static void jphase_core(double lambda, double omega, double theta_deg, int pol,
                        jet n0, jet ns, const jet *layerN, const double *thick, int N,
                        const jet *sinJet, double *out10) {
    jet wavelength = jwavelength(lambda, omega);
    jangle angle;
    jincidence(theta_deg, sinJet, &angle);
    jet incidentEta = jadmittance(n0, jincident_cos(n0, angle.sine, angle.cosine), pol);
    jet substrateCosine = jsnell_cos(n0, angle.sine, angle.cosine, ns);
    jet substrateEta = jadmittance(ns, substrateCosine, pol);

    jmat2 M = jidentity();
    double logScale = 0.0, transmissionLogScale = 0.0;
    for (int k = 0; k < N; k++) {
        if (!(thick[k] > 0.0)) continue;
        jet cosine = jsnell_cos(n0, angle.sine, angle.cosine, layerN[k]);
        double excess;
        M = jmatmul(M, jlayer_matrix(layerN[k], thick[k], wavelength, cosine, pol, &excess));
        logScale += jrescale(&M, MATRIX_RESCALE_THRESHOLD);
        transmissionLogScale += excess;
    }
    jcoef coefficients = jcoef_from_matrix(M, incidentEta, substrateEta, exp(-logScale));
    jphase(coefficients.reflection, &out10[0]);
    jphase(coefficients.transmission, &out10[5]);
    /* The phase of t is the held one's; |t|² takes what the held layers leave
     * out, as the logScale of coefficientPhaseDispersion in phase.js. */
    if (transmissionLogScale > 0.0) out10[9] *= exp(-2.0 * transmissionLogScale);
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
    if (!layerN) { fill_nan(out, 10); return; }
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
    if (!layerN) { fill_nan(out, 10L * nLam); return; }
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
 *   deriv : 10 × N    [side][quantity][layer], side 0 = r, 1 = t,
 *                     quantity 0 = dPhaseDeg, 1 = dGD, 2 = dGDD, 3 = dTOD,
 *                     4 = d(ln |coefficient|²)/dd, the relative intensity
 *                     derivative, which the amplitude ratio in ellipsometry
 *                     needs and which costs nothing here: it is the real part
 *                     of the same logarithmic derivative the phase is the
 *                     imaginary part of
 *
 * The prefix and suffix products carry binary exponents, as in the derivative
 * kernels above, so opaque stacks whose products leave double range still get
 * their derivatives. */

/* Scratch for one wavelength of the phase Jacobian: per-layer matrices with
 * their thickness derivatives, and the prefix/suffix products with their
 * exponents. Allocated once per call by the entry points below, so the batched
 * one reuses it across the whole grid. */
typedef struct {
    jmat2 *layerM, *layerDM, *prefix, *suffix;
    int *prefixExp, *suffixExp;
} jphase_jacobian_scratch;

static void jphase_jacobian_free(jphase_jacobian_scratch s) {
    free(s.layerM); free(s.layerDM); free(s.prefix); free(s.suffix);
    free(s.prefixExp); free(s.suffixExp);
}

/* Returns 1, or 0 with nothing left allocated. */
static int jphase_jacobian_alloc(int N, jphase_jacobian_scratch *out) {
    int M = (N > 0 ? N : 1);
    jphase_jacobian_scratch s;
    s.layerM    = (jmat2 *)malloc(sizeof(jmat2) * M);
    s.layerDM   = (jmat2 *)malloc(sizeof(jmat2) * M);
    s.prefix    = (jmat2 *)malloc(sizeof(jmat2) * (N + 1));
    s.suffix    = (jmat2 *)malloc(sizeof(jmat2) * (N + 1));
    s.prefixExp = (int *)  malloc(sizeof(int) * (N + 1));
    s.suffixExp = (int *)  malloc(sizeof(int) * (N + 1));
    void *const blocks[] = { s.layerM, s.layerDM, s.prefix, s.suffix, s.prefixExp, s.suffixExp };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        jphase_jacobian_free(s);
        return 0;
    }
    *out = s;
    return 1;
}

/* One wavelength of the phase Jacobian; see tmm_phase_jacobian for the layout
 * of `out` and `deriv`. */
static void jphase_jacobian_core(double lambda, double omega, double theta_deg, int pol,
                                 jet n0, jet ns, const jet *layerN, const double *thick, int N,
                                 const jet *sinePtr, double *out, double *deriv,
                                 jphase_jacobian_scratch s) {
    jmat2 *layerM = s.layerM, *layerDM = s.layerDM, *prefix = s.prefix, *suffix = s.suffix;
    int *prefixExp = s.prefixExp, *suffixExp = s.suffixExp;
    jet wavelength = jwavelength(lambda, omega);
    jangle angle;
    jincidence(theta_deg, sinePtr, &angle);
    jet incidentEta = jadmittance(n0, jincident_cos(n0, angle.sine, angle.cosine), pol);
    jet substrateCosine = jsnell_cos(n0, angle.sine, angle.cosine, ns);
    jet substrateEta = jadmittance(ns, substrateCosine, pol);

    double transmissionLogScale = 0.0;
    for (int k = 0; k < N; k++) {
        /* Match jphase_core's skip rule for invalid negative/NaN thicknesses,
         * while retaining the useful derivative of a zero-thickness layer. */
        if (!(thick[k] >= 0.0)) {
            layerM[k] = jidentity();
            layerDM[k] = jzero();
            continue;
        }
        jet cosine = jsnell_cos(n0, angle.sine, angle.cosine, layerN[k]);
        double excess;
        jlayer_matrix_dd(layerN[k], thick[k], wavelength, cosine, pol,
                         &layerM[k], &layerDM[k], &excess);
        transmissionLogScale += excess;
    }

    prefix[0] = jidentity();
    prefixExp[0] = 0;
    for (int k = 0; k < N; k++) {
        prefix[k + 1] = jmatmul(prefix[k], layerM[k]);
        prefixExp[k + 1] = prefixExp[k] + jrescale_binary(&prefix[k + 1]);
    }
    suffix[N] = jidentity();
    suffixExp[N] = 0;
    for (int k = N - 1; k >= 0; k--) {
        suffix[k] = jmatmul(layerM[k], suffix[k + 1]);
        suffixExp[k] = suffixExp[k + 1] + jrescale_binary(&suffix[k]);
    }

    int totalExp = prefixExp[N];
    jcoef coefficients = jcoef_from_matrix(prefix[N], incidentEta, substrateEta,
                                           ldexp(1.0, -totalExp));
    jphase(coefficients.reflection, &out[0]);
    jphase(coefficients.transmission, &out[5]);
    if (transmissionLogScale > 0.0) out[9] *= exp(-2.0 * transmissionLogScale);

    for (int k = 0; k < N; k++) {
        jmat2 matrixDerivative = jmatmul(jmatmul(prefix[k], layerDM[k]), suffix[k + 1]);
        jmat_scale(&matrixDerivative, ldexp(1.0, prefixExp[k] + suffixExp[k + 1] - totalExp));
        jet dReflection, dTransmission;
        jcoef_thickness(matrixDerivative, coefficients, incidentEta, substrateEta,
                        &dReflection, &dTransmission);
        const jet sides[2] = { dReflection, dTransmission };
        const jet base[2] = { coefficients.reflection, coefficients.transmission };
        for (int s2 = 0; s2 < 2; s2++) {
            cx dd[4];
            /* d(ln coefficient)/dd: the imaginary part is the phase derivative
             * (negated for Macleod's sign), the real part is d(ln |c|)/dd, so
             * twice it is the relative derivative of |c|². */
            jderivs(jdiv(sides[s2], base[s2]), dd);
            deriv[((long)s2 * 5 + 0) * N + k] = -dd[0].im * 180.0 / PI;
            deriv[((long)s2 * 5 + 1) * N + k] = dd[1].im;
            deriv[((long)s2 * 5 + 2) * N + k] = dd[2].im;
            deriv[((long)s2 * 5 + 3) * N + k] = dd[3].im;
            deriv[((long)s2 * 5 + 4) * N + k] = 2.0 * dd[0].re;
        }
    }
}

TMM_EXPORT("tmm_phase_jacobian")
void tmm_phase_jacobian(double lambda, double omega, double theta_deg, int pol,
                        const double *n0jet, const double *nsjet,
                        const double *layerJets, const double *thick, int N,
                        const double *sinJet, double *out, double *deriv) {
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    jphase_jacobian_scratch scratch;
    if (!layerN || !jphase_jacobian_alloc(N, &scratch)) {
        free(layerN);
        fill_nan(out, 10); fill_nan(deriv, 10L * N);
        return;
    }
    for (int k = 0; k < N; k++) layerN[k] = jread(&layerJets[8 * k]);
    jet sine;
    if (sinJet) sine = jread(sinJet);
    jphase_jacobian_core(lambda, omega, theta_deg, pol, jread(n0jet), jread(nsjet),
                         layerN, thick, N, sinJet ? &sine : NULL, out, deriv, scratch);
    jphase_jacobian_free(scratch);
    free(layerN);
}

/* ── Exported: batched phase Jacobian over a wavelength grid ─────────────────
 * tmm_phase_jacobian at every wavelength of a grid in one call, for a caller
 * fitting a whole measured spectrum of phase-derived quantities: ellipsometric
 * Ψ and Δ, or group delay. Inputs are laid out as for tmm_phase_spectrum.
 *
 *   out   : nLam × 10        as tmm_phase_spectrum
 *   deriv : nLam × 10 × N    per λ, the 10 × N block of tmm_phase_jacobian */

TMM_EXPORT("tmm_phase_jacobian_spectrum")
void tmm_phase_jacobian_spectrum(const double *lambdas, const double *omegas, int nLam,
                                 const double *n0jets, const double *nsjets,
                                 const double *matJets, const double *thick, int N,
                                 double theta_deg, int pol,
                                 const double *sinJets, double *out, double *deriv) {
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    jphase_jacobian_scratch scratch;
    if (!layerN || !jphase_jacobian_alloc(N, &scratch)) {
        free(layerN);
        fill_nan(out, 10L * nLam); fill_nan(deriv, 10L * N * nLam);
        return;
    }
    for (int li = 0; li < nLam; li++) {
        for (int k = 0; k < N; k++) {
            long base = ((long)k * nLam + li) * 8;
            layerN[k] = jread(&matJets[base]);
        }
        jet sine;
        if (sinJets) sine = jread(&sinJets[8 * (long)li]);
        jphase_jacobian_core(lambdas[li], omegas[li], theta_deg, pol,
                             jread(&n0jets[8 * (long)li]), jread(&nsjets[8 * (long)li]),
                             layerN, thick, N, sinJets ? &sine : NULL,
                             &out[10 * (long)li], &deriv[10L * N * li], scratch);
    }
    jphase_jacobian_free(scratch);
    free(layerN);
}
