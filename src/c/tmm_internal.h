/*
 * tmm_internal.h : what the sources of the tmmcore C kernel share.
 *
 * The kernel computes reflectance, transmittance and absorptance for a stack
 * of absorbing, dispersive layers at arbitrary angle of incidence, in s and p
 * polarization, together with exact analytic derivatives: the thickness
 * Jacobian, the thickness Hessian, and the needle-insertion P-function.
 *
 * Also computes phase, group delay, group delay dispersion and third-order
 * dispersion, by carrying the same characteristic matrix in third-order Taylor
 * arithmetic with angular frequency as the differentiation variable, and
 * exports that Taylor arithmetic so a caller can build dispersive index jets.
 *
 * Its sources, compiled together:
 *   rta.c          R, T and A at one wavelength and over a grid
 *   growing.c      a growing top layer: the monitoring curve, the spectra
 *                  after each deposited layer, a growing stack over a grid
 *   derivatives.c  the thickness Jacobian, the thickness Hessian and the
 *                  needle scan
 *   phase.c        phase, GD, GDD and TOD, and their thickness derivatives
 *   jet.c          the Taylor-jet arithmetic, over a grid of jets
 *
 * This header holds what they share: complex numbers, 2×2 matrices, Snell's
 * law, the layer matrix, the boundaries, the layout of a grid call, and the
 * Taylor jets.
 *
 * Conventions:
 *   ñ = n + ik              k > 0 for absorbing media
 *   time factor exp(−iωt)   a wave exp(i(kz − ωt)) decays for k > 0
 *   off-diagonals of the characteristic matrix carry −i
 *   p admittance n/cosθ     so r_p = r_s at normal incidence (Macleod §9.7.1)
 *   wavelengths and thicknesses in nm, angles in degrees from normal
 *   pol 0 = s, 1 = p
 *
 * Every exported call takes its arguments in one order: wavelength or
 * wavelength grid, angle, polarization, incident medium, substrate, layers,
 * then any further inputs, then the outputs. A grid call takes the angle as
 * `theta_deg` with `nTheta` values: 1 for one angle across the grid, or one
 * per wavelength.
 *
 * All arithmetic is double precision. The JavaScript implementation in src/js
 * is the behavioural reference, and the comments here name the file of it
 * that each piece mirrors; the two agree to float64 round-off (see tests/).
 *
 * Builds with any C99 compiler; no dependencies beyond libm.
 *
 * References:
 *   • Macleod, Thin-Film Optical Filters 5th ed., §2.4, Eqs. 2.111, 2.123–2.125;
 *     Eq. 2.83; §10.2, Eqs. 10.11-10.13 (tilted admittances)
 *   • Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6)
 *   • Tikhonravov, Trubetskov & DeBell, Appl. Opt. 35, 5493 (1996)
 *   • Birge & Kärtner, Appl. Opt. 45, 1478 (2006)   [phase dispersion]
 */

#ifndef TMM_INTERNAL_H
#define TMM_INTERNAL_H

#include <math.h>
#include <stdlib.h>

/* WebAssembly builds export by name; a Windows DLL exports what is marked
 * dllexport; elsewhere every external function of a shared library is
 * visible already. */
#if defined(__wasm__) || defined(__EMSCRIPTEN__)
#define TMM_EXPORT(name) __attribute__((export_name(name)))
#elif defined(_WIN32)
#define TMM_EXPORT(name) __declspec(dllexport)
#else
#define TMM_EXPORT(name)
#endif

#define PI 3.14159265358979323846
#define MAX_IM_DELTA 50.0
#define MATRIX_RESCALE_THRESHOLD 1e100
#define BINARY_RESCALE_LIMIT 18446744073709551616.0   /* 2^64 */

/* ── Complex number arithmetic (cx = {re, im}) ───────────────────────────────
 * Mirrors cadd/csub/cmul/cdiv/cabs2/cconj/csqrt/ccos/csin in primitives.js. */

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
static inline void fill_nan(double *out, long count) {
    for (long i = 0; i < count; i++) out[i] = NAN;
}

/* 1 when every one of `count` allocations succeeded. */
static inline int all_allocated(void *const *blocks, int count) {
    for (int i = 0; i < count; i++) if (!blocks[i]) return 0;
    return 1;
}

/* The angle of incidence at wavelength index `li` of a grid: one angle for the
 * whole grid when nTheta is 1, else one per wavelength. */
static inline double theta_at(const double *theta_deg, int nTheta, int li) {
    return theta_deg[nTheta > 1 ? li : 0];
}

/* The inputs of a grid call, laid out as tmm_rta_spectrum documents. */
typedef struct {
    const double *lambdas; int nLam;
    const double *theta_deg; int nTheta;
    const double *n0arr, *nsarr;
    const double *matNK, *thick; int N;
} grid_layout;

/* The layers of a grid call at wavelength index `li`, as the N triples
 * [n_re, n_im, d] a one-wavelength call takes: matNK is N × nLam × 2, layout
 * [layer][λ][re, im]. */
static inline void layers_at(const grid_layout *g, int li, double *layers) {
    for (int k = 0; k < g->N; k++) {
        long base = ((long)k * g->nLam + li) * 2;
        layers[3 * k + 0] = g->matNK[base + 0];
        layers[3 * k + 1] = g->matNK[base + 1];
        layers[3 * k + 2] = g->thick[k];
    }
}

/* ── Snell's law: cosθ_j from incident (n0, θ0) into medium nj ──────────────
 * The transverse invariant is Re(n0) sinθ0. At oblique incidence Macleod keeps
 * the incident medium free of absorption (5th ed., §2.3.5 and §10.2); for an
 * absorbing one tmmcore takes the real part, so that medium carries a wave
 * whose amplitude falls along the normal only, and every medium's cosθ follows
 * from that one real invariant. Mirrors incidence / snellCosTheta /
 * incidentCosTheta in primitives.js, which gives the reason. */

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
 * layerMatrix in primitives.js. */

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
 * name in primitives.js. */

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
 * vanish. Mirrors substrateVector / transmittedFlux in primitives.js. */

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
 * layerAbsorptance in primitives.js. */
static inline double layerAbsorptance(double R, double T, cx eta0, cx r) {
    return 1.0 - R - T + 2.0 * (eta0.im / eta0.re) * r.im;
}

/* ── Third-order Taylor jets ─────────────────────────────────────────────────
 * A jet holds [f, f', f''/2!, f'''/3!], each entry complex. Ordinary power-
 * series algebra on these differentiates a function exactly, with no finite
 * differences. Port of src/js/jet.js; where the JS multiplies by a reciprocal
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

static inline void jsincos(jet a, jet *sine, jet *cosine) {
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

/* The real part of a jet, order by order. */
static inline jet jreal(jet a) {
    jet o;
    for (int i = 0; i < JET_N; i++) o.c[i] = cmk(a.c[i].re, 0.0);
    return o;
}

#endif
