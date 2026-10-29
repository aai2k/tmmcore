/*
 * phase.c : phase, group delay, GDD and TOD of r and t, and their thickness
 * derivatives, at one wavelength and over a grid, from the characteristic
 * matrix carried in Taylor jets of angular frequency. Conventions are in
 * tmm_internal.h.
 */

#include "tmm_internal.h"

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
 * Mirrors the phase call in phase.js. Index jets are 8 doubles each, [re,im]
 * per order. `omega` is the angular frequency at `lambda` and sets the time
 * unit of the result; `sinJet`, the incident sine as a jet, may be NULL. */

TMM_EXPORT("tmm_phase")
void tmm_phase(double lambda, double theta_deg, int pol,
               const double *n0jet, const double *nsjet,
               const double *layerJets, const double *thick, int N,
               double omega, const double *sinJet, double *out) {
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
 * way tmm_rta_spectrum does. Polarization is an argument rather than both-at-once,
 * because this kernel is an order of magnitude dearer per sample than the plain
 * spectrum and callers at normal incidence would pay twice for nothing.
 *
 * Memory layout (all f64, caller-owned):
 *   lambdas   : nLam
 *   theta_deg : nTheta, 1 or nLam
 *   n0jets    : nLam × 8               incident-medium index jet per λ
 *   nsjets    : nLam × 8               substrate index jet per λ
 *   matJets   : N × nLam × 8           per-layer index jet, layout [layer][λ]
 *   thick     : N
 *   omegas    : nLam                   angular frequency per λ; sets the unit
 *   sinJets   : nLam × 8, or NULL      incident sine as a jet; replaces θ
 *   out       : nLam × 10              [r: phaseRad,GD,GDD,TOD,|r|²][t: same] */

TMM_EXPORT("tmm_phase_spectrum")
void tmm_phase_spectrum(const double *lambdas, int nLam,
                        const double *theta_deg, int nTheta, int pol,
                        const double *n0jets, const double *nsjets,
                        const double *matJets, const double *thick, int N,
                        const double *omegas, const double *sinJets, double *out) {
    jet *layerN = (jet *)malloc(sizeof(jet) * (N > 0 ? N : 1));
    if (!layerN) { fill_nan(out, 10L * nLam); return; }
    for (int li = 0; li < nLam; li++) {
        for (int k = 0; k < N; k++) {
            long base = ((long)k * nLam + li) * 8;
            layerN[k] = jread(&matJets[base]);
        }
        jet sine;
        if (sinJets) sine = jread(&sinJets[8 * (long)li]);
        jphase_core(lambdas[li], omegas[li], theta_at(theta_deg, nTheta, li), pol,
                    jread(&n0jets[8 * (long)li]), jread(&nsjets[8 * (long)li]),
                    layerN, thick, N,
                    sinJets ? &sine : NULL, &out[10 * (long)li]);
    }
    free(layerN);
}

/* ── Exported: phase dispersion plus exact thickness derivatives ──────────────
 * Mirrors the phase Jacobian in phase.js, with the arguments of tmm_phase.
 * Frequency stays the Taylor variable, so each thickness derivative is itself
 * a third-order frequency jet.
 * Zero-thickness layers are retained so derivative indices line up with the
 * caller's design array. Negative and NaN thicknesses are skipped with a zero
 * derivative, matching the point evaluator's base result.
 *
 *   out   : 10        as tmm_phase
 *  deriv : 10 × N    [side][quantity][layer], side 0 = r, 1 = t,
 *                     quantity 0 = dPhaseDeg, 1 = dGD, 2 = dGDD, 3 = dTOD,
 *                     4 = d(ln |coefficient|²)/dd, the relative intensity
 *                     derivative, which the amplitude ratio in ellipsometry
 *                     needs and which costs nothing here: it is the real part
 *                     of the same logarithmic derivative the phase is the
 *                     imaginary part of
 *
 * The prefix and suffix products carry binary exponents, as in the derivative
 * kernels of derivatives.c, so opaque stacks whose products leave double range
 * still get their derivatives. */

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
void tmm_phase_jacobian(double lambda, double theta_deg, int pol,
                        const double *n0jet, const double *nsjet,
                        const double *layerJets, const double *thick, int N,
                        double omega, const double *sinJet, double *out, double *deriv) {
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
void tmm_phase_jacobian_spectrum(const double *lambdas, int nLam,
                                 const double *theta_deg, int nTheta, int pol,
                                 const double *n0jets, const double *nsjets,
                                 const double *matJets, const double *thick, int N,
                                 const double *omegas, const double *sinJets,
                                 double *out, double *deriv) {
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
        jphase_jacobian_core(lambdas[li], omegas[li], theta_at(theta_deg, nTheta, li), pol,
                             jread(&n0jets[8 * (long)li]), jread(&nsjets[8 * (long)li]),
                             layerN, thick, N, sinJets ? &sine : NULL,
                             &out[10 * (long)li], &deriv[10L * N * li], scratch);
    }
    jphase_jacobian_free(scratch);
    free(layerN);
}

/* ── Exported: angular frequency of a vacuum wavelength ──────────────────────
 * ω = 2πc/λ in rad/fs for λ in nm, the default frequency of the phase calls,
 * which gives GD, GDD and TOD in fs, fs² and fs³. */

#define C_NM_PER_FS 299.792458

TMM_EXPORT("tmm_omega_from_lambda_nm")
void tmm_omega_from_lambda_nm(const double *lambda_nm, int count, double *out) {
    for (int i = 0; i < count; i++) out[i] = 2.0 * PI * C_NM_PER_FS / lambda_nm[i];
}
