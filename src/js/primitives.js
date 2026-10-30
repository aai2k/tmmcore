/**
 * primitives.js : the pieces of the transfer-matrix method that every kernel
 * of the JavaScript reference implementation is built from: complex numbers,
 * 2×2 matrices, Snell's law, the layer matrix and the boundaries.
 *
 * System model:
 *   incident medium (n0, θ0) → layer1 → … → layerN → substrate
 *
 * Conventions:
 *   ñ = n + ik              k ≥ 0 for absorbing media
 *   time factor exp(−iωt)   a wave exp(i(kz − ωt)) decays for k > 0
 *   off-diagonals of the characteristic matrix carry −i
 *   p admittance n/cosθ     so r_p = r_s at normal incidence (Macleod §9.7.1)
 *   wavelengths and thicknesses in nm, angles in degrees from normal
 *
 * This is the complex conjugate of Macleod's convention (ñ = n − ik, exp(+iωt),
 * +i on the off-diagonals). R, T and A are identical under conjugation; phase-
 * sensitive quantities negate the raw TMM phase to recover Macleod's sign.
 *
 * Complex numbers are [re, im] pairs. All arithmetic is double precision. The C
 * kernel in src/c is a line-by-line port of the JavaScript in this folder and
 * agrees with it to float64 round-off (see tests/).
 *
 * References:
 *   • Macleod, Thin-Film Optical Filters 5th ed., §2.4, Eqs. 2.111, 2.123–2.125;
 *     Eq. 2.83 (transmittance out of an absorbing incident medium);
 *     Eq. 2.85 and Eqs. 10.5-10.7 (tilted admittances from the real invariant n0 sinθ0)
 */

// ── Complex number arithmetic ─────────────────────────────────────────────────
// All complex numbers are [re, im] arrays.

function cadd([ar, ai], [br, bi]) { return [ar + br, ai + bi]; }
function csub([ar, ai], [br, bi]) { return [ar - br, ai - bi]; }
function cmul([ar, ai], [br, bi]) { return [ar * br - ai * bi, ar * bi + ai * br]; }
function cdiv([ar, ai], [br, bi]) {
    const d = br * br + bi * bi;
    return [(ar * br + ai * bi) / d, (ai * br - ar * bi) / d];
}
function cabs2([ar, ai]) { return ar * ar + ai * ai; }
function cconj([ar, ai]) { return [ar, -ai]; }
function csqrt([ar, ai]) {
    const r = Math.sqrt(Math.sqrt(ar * ar + ai * ai));
    const theta = Math.atan2(ai, ar) / 2;
    return [r * Math.cos(theta), r * Math.sin(theta)];
}
function ccos([ar, ai]) {
    return [Math.cos(ar) * Math.cosh(ai), -Math.sin(ar) * Math.sinh(ai)];
}
function csin([ar, ai]) {
    return [Math.sin(ar) * Math.cosh(ai), Math.cos(ar) * Math.sinh(ai)];
}
function creal([ar]) { return ar; }
function cimag([, ai]) { return ai; }
function cscale([re, im], s) { return [re * s, im * s]; }

// ── 2×2 complex matrix multiply ───────────────────────────────────────────────

function matmul(A, B) {
    return [
        [
            cadd(cmul(A[0][0], B[0][0]), cmul(A[0][1], B[1][0])),
            cadd(cmul(A[0][0], B[0][1]), cmul(A[0][1], B[1][1]))
        ],
        [
            cadd(cmul(A[1][0], B[0][0]), cmul(A[1][1], B[1][0])),
            cadd(cmul(A[1][0], B[0][1]), cmul(A[1][1], B[1][1]))
        ]
    ];
}

// M·v, a 2×2 complex matrix times a 2-vector.
function cmatvec(M, v) {
    return [
        cadd(cmul(M[0][0], v[0]), cmul(M[0][1], v[1])),
        cadd(cmul(M[1][0], v[0]), cmul(M[1][1], v[1])),
    ];
}

const MATRIX_RESCALE_THRESHOLD = 1e100;

// A common real scale factor cancels from reflectance (Macleod 5th ed.,
// Eq. 2.123, p. 45). Callers retain it for transmittance (Eq. 2.125) while the
// bounded characteristic-matrix product prevents opaque-stack overflow.
function rescaleMatrix(M) {
    let scale = 0;
    for (const row of M) {
        for (const [re, im] of row) {
            scale = Math.max(scale, Math.abs(re), Math.abs(im));
        }
    }
    if (scale <= MATRIX_RESCALE_THRESHOLD) return 0;
    const inverse = 1 / scale;
    for (const row of M) {
        for (const value of row) {
            value[0] *= inverse;
            value[1] *= inverse;
        }
    }
    return Math.log(scale);
}

// ── Snell's law ───────────────────────────────────────────────────────────────
//
// Every wave in the stack shares one transverse wavevector, n0 sinθ0 with n0
// the REAL part of the incident index. Macleod writes Snell's law and the
// tilted admittances for a transparent incident medium, so n0 sinθ0 is real
// (5th ed., Eq. 2.85 and Eqs. 10.5-10.7). Essential Macleod takes the angles
// of an absorbing incident medium from the real part of its index too. An
// absorbing incident medium then carries an inhomogeneous wave whose amplitude
// falls along the normal only, and each medium's cosθ follows from that one
// real invariant. Putting the complex index into the invariant would make the
// amplitude vary along the interface, so energy would flow sideways inside
// lossless layers and R + T would be wrong by a term linear in the incident k
// at oblique incidence.

// sinθ0 and cosθ0 of the angle of incidence, each as [value, 0]. The cosine is
// taken directly: near grazing incidence sqrt(1 − sin²θ0) keeps three correct
// digits of it at 89.99999° and none past 89.9999999°.
function incidence(theta_deg) {
    const rad = theta_deg * Math.PI / 180;
    return { sinTheta0: [Math.sin(rad), 0], cosTheta0: [Math.cos(rad), 0] };
}

// cosθ in medium nj. With a = n0 sinθ0 the real invariant,
//     cos²θ = 1 − (a/nj)²,
// and since n0² = a² + (n0 cosθ0)² for a real angle, the same number is
//     ((nj − n0)(nj + n0) + (n0 cosθ0)²) / nj².
// The first form cancels where (a/nj)² is near 1, and at grazing incidence
// into a medium of the incident index nothing of cos²θ survives it. The
// second cancels where nj² − n0² and (n0 cosθ0)² nearly offset. Their
// rounding errors scale as a² and as |nj² − n0²| + (n0 cosθ0)², so the second
// is taken where |nj² − n0²| + 2 (n0 cosθ0)² < n0². That never holds at
// normal incidence, where cosθ = 1 stays exact. `cosTheta0`, the cosine of the
// real angle of incidence as [cosθ0, 0], is optional; without it the first
// form is used.
function snellCosTheta(n0, sinTheta0, nj, cosTheta0) {
    if (cosTheta0) {
        const nr = n0[0];
        const nc = nr * cosTheta0[0];
        const q = nc * nc;
        // n0² − 2 (n0 cosθ0)², positive only past 45°. Since
        // |nj² − n0²| ≥ ||nj|² − n0²|, a medium far from the incident index
        // fails the test before (nj − n0)(nj + n0) is formed.
        const room = nr * nr - 2 * q;
        const near = room > 0 && Math.abs(cabs2(nj) - nr * nr) < room;
        const m = near ? cmul(csub(nj, [nr, 0]), cadd(nj, [nr, 0])) : null;
        if (m && cabs2(m) < room * room) return csqrt(cdiv(cadd(m, [q, 0]), cmul(nj, nj)));
    }
    // sinThetaJ = Re(n0) * sinTheta0 / nj   (complex)
    const sinThetaJ = cdiv(cmul([n0[0], 0], sinTheta0), nj);
    // cosTheta = sqrt(1 - sin²θ)
    return csqrt(csub([1, 0], cmul(sinThetaJ, sinThetaJ)));
}

// cosθ0 of the incident medium. A transparent medium keeps its own cosine,
// cosTheta0 when given, else sqrt(1 − sin²θ0). An absorbing one takes its
// cosine from the same real invariant as the layers, so its tilted admittance
// is sqrt(N0² − n0² sin²θ0) and the incident wave belongs to the same boundary
// problem as the rest of the stack. R + T then departs from 1 for lossless
// layers only by the interference of the incident and reflected waves in the
// absorbing medium, the term Macleod's Eq. 2.83 carries, which is second order
// in k0 for a bare interface. Essential Macleod reports R and T at the real
// part of this admittance, which leaves that term out; tests/macleod.mjs
// derives its values from this r through the coating's surface admittance.
function incidentCosTheta(n0, sinTheta0, cosTheta0) {
    if (n0[1] !== 0) return snellCosTheta(n0, sinTheta0, n0, cosTheta0);
    return cosTheta0
        ? [cosTheta0[0], cosTheta0[1]]
        : csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));
}

// ── Layer characteristic matrix ───────────────────────────────────────────────

// Overflow guard on the phase thickness. ccos/csin use cosh/sinh(Im δ), and
// cosh(710) = Inf would turn the whole product into NaN. So |Im δ| is held at
// MAX_IM_DELTA, where cosh(50) ≈ 2.6e21 stays finite under products. Past the
// bound cosδ and sinδ are each one growing exponential to within e^{−100}
// relative, so the held matrix is the true one divided by e^{|Im δ| − 50}: the
// same shape, only smaller. That factor cancels from r and from every
// admittance, and is carried alongside as a log scale (layerLogScale) for t,
// T and fields, which it does not cancel from. The bound tests Im δ whatever
// k is: a lossless layer past the critical angle carries an evanescent wave
// with imaginary δ, and a thick enough one reaches the bound too.
const MAX_IM_DELTA = 50;

// Phase thickness δ = (2π/λ) n d cosθ (Macleod 5th ed., Eq. 9.2), complex,
// with Im δ held to ±MAX_IM_DELTA. `excess` is |Im δ| − MAX_IM_DELTA where
// the bound held it and 0 elsewhere: the natural log of the factor the held
// matrix falls short of the true one by.
function layerPhase(nj, dj_nm, lambda_nm, cosTheta_j) {
    const k0 = (2 * Math.PI) / lambda_nm;
    const delta = cmul(cmul(nj, [k0 * dj_nm, 0]), cosTheta_j);
    const excess = Math.abs(delta[1]) > MAX_IM_DELTA ? Math.abs(delta[1]) - MAX_IM_DELTA : 0;
    if (delta[1] > MAX_IM_DELTA) delta[1] = MAX_IM_DELTA;
    else if (delta[1] < -MAX_IM_DELTA) delta[1] = -MAX_IM_DELTA;
    return { delta, excess };
}

// The log scale layerMatrix leaves out for the same arguments: its matrix
// times e^{layerLogScale} is the layer's true matrix. A caller building a
// product adds it to the log scale it carries for t, exactly as it adds the
// return of rescaleMatrix; T then takes e^{−2·logScale}. Zero for every layer
// short of the bound, so an ordinary stack is unaffected bit for bit.
function layerLogScale(nj, dj_nm, lambda_nm, cosTheta_j) {
    return layerPhase(nj, dj_nm, lambda_nm, cosTheta_j).excess;
}

// Tilted admittance in units of the free-space admittance: n cosθ for s,
// n / cosθ for p (Macleod 5th ed., §9.2).
function layerAdmittance(nj, cosTheta_j, pol) {
    return pol === 's' ? cmul(nj, cosTheta_j) : cdiv(nj, cosTheta_j);
}

// M = [[cosδ, −i sinδ/η], [−i η sinδ, cosδ]]. Two adjacent layers of one
// material act as one layer of their summed thickness, so matrices of one
// admittance compose by adding phases: M(a)·M(b) = M(a + b).
function phaseMatrix(delta, eta) {
    const cosD = ccos(delta);
    const sinD = csin(delta);
    const iSinD_div_eta = cmul([0, -1], cdiv(sinD, eta));
    const iEta_sinD     = cmul([0, -1], cmul(eta, sinD));
    return [
        [cosD,        iSinD_div_eta],
        [iEta_sinD,   cosD         ]
    ];
}

function layerMatrix(nj, dj_nm, lambda_nm, cosTheta_j, pol) {
    if (atCriticalAngle(cosTheta_j)) {
        return criticalMatrix(criticalGenerator(nj, (2 * Math.PI) / lambda_nm, pol), dj_nm);
    }
    return phaseMatrix(layerPhase(nj, dj_nm, lambda_nm, cosTheta_j).delta,
                       layerAdmittance(nj, cosTheta_j, pol));
}

// ── The critical angle ───────────────────────────────────────────────────────
//
// Where cosθ in a medium is exactly zero the wave in it runs along the
// interfaces. δ and the s admittance n cosθ vanish together and the p
// admittance n/cosθ diverges, so phaseMatrix is 0/0 there while its limit is
// finite and linear in the thickness:
//     M = I + d·G,   G = [[0, −iP], [−iS, 0]],   P = Q/η,  S = Qη,
// with Q = dδ/dd = (2π/λ) n cosθ. Q/η_s = 2π/λ and Q·η_p = (2π/λ) n² hold at
// every angle, and the other two products carry cos²θ, so the pair (P, S) is
// (2π/λ, 0) for s and (0, (2π/λ) n²) for p. G is also the thickness
// derivative of any layer matrix at δ = 0, the needle's matrix derivative.

const atCriticalAngle = cosTheta => cosTheta[0] === 0 && cosTheta[1] === 0;

function criticalGenerator(nj, k0, pol) {
    return pol === 's'
        ? { P: [k0, 0], S: [0, 0] }
        : { P: [0, 0], S: cscale(cmul(nj, nj), k0) };
}

// [[0, −iP], [−iS, 0]]
function generatorMatrix({ P, S }) {
    return [[[0, 0], cmul([0, -1], P)], [cmul([0, -1], S), [0, 0]]];
}

// I + d·G
function criticalMatrix({ P, S }, dj_nm) {
    return [[[1, 0], cmul([0, -dj_nm], P)], [cmul([0, -dj_nm], S), [1, 0]]];
}

// ── Boundaries ───────────────────────────────────────────────────────────────

// The substrate closes the product as [B, C] = M·[1, ηs] (Macleod 5th ed.,
// Eq. 2.113). At its critical angle ηs diverges in p. r depends only on the
// direction of [B, C], so the vector becomes [0, 1]; t = 2η0 bs/(η0B + C) and
// the transmitted flux Re(cs·conj(bs)) both vanish with bs, and no energy
// crosses into the substrate.
function substrateVector(ns, cosThetaS, pol) {
    return pol === 'p' && atCriticalAngle(cosThetaS)
        ? [[0, 0], [1, 0]]
        : [[1, 0], layerAdmittance(ns, cosThetaS, pol)];
}

// Re(cs·conj(bs)): Re(ηs) for [1, ηs] and 0 for [0, 1].
const transmittedFlux = ([bs, cs]) => creal(cmul(cs, cconj(bs)));

// Absorptance of the layers: the net irradiance entering the front surface,
// ½Re(E H*) there (Macleod 5th ed., Eq. 2.120), less the transmitted, per unit
// incident irradiance ½Re(η0)|E⁺|². With E = E⁺(1 + r) and H = η0 E⁺(1 − r),
//     A = 1 − R − T + 2 (Im η0 / Re η0) Im r.
// The last term is zero for a transparent incident medium. In an absorbing
// one the incident and reflected waves interfere in the irradiance (the mixed
// Poynting vector, Macleod 5th ed., §2.13, Eq. 2.163), R + T + A differs from 1
// by exactly that term, and 1 − R − T is not the absorptance.
function layerAbsorptance(R, T, eta0, r) {
    return 1 - R - T + 2 * (eta0[1] / eta0[0]) * r[1];
}

// Shared by the kernels of this folder. The package index exports the ones a
// caller building transfer-matrix variants on the same conventions needs;
// they are lower level than the calls, and correspondingly less stable across
// versions.
export {
    cadd, csub, cmul, cdiv, cabs2, cconj, csqrt, ccos, csin, creal, cimag, cscale,
    matmul, cmatvec, rescaleMatrix, incidence, snellCosTheta, incidentCosTheta,
    layerPhase, layerLogScale, layerAdmittance, phaseMatrix, layerMatrix,
    atCriticalAngle, criticalGenerator, generatorMatrix, criticalMatrix,
    substrateVector, transmittedFlux, layerAbsorptance,
};
