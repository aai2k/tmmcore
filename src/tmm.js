/**
 * Transfer-matrix method for multilayer thin films: JavaScript reference
 * implementation.
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
 * kernel shipped alongside is a line-by-line port of this file and agrees with
 * it to float64 round-off (see tests/).
 *
 * References:
 *   • Macleod, Thin-Film Optical Filters 5th ed., §2.4, Eqs. 2.111, 2.123–2.125;
 *     Eq. 2.83 (transmittance out of an absorbing incident medium);
 *     §10.2, Eqs. 10.11–10.13 (tilted admittances from the real invariant n0 sinθ0)
 *   • Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6)
 *   • Tikhonravov, Trubetskov & DeBell, Appl. Opt. 35, 5493 (1996)
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
// the REAL part of the incident index (Macleod 5th ed., §10.2, Eqs. 10.11 to
// 10.13). An absorbing incident medium then carries an inhomogeneous wave whose
// amplitude falls along the normal only, and each medium's cosθ follows from
// that one real invariant. Putting the complex index into the invariant would
// make the amplitude vary along the interface, so energy would flow sideways
// inside lossless layers and R + T would be wrong by a term linear in the
// incident k at oblique incidence.

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
// in k0 for a bare interface.
function incidentCosTheta(n0, sinTheta0, cosTheta0) {
    if (n0[1] !== 0) return snellCosTheta(n0, sinTheta0, n0, cosTheta0);
    return cosTheta0
        ? [cosTheta0[0], cosTheta0[1]]
        : csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));
}

// ── Layer characteristic matrix ───────────────────────────────────────────────

// Overflow guard on the phase thickness. ccos/csin use cosh/sinh(Im δ), and
// cosh(710) = Inf would turn the whole product into NaN. By |Im δ| of a few
// tens the layer is optically opaque, its single-pass transmittance e^{−2 Im δ}
// below 4e−44 at the bound, and the surface reflectance has converged, so
// holding Im δ there is exact to machine precision while cosh(50) ≈ 2.6e21
// stays finite under products. The bound tests Im δ whatever k is: a lossless
// layer past the critical angle carries an evanescent wave with imaginary δ,
// and a thick enough one reaches the bound too.
const MAX_IM_DELTA = 50;

// Phase thickness δ = (2π/λ) n d cosθ (Macleod 5th ed., Eq. 9.2), complex,
// with Im δ held to ±MAX_IM_DELTA. `clamped` reports whether the bound held
// it, in which case δ depends on d through its real part only.
function layerPhase(nj, dj_nm, lambda_nm, cosTheta_j) {
    const k0 = (2 * Math.PI) / lambda_nm;
    const delta = cmul(cmul(nj, [k0 * dj_nm, 0]), cosTheta_j);
    const clamped = Math.abs(delta[1]) > MAX_IM_DELTA;
    if (delta[1] > MAX_IM_DELTA) delta[1] = MAX_IM_DELTA;
    else if (delta[1] < -MAX_IM_DELTA) delta[1] = -MAX_IM_DELTA;
    return { delta, clamped };
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

// ── Core TMM for one wavelength ───────────────────────────────────────────────

/**
 * @param {number}   lambda_nm  wavelength in nm
 * @param {number}   theta_deg  angle of incidence in degrees
 * @param {string}   pol        's' or 'p'
 * @param {[number,number]} n0  complex n of incident medium
 * @param {[number,number]} ns  complex n of substrate (exit medium)
 * @param {{ n:[number,number], d:number }[]} layers  each layer { n: [re,im], d: thickness_nm }
 * @returns {{ R:number, T:number, A:number }}
 */
export function tmm(lambda_nm, theta_deg, pol, n0, ns, layers) {
    const { sinTheta0, cosTheta0 } = incidence(theta_deg);

    // Admittance of incident medium
    const eta0 = layerAdmittance(n0, incidentCosTheta(n0, sinTheta0, cosTheta0), pol);

    // [1, ηs], or [0, 1] at the substrate's critical angle in p
    const substrate = substrateVector(ns, snellCosTheta(n0, sinTheta0, ns, cosTheta0), pol);

    // Build total transfer matrix M = M1 × M2 × ... × MN
    let M = [[  [1, 0], [0, 0]  ], [  [0, 0], [1, 0]  ]]; // identity
    let logScale = 0;

    for (const { n, d } of layers) {
        // Not positive, NaN included: the layer is absent.
        if (!(d > 0)) continue;
        const cosThetaJ = snellCosTheta(n0, sinTheta0, n, cosTheta0);
        const Mj = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
        M = matmul(M, Mj);
        logScale += rescaleMatrix(M);
    }

    // [B, C]^T = M × substrate
    const [B, C] = cmatvec(M, substrate);

    // r = (η0 B - C) / (η0 B + C)
    const eta0B = cmul(eta0, B);
    const r = cdiv(csub(eta0B, C), cadd(eta0B, C));

    // t = 2 η0 bs / (η0 B + C)
    const t = cmul(cdiv(cmul([2, 0], eta0), cadd(eta0B, C)), substrate[0]);

    const R = cabs2(r);

    // T = Re(etaS) / Re(eta0) * |t|²
    const T = Math.max(0, transmittedFlux(substrate) / creal(eta0) * cabs2(t) * Math.exp(-2 * logScale));

    const A = Math.max(0, layerAbsorptance(R, T, eta0, r));

    return { R, T, A };
}

// ── Shared pieces of the derivative kernels ──────────────────────────────────
//
// Each derivative kernel splits the stack at every layer, the pre/post form of
// Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6):
//   Pre[j]  = M_0·…·M_{j−1}              a 2×2 matrix
//   Post[j] = M_j·…·M_{N−1}·[1, ηs]      a 2-vector, so [B, C] = Post[0]
// ([0, 1] in place of [1, ηs] at the substrate's critical angle in p).
// Through opaque or evanescent layers both grow without bound; past 1e77 a
// complex division by den², den = η0 B + C, forms |den|⁴ and overflows, and
// the whole wavelength's derivatives would come out zero. So each is stored as a
// mantissa times 2^e, with the integer e kept alongside. A power of two scales
// a double exactly, so the mantissas hold the same significant bits the
// unscaled products would. Every derivative is one prefix times one suffix
// divided by den, so it takes the exponent e_pre + e_post − e_0 and the scales
// leave no trace in the result.

const BINARY_RESCALE_LIMIT = 2 ** 64;

function cscale([re, im], s) { return [re * s, im * s]; }

// Divides the complex values in place by a power of two when the largest part
// exceeds BINARY_RESCALE_LIMIT, and returns that power's exponent.
function rescaleBinary(values) {
    let largest = 0;
    for (const [re, im] of values) largest = Math.max(largest, Math.abs(re), Math.abs(im));
    if (!(largest > BINARY_RESCALE_LIMIT) || largest === Infinity) return 0;
    const exponent = Math.floor(Math.log2(largest));
    const factor = 2 ** -exponent;
    for (const value of values) {
        value[0] *= factor;
        value[1] *= factor;
    }
    return exponent;
}

const matrixEntries = M => [M[0][0], M[0][1], M[1][0], M[1][1]];

// The pair (P, S) = (Q/η, Qη) of a medium at cosθ, with Q = dδ/dd as the
// matrix in the product moves: dM/dd = [[−Q sinδ, −iP cosδ], [−iS cosδ,
// −Q sinδ]], and at δ = 0 that is the generator G. Neither P nor S divides by
// an admittance that vanishes or diverges, since at the critical angle the
// pair takes its limit.
function generator(nj, cosTheta, Q, k0, pol) {
    if (atCriticalAngle(cosTheta)) return criticalGenerator(nj, k0, pol);
    const eta = layerAdmittance(nj, cosTheta, pol);
    return { P: cdiv(Q, eta), S: cmul(Q, eta) };
}

// One layer as the derivative kernels use it. A thickness that is not a number
// ≥ 0 leaves the layer out, as tmm() does: identity matrix and zero
// derivative. Zero thickness keeps its derivative, which is what a candidate
// layer grows from.
function derivativeLayer(n, d, { n0, sinTheta0, cosTheta0, k0, lambda_nm, pol }) {
    const present = d >= 0;
    const thickness = present ? d : 0;
    const cosTheta = snellCosTheta(n0, sinTheta0, n, cosTheta0);
    const { delta, clamped } = layerPhase(n, thickness, lambda_nm, cosTheta);
    const eta = layerAdmittance(n, cosTheta, pol);
    // Q = dδ/dd = (2π/λ) n cosθ. Where the bound holds Im δ, the matrix in
    // the product moves with d through Re δ alone, so Q loses its imaginary
    // part and every derivative is taken of that same matrix.
    const rawQ = cmul(cmul(n, [k0, 0]), cosTheta);
    const Q = !present ? [0, 0] : clamped ? [rawQ[0], 0] : rawQ;
    const { P, S } = present
        ? generator(n, cosTheta, Q, k0, pol)
        : { P: [0, 0], S: [0, 0] };
    const layer = { n, thickness, cosTheta, delta, eta, critical: atCriticalAngle(cosTheta), Q, P, S };
    layer.M = partMatrix(layer, delta, thickness);
    return layer;
}

// The matrix of `thickness` of a layer whose phase over it is `delta`.
function partMatrix(layer, delta, thickness) {
    return layer.critical ? criticalMatrix(layer, thickness) : phaseMatrix(delta, layer.eta);
}

// Everything the derivative kernels share at one (λ, θ, pol), from their
// common arguments. Layers are used as given, with no d > 0 filter, so
// derivative indices match the caller's design.
function derivativeStack({ lambda_nm, theta_deg, pol, n0, ns, layers }) {
    const { sinTheta0, cosTheta0 } = incidence(theta_deg);
    const eta0 = layerAdmittance(n0, incidentCosTheta(n0, sinTheta0, cosTheta0), pol);
    const substrate = substrateVector(ns, snellCosTheta(n0, sinTheta0, ns, cosTheta0), pol);
    const k0 = (2 * Math.PI) / lambda_nm;
    const N = layers.length;
    const geometry = { n0, sinTheta0, cosTheta0, k0, lambda_nm, pol };

    const data = layers.map(({ n, d }) => derivativeLayer(n, d, geometry));

    const Pre = new Array(N + 1);
    const preExp = new Array(N + 1);
    Pre[0] = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
    preExp[0] = 0;
    for (let j = 0; j < N; j++) {
        Pre[j + 1] = matmul(Pre[j], data[j].M);
        preExp[j + 1] = preExp[j] + rescaleBinary(matrixEntries(Pre[j + 1]));
    }
    const Post = new Array(N + 1);
    const postExp = new Array(N + 1);
    Post[N] = substrate;
    postExp[N] = 0;
    for (let j = N - 1; j >= 0; j--) {
        Post[j] = cmatvec(data[j].M, Post[j + 1]);
        postExp[j] = postExp[j + 1] + rescaleBinary(Post[j]);
    }

    return { geometry, eta0, substrate, N, data, Pre, preExp, Post, postExp };
}

// R, T, A of the whole stack from [B, C] = Post[0] (Macleod 5th ed., Eqs.
// 2.123–2.125), with the pieces the chain rule below needs: r, the true t,
// b = B/den and c = C/den, and g = 2 Im η0 / Re η0, the weight of Im r in the
// absorptance. Only t carries the exponent of [B, C].
function stackResponse({ eta0, substrate, Post, postExp }) {
    const [B, C] = Post[0];
    const den = cadd(cmul(eta0, B), C);
    const inverseDen = cdiv([1, 0], den);
    const r = cdiv(csub(cmul(eta0, B), C), den);
    const t = cscale(cmul(cdiv(cmul([2, 0], eta0), den), substrate[0]), 2 ** -postExp[0]);
    const Tfac = transmittedFlux(substrate) / creal(eta0);
    const R = cabs2(r);
    const T = Math.max(0, Tfac * cabs2(t));
    const A = Math.max(0, layerAbsorptance(R, T, eta0, r));
    return {
        R, T, A, r, t, Tfac, eta0, inverseDen, g: 2 * (eta0[1] / eta0[0]),
        b: cmul(B, inverseDen), c: cmul(C, inverseDen), exponent: postExp[0],
    };
}

// A derivative of [B, C], given as a mantissa vector and its exponent, divided
// by den: (u, w) = (dB/den, dC/den), free of every scale.
function perDen(response, [dB, dC], exponent) {
    const scale = 2 ** (exponent - response.exponent);
    const { inverseDen } = response;
    return [cscale(cmul(dB, inverseDen), scale), cscale(cmul(dC, inverseDen), scale)];
}

// {dR, dT, dA} from (u, w) = d[B, C]/den:
//   dr = 2η0 (c·u − b·w),   dt = −t (η0·u + w)
// the chain rule dr = f (C dB − B dC), dt = −f (η0 dB + dC) with f = 2η0/den²,
// written so that den is divided out once rather than squared. dA follows
// layerAbsorptance: −(dR + dT) + g Im dr.
function responseDerivative(response, [u, w]) {
    const { eta0, r, t, b, c, Tfac, g } = response;
    const dr = cmul(cmul([2, 0], eta0), csub(cmul(c, u), cmul(b, w)));
    const dR = 2 * creal(cmul(cconj(r), dr));
    const dt = cmul([-1, 0], cmul(t, cadd(cmul(eta0, u), w)));
    const dT = Tfac * 2 * creal(cmul(cconj(t), dt));
    return { dR, dT, dA: -(dR + dT) + g * dr[1] };
}

// ── Analytic needle P-function kernel ─────────────────────────────────────────
//
// Returns the ANALYTIC merit-function gradient dF/dd of inserting an
// infinitesimally thin needle, for every insertion position × candidate
// material, at one (λ, θ, pol).  This is the d→0 limit of Sullivan's
// numerical pre/post method, i.e. Tikhonravov's analytic P-function.
//
// Derivation (citations):
//   • Characteristic matrix & [B,C], r, t:  Macleod, Thin-Film Optical
//     Filters 5th ed., §2.4 Eqs. 2.111, 2.123–2.125 (JS sign convention:
//     off-diagonals carry −i, see layerMatrix above).
//   • Pre/post decomposition  M = M_pre · M_k · M_post  and needle
//     insertion:  Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996),
//     Eqs. (3)–(6).
//   • Needle series  dF = P₁ d + P₂ d² + …, insert where P₁<0:
//     Tikhonravov, Trubetskov & DeBell, Appl. Opt. 35, 5493 (1996),
//     Eqs. (1)–(2).
//
// Needle matrix of index nₐ, thickness d:
//   M_n(δ) = [[cosδ, −i sinδ/ηₐ], [−i ηₐ sinδ, cosδ]],  δ = (2π/λ) nₐ d cosθₐ
// As d→0:  M_n = I + A·d + O(d²),  with
//   A = [[0, −i Q/ηₐ], [−i ηₐ Q, 0]],   Q = (2π/λ) nₐ cosθₐ.
// Insertion at gap `pos`:  [B,C] = Pre·Post, and to first order
//   d[B,C]/dd = Pre · A · Post.
// Then with den = η₀B + C:
//   dr/dd = (2η₀/den²)·(C·dB − B·dC),   dR/dd = 2 Re[ r̄ · dr/dd ]
//   dt/dd = −(2η₀/den²)·(η₀·dB + dC),   dT/dd = (Re ηs/Re η₀)·2 Re[ t̄ · dt/dd ]
//   dA/dd = −(dR/dd + dT/dd) + (2 Im η₀/Re η₀)·Im(dr/dd)
// The host layer cancels automatically through Pre/Post (a needle of the
// host index at an interior point gives ~0), so no nₐ²−n_host² term is
// needed, exactly as in Sullivan's scheme.
//
// Returns { R, T, A, gaps, intra } where
//   gaps[pos]            = [{dR,dT,dA} per candidate]   pos = 0..N
//   intra[k][fi]         = { frac, perCand:[{dR,dT,dA}] }   (host-split)
function cmatvec(M, v) {
    return [
        cadd(cmul(M[0][0], v[0]), cmul(M[0][1], v[1])),
        cadd(cmul(M[1][0], v[0]), cmul(M[1][1], v[1])),
    ];
}

export function tmmNeedleScan(lambda_nm, theta_deg, pol, n0, ns, layers,
                              candidateNs, intraFracs = []) {
    // Layers are used as given, so gap and intra indices match the caller's
    // design; a zero-thickness layer is the identity matrix and harmless.
    const stack = derivativeStack({ lambda_nm, theta_deg, pol, n0, ns, layers });
    const { N, data, Pre, preExp, Post, postExp } = stack;
    const { sinTheta0, cosTheta0, k0 } = stack.geometry;
    const response = stackResponse(stack);

    // Needle derivative-matrix A = G for a candidate index:
    // [[0, −i Q/ηₐ], [−i ηₐ Q, 0]] with Q = (2π/λ) nₐ cosθₐ.
    function needleA(nA) {
        const cthA = snellCosTheta(n0, sinTheta0, nA, cosTheta0);
        const Q = cmul(cmul(nA, [k0, 0]), cthA);
        return generatorMatrix(generator(nA, cthA, Q, k0, pol));
    }
    const Acache = candidateNs.map(needleA);

    // {dR,dT,dA} per candidate for a needle between a prefix and a suffix.
    const insertAt = (pre, post, exponent) => Acache.map(Amat =>
        responseDerivative(response, perDen(response, cmatvec(pre, cmatvec(Amat, post)), exponent)));

    const gaps = new Array(N + 1);
    for (let pos = 0; pos <= N; pos++) {
        gaps[pos] = insertAt(Pre[pos], Post[pos], preExp[pos] + postExp[pos]);
    }

    const intra = [];
    if (intraFracs.length) {
        for (let k = 0; k < N; k++) {
            const host = data[k];
            const { n, thickness, cosTheta, delta } = host;
            intra.push(intraFracs.map(frac => {
                // The host splits into halves whose product is M_k itself, bound
                // included, so the needle sits in the stack R and T came from.
                // The front half keeps its own phase, bounded only if it alone
                // passes the bound, so the field at the needle's depth is the
                // true one; the back half takes the remainder of δ.
                const front = layerPhase(n, frac * thickness, lambda_nm, cosTheta).delta;
                const preIn = matmul(Pre[k], partMatrix(host, front, frac * thickness));
                const postIn = cmatvec(partMatrix(host, csub(delta, front), (1 - frac) * thickness),
                    Post[k + 1]);
                return { frac, perCand: insertAt(preIn, postIn, preExp[k] + postExp[k + 1]) };
            }));
        }
    }

    return { R: response.R, T: response.T, A: response.A, gaps, intra, N };
}

// ── Analytic thickness-Jacobian kernel ────────────────────────────────────────
//
// Returns the EXACT analytic derivatives dR/dd_k, dT/dd_k, dA/dd_k of every
// existing layer's thickness, at one (λ, θ, pol).  Replaces the central-
// difference Jacobian in the DLS refiner (2·N fewer TMM evals per step).
//
// Derivation (citations):
//   • Characteristic matrix Eq. 2.111 and product form Eq. 2.113 with
//       δ_r = 2π N_r d_r cosθ_r / λ
//     Macleod, Thin-Film Optical Filters 5th ed., §2.4 (verified verbatim).
//     This module's documented sign convention puts −i on the off-diagonals
//     (see file header & layerMatrix); the derivative below is taken of THAT
//     matrix, not Macleod's +i form, so it stays byte-consistent with tmm().
//   • Pre/post decomposition  [B,C] = Pre·M_k·Post  and the parametric
//     derivative  ∂[B,C]/∂p = Pre·(∂M_k/∂p)·Post:
//     Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6).
//
// Only δ depends on d_k (η = n cosθ for s, n/cosθ for p, and cosθ depend on
// n, θ, λ only).  With Q ≡ dδ/dd = (2π/λ) n cosθ:
//
//   dM_k/dd_k = Q · [[ −sinδ,      −i cosδ / η ],
//                    [ −i η cosδ,  −sinδ       ]]
//
// As δ→0 this collapses to [[0,−iQ/η],[−iQη,0]], exactly the needle
// A-matrix in tmmNeedleScan (needleA), i.e. the needle kernel is the δ=0
// special case of this; a strong internal-consistency check.
//
// δ and Q are those of the matrix in the product: past the imaginary-phase
// bound, the held δ and the real part of Q (see derivativeLayer). The
// off-diagonals are written with P = Q/η and S = Qη, which stay finite at the
// critical angle.
function layerMatrixDerivative({ delta, Q, P, S }) {
    const cD = ccos(delta), sD = csin(delta);
    return [
        [ cmul([-1, 0], cmul(Q, sD)),  cmul([0, -1], cmul(P, cD)) ],
        [ cmul([0, -1], cmul(S, cD)),  cmul([-1, 0], cmul(Q, sD)) ],
    ];
}

export function tmmThicknessJacobian(lambda_nm, theta_deg, pol, n0, ns, layers) {
    // Layers are used as given (index parity, see tmmNeedleScan).
    const stack = derivativeStack({ lambda_nm, theta_deg, pol, n0, ns, layers });
    const { N, data, Pre, preExp, Post, postExp } = stack;
    const response = stackResponse(stack);

    const dRdd = new Array(N), dTdd = new Array(N), dAdd = new Array(N);
    for (let k = 0; k < N; k++) {
        const dV = cmatvec(Pre[k], cmatvec(layerMatrixDerivative(data[k]), Post[k + 1]));
        const m = responseDerivative(response, perDen(response, dV, preExp[k] + postExp[k + 1]));
        dRdd[k] = m.dR; dTdd[k] = m.dT; dAdd[k] = m.dA;
    }

    return { R: response.R, T: response.T, A: response.A, dRdd, dTdd, dAdd, N };
}

// ── Analytic thickness-Hessian kernel ─────────────────────────────────────────
//
// Returns the EXACT analytic SECOND derivatives ∂²R/∂dᵢ∂dⱼ, ∂²T/∂dᵢ∂dⱼ,
// ∂²A/∂dᵢ∂dⱼ (full N×N symmetric matrices) plus the first derivatives, at one
// (λ, θ, pol). This is the second-order extension of tmmThicknessJacobian and
// enables true Newton refinement (Tikhonov–Tikhonravov–Trubetskov, "Second
// order optimization methods in the synthesis of multilayer coatings," Comp.
// Maths. Math. Phys. 33, 1339 (1993)).
//
// Derivation (same Abelès matrix calculus as the Jacobian, Macleod Eq.
// 2.111/2.113; pre/post decomposition Sullivan & Dobrowolski 1996):
//   [B,C] = M₀···M_{N-1}·[1,ηs];  ∂[B,C]/∂dₖ = Pre[k]·(dMₖ)·Post[k+1].
//   Mixed second partials (i < j, position-ordered):
//     ∂²[B,C]/∂dᵢ∂dⱼ = Pre[i]·dMᵢ·(M_{i+1}···M_{j-1})·dMⱼ·Post[j+1]
//   Diagonal (i = j):
//     ∂²[B,C]/∂dᵢ² = Pre[i]·(d²Mᵢ/ddᵢ²)·Post[i+1],
//     d²Mₖ/ddₖ² = Q²·[[ −cosδ,  i sinδ/η ], [ i η sinδ,  −cosδ ]],  Q ≡ (2π/λ)n cosθ
//   (d²Mₖ is the d-derivative of dMₖ = Q[[−sinδ,−i cosδ/η],[−iη cosδ,−sinδ]];
//    as δ→0 it → Q²·[[−1,0],[0,−1]], the curvature of an emerging needle.)
//
// Second-order chain rule R = |r|², r = (η₀B−C)/den, den = η₀B+C, f = 2η₀/den²:
//   drₖ = f(C dBₖ − B dCₖ),  ddenₖ = η₀ dBₖ + dCₖ
//   d²r_ij = f(dCᵢdBⱼ + C d²B_ij − dBᵢdCⱼ − B d²C_ij) − 2·drⱼ·ddenᵢ/den
//   d²R_ij = 2 Re( conj(drᵢ)drⱼ + conj(r) d²r_ij )
//   t = 2η₀/den, dtₖ = −f·ddenₖ
//   d²t_ij = −2η₀ d²den_ij/den² + 4η₀ ddenᵢddenⱼ/den³,  d²den_ij = η₀ d²B_ij + d²C_ij
//   d²T_ij = Tfac·2 Re( conj(dtᵢ)dtⱼ + conj(t) d²t_ij ),
//   d²A_ij = −(d²R + d²T) + (2 Im η₀/Re η₀)·Im d²r_ij
// evaluated, as in responseDerivative, on every derivative of [B,C] divided by
// den once; see responseSecondDerivative.
//
// Cost: O(N²) small-matrix ops per (λ,θ,pol) via cached Pre/Post + an
// incrementally-built middle product. Checked against finite differences of
// tmm() in tests/derivatives_fd.mjs.

// d²Mₖ/ddₖ² = Q²·[[ −cosδ, i sinδ/η ], [ i η sinδ, −cosδ ]], at the δ and Q of
// the matrix in the product, with Q²/η = QP and Q²η = QS as in
// layerMatrixDerivative.
function layerMatrixSecondDerivative({ delta, Q, P, S }) {
    const cD = ccos(delta), sD = csin(delta);
    const Q2 = cmul(Q, Q);
    return [
        [ cmul([-1, 0], cmul(Q2, cD)),          cmul([0, 1], cmul(cmul(Q, P), sD)) ],
        [ cmul([0, 1], cmul(cmul(Q, S), sD)),   cmul([-1, 0], cmul(Q2, cD))        ],
    ];
}

// d²R, d²T and d²A from the first derivatives (uᵢ, wᵢ), (uⱼ, wⱼ) and the mixed
// partial (U, W), each a derivative of [B,C] divided by den. With b = B/den,
// c = C/den, gₖ = η₀uₖ + wₖ = ddenₖ/den and G = η₀U + W, the chain rule above
// reads
//   d²r = 2η₀ (wᵢuⱼ + cU − uᵢwⱼ − bW) − 2 drⱼ gᵢ,   d²t = −t G + 2t gᵢgⱼ
function responseSecondDerivative(response, [ui, wi], [uj, wj], [U, W]) {
    const { eta0, r, t, b, c, Tfac, g } = response;
    const twoEta0 = cmul([2, 0], eta0);
    const dr_i = cmul(twoEta0, csub(cmul(c, ui), cmul(b, wi)));
    const dr_j = cmul(twoEta0, csub(cmul(c, uj), cmul(b, wj)));
    const g_i = cadd(cmul(eta0, ui), wi);
    const g_j = cadd(cmul(eta0, uj), wj);
    const innerR = csub(cadd(cmul(wi, uj), cmul(c, U)), cadd(cmul(ui, wj), cmul(b, W)));
    const d2r = csub(cmul(twoEta0, innerR), cmul(cmul([2, 0], dr_j), g_i));
    const d2Rij = 2 * (creal(cmul(cconj(dr_i), dr_j)) + creal(cmul(cconj(r), d2r)));
    const dt_i = cmul([-1, 0], cmul(t, g_i));
    const dt_j = cmul([-1, 0], cmul(t, g_j));
    const G = cadd(cmul(eta0, U), W);
    const d2t = cadd(cmul([-1, 0], cmul(t, G)), cmul(cmul([2, 0], t), cmul(g_i, g_j)));
    const d2Tij = Tfac * 2 * (creal(cmul(cconj(dt_i), dt_j)) + creal(cmul(cconj(t), d2t)));
    return { d2Rij, d2Tij, d2Aij: -(d2Rij + d2Tij) + g * d2r[1] };
}

export function tmmThicknessHessian(lambda_nm, theta_deg, pol, n0, ns, layers) {
    const stack = derivativeStack({ lambda_nm, theta_deg, pol, n0, ns, layers });
    const { N, data, Pre, preExp, Post, postExp } = stack;
    const response = stackResponse(stack);

    // Per-layer first-derivative pieces: dM[k], its right-applied vector
    // v[k] = dMₖ·Post[k+1] (exponent postExp[k+1]), and d[B,C]/den.
    const dM = data.map(layerMatrixDerivative);
    const v = new Array(N), first = new Array(N);
    const dRdd = new Array(N), dTdd = new Array(N), dAdd = new Array(N);
    for (let k = 0; k < N; k++) {
        v[k] = cmatvec(dM[k], Post[k + 1]);
        first[k] = perDen(response, cmatvec(Pre[k], v[k]), preExp[k] + postExp[k + 1]);
        const m = responseDerivative(response, first[k]);
        dRdd[k] = m.dR; dTdd[k] = m.dT; dAdd[k] = m.dA;
    }

    const d2Rdd = Array.from({ length: N }, () => new Array(N).fill(0));
    const d2Tdd = Array.from({ length: N }, () => new Array(N).fill(0));
    const d2Add = Array.from({ length: N }, () => new Array(N).fill(0));
    for (let i = 0; i < N; i++) {
        const Wmat_i = matmul(Pre[i], dM[i]);   // Pre[i]·dMᵢ  (used for j>i)
        // M_{i+1}···M_{j-1} with its own exponent, starting empty at j = i+1.
        let Cmid = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
        let midExp = 0;
        for (let j = i; j < N; j++) {
            const second = j === i
                ? perDen(response, cmatvec(Pre[i],
                    cmatvec(layerMatrixSecondDerivative(data[i]), Post[i + 1])), preExp[i] + postExp[i + 1])
                // ∂²[B,C]/∂dᵢ∂dⱼ = (Pre[i]·dMᵢ)·(M_{i+1}···M_{j-1})·(dMⱼ·Post[j+1])
                : perDen(response, cmatvec(Wmat_i, cmatvec(Cmid, v[j])),
                    preExp[i] + midExp + postExp[j + 1]);
            const { d2Rij, d2Tij, d2Aij } = responseSecondDerivative(response, first[i], first[j], second);
            d2Rdd[i][j] = d2Rdd[j][i] = d2Rij;
            d2Tdd[i][j] = d2Tdd[j][i] = d2Tij;
            d2Add[i][j] = d2Add[j][i] = d2Aij;
            if (j >= i + 1) {                    // advance middle: include M_j
                Cmid = matmul(Cmid, data[j].M);
                midExp += rescaleBinary(matrixEntries(Cmid));
            }
        }
    }

    return {
        R: response.R, T: response.T, A: response.A,
        dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N,
    };
}

// ── Low-level primitives ─────────────────────────────────────────────────────
// Exported for callers building transfer-matrix variants on the same
// conventions. Lower level than the four functions above, and correspondingly
// less stable across versions.

export {
    cadd, csub, cmul, cdiv, cabs2, cconj, csqrt, ccos, csin, creal, cimag,
    matmul, rescaleMatrix, snellCosTheta, incidentCosTheta, layerMatrix, cmatvec
};
