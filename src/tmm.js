/**
 * Transfer-matrix method for multilayer thin films — JavaScript reference
 * implementation.
 *
 * System model:
 *   incident medium (n0, θ0) → layer1 → … → layerN → substrate
 *
 * Conventions:
 *   ñ = n + ik              k ≥ 0 for absorbing media
 *   time factor exp(−iωt)   a wave exp(i(kz − ωt)) decays for k > 0
 *   off-diagonals of the characteristic matrix carry −i
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
 *   • Macleod, Thin-Film Optical Filters 5th ed., §2.4, Eqs. 2.111, 2.123–2.125
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

function snellCosTheta(n0, sinTheta0, nj) {
    // sinThetaJ = n0 * sinTheta0 / nj   (complex)
    const sinThetaJ = cdiv(cmul(n0, sinTheta0), nj);
    // cosTheta = sqrt(1 - sin²θ)
    return csqrt(csub([1, 0], cmul(sinThetaJ, sinThetaJ)));
}

// ── Layer characteristic matrix ───────────────────────────────────────────────

function layerMatrix(nj, dj_nm, lambda_nm, cosTheta_j, pol) {
    // Phase thickness: delta = (2π/λ) * n * d * cosθ  (complex)
    const k0 = (2 * Math.PI) / lambda_nm;
    const delta = cmul(cmul(nj, [k0 * dj_nm, 0]), cosTheta_j);

    // Numerical-overflow guard for very thick ABSORBING layers (k>0). ccos/csin
    // use cosh/sinh(Im δ), and cosh(710)=Inf → the whole TMM returns NaN. But by
    // |Im δ| ≳ a few tens the layer is already optically opaque (single-pass
    // transmittance e^{−2 Im δ} ≈ 0) AND the surface reflectance has fully
    // converged, so clamping Im δ here is exact to machine precision while
    // keeping the characteristic matrix finite. For non-absorbing / thin layers
    // (|Im δ| < MAX_IM_DELTA) this is a no-op → results are bit-identical.
    const MAX_IM_DELTA = 50; // e^{−100} ≈ 4e−44 ; cosh(50) ≈ 2.6e21 (safe under products)
    if (delta[1] > MAX_IM_DELTA) delta[1] = MAX_IM_DELTA;
    else if (delta[1] < -MAX_IM_DELTA) delta[1] = -MAX_IM_DELTA;

    const cosD = ccos(delta);
    const sinD = csin(delta);

    // Admittance eta
    let eta;
    if (pol === 's') {
        eta = cmul(nj, cosTheta_j);          // n cosθ
    } else {
        eta = cdiv(nj, cosTheta_j);          // n / cosθ
    }

    // M = [[cosD, -i sinD / eta], [-i eta sinD, cosD]]
    const iSinD_div_eta = cmul([0, -1], cdiv(sinD, eta));
    const iEta_sinD     = cmul([0, -1], cmul(eta, sinD));

    return [
        [cosD,        iSinD_div_eta],
        [iEta_sinD,   cosD         ]
    ];
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
    const sinTheta0 = [Math.sin(theta_deg * Math.PI / 180), 0];
    const cosTheta0 = csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));

    // Admittance of incident medium
    const eta0 = pol === 's'
        ? cmul(n0, cosTheta0)
        : cdiv(n0, cosTheta0);

    // Admittance of substrate
    const cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    const etaS = pol === 's'
        ? cmul(ns, cosThetaS)
        : cdiv(ns, cosThetaS);

    // Build total transfer matrix M = M1 × M2 × ... × MN
    let M = [[  [1, 0], [0, 0]  ], [  [0, 0], [1, 0]  ]]; // identity
    let logScale = 0;

    for (const { n, d } of layers) {
        if (d <= 0) continue;
        const cosThetaJ = snellCosTheta(n0, sinTheta0, n);
        const Mj = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
        M = matmul(M, Mj);
        logScale += rescaleMatrix(M);
    }

    // [B, C]^T = M × [1, eta_s]^T
    const B = cadd(M[0][0], cmul(M[0][1], etaS));
    const C = cadd(M[1][0], cmul(M[1][1], etaS));

    // r = (η0 B - C) / (η0 B + C)
    const eta0B = cmul(eta0, B);
    const r = cdiv(csub(eta0B, C), cadd(eta0B, C));

    // t = 2 η0 / (η0 B + C)
    const t = cdiv(cmul([2, 0], eta0), cadd(eta0B, C));

    const R = cabs2(r);

    // T = Re(etaS) / Re(eta0) * |t|²
    const T = Math.max(0, creal(etaS) / creal(eta0) * cabs2(t) * Math.exp(-2 * logScale));

    const A = Math.max(0, 1 - R - T);

    return { R, T, A };
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
//   dA/dd = −(dR/dd + dT/dd)
// The host layer cancels automatically through Pre/Post (a needle of the
// host index at an interior point gives ~0), so no nₐ²−n_host² term is
// needed — exactly as in Sullivan's scheme.
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
    const sinTheta0 = [Math.sin(theta_deg * Math.PI / 180), 0];
    const cosTheta0 = csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));
    const eta0 = pol === 's' ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    const cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    const etaS = pol === 's' ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    // NOTE: layers are used as-is (no d>0 filter) so gap/intra indices match
    // the caller's design.frontLayers exactly. A zero-thickness layer yields
    // an identity characteristic matrix (δ=0), which is harmless.
    const valid = layers;
    const N = valid.length;
    const cosThJ = valid.map(({ n }) => snellCosTheta(n0, sinTheta0, n));
    const Ms = valid.map(({ n, d }, k) => layerMatrix(n, Math.max(d, 0), lambda_nm, cosThJ[k], pol));

    const I = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
    // Pre[j] = M_0·…·M_{j-1};  Post[j] = M_j·…·M_{N-1}·[1,ηs]
    const Pre = new Array(N + 1);
    Pre[0] = I;
    for (let j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    const Post = new Array(N + 1);
    Post[N] = [[1, 0], etaS];
    for (let j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    // Base spectral quantities from the full matrix
    const Bv = Post[0][0], Cv = Post[0][1];
    const den = cadd(cmul(eta0, Bv), Cv);
    const den2 = cmul(den, den);
    const r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    const t = cdiv(cmul([2, 0], eta0), den);
    const R = cabs2(r);
    const Tfac = creal(etaS) / creal(eta0);
    const T = Math.max(0, Tfac * cabs2(t));
    const A = Math.max(0, 1 - R - T);

    const k0 = (2 * Math.PI) / lambda_nm;

    // Needle derivative-matrix A for a candidate index, at a given cosθ.
    function needleA(nA) {
        const cthA = snellCosTheta(n0, sinTheta0, nA);
        const etaA = pol === 's' ? cmul(nA, cthA) : cdiv(nA, cthA);
        const Q = cmul(cmul(nA, [k0, 0]), cthA);          // (2π/λ) nₐ cosθₐ
        return [
            [[0, 0], cmul([0, -1], cdiv(Q, etaA))],       // −i Q/ηₐ
            [cmul([0, -1], cmul(etaA, Q)), [0, 0]],        // −i ηₐ Q
        ];
    }

    // Given d[B,C]/dd = (dB,dC), produce {dR,dT,dA}.
    function metrics(dB, dC) {
        // dr/dd = (2η₀/den²)(C·dB − B·dC)
        const f = cdiv(cmul([2, 0], eta0), den2);
        const dr = cmul(f, csub(cmul(Cv, dB), cmul(Bv, dC)));
        const dR = 2 * creal(cmul(cconj(r), dr));
        // dt/dd = −(2η₀/den²)(η₀·dB + dC)
        const dt = cmul([-1, 0], cmul(f, cadd(cmul(eta0, dB), dC)));
        const dT = Tfac * 2 * creal(cmul(cconj(t), dt));
        return { dR, dT, dA: -(dR + dT) };
    }

    // Precompute A·Post[pos] is position-dependent; do per (pos, cand).
    const Acache = candidateNs.map(needleA);

    const gaps = new Array(N + 1);
    for (let pos = 0; pos <= N; pos++) {
        const pre = Pre[pos], post = Post[pos];
        gaps[pos] = Acache.map(Amat => {
            const dV = cmatvec(pre, cmatvec(Amat, post));
            return metrics(dV[0], dV[1]);
        });
    }

    const intra = [];
    if (intraFracs.length) {
        for (let k = 0; k < N; k++) {
            const { n, d } = valid[k];
            const cth = cosThJ[k];
            const rowK = [];
            for (const frac of intraFracs) {
                const Mleft  = layerMatrix(n, Math.max(frac * d, 1e-9),       lambda_nm, cth, pol);
                const Mright = layerMatrix(n, Math.max((1 - frac) * d, 1e-9), lambda_nm, cth, pol);
                const preIn  = matmul(Pre[k], Mleft);
                const postIn = cmatvec(Mright, Post[k + 1]);
                rowK.push({
                    frac,
                    perCand: Acache.map(Amat => {
                        const dV = cmatvec(preIn, cmatvec(Amat, postIn));
                        return metrics(dV[0], dV[1]);
                    }),
                });
            }
            intra.push(rowK);
        }
    }

    return { R, T, A, gaps, intra, N };
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
// As δ→0 this collapses to [[0,−iQ/η],[−iQη,0]] — exactly the needle
// A-matrix in tmmNeedleScan (needleA), i.e. the needle kernel is the δ=0
// special case of this; a strong internal-consistency check.
//
// The {dR,dT,dA} chain rule below is identical (verbatim) to the validated
// `metrics()` in tmmNeedleScan.
export function tmmThicknessJacobian(lambda_nm, theta_deg, pol, n0, ns, layers) {
    const sinTheta0 = [Math.sin(theta_deg * Math.PI / 180), 0];
    const cosTheta0 = csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));
    const eta0 = pol === 's' ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    const cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    const etaS = pol === 's' ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    const valid = layers;                    // used as-is (index parity, see needle)
    const N = valid.length;
    const cosThJ = valid.map(({ n }) => snellCosTheta(n0, sinTheta0, n));
    const Ms = valid.map(({ n, d }, k) =>
        layerMatrix(n, Math.max(d, 0), lambda_nm, cosThJ[k], pol));

    const I = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
    // Pre[j] = M_0·…·M_{j-1};  Post[j] = M_j·…·M_{N-1}·[1,ηs]
    const Pre = new Array(N + 1);
    Pre[0] = I;
    for (let j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    const Post = new Array(N + 1);
    Post[N] = [[1, 0], etaS];
    for (let j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    // Base spectral quantities from the full matrix (same as tmmNeedleScan).
    const Bv = Post[0][0], Cv = Post[0][1];
    const den = cadd(cmul(eta0, Bv), Cv);
    const den2 = cmul(den, den);
    const r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    const t = cdiv(cmul([2, 0], eta0), den);
    const R = cabs2(r);
    const Tfac = creal(etaS) / creal(eta0);
    const T = Math.max(0, Tfac * cabs2(t));
    const A = Math.max(0, 1 - R - T);

    // [B,C] → {dR,dT,dA}  (verbatim from validated tmmNeedleScan.metrics)
    const f = cdiv(cmul([2, 0], eta0), den2);
    function metrics(dB, dC) {
        const dr = cmul(f, csub(cmul(Cv, dB), cmul(Bv, dC)));
        const dR = 2 * creal(cmul(cconj(r), dr));
        const dt = cmul([-1, 0], cmul(f, cadd(cmul(eta0, dB), dC)));
        const dT = Tfac * 2 * creal(cmul(cconj(t), dt));
        return { dR, dT, dA: -(dR + dT) };
    }

    const k0 = (2 * Math.PI) / lambda_nm;
    const dRdd = new Array(N), dTdd = new Array(N), dAdd = new Array(N);
    for (let k = 0; k < N; k++) {
        const { n, d } = valid[k];
        const cth  = cosThJ[k];
        const etaK = pol === 's' ? cmul(n, cth) : cdiv(n, cth);
        const Q    = cmul(cmul(n, [k0, 0]), cth);             // (2π/λ) n cosθ
        const delta = cmul(cmul(n, [k0 * Math.max(d, 0), 0]), cth);
        const cD = ccos(delta), sD = csin(delta);
        // dM_k/dd_k = Q · [[ −sinδ, −i cosδ/η ], [ −i η cosδ, −sinδ ]]
        const dMk = [
            [ cmul(Q, cmul([-1, 0], sD)),                cmul(Q, cmul([0, -1], cdiv(cD, etaK))) ],
            [ cmul(Q, cmul([0, -1], cmul(etaK, cD))),    cmul(Q, cmul([-1, 0], sD))             ],
        ];
        const dV = cmatvec(Pre[k], cmatvec(dMk, Post[k + 1]));
        const m  = metrics(dV[0], dV[1]);
        dRdd[k] = m.dR; dTdd[k] = m.dT; dAdd[k] = m.dA;
    }

    return { R, T, A, dRdd, dTdd, dAdd, N };
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
// Derivation (same Abelès matrix calculus as the Jacobian — Macleod Eq.
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
//   d²T_ij = Tfac·2 Re( conj(dtᵢ)dtⱼ + conj(t) d²t_ij ),  d²A = −(d²R + d²T)
//
// Cost: O(N²) small-matrix ops per (λ,θ,pol) via cached Pre/Post + an
// incrementally-built middle product. *Must be FD-validated before trust
// (tests/hessian_fd_validation.mjs).*
export function tmmThicknessHessian(lambda_nm, theta_deg, pol, n0, ns, layers) {
    const sinTheta0 = [Math.sin(theta_deg * Math.PI / 180), 0];
    const cosTheta0 = csqrt(csub([1, 0], cmul(sinTheta0, sinTheta0)));
    const eta0 = pol === 's' ? cmul(n0, cosTheta0) : cdiv(n0, cosTheta0);
    const cosThetaS = snellCosTheta(n0, sinTheta0, ns);
    const etaS = pol === 's' ? cmul(ns, cosThetaS) : cdiv(ns, cosThetaS);

    const valid = layers;
    const N = valid.length;
    const cosThJ = valid.map(({ n }) => snellCosTheta(n0, sinTheta0, n));
    const Ms = valid.map(({ n, d }, k) =>
        layerMatrix(n, Math.max(d, 0), lambda_nm, cosThJ[k], pol));

    const I = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
    const Pre = new Array(N + 1);
    Pre[0] = I;
    for (let j = 0; j < N; j++) Pre[j + 1] = matmul(Pre[j], Ms[j]);
    const Post = new Array(N + 1);
    Post[N] = [[1, 0], etaS];
    for (let j = N - 1; j >= 0; j--) Post[j] = cmatvec(Ms[j], Post[j + 1]);

    const Bv = Post[0][0], Cv = Post[0][1];
    const den  = cadd(cmul(eta0, Bv), Cv);
    const den2 = cmul(den, den);
    const den3 = cmul(den2, den);
    const r = cdiv(csub(cmul(eta0, Bv), Cv), den);
    const t = cdiv(cmul([2, 0], eta0), den);
    const R = cabs2(r);
    const Tfac = creal(etaS) / creal(eta0);
    const T = Math.max(0, Tfac * cabs2(t));
    const A = Math.max(0, 1 - R - T);
    const f = cdiv(cmul([2, 0], eta0), den2);

    // First-derivative metrics (verbatim from tmmThicknessJacobian).
    function metrics(dB, dC) {
        const fmet = cdiv(cmul([2, 0], eta0), den2);
        const dr = cmul(fmet, csub(cmul(Cv, dB), cmul(Bv, dC)));
        const dR = 2 * creal(cmul(cconj(r), dr));
        const dt = cmul([-1, 0], cmul(fmet, cadd(cmul(eta0, dB), dC)));
        const dT = Tfac * 2 * creal(cmul(cconj(t), dt));
        return { dR, dT, dA: -(dR + dT) };
    }

    const k0 = (2 * Math.PI) / lambda_nm;
    // Per-layer first-derivative pieces: dM[k], its right-applied vector v[k] =
    // dMₖ·Post[k+1], the [dB,dC] vector, and the diagonal second-derivative
    // matrix d2M[k].
    const dM = new Array(N), d2M = new Array(N), v = new Array(N);
    const dB = new Array(N), dC = new Array(N);
    const dRdd = new Array(N), dTdd = new Array(N), dAdd = new Array(N);
    for (let k = 0; k < N; k++) {
        const { n, d } = valid[k];
        const cth  = cosThJ[k];
        const etaK = pol === 's' ? cmul(n, cth) : cdiv(n, cth);
        const Q    = cmul(cmul(n, [k0, 0]), cth);              // (2π/λ) n cosθ
        const Q2   = cmul(Q, Q);
        const delta = cmul(cmul(n, [k0 * Math.max(d, 0), 0]), cth);
        const cD = ccos(delta), sD = csin(delta);
        // dMₖ/ddₖ = Q·[[ −sinδ, −i cosδ/η ], [ −i η cosδ, −sinδ ]]
        dM[k] = [
            [ cmul(Q, cmul([-1, 0], sD)),             cmul(Q, cmul([0, -1], cdiv(cD, etaK))) ],
            [ cmul(Q, cmul([0, -1], cmul(etaK, cD))), cmul(Q, cmul([-1, 0], sD))             ],
        ];
        // d²Mₖ/ddₖ² = Q²·[[ −cosδ, i sinδ/η ], [ i η sinδ, −cosδ ]]
        d2M[k] = [
            [ cmul(Q2, cmul([-1, 0], cD)),            cmul(Q2, cmul([0, 1], cdiv(sD, etaK))) ],
            [ cmul(Q2, cmul([0, 1], cmul(etaK, sD))), cmul(Q2, cmul([-1, 0], cD))            ],
        ];
        v[k] = cmatvec(dM[k], Post[k + 1]);
        const dVk = cmatvec(Pre[k], v[k]);
        dB[k] = dVk[0]; dC[k] = dVk[1];
        const m = metrics(dB[k], dC[k]);
        dRdd[k] = m.dR; dTdd[k] = m.dT; dAdd[k] = m.dA;
    }

    // Second-order metrics for a pair (i,j) given the mixed partial [d2B,d2C].
    function hessMetrics(i, j, d2Bv, d2Cv) {
        const dBi = dB[i], dCi = dC[i], dBj = dB[j], dCj = dC[j];
        const dr_i = cmul(f, csub(cmul(Cv, dBi), cmul(Bv, dCi)));
        const dr_j = cmul(f, csub(cmul(Cv, dBj), cmul(Bv, dCj)));
        const dden_i = cadd(cmul(eta0, dBi), dCi);
        const dden_j = cadd(cmul(eta0, dBj), dCj);
        // d²r_ij = f(dCᵢdBⱼ + C d²B − dBᵢdCⱼ − B d²C) − 2 drⱼ ddenᵢ/den
        const innerR = csub(
            cadd(cmul(dCi, dBj), cmul(Cv, d2Bv)),
            cadd(cmul(dBi, dCj), cmul(Bv, d2Cv))
        );
        const d2r = csub(cmul(f, innerR),
                         cdiv(cmul(cmul([2, 0], dr_j), dden_i), den));
        const d2Rij = 2 * (creal(cmul(cconj(dr_i), dr_j)) + creal(cmul(cconj(r), d2r)));
        // T
        const dt_i = cmul([-1, 0], cmul(f, dden_i));
        const dt_j = cmul([-1, 0], cmul(f, dden_j));
        const d2den = cadd(cmul(eta0, d2Bv), d2Cv);
        const d2t = cadd(
            cmul(cmul([-2, 0], eta0), cdiv(d2den, den2)),
            cmul(cmul([4, 0], eta0), cdiv(cmul(dden_i, dden_j), den3))
        );
        const d2Tij = Tfac * 2 * (creal(cmul(cconj(dt_i), dt_j)) + creal(cmul(cconj(t), d2t)));
        return { d2Rij, d2Tij };
    }

    const d2Rdd = Array.from({ length: N }, () => new Array(N).fill(0));
    const d2Tdd = Array.from({ length: N }, () => new Array(N).fill(0));
    const d2Add = Array.from({ length: N }, () => new Array(N).fill(0));
    for (let i = 0; i < N; i++) {
        const Wmat_i = matmul(Pre[i], dM[i]);   // Pre[i]·dMᵢ  (used for j>i)
        let Cmid = I;                            // M_{i+1}···M_{j-1}, starts empty at j=i+1
        for (let j = i; j < N; j++) {
            let d2Bv, d2Cv;
            if (j === i) {
                const w = cmatvec(Pre[i], cmatvec(d2M[i], Post[i + 1]));
                d2Bv = w[0]; d2Cv = w[1];
            } else {
                // ∂²[B,C]/∂dᵢ∂dⱼ = (Pre[i]·dMᵢ)·(M_{i+1}···M_{j-1})·(dMⱼ·Post[j+1])
                const w = cmatvec(Wmat_i, cmatvec(Cmid, v[j]));
                d2Bv = w[0]; d2Cv = w[1];
            }
            const { d2Rij, d2Tij } = hessMetrics(i, j, d2Bv, d2Cv);
            d2Rdd[i][j] = d2Rdd[j][i] = d2Rij;
            d2Tdd[i][j] = d2Tdd[j][i] = d2Tij;
            d2Add[i][j] = d2Add[j][i] = -(d2Rij + d2Tij);
            if (j >= i + 1) Cmid = matmul(Cmid, Ms[j]); // advance middle: include M_j
        }
    }

    return { R, T, A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N };
}

// ── Low-level primitives ─────────────────────────────────────────────────────
// Exported for callers building transfer-matrix variants on the same
// conventions. Lower level than the four functions above, and correspondingly
// less stable across versions.

export {
    cadd, csub, cmul, cdiv, cabs2, cconj, csqrt, ccos, csin, creal, cimag,
    matmul, rescaleMatrix, snellCosTheta, layerMatrix, cmatvec
};
