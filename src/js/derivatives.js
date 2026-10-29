/**
 * derivatives.js : exact analytic derivatives of R, T and A with respect to
 * the layer thicknesses, at one wavelength and over a wavelength grid: the
 * thickness Jacobian, the thickness Hessian, and the needle scan, the
 * derivative of inserting a thin layer at any depth. Conventions are in
 * primitives.js.
 *
 * References:
 *   • Sullivan & Dobrowolski, Appl. Opt. 35, 5484 (1996), Eqs. (3)–(6)
 *   • Tikhonravov, Trubetskov & DeBell, Appl. Opt. 35, 5493 (1996)
 */

import {
    cadd, csub, cmul, cdiv, cabs2, cconj, ccos, csin, creal, cscale, matmul, cmatvec, incidence,
    snellCosTheta, incidentCosTheta, layerPhase, layerAdmittance, phaseMatrix, atCriticalAngle,
    criticalGenerator, generatorMatrix, criticalMatrix, substrateVector, transmittedFlux,
    layerAbsorptance,
} from './primitives.js';
import { checkLayerKeys, pointArgs } from './stack.js';

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
    const { delta, excess } = layerPhase(n, thickness, lambda_nm, cosTheta);
    const eta = layerAdmittance(n, cosTheta, pol);
    // Q = dδ/dd = (2π/λ) n cosθ. Past the bound the layer is the held matrix
    // times e^{excess}, and together they move with d as the true matrix does:
    // Q in full, with the held δ in the sines and cosines, is their derivative
    // to within e^{−100} relative. Its imaginary part is what makes T fall
    // with d; it cancels from r.
    const Q = present ? cmul(cmul(n, [k0, 0]), cosTheta) : [0, 0];
    const { P, S } = present
        ? generator(n, cosTheta, Q, k0, pol)
        : { P: [0, 0], S: [0, 0] };
    const layer = {
        n, thickness, cosTheta, delta, excess, eta, critical: atCriticalAngle(cosTheta), Q, P, S,
    };
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
    checkLayerKeys(layers, 'n', 'd');
    const { sinTheta0, cosTheta0 } = incidence(theta_deg);
    const eta0 = layerAdmittance(n0, incidentCosTheta(n0, sinTheta0, cosTheta0), pol);
    const substrate = substrateVector(ns, snellCosTheta(n0, sinTheta0, ns, cosTheta0), pol);
    const k0 = (2 * Math.PI) / lambda_nm;
    const N = layers.length;
    const geometry = { n0, sinTheta0, cosTheta0, k0, lambda_nm, pol };

    const data = layers.map(({ n, d }) => derivativeLayer(n, d, geometry));
    // What the held layers leave out of the product, for t.
    const logScale = data.reduce((sum, layer) => sum + layer.excess, 0);

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

    return { geometry, eta0, substrate, N, data, Pre, preExp, Post, postExp, logScale };
}

// R, T, A of the whole stack from [B, C] = Post[0] (Macleod 5th ed., Eqs.
// 2.123–2.125), with the pieces the chain rule below needs: r, the true t,
// b = B/den and c = C/den, and g = 2 Im η0 / Re η0, the weight of Im r in the
// absorptance. Only t carries the exponent of [B, C] and the log scale of the
// held layers.
function stackResponse({ eta0, substrate, Post, postExp, logScale }) {
    const [B, C] = Post[0];
    const den = cadd(cmul(eta0, B), C);
    const inverseDen = cdiv([1, 0], den);
    const r = cdiv(csub(cmul(eta0, B), C), den);
    let t = cscale(cmul(cdiv(cmul([2, 0], eta0), den), substrate[0]), 2 ** -postExp[0]);
    if (logScale > 0) t = cscale(t, Math.exp(-logScale));
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
//     off-diagonals carry −i, see layerMatrix in primitives.js).
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
    // `factor` is a log-scale correction to the pair, 1 unless it straddles a
    // held layer.
    const insertAt = (pre, post, exponent, factor = 1) => Acache.map(Amat => {
        let v = cmatvec(pre, cmatvec(Amat, post));
        if (factor !== 1) v = [cscale(v[0], factor), cscale(v[1], factor)];
        return responseDerivative(response, perDen(response, v, exponent));
    });

    const gaps = new Array(N + 1);
    for (let pos = 0; pos <= N; pos++) {
        gaps[pos] = insertAt(Pre[pos], Post[pos], preExp[pos] + postExp[pos]);
    }

    const intra = [];
    if (intraFracs.length) {
        for (let k = 0; k < N; k++) {
            const host = data[k];
            const { n, thickness, cosTheta, delta } = host;
            intra.push(Array.from(intraFracs, frac => {
                const exponent = preExp[k] + postExp[k + 1];
                const front = layerPhase(n, frac * thickness, lambda_nm, cosTheta);
                const preIn = matmul(Pre[k], partMatrix(host, front.delta, frac * thickness));
                if (host.excess === 0) {
                    // The halves multiply to M_k itself, so the needle sits in
                    // the stack R and T came from.
                    const postIn = cmatvec(partMatrix(host, csub(delta, front.delta),
                        (1 - frac) * thickness), Post[k + 1]);
                    return { frac, perCand: insertAt(preIn, postIn, exponent) };
                }
                // Past the bound each half is held on its own. The halves times
                // e^{front.excess + back.excess} are the true host, and the
                // product carries the host as its held matrix, the true one over
                // e^{host.excess}, so the pair takes the ratio: between e^{−50}
                // and 1. The field at the needle's depth is then the true one.
                const back = layerPhase(n, (1 - frac) * thickness, lambda_nm, cosTheta);
                const postIn = cmatvec(partMatrix(host, back.delta, (1 - frac) * thickness), Post[k + 1]);
                const factor = Math.exp(front.excess + back.excess - host.excess);
                return { frac, perCand: insertAt(preIn, postIn, exponent, factor) };
            }));
        }
    }

    return { R: response.R, T: response.T, A: response.A, gaps, intra, N };
}

// ── Analytic thickness-Jacobian kernel ────────────────────────────────────────
//
// Returns the EXACT analytic derivatives dR/dd_k, dT/dd_k, dA/dd_k of every
// existing layer's thickness, at one (λ, θ, pol), from one set of prefix and
// suffix products rather than the 2·N stack evaluations of central differences.
//
// Derivation (citations):
//   • Characteristic matrix Eq. 2.111 and product form Eq. 2.113 with
//       δ_r = 2π N_r d_r cosθ_r / λ
//     Macleod, Thin-Film Optical Filters 5th ed., §2.4 (verified verbatim).
//     This package's documented sign convention puts −i on the off-diagonals
//     (see layerMatrix in primitives.js); the derivative below is taken of THAT
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
// bound, the held δ and Q in full (see derivativeLayer). The off-diagonals are
// written with P = Q/η and S = Qη, which stay finite at the critical angle.
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

// ── Over a wavelength grid ───────────────────────────────────────────────────

const FIRST = ['dRdd', 'dTdd', 'dAdd'];
const SECOND = ['d2Rdd', 'd2Tdd', 'd2Add'];

// R, T, A per wavelength and the first derivatives, nLam × N, from one kernel
// call per wavelength; with `second` also the second derivatives, nLam × N × N.
function derivativeGrid(kernel, second, grid, pol) {
    const nLam = grid.lambdas.length;
    const N = grid.thick.length;
    const out = { R: new Float64Array(nLam), T: new Float64Array(nLam), A: new Float64Array(nLam), N };
    for (const key of FIRST) out[key] = new Float64Array(nLam * N);
    if (second) for (const key of SECOND) out[key] = new Float64Array(nLam * N * N);
    for (let i = 0; i < nLam; i++) {
        const one = kernel(...pointArgs(grid, pol, i));
        out.R[i] = one.R; out.T[i] = one.T; out.A[i] = one.A;
        for (const key of FIRST) out[key].set(one[key], i * N);
        if (!second) continue;
        for (const key of SECOND) one[key].forEach((row, j) => out[key].set(row, (i * N + j) * N));
    }
    return out;
}

/**
 * The thickness Jacobian at every wavelength, for one polarization:
 * { R, T, A } per wavelength and dRdd, dTdd, dAdd, nLam × N.
 */
export function jacobianSpectrumKernel(grid, pol) {
    return derivativeGrid(tmmThicknessJacobian, false, grid, pol);
}

/**
 * The thickness Hessian at every wavelength, for one polarization: the
 * Jacobian's fields and d2Rdd, d2Tdd, d2Add, nLam × N × N, the value for
 * wavelength i and layers j, k at [(i·N + j)·N + k].
 */
export function hessianSpectrumKernel(grid, pol) {
    return derivativeGrid(tmmThicknessHessian, true, grid, pol);
}

/**
 * The needle scan at every wavelength, for one polarization, laid out as
 * tmm_needle_scan_spectrum writes it: base nLam × 3 [R, T, A], gaps
 * nLam × (N + 1) × nCand × 3 and intra nLam × N × nFrac × nCand × 3, the
 * innermost triple [dR, dT, dA]. candNK is [candidate][λ] ñ; fracs the split
 * fractions, the same at every wavelength.
 */
export function needleScanSpectrumKernel(grid, pol, { candNK, fracs }) {
    const nLam = grid.lambdas.length;
    const N = grid.thick.length;
    const nCand = candNK.length;
    const nFrac = fracs.length;
    const gapBlock = (N + 1) * nCand * 3;
    const intraBlock = N * nFrac * nCand * 3;
    const out = {
        base: new Float64Array(nLam * 3),
        gaps: new Float64Array(nLam * gapBlock),
        intra: new Float64Array(nLam * intraBlock),
    };
    const put = (array, at, { dR, dT, dA }) => { array[at] = dR; array[at + 1] = dT; array[at + 2] = dA; };
    for (let i = 0; i < nLam; i++) {
        const scan = tmmNeedleScan(...pointArgs(grid, pol, i), candNK.map(row => row[i]), fracs);
        out.base.set([scan.R, scan.T, scan.A], 3 * i);
        scan.gaps.forEach((row, pos) => row.forEach((entry, c) =>
            put(out.gaps, i * gapBlock + (pos * nCand + c) * 3, entry)));
        scan.intra.forEach((row, k) => row.forEach(({ perCand }, fi) => perCand.forEach((entry, c) =>
            put(out.intra, i * intraBlock + ((k * nFrac + fi) * nCand + c) * 3, entry))));
    }
    return out;
}
