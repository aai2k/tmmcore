/**
 * A growing coating in JavaScript: the monitoring curve of one layer, the
 * spectra after each deposited layer, and a growing top layer over a
 * wavelength grid. Line-by-line ports of tmm_monitor_curve,
 * tmm_deposition_spectra and tmm_growing_stack_* in src/c/growing.c.
 *
 * The incremental control algorithm of Tikhonravov & Trubetskov, Appl. Opt.
 * 44, 6877 (2005): the finished layers beneath the growing one do not change
 * while it grows, so their characteristic-matrix product is built once and
 * each sample costs one layer matrix, one 2×2 product and the tails. The
 * growing layer faces the incident medium, so its matrix multiplies the
 * finished product from the left.
 *
 * Each call returns, per polarization, the forward R and T of the coated
 * surface (Rs, Ts, Rp, Tp) and its reflectance from the substrate side (Rrs,
 * Rrp), for a caller that treats the witness as an incoherent slab with a bare
 * back face. A is not formed here: such a caller forms it after the slab sum,
 * and one with a semi-infinite substrate takes 1 − R − T.
 *
 * Grids take the layout of stack.js: lambdas, theta_deg (one angle or one per
 * wavelength), n0List, nsList, layerNK ([layer][λ] ñ) and thick (nm).
 */

import {
    cadd, csub, cmul, cdiv, cabs2, creal, matmul, rescaleMatrix, cmatvec,
    incidence, snellCosTheta, incidentCosTheta, layerAdmittance, substrateVector, transmittedFlux,
    layerMatrix, layerLogScale,
} from './primitives.js';
import { thetaAt } from './stack.js';

const POLS = ['s', 'p'];
const OUTPUTS = ['Rs', 'Ts', 'Rp', 'Tp', 'Rrs', 'Rrp'];

const identity = () => [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];

function outputs(length) {
    const out = {};
    for (const key of OUTPUTS) out[key] = new Float64Array(length);
    return out;
}

// The output arrays of one polarization: forward R, forward T, reverse R.
const sideOf = (out, pol) => (pol === 's' ? [out.Rs, out.Ts, out.Rrs] : [out.Rp, out.Tp, out.Rrp]);

// Forward R and T of a built product: the [B, C] → r, t → R, T block of tmm().
function forwardTail(M, eta0, substrate, logScale) {
    const [B, C] = cmatvec(M, substrate);
    const eta0B = cmul(eta0, B);
    const r = cdiv(csub(eta0B, C), cadd(eta0B, C));
    const t = cmul(cdiv(cmul([2, 0], eta0), cadd(eta0B, C)), substrate[0]);
    const R = cabs2(r);
    const T = Math.max(0, transmittedFlux(substrate) / creal(eta0) * cabs2(t) * Math.exp(-2 * logScale));
    return { R, T };
}

// Reflectance of the same product seen from the substrate side. The reversed
// stack's matrix is the anti-transpose of the forward one, since every layer
// matrix is unchanged by anti-transposition, so this pass costs no second
// product; the common rescale factor cancels in the amplitude ratio. The
// substrate admittance cs/bs is the incident one of this pass.
function reverseTail(M, eta0, substrate) {
    const B = cadd(M[1][1], cmul(M[0][1], eta0));
    const C = cadd(M[1][0], cmul(M[0][0], eta0));
    const etaSB = cmul(substrate[1], B);
    const bsC = cmul(substrate[0], C);
    return cabs2(cdiv(csub(etaSB, bsC), cadd(etaSB, bsC)));
}

// The incident admittance and the substrate vector of one polarization.
function boundaries(n0, ns, angle, cosTheta0, pol) {
    return {
        eta0: layerAdmittance(n0, cosTheta0, pol),
        substrate: substrateVector(ns, snellCosTheta(n0, angle.sinTheta0, ns, angle.cosTheta0), pol),
    };
}

// One layer of thickness d folded into a product: on the right for a stack
// listed outermost first, on the left for a layer deposited on top.
function fold(product, n, d, at, onTop) {
    const { n0, angle, lambda_nm, pol } = at;
    const cosThetaJ = snellCosTheta(n0, angle.sinTheta0, n, angle.cosTheta0);
    const layer = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
    const M = onTop ? matmul(layer, product.M) : matmul(product.M, layer);
    const logScale = product.logScale + layerLogScale(n, d, lambda_nm, cosThetaJ) + rescaleMatrix(M);
    return { M, logScale };
}

// Writes the forward and reverse tails of a product to entry `at` of one side.
function writeTails([R, T, Rr], at, product, { eta0, substrate }) {
    const forward = forwardTail(product.M, eta0, substrate, product.logScale);
    R[at] = forward.R;
    T[at] = forward.T;
    Rr[at] = reverseTail(product.M, eta0, substrate);
}

/**
 * The signal of one layer as it grows on a finished stack, at one wavelength.
 * `base` is the finished stack as { n: [re, im], d } outermost first, zero
 * thicknesses skipped; `top` the growing layer's ñ; `thicknesses` its sample
 * thicknesses (nm), where d ≤ 0 gives the bare finished stack.
 * Returns { Rs, Ts, Rp, Tp, Rrs, Rrp }, one value per sample thickness.
 */
export function monitorCurveKernel({ lambda_nm, theta_deg, n0, ns, base, top, thicknesses }) {
    const out = outputs(thicknesses.length);
    const angle = incidence(theta_deg);
    const cosTheta0 = incidentCosTheta(n0, angle.sinTheta0, angle.cosTheta0);
    for (const pol of POLS) {
        const at = { n0, angle, lambda_nm, pol };
        const ends = boundaries(n0, ns, angle, cosTheta0, pol);
        let finished = { M: identity(), logScale: 0 };
        for (const { n, d } of base) {
            if (d > 0) finished = fold(finished, n, d, at, false);
        }
        const side = sideOf(out, pol);
        thicknesses.forEach((d, k) => {
            writeTails(side, k, d > 0 ? fold(finished, top, d, at, true) : finished, ends);
        });
    }
    return out;
}

/**
 * The coated surface after each deposited layer, over a wavelength grid. The
 * layers arrive in deposition order, first deposited (substrate side) first;
 * each folds into the running product from the left, and the tails after
 * every fold are written out. A step of zero thickness repeats the previous
 * one. Returns { Rs, Ts, Rp, Tp, Rrs, Rrp }, each N × nLam, step-major: step
 * k at wavelength i at [k·nLam + i].
 */
export function depositionSpectraKernel(grid) {
    const { lambdas, n0List, nsList, layerNK, thick } = grid;
    const nLam = lambdas.length;
    const out = outputs(thick.length * nLam);
    const states = Array.from(lambdas, (lambda_nm, i) => {
        const n0 = n0List[i];
        const angle = incidence(thetaAt(grid.theta_deg, i));
        const cosTheta0 = incidentCosTheta(n0, angle.sinTheta0, angle.cosTheta0);
        return POLS.map(pol => ({
            at: { n0, angle, lambda_nm, pol },
            ends: boundaries(n0, nsList[i], angle, cosTheta0, pol),
            product: { M: identity(), logScale: 0 },
        }));
    });
    thick.forEach((d, k) => {
        for (let i = 0; i < nLam; i++) {
            for (const state of states[i]) {
                if (d > 0) state.product = fold(state.product, layerNK[k][i], d, state.at, true);
                writeTails(sideOf(out, state.at.pol), k * nLam + i, state.product, state.ends);
            }
        }
    });
    return out;
}

/**
 * A growing top layer over a wavelength grid: the finished stack, outermost
 * first, folds once per wavelength and polarization; setTop(ngList) declares
 * the growing layer's ñ per wavelength; sample(d, out?) returns
 * { Rs, Ts, Rp, Tp, Rrs, Rrp } over the grid for thickness d of it, written
 * into `out` when given. d ≤ 0 samples the finished stack alone and needs no
 * top. After free(), setTop and sample throw; free() again does nothing.
 */
export function growingStackKernel(grid) {
    const { lambdas, n0List, nsList, layerNK, thick } = grid;
    const nLam = lambdas.length;
    const states = Array.from(lambdas, (lambda_nm, i) => {
        const n0 = n0List[i];
        const angle = incidence(thetaAt(grid.theta_deg, i));
        const cosTheta0 = incidentCosTheta(n0, angle.sinTheta0, angle.cosTheta0);
        return POLS.map(pol => {
            const at = { n0, angle, lambda_nm, pol };
            let finished = { M: identity(), logScale: 0 };
            thick.forEach((d, k) => {
                if (d > 0) finished = fold(finished, layerNK[k][i], d, at, false);
            });
            return { at, ends: boundaries(n0, nsList[i], angle, cosTheta0, pol), finished };
        });
    });
    let top = null;
    let freed = false;
    const live = () => { if (freed) throw new Error('tmmcore: growing stack already freed'); };
    return {
        setTop(ngList) {
            live();
            top = Array.from(ngList);
        },
        sample(d, out = null) {
            live();
            if (d > 0 && !top) throw new Error('tmmcore: growing stack sampled before setTop');
            const result = out || outputs(nLam);
            for (let i = 0; i < nLam; i++) {
                for (const { at, ends, finished } of states[i]) {
                    const product = d > 0 ? fold(finished, top[i], d, at, true) : finished;
                    writeTails(sideOf(result, at.pol), i, product, ends);
                }
            }
            return result;
        },
        free() { freed = true; },
    };
}
