/**
 * derivatives.js : the thickness Jacobian, the thickness Hessian and the
 * needle scan on a WebAssembly instance, at one wavelength and over a
 * wavelength grid, with the methods of 0.4 that run them.
 */

import { checkLayerKeys } from '../js/stack.js';
import { warnOnce } from '../js/deprecated.js';
import { onGrid, named, PAIRS } from './memory.js';

const FIRST = ['dRdd', 'dTdd', 'dAdd'];
const SECOND = ['d2Rdd', 'd2Tdd', 'd2Add'];

export const derivativeMethods = {
    /**
     * Analytic thickness Jacobian for one (λ, θ, pol) : mirrors
     * tmmThicknessJacobian(). layers used AS-IS (index parity).
     * @returns {{R,T,A, dRdd:Float64Array, dTdd, dAdd, N}}
     */
    _jacobian(lambda_nm, theta_deg, polCode, { n0, ns, layers }) {
        checkLayerKeys(layers, 'n', 'd');
        const N = layers.length;
        const M = Math.max(1, N);
        // arena layout: layers[3N] | dRdd[M] | dTdd[M] | dAdd[M] | base[3]
        const oLay = 0, oDR = 3 * N, oDT = oDR + M, oDA = oDT + M, oBase = oDA + M;
        const need = oBase + 3;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            buf[3 * i + 0] = layers[i].n[0];
            buf[3 * i + 1] = layers[i].n[1];
            buf[3 * i + 2] = layers[i].d;
        }
        const P = (off) => ptr + off * 8;
        this._tmm_jacobian(lambda_nm, theta_deg, polCode | 0,
            n0[0], n0[1], ns[0], ns[1], P(oLay), N, P(oDR), P(oDT), P(oDA), P(oBase));
        // Re-create the view AFTER the kernel call: under ALLOW_MEMORY_GROWTH
        // the kernel may grow wasm memory, which detaches the ArrayBuffer `buf`
        // was created over : reading the stale `buf` then yields garbage or
        // throws. `out` is a fresh view over the current buffer.
        const out = this._view(ptr, need);
        return {
            R: out[oBase], T: out[oBase + 1], A: out[oBase + 2], N,
            dRdd: out.slice(oDR, oDR + N),
            dTdd: out.slice(oDT, oDT + N),
            dAdd: out.slice(oDA, oDA + N),
        };
    },

    /**
     * Deprecated, removed in 0.6.0: use jacobian(stack, wavelength,
     * { angle, pol }). Arguments as tmmOne.
     * @returns {{R,T,A, dRdd:Float64Array, dTdd, dAdd, N}}
     */
    tmmJacobian(lambda_nm, theta_deg, polCode, n0, ns, layers) {
        warnOnce('tmmJacobian', 'jacobian');
        return this._jacobian(lambda_nm, theta_deg, polCode, { n0, ns, layers });
    },

    /** True: every build since 0.5.0 carries the Hessian kernel. */
    hasHessian() { return true; },

    /**
     * Analytic thickness Hessian for one (λ, θ, pol) : mirrors
     * tmmThicknessHessian(). Returns first AND second derivatives; the N×N
     * second-derivative blocks are reshaped into nested arrays (one Float64Array
     * row per layer, FULL symmetric) so the shape matches the JS oracle exactly.
     * @returns {{R,T,A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N}}
     */
    _hessian(lambda_nm, theta_deg, polCode, { n0, ns, layers }) {
        checkLayerKeys(layers, 'n', 'd');
        const N = layers.length;
        const M = Math.max(1, N);
        const NN = Math.max(1, N * N);
        // arena: layers[3N] | dRdd[M] | dTdd[M] | dAdd[M] | d2R[NN] | d2T[NN] | d2A[NN] | base[3]
        const oLay = 0, oDR = 3 * N, oDT = oDR + M, oDA = oDT + M,
              oR2 = oDA + M, oT2 = oR2 + NN, oA2 = oT2 + NN, oBase = oA2 + NN;
        const need = oBase + 3;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            buf[3 * i + 0] = layers[i].n[0];
            buf[3 * i + 1] = layers[i].n[1];
            buf[3 * i + 2] = layers[i].d;
        }
        const P = (off) => ptr + off * 8;
        this._tmm_hessian(lambda_nm, theta_deg, polCode | 0,
            n0[0], n0[1], ns[0], ns[1], P(oLay), N,
            P(oDR), P(oDT), P(oDA), P(oR2), P(oT2), P(oA2), P(oBase));
        // Fresh view after the call (memory may have grown → buf detached).
        const out = this._view(ptr, need);
        const reshape = (off) => {
            const rows = new Array(N);
            for (let i = 0; i < N; i++) rows[i] = out.slice(off + i * N, off + i * N + N);
            return rows;
        };
        return {
            R: out[oBase], T: out[oBase + 1], A: out[oBase + 2], N,
            dRdd: out.slice(oDR, oDR + N),
            dTdd: out.slice(oDT, oDT + N),
            dAdd: out.slice(oDA, oDA + N),
            d2Rdd: reshape(oR2), d2Tdd: reshape(oT2), d2Add: reshape(oA2),
        };
    },

    /**
     * Deprecated, removed in 0.6.0: use hessian(stack, wavelength,
     * { angle, pol }). Arguments as tmmOne.
     * @returns {{R,T,A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N}}
     */
    tmmHessian(lambda_nm, theta_deg, polCode, n0, ns, layers) {
        warnOnce('tmmHessian', 'hessian');
        return this._hessian(lambda_nm, theta_deg, polCode, { n0, ns, layers });
    },

    /**
     * Analytic needle P-function scan : mirrors tmmNeedleScan() in src/js/derivatives.js,
     * reshaping the flat kernel output into the same nested structure.
     * @param {{n:[number,number],d:number}[]} layers   used AS-IS (index parity)
     * @param {[number,number][]} candidateNs           candidate ñ
     * @param {number[]} intraFracs                      intra-layer split fractions
     * @returns {{R,T,A,N, gaps:Array, intra:Array}}
     *   gaps[pos][ci] = {dR,dT,dA}  (pos = 0..N)
     *   intra[k][fi]  = {frac, perCand:[{dR,dT,dA}]}
     */
    tmmNeedleScan(lambda_nm, theta_deg, polCode, n0, ns, layers, candidateNs, intraFracs = []) {
        checkLayerKeys(layers, 'n', 'd');
        const N = layers.length;
        const nCand = candidateNs.length;
        const nFrac = intraFracs.length;
        const nGap = (N + 1) * nCand * 3;
        const nIntra = Math.max(1, N * nFrac * nCand * 3);

        const layPtr = this._alloc(Math.max(1, 3 * N));
        const candPtr = this._alloc(Math.max(1, 2 * nCand));
        const fracPtr = this._alloc(Math.max(1, nFrac));
        const basePtr = this._alloc(3);
        const gapPtr = this._alloc(Math.max(1, nGap));
        const intraPtr = this._alloc(nIntra);

        const lay = this._view(layPtr, Math.max(1, 3 * N));
        for (let i = 0; i < N; i++) {
            lay[3 * i + 0] = layers[i].n[0];
            lay[3 * i + 1] = layers[i].n[1];
            lay[3 * i + 2] = layers[i].d;
        }
        const cand = this._view(candPtr, Math.max(1, 2 * nCand));
        for (let c = 0; c < nCand; c++) { cand[2 * c] = candidateNs[c][0]; cand[2 * c + 1] = candidateNs[c][1]; }
        const frac = this._view(fracPtr, Math.max(1, nFrac));
        for (let i = 0; i < nFrac; i++) frac[i] = intraFracs[i];

        this._tmm_needle_scan(lambda_nm, theta_deg, polCode | 0,
            n0[0], n0[1], ns[0], ns[1], layPtr, N, candPtr, nCand, fracPtr, nFrac,
            basePtr, gapPtr, intraPtr);

        // Copy outputs out before freeing, reshaping to the JS nested layout.
        const base = this._view(basePtr, 3);
        const R = base[0], T = base[1], A = base[2];
        const gapV = this._view(gapPtr, Math.max(1, nGap));
        const gaps = new Array(N + 1);
        for (let pos = 0; pos <= N; pos++) {
            const row = new Array(nCand);
            for (let c = 0; c < nCand; c++) {
                const o = (pos * nCand + c) * 3;
                row[c] = { dR: gapV[o], dT: gapV[o + 1], dA: gapV[o + 2] };
            }
            gaps[pos] = row;
        }
        const intra = [];
        if (nFrac > 0) {
            const intraV = this._view(intraPtr, nIntra);
            for (let k = 0; k < N; k++) {
                const rowK = [];
                for (let fi = 0; fi < nFrac; fi++) {
                    const perCand = new Array(nCand);
                    for (let c = 0; c < nCand; c++) {
                        const o = ((k * nFrac + fi) * nCand + c) * 3;
                        perCand[c] = { dR: intraV[o], dT: intraV[o + 1], dA: intraV[o + 2] };
                    }
                    rowK.push({ frac: intraFracs[fi], perCand });
                }
                intra.push(rowK);
            }
        }

        for (const p of [layPtr, candPtr, fracPtr, basePtr, gapPtr, intraPtr]) this.free(p);
        return { R, T, A, gaps, intra, N };
    },
};

// R, T and A per wavelength, from the kernel's [R, T, A] triple per wavelength.
function splitBase(values, nLam) {
    const [R, T, A] = [0, 1, 2].map(q => Float64Array.from({ length: nLam }, (_, i) => values[3 * i + q]));
    return { R, T, A };
}

function derivativesWasm(wasm, grid, polCode, second) {
    const nLam = grid.lambdas.length;
    const N = grid.thick.length;
    const keys = second ? [...FIRST, ...SECOND] : FIRST;
    const sizes = [...keys.map((_, q) => (q < 3 ? nLam * N : nLam * N * N)), 3 * nLam];
    const kernel = second ? wasm._tmm_hessian_spectrum : wasm._tmm_jacobian_spectrum;
    return onGrid(wasm, grid, PAIRS, ({ lead, out, read }) => {
        const outs = sizes.map(size => out(size));
        kernel(...lead(polCode), ...outs);
        const values = outs.map((ptr, i) => read(ptr, sizes[i]));
        return { ...splitBase(values.pop(), nLam), ...named(keys, values), N };
    });
}

/** The thickness Jacobian over the grid for one polarization, as jacobianSpectrumKernel. */
export function jacobianSpectrumWasm(wasm, grid, polCode) {
    return derivativesWasm(wasm, grid, polCode, false);
}

/** The thickness Hessian over the grid for one polarization, as hessianSpectrumKernel. */
export function hessianSpectrumWasm(wasm, grid, polCode) {
    return derivativesWasm(wasm, grid, polCode, true);
}

/**
 * The needle scan over the grid for one polarization, as
 * needleScanSpectrumKernel: { base, gaps, intra } in the kernel's layout.
 * candNK is [candidate][λ] ñ; fracs the split fractions.
 */
export function needleScanSpectrumWasm(wasm, grid, polCode, { candNK, fracs }) {
    const nLam = grid.lambdas.length;
    const N = grid.thick.length;
    const nCand = candNK.length;
    const sizes = [3 * nLam, nLam * (N + 1) * nCand * 3, nLam * N * fracs.length * nCand * 3];
    return onGrid(wasm, grid, PAIRS, ({ lead, put, putTable, out, read }) => {
        const outs = sizes.map(size => out(size));
        wasm._tmm_needle_scan_spectrum(...lead(polCode), putTable(candNK, nCand, PAIRS), nCand,
            put(fracs), fracs.length, ...outs);
        return named(['base', 'gaps', 'intra'], outs.map((ptr, i) => read(ptr, sizes[i])));
    });
}
