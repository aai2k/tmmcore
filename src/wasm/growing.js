/**
 * growing.js : a growing top layer on a WebAssembly instance: the monitoring
 * curve of one layer, the spectra after each deposited layer, and a growing
 * stack over a wavelength grid. The JavaScript ports of the same kernels are
 * in src/js/growing.js.
 */

import { checkLayerKeys } from '../js/stack.js';
import { onGrid, named, PAIRS } from './memory.js';

const GROWING_KEYS = ['Rs', 'Ts', 'Rp', 'Tp', 'Rrs', 'Rrp'];

export const growingMethods = {
    /** True: every build since 0.5.0 carries the growing-stack kernels. */
    hasGrowingKernels() { return true; },

    /** True: every build since 0.5.0 carries the persistent growing-layer evaluator. */
    hasGrowingEval() { return true; },

    /**
     * Monitor curve of one growing layer at one wavelength: the completed
     * stack's matrix is built once, then every sample thickness costs one 2×2
     * multiply. Returns the forward and substrate-side passes of the coated
     * surface, both polarizations; the caller does the incoherent slab
     * combination (or forms A = 1−R−T for a semi-infinite substrate).
     * @param {[number,number]} n0   incident ñ
     * @param {[number,number]} ns   substrate ñ
     * @param {{n:[number,number],d:number}[]} baseLayers  completed stack,
     *        outermost first
     * @param {[number,number]} ngNK growing layer ñ
     * @param {ArrayLike<number>} dArr  sample thicknesses (nm)
     * @returns {{Rs,Ts,Rp,Tp,Rrs,Rrp}} each a Float64Array(dArr.length)
     */
    monitorCurve(lambda_nm, theta_deg, n0, ns, baseLayers, ngNK, dArr) {
        checkLayerKeys(baseLayers, 'n', 'd');
        const NB = baseLayers.length;
        const nD = dArr.length;
        // arena: base[3NB] | d[nD] | Rs,Ts,Rp,Tp,Rrs,Rrp [6·nD]
        const oBase = 0, oD = 3 * NB, oOut = oD + nD;
        const need = oOut + 6 * nD;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < NB; i++) {
            buf[3 * i + 0] = baseLayers[i].n[0];
            buf[3 * i + 1] = baseLayers[i].n[1];
            buf[3 * i + 2] = baseLayers[i].d;
        }
        for (let k = 0; k < nD; k++) buf[oD + k] = dArr[k];
        const P = (off) => ptr + off * 8;
        this._tmm_monitor_curve(lambda_nm, theta_deg,
            n0[0], n0[1], ns[0], ns[1], P(oBase), NB, ngNK[0], ngNK[1],
            P(oD), nD,
            P(oOut), P(oOut + nD), P(oOut + 2 * nD), P(oOut + 3 * nD),
            P(oOut + 4 * nD), P(oOut + 5 * nD));
        // Fresh view after the call: memory growth detaches the old buffer.
        const out = this._view(ptr, need);
        return {
            Rs:  out.slice(oOut, oOut + nD),
            Ts:  out.slice(oOut + nD, oOut + 2 * nD),
            Rp:  out.slice(oOut + 2 * nD, oOut + 3 * nD),
            Tp:  out.slice(oOut + 3 * nD, oOut + 4 * nD),
            Rrs: out.slice(oOut + 4 * nD, oOut + 5 * nD),
            Rrp: out.slice(oOut + 5 * nD, oOut + 6 * nD),
        };
    },

    /**
     * Per-step spectra of a growing stack: one call returns the forward and
     * substrate-side passes of the coated surface after every deposited layer,
     * so all N step spectra cost about what the final one costs alone. Layers
     * in DEPOSITION order (first deposited first).
     * @param {number[]} lambdas
     * @param {[number,number][]} n0List  incident ñ per λ
     * @param {[number,number][]} nsList  substrate ñ per λ
     * @param {[number,number][][]} layerNK  [layer][λ] = ñ, deposition order
     * @param {number[]} thick   thicknesses (nm); 0 repeats the previous step
     * @param {number}   theta_deg
     * @returns {{Rs,Ts,Rp,Tp,Rrs,Rrp}} each a Float64Array(N · nLam),
     *          step-major: value for step k at wavelength i is at [k·nLam + i]
     */
    depositionSpectra(lambdas, n0List, nsList, layerNK, thick, theta_deg) {
        return depositionSpectraWasm(this, { lambdas, n0List, nsList, layerNK, thick, theta_deg: [Number(theta_deg)] });
    },

    /**
     * Persistent growing-layer evaluator: the completed stack's products for
     * the whole wavelength grid are folded once and kept in kernel memory;
     * each sample() then answers one thickness of the growing layer across
     * the grid. This is the shape a broadband monitor scan needs (a spectrum
     * per scan while the layers beneath stay fixed), where monitorCurve is
     * the single-λ, many-thicknesses shape.
     *
     * The completed stack arrives OUTERMOST FIRST, layer-major like
     * depositionSpectra's layerNK; zero-thickness entries are skipped.
     *
     * The returned handle owns kernel memory: call free() when the layer is
     * done. A dropped handle is reclaimed by a finalizer eventually, but
     * deterministic free() is what keeps a long run's footprint flat.
     *
     * @param {number[]} lambdas
     * @param {[number,number][]} n0List  incident ñ per λ
     * @param {[number,number][]} nsList  substrate ñ per λ
     * @param {[number,number][][]} layerNK  [layer][λ] = ñ, outermost first
     * @param {number[]} thick   completed thicknesses (nm)
     * @param {number}   theta_deg
     * @returns {{setTop, sample, free}} setTop(ngList) declares the growing
     *          layer's ñ per λ; sample(d, out?) returns {Rs,Ts,Rp,Tp,Rrs,Rrp}
     *          (each Float64Array(nLam), written into `out` when given, so a
     *          scan loop can reuse one set of buffers).
     */
    growingEval(lambdas, n0List, nsList, layerNK, thick, theta_deg) {
        return growingStackWasm(this, { lambdas, n0List, nsList, layerNK, thick, theta_deg: [Number(theta_deg)] });
    },
};

/**
 * The coated surface after each deposited layer, as depositionSpectraKernel:
 * { Rs, Ts, Rp, Tp, Rrs, Rrp }, each N × nLam, step-major. The layers arrive
 * in deposition order, first deposited first.
 */
export function depositionSpectraWasm(wasm, grid) {
    const size = grid.thick.length * grid.lambdas.length;
    return onGrid(wasm, grid, PAIRS, ({ lead, out, read }) => {
        const outs = GROWING_KEYS.map(() => out(size));
        if (!wasm._tmm_deposition_spectra(...lead(), ...outs)) {
            // The outputs were never written; returning them would hand the
            // caller uninitialized memory as spectra.
            throw new Error('tmmWasm: deposition-spectra state allocation failed');
        }
        return named(GROWING_KEYS, outs.map(ptr => read(ptr, size)));
    });
}

// Backstop for growing-stack handles that are dropped without free():
// reclaims their kernel memory when the JS handle is collected. Deterministic
// free() remains the contract; this only keeps a leak from being permanent.
const growingStackFinalizer = typeof FinalizationRegistry !== 'undefined'
    ? new FinalizationRegistry(({ wasm, ptr }) => {
        try { wasm._tmm_growing_stack_free(ptr); } catch (_) { /* instance gone */ }
    })
    : null;

/**
 * A growing top layer over the grid, as growingStackKernel: the finished
 * stack, outermost first, is folded once and kept in kernel memory.
 * setTop(ngList) declares the growing layer's ñ per wavelength; sample(d, out?)
 * returns { Rs, Ts, Rp, Tp, Rrs, Rrp }, each Float64Array(nLam), written into
 * `out` when given; free() releases the kernel memory.
 */
export function growingStackWasm(wasm, grid) {
    const nLam = grid.lambdas.length;
    const handlePtr = onGrid(wasm, grid, PAIRS, ({ lead }) => wasm._tmm_growing_stack_create(...lead()));
    if (!handlePtr) throw new Error('tmmWasm: growing evaluator allocation failed');

    const handle = {
        ptr: handlePtr,
        setTop(ngList) {
            if (!this.ptr) throw new Error('tmmWasm: growing evaluator already freed');
            const p = wasm._scratch(2 * nLam);
            const buf = wasm._view(p, 2 * nLam);
            for (let i = 0; i < nLam; i++) {
                buf[2 * i] = ngList[i][0]; buf[2 * i + 1] = ngList[i][1];
            }
            wasm._tmm_growing_stack_set_top(this.ptr, p);
        },
        sample(d, out = null) {
            if (!this.ptr) throw new Error('tmmWasm: growing evaluator already freed');
            const p = wasm._scratch(6 * nLam);
            const P = (block) => p + block * nLam * 8;
            const ok = wasm._tmm_growing_stack_sample(this.ptr, d,
                P(0), P(1), P(2), P(3), P(4), P(5));
            // The kernel writes nothing without a declared growing layer;
            // copying the arena anyway would hand back stale memory.
            if (!ok) throw new Error('tmmWasm: growing evaluator sampled before setTop');
            // Fresh view AFTER the call (memory may have grown → detach).
            const buf = wasm._view(p, 6 * nLam);
            const result = out || named(GROWING_KEYS, GROWING_KEYS.map(() => new Float64Array(nLam)));
            GROWING_KEYS.forEach((key, q) => result[key].set(buf.subarray(q * nLam, (q + 1) * nLam)));
            return result;
        },
        free() {
            if (!this.ptr) return;
            growingStackFinalizer?.unregister(this);
            wasm._tmm_growing_stack_free(this.ptr);
            this.ptr = 0;
        },
    };
    growingStackFinalizer?.register(handle, { wasm, ptr: handlePtr }, handle);
    return handle;
}
