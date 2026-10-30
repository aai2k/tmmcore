/**
 * growing.js : a growing top layer on a WebAssembly instance: the monitoring
 * curve of one layer, the spectra after each deposited layer, and a growing
 * stack over a wavelength grid, with the forms of 0.4 that run them. The
 * JavaScript ports of the same kernels are in src/js/growing.js.
 */

import { checkLayerKeys, isNumbers } from '../js/stack.js';
import { warnOnce } from '../js/deprecated.js';
import { onGrid, named, PAIRS, writePairLayers } from './memory.js';

const GROWING_KEYS = ['Rs', 'Ts', 'Rp', 'Tp', 'Rrs', 'Rrp'];

export const growingMethods = {
    /** True: every build since 0.5.0 carries the growing-stack kernels. */
    hasGrowingKernels() { return true; },

    /** True: every build since 0.5.0 carries the persistent growing-layer evaluator. */
    hasGrowingEval() { return true; },

    /**
     * Deprecated, removed in 0.6.0: use growingStack(stack, wavelengths,
     * { angle }). The arguments of depositionSpectra in 0.4, with layerNK and
     * thick the finished stack, the layer facing the incident medium first.
     * Returns the handle of growingStackWasm, whose setTop takes the growing
     * layer's ñ per wavelength as [re, im].
     */
    growingEval(lambdas, n0List, nsList, layerNK, thick, theta_deg) {
        warnOnce('growingEval', 'growingStack');
        return growingStackWasm(this, { lambdas, n0List, nsList, layerNK, thick, theta_deg: [Number(theta_deg)] },
            PAIRS);
    },
};

/**
 * The monitoring curve at one wavelength, as monitorCurveKernel: { Rs, Ts, Rp,
 * Tp, Rrs, Rrp }, one value per sample thickness. The finished layers `base`,
 * outermost first, are written by `writeLayers` (memory.js) and their keys are
 * checked by the caller; the media and the growing layer `top` are [re, im].
 */
export function monitorCurveWasm(wasm, { lambda_nm, theta_deg, n0, ns, base, top, thicknesses }, writeLayers) {
    const NB = base.length;
    const nD = thicknesses.length;
    // arena: base[3NB] | d[nD] | Rs,Ts,Rp,Tp,Rrs,Rrp [6·nD]
    const oD = 3 * NB, oOut = oD + nD;
    const need = oOut + 6 * nD;
    const ptr = wasm._scratch(need);
    const buf = wasm._view(ptr, need);
    writeLayers(buf, base);
    for (let k = 0; k < nD; k++) buf[oD + k] = thicknesses[k];
    const P = (off) => ptr + off * 8;
    const outs = GROWING_KEYS.map((_, q) => oOut + q * nD);
    wasm._tmm_monitor_curve(lambda_nm, theta_deg, n0[0], n0[1], ns[0], ns[1], P(0), NB, top[0], top[1],
        P(oD), nD, ...outs.map(P));
    // Fresh view after the call: memory growth detaches the old buffer.
    const out = wasm._view(ptr, need);
    return named(GROWING_KEYS, outs.map(at => out.slice(at, at + nD)));
}

/**
 * The coated surface after each deposited layer, as depositionSpectraKernel:
 * { Rs, Ts, Rp, Tp, Rrs, Rrp }, each N × nLam, step-major. The layers arrive
 * in deposition order, first deposited first. `indices` says how the grid
 * holds its indices (memory.js), PAIRS or STACK.
 */
export function depositionSpectraWasm(wasm, grid, indices) {
    const size = grid.thick.length * grid.lambdas.length;
    return onGrid(wasm, grid, indices, ({ lead, out, read }) => {
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
 * stack, outermost first, is folded once and kept in kernel memory. `indices`
 * says how the grid holds its indices (memory.js), PAIRS or STACK.
 * setTop(ngList) declares the growing layer's ñ per wavelength as [re, im];
 * sample(d, out?) returns { Rs, Ts, Rp, Tp, Rrs, Rrp }, each
 * Float64Array(nLam), written into `out` when given; free() releases the
 * kernel memory.
 */
export function growingStackWasm(wasm, grid, indices) {
    const nLam = grid.lambdas.length;
    const handlePtr = onGrid(wasm, grid, indices, ({ lead }) => wasm._tmm_growing_stack_create(...lead()));
    if (!handlePtr) throw new Error('tmmWasm: growing-stack state allocation failed');

    const handle = {
        ptr: handlePtr,
        setTop(ngList) {
            if (!this.ptr) throw new Error('tmmcore: growing stack already freed');
            const p = wasm._scratch(2 * nLam);
            const buf = wasm._view(p, 2 * nLam);
            for (let i = 0; i < nLam; i++) {
                buf[2 * i] = ngList[i][0]; buf[2 * i + 1] = ngList[i][1];
            }
            wasm._tmm_growing_stack_set_top(this.ptr, p);
        },
        sample(d, out = null) {
            if (!this.ptr) throw new Error('tmmcore: growing stack already freed');
            const p = wasm._scratch(6 * nLam);
            const P = (block) => p + block * nLam * 8;
            const ok = wasm._tmm_growing_stack_sample(this.ptr, d,
                P(0), P(1), P(2), P(3), P(4), P(5));
            // The kernel writes nothing without a declared growing layer;
            // copying the arena anyway would hand back stale memory.
            if (!ok) throw new Error('tmmcore: growing stack sampled before setTop');
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

// monitorCurve with its arguments of 0.4: the media and the growing layer as
// [re, im], the finished layers as { n: [re, im], d }, outermost first.
function monitorCurveOld(wasm, [lambda_nm, theta_deg, n0, ns, baseLayers, ngNK, dArr]) {
    checkLayerKeys(baseLayers, 'n', 'd');
    return monitorCurveWasm(wasm,
        { lambda_nm, theta_deg, n0, ns, base: baseLayers, top: ngNK, thicknesses: dArr }, writePairLayers);
}

// depositionSpectra with its arguments of 0.4: the lists of the grid one by
// one, each ñ as [re, im], the layers in deposition order, and one angle.
function depositionSpectraOld(wasm, [lambdas, n0List, nsList, layerNK, thick, theta_deg]) {
    const grid = { lambdas, n0List, nsList, layerNK, thick, theta_deg: [Number(theta_deg)] };
    return depositionSpectraWasm(wasm, grid, PAIRS);
}

/**
 * monitorCurve and depositionSpectra of an instance, built on its calls of
 * 0.5, `calls`. Through 0.5.x each also takes its arguments of 0.4, told apart
 * by the first: a wavelength or a list of them in 0.4, and anything else read
 * as the stack of 0.5, which refuses what is not one. The 0.4 form prints one
 * warning on its first use.
 */
export function growingCalls(wasm, calls) {
    return {
        monitorCurve(stack, wavelength, options = {}) {
            if (typeof stack !== 'number') return calls.monitorCurve(stack, wavelength, options);
            warnOnce('monitorCurve(lambda_nm, theta_deg, n0, ns, baseLayers, ngNK, dArr)',
                'monitorCurve(stack, wavelength, { angle, top, thicknesses })');
            return monitorCurveOld(wasm, arguments);
        },
        depositionSpectra(stack, wavelengths, options = {}) {
            if (!isNumbers(stack, stack?.length)) return calls.depositionSpectra(stack, wavelengths, options);
            warnOnce('depositionSpectra(lambdas, n0List, nsList, layerNK, thick, theta_deg)',
                'depositionSpectra(stack, wavelengths, { angle })');
            return depositionSpectraOld(wasm, arguments);
        },
    };
}
