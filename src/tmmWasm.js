/**
 * tmmWasm.js : loader and ergonomic wrappers for the WebAssembly TMM kernel.
 *
 * The kernel (`tmm_kernel.c`, built to `tmm_kernel.wasm`) is a line-by-line port
 * of the JavaScript TMM in `tmm.js`. This module instantiates it : in a browser
 * main thread, a Web Worker, or Node : and exposes wrappers whose signatures
 * mirror the JS functions.
 *
 * Acceleration is opt-in and falls back to JavaScript: if the `.wasm` is
 * unavailable, instantiation fails, or the feature flag is off, the wrappers
 * return `null` and callers use the JS path. Results are identical either way
 * to float64 round-off.
 *
 * Instances are not shared across threads : there is no shared memory, so each
 * context instantiates its own from the same bytes. Use `instantiateTmmWasm()`
 * where you already hold the bytes (a worker receives them in its init message)
 * and `initTmmWasmFromUrl()` where the artifact is fetchable.
 */

import { omegaFromLambdaNm } from './phase.js';

let _instance = null;     // TmmWasmInstance | null
let _enabled = false;     // feature flag (default OFF)
let _initPromise = null;  // de-dupe concurrent init

// ── Jet marshalling ──────────────────────────────────────────────────────────
// A jet crosses the boundary as 8 doubles, [re, im] per order.

function writeJet(buf, offset, jet) {
    for (let i = 0; i < 4; i++) {
        buf[offset + 2 * i] = jet[i][0];
        buf[offset + 2 * i + 1] = jet[i][1];
    }
}

// The kernel writes NaN where the coefficient is exactly zero and the phase is
// undefined; the JS reference returns null there, so agree with it.
function readPhase(buf, offset) {
    const magnitudeSquared = buf[offset + 4];
    if (Number.isNaN(magnitudeSquared)) return null;
    const phaseRad = buf[offset];
    return {
        phaseRad,
        phaseDeg: phaseRad * 180 / Math.PI,
        gd: buf[offset + 1],
        gdd: buf[offset + 2],
        tod: buf[offset + 3],
        magnitudeSquared,
    };
}

// Permissive imports: a STANDALONE_WASM build of pure-math C usually needs no
// imports, but ALLOW_MEMORY_GROWTH may emit `emscripten_notify_memory_growth`,
// and some toolchains emit WASI stubs. Cover them so instantiation never throws.
function wasmImports() {
    return {
        env: { emscripten_notify_memory_growth: () => {} },
        wasi_snapshot_preview1: new Proxy({}, { get: () => () => 0 }),
    };
}

// Backstop for growing-evaluator handles that are dropped without free():
// reclaims their kernel memory when the JS handle is collected. Deterministic
// free() remains the contract; this only keeps a leak from being permanent.
const growingEvalFinalizer = typeof FinalizationRegistry !== 'undefined'
    ? new FinalizationRegistry(({ wasm, ptr }) => {
        try { wasm._growing_eval_free(ptr); } catch (_) { /* instance gone */ }
    })
    : null;

export class TmmWasmInstance {
    constructor(instance) {
        const ex = instance.exports;
        this.exports = ex;
        // STANDALONE_WASM reactor modules expose an initializer that must run
        // before malloc (sets up the allocator + any static ctors). Call it once.
        if (typeof ex._initialize === 'function') ex._initialize();
        else if (typeof ex.__wasm_call_ctors === 'function') ex.__wasm_call_ctors();
        this.memory = ex.memory;
        this.malloc = ex.malloc || ex._malloc;
        this.free = ex.free || ex._free;
        this._tmm_one = ex.tmm_one || ex._tmm_one;
        this._tmm_spectrum = ex.tmm_spectrum || ex._tmm_spectrum;
        this._tmm_jacobian = ex.tmm_jacobian || ex._tmm_jacobian;
        this._tmm_needle_scan = ex.tmm_needle_scan || ex._tmm_needle_scan;
        // Optional (added later for SQP/Newton accel): a .wasm built before the
        // Hessian kernel existed simply lacks it → callers fall back to JS.
        this._tmm_hessian = ex.tmm_hessian || ex._tmm_hessian || null;
        // Optional, same reason: the phase-dispersion kernel arrived after the
        // spectral one, so an older artifact lacks these three.
        this._tmm_phase_one = ex.tmm_phase_one || ex._tmm_phase_one || null;
        this._tmm_phase_spectrum = ex.tmm_phase_spectrum || ex._tmm_phase_spectrum || null;
        this._tmm_phase_jacobian = ex.tmm_phase_jacobian || ex._tmm_phase_jacobian || null;
        // Optional, same reason: the batched phase Jacobian arrived after the
        // three above, for fitting a whole measured Ψ/Δ or group-delay spectrum.
        this._tmm_phase_jacobian_spectrum = ex.tmm_phase_jacobian_spectrum
            || ex._tmm_phase_jacobian_spectrum || null;
        // Optional, same reason: the growing-stack kernels (monitor curve and
        // per-step deposition spectra) arrived after all of the above.
        this._tmm_monitor_curve = ex.tmm_monitor_curve || ex._tmm_monitor_curve || null;
        this._tmm_deposition_spectra = ex.tmm_deposition_spectra || ex._tmm_deposition_spectra || null;
        // Optional, same reason: the persistent growing-layer evaluator
        // (wavelength-grid-per-call) arrived after the two kernels above.
        this._growing_eval_create = ex.tmm_growing_eval_create || ex._tmm_growing_eval_create || null;
        this._growing_eval_set_top = ex.tmm_growing_eval_set_top || ex._tmm_growing_eval_set_top || null;
        this._growing_eval_sample = ex.tmm_growing_eval_sample || ex._tmm_growing_eval_sample || null;
        this._growing_eval_free = ex.tmm_growing_eval_free || ex._tmm_growing_eval_free || null;
        const missingExports = !this.malloc || !this.free || !this._tmm_one ||
            !this._tmm_spectrum || !this._tmm_jacobian || !this._tmm_needle_scan;
        if (missingExports) {
            throw new Error('tmmWasm: required exports missing from module');
        }
        this._scratchPtr = 0;   // persistent per-call scratch arena (lazy)
        this._scratchN = 0;
    }

    _alloc(nDoubles) {
        const ptr = this.malloc(nDoubles * 8);
        if (!ptr) throw new Error('tmmWasm: malloc failed');
        return ptr;
    }
    // Fresh view : memory.buffer is detached after any growth, so re-create
    // views AFTER all mallocs for a call are done.
    _view(ptr, nDoubles) {
        return new Float64Array(this.memory.buffer, ptr, nDoubles);
    }

    // Persistent scratch arena for the per-call hot paths (tmmOne/tmmJacobian).
    // These are invoked thousands of times per optimization run; malloc/free +
    // typed-array churn per call would dominate the JS↔WASM boundary cost and
    // can make per-call WASM SLOWER than JS. Reusing one buffer (grown on demand,
    // never freed between calls) makes each call just write-args / read-result.
    _scratch(nDoubles) {
        if (this._scratchN < nDoubles) {
            if (this._scratchPtr) this.free(this._scratchPtr);
            this._scratchPtr = this._alloc(nDoubles);
            this._scratchN = nDoubles;
        }
        return this._scratchPtr;
    }

    /**
     * Single (λ, θ, pol) : mirrors tmm() in thinFilmMath.js.
     * @returns {{R:number,T:number,A:number}}
     */
    tmmOne(lambda_nm, theta_deg, polCode /* 0=s,1=p */, n0, ns, layers) {
        const N = layers.length;
        const need = 3 * N + 3;                  // layers [0..3N) + out [3N..3N+3)
        const ptr = this._scratch(need);         // may grow→detach; view created AFTER
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            buf[3 * i + 0] = layers[i].n[0];
            buf[3 * i + 1] = layers[i].n[1];
            buf[3 * i + 2] = layers[i].d;
        }
        const outPtr = ptr + 3 * N * 8;
        this._tmm_one(lambda_nm, theta_deg, polCode | 0,
            n0[0], n0[1], ns[0], ns[1], ptr, N, outPtr);
        return { R: buf[3 * N], T: buf[3 * N + 1], A: buf[3 * N + 2] };
    }

    /**
     * Batched spectrum over a λ grid for BOTH polarizations : the boundary-
     * amortizing path behind evaluateSpectrum().
     * @param {number[]} lambdas
     * @param {[number,number][]} n0List  incident ñ per λ
     * @param {[number,number][]} nsList  substrate ñ per λ
     * @param {[number,number][][]} layerNK  [layer][λ] = ñ
     * @param {number[]} thick   layer thicknesses (nm), length N
     * @param {number}   theta_deg
     * @returns {{Rs,Ts,As,Rp,Tp,Ap}} each a Float64Array(nLam)
     */
    tmmSpectrum(lambdas, n0List, nsList, layerNK, thick, theta_deg) {
        const nLam = lambdas.length;
        const N = thick.length;

        const lamPtr = this._alloc(nLam);
        const n0Ptr  = this._alloc(2 * nLam);
        const nsPtr  = this._alloc(2 * nLam);
        const mPtr   = this._alloc(Math.max(1, 2 * N * nLam));
        const thPtr  = this._alloc(Math.max(1, N));
        const rsPtr = this._alloc(nLam), tsPtr = this._alloc(nLam), asPtr = this._alloc(nLam);
        const rpPtr = this._alloc(nLam), tpPtr = this._alloc(nLam), apPtr = this._alloc(nLam);

        // Views created after all mallocs (buffer may have grown/detached).
        const lam = this._view(lamPtr, nLam);
        const n0v = this._view(n0Ptr, 2 * nLam);
        const nsv = this._view(nsPtr, 2 * nLam);
        const mv  = this._view(mPtr, Math.max(1, 2 * N * nLam));
        const thv = this._view(thPtr, Math.max(1, N));
        for (let i = 0; i < nLam; i++) {
            lam[i] = lambdas[i];
            n0v[2 * i] = n0List[i][0]; n0v[2 * i + 1] = n0List[i][1];
            nsv[2 * i] = nsList[i][0]; nsv[2 * i + 1] = nsList[i][1];
        }
        for (let k = 0; k < N; k++) {
            thv[k] = thick[k];
            const row = layerNK[k];
            const base = k * nLam * 2;
            for (let i = 0; i < nLam; i++) {
                mv[base + 2 * i]     = row[i][0];
                mv[base + 2 * i + 1] = row[i][1];
            }
        }

        this._tmm_spectrum(lamPtr, nLam, n0Ptr, nsPtr, mPtr, thPtr, N, theta_deg,
            rsPtr, tsPtr, asPtr, rpPtr, tpPtr, apPtr);

        // Copy outputs out of wasm memory before freeing.
        const cp = (p) => Float64Array.from(this._view(p, nLam));
        const res = { Rs: cp(rsPtr), Ts: cp(tsPtr), As: cp(asPtr),
                      Rp: cp(rpPtr), Tp: cp(tpPtr), Ap: cp(apPtr) };
        for (const p of [lamPtr, n0Ptr, nsPtr, mPtr, thPtr,
                         rsPtr, tsPtr, asPtr, rpPtr, tpPtr, apPtr]) this.free(p);
        return res;
    }

    /** True if the loaded module carries the growing-stack kernels. */
    hasGrowingKernels() {
        return !!(this._tmm_monitor_curve && this._tmm_deposition_spectra);
    }

    /** True if the loaded module carries the persistent growing-layer evaluator. */
    hasGrowingEval() {
        return !!(this._growing_eval_create && this._growing_eval_set_top
            && this._growing_eval_sample && this._growing_eval_free);
    }

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
    }

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
        const nLam = lambdas.length;
        const N = thick.length;
        const nOut = Math.max(1, N * nLam);

        const lamPtr = this._alloc(nLam);
        const n0Ptr  = this._alloc(2 * nLam);
        const nsPtr  = this._alloc(2 * nLam);
        const mPtr   = this._alloc(Math.max(1, 2 * N * nLam));
        const thPtr  = this._alloc(Math.max(1, N));
        const outPtrs = Array.from({ length: 6 }, () => this._alloc(nOut));

        const lam = this._view(lamPtr, nLam);
        const n0v = this._view(n0Ptr, 2 * nLam);
        const nsv = this._view(nsPtr, 2 * nLam);
        const mv  = this._view(mPtr, Math.max(1, 2 * N * nLam));
        const thv = this._view(thPtr, Math.max(1, N));
        for (let i = 0; i < nLam; i++) {
            lam[i] = lambdas[i];
            n0v[2 * i] = n0List[i][0]; n0v[2 * i + 1] = n0List[i][1];
            nsv[2 * i] = nsList[i][0]; nsv[2 * i + 1] = nsList[i][1];
        }
        for (let k = 0; k < N; k++) {
            thv[k] = thick[k];
            const row = layerNK[k];
            const base = k * nLam * 2;
            for (let i = 0; i < nLam; i++) {
                mv[base + 2 * i]     = row[i][0];
                mv[base + 2 * i + 1] = row[i][1];
            }
        }

        const ok = this._tmm_deposition_spectra(lamPtr, nLam, n0Ptr, nsPtr, mPtr, thPtr, N,
            theta_deg, ...outPtrs);
        if (!ok) {
            for (const p of [lamPtr, n0Ptr, nsPtr, mPtr, thPtr, ...outPtrs]) this.free(p);
            // The outputs were never written; returning them would hand the
            // caller uninitialized memory as spectra.
            throw new Error('tmmWasm: deposition-spectra state allocation failed');
        }

        const cp = (p) => Float64Array.from(this._view(p, nOut));
        const [Rs, Ts, Rp, Tp, Rrs, Rrp] = outPtrs.map(cp);
        for (const p of [lamPtr, n0Ptr, nsPtr, mPtr, thPtr, ...outPtrs]) this.free(p);
        return { Rs, Ts, Rp, Tp, Rrs, Rrp };
    }

    /**
     * Persistent growing-layer evaluator: the completed stack's products for
     * the whole wavelength grid are folded once and kept in kernel memory;
     * each sample() then answers one thickness of the growing layer across
     * the grid. This is the shape a broadband monitor scan needs (a spectrum
     * per scan while the layers beneath stay fixed), where monitorCurve is
     * the single-λ, many-thicknesses shape.
     *
     * The completed stack arrives OUTERMOST FIRST, layer-major like
     * depositionSpectra's matNK; zero-thickness entries are skipped.
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
        if (!this.hasGrowingEval()) throw new Error('tmmWasm: growing evaluator kernel not in this build');
        const nLam = lambdas.length;
        const NB = thick.length;

        const lamPtr = this._alloc(nLam);
        const n0Ptr  = this._alloc(2 * nLam);
        const nsPtr  = this._alloc(2 * nLam);
        const mPtr   = this._alloc(Math.max(1, 2 * NB * nLam));
        const thPtr  = this._alloc(Math.max(1, NB));
        const lam = this._view(lamPtr, nLam);
        const n0v = this._view(n0Ptr, 2 * nLam);
        const nsv = this._view(nsPtr, 2 * nLam);
        const mv  = this._view(mPtr, Math.max(1, 2 * NB * nLam));
        const thv = this._view(thPtr, Math.max(1, NB));
        for (let i = 0; i < nLam; i++) {
            lam[i] = lambdas[i];
            n0v[2 * i] = n0List[i][0]; n0v[2 * i + 1] = n0List[i][1];
            nsv[2 * i] = nsList[i][0]; nsv[2 * i + 1] = nsList[i][1];
        }
        for (let k = 0; k < NB; k++) {
            thv[k] = thick[k];
            const row = layerNK[k];
            const base = k * nLam * 2;
            for (let i = 0; i < nLam; i++) {
                mv[base + 2 * i]     = row[i][0];
                mv[base + 2 * i + 1] = row[i][1];
            }
        }
        const handlePtr = this._growing_eval_create(lamPtr, nLam, theta_deg,
            n0Ptr, nsPtr, mPtr, thPtr, NB);
        for (const p of [lamPtr, n0Ptr, nsPtr, mPtr, thPtr]) this.free(p);
        if (!handlePtr) throw new Error('tmmWasm: growing evaluator allocation failed');

        const self = this;
        const handle = {
            ptr: handlePtr,
            setTop(ngList) {
                if (!this.ptr) throw new Error('tmmWasm: growing evaluator already freed');
                const p = self._scratch(2 * nLam);
                const buf = self._view(p, 2 * nLam);
                for (let i = 0; i < nLam; i++) {
                    buf[2 * i] = ngList[i][0]; buf[2 * i + 1] = ngList[i][1];
                }
                self._growing_eval_set_top(this.ptr, p);
            },
            sample(d, out = null) {
                if (!this.ptr) throw new Error('tmmWasm: growing evaluator already freed');
                const p = self._scratch(6 * nLam);
                const P = (block) => p + block * nLam * 8;
                const ok = self._growing_eval_sample(this.ptr, d,
                    P(0), P(1), P(2), P(3), P(4), P(5));
                // The kernel writes nothing without a declared growing layer;
                // copying the arena anyway would hand back stale memory.
                if (!ok) throw new Error('tmmWasm: growing evaluator sampled before setTop');
                // Fresh view AFTER the call (memory may have grown → detach).
                const buf = self._view(p, 6 * nLam);
                if (!out) {
                    out = {
                        Rs: new Float64Array(nLam), Ts: new Float64Array(nLam),
                        Rp: new Float64Array(nLam), Tp: new Float64Array(nLam),
                        Rrs: new Float64Array(nLam), Rrp: new Float64Array(nLam),
                    };
                }
                out.Rs.set(buf.subarray(0, nLam));
                out.Ts.set(buf.subarray(nLam, 2 * nLam));
                out.Rp.set(buf.subarray(2 * nLam, 3 * nLam));
                out.Tp.set(buf.subarray(3 * nLam, 4 * nLam));
                out.Rrs.set(buf.subarray(4 * nLam, 5 * nLam));
                out.Rrp.set(buf.subarray(5 * nLam, 6 * nLam));
                return out;
            },
            free() {
                if (!this.ptr) return;
                growingEvalFinalizer?.unregister(this);
                self._growing_eval_free(this.ptr);
                this.ptr = 0;
            },
        };
        growingEvalFinalizer?.register(handle,
            { wasm: this, ptr: handlePtr }, handle);
        return handle;
    }

    /**
     * Analytic thickness Jacobian for one (λ, θ, pol) : mirrors
     * tmmThicknessJacobian(). layers used AS-IS (index parity).
     * @returns {{R,T,A, dRdd:Float64Array, dTdd, dAdd, N}}
     */
    tmmJacobian(lambda_nm, theta_deg, polCode, n0, ns, layers) {
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
        // Re-create the view AFTER the kernel call (like tmmSpectrum): under
        // ALLOW_MEMORY_GROWTH the kernel may grow wasm memory, which detaches the
        // ArrayBuffer `buf` was created over : reading the stale `buf` then yields
        // garbage / throws. `out` is a fresh view over the current buffer.
        const out = this._view(ptr, need);
        return {
            R: out[oBase], T: out[oBase + 1], A: out[oBase + 2], N,
            dRdd: out.slice(oDR, oDR + N),
            dTdd: out.slice(oDT, oDT + N),
            dAdd: out.slice(oDA, oDA + N),
        };
    }

    /** True if the loaded module carries the Hessian kernel (newer build). */
    hasHessian() { return !!this._tmm_hessian; }

    /**
     * Analytic thickness Hessian for one (λ, θ, pol) : mirrors
     * tmmThicknessHessian(). Returns first AND second derivatives; the N×N
     * second-derivative blocks are reshaped into nested arrays (one Float64Array
     * row per layer, FULL symmetric) so the shape matches the JS oracle exactly.
     * @returns {{R,T,A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N}}
     */
    tmmHessian(lambda_nm, theta_deg, polCode, n0, ns, layers) {
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
    }

    /** True if the loaded module carries the phase-dispersion kernel. */
    hasPhase() { return !!this._tmm_phase_one; }

    /**
     * Phase, group delay, GDD and TOD at one wavelength : mirrors
     * tmmPhaseDispersion() in phase.js.
     *
     * @param {number[][]} n0Jet  incident-medium index jet, 4 × [re, im]
     * @param {number[][]} nsJet  substrate index jet
     * @param {{nJet:number[][], d:number}[]} layers
     * @param {{omega?:number, sinTheta0Jet?:number[][]}} [options]
     * @returns {{r: object|null, t: object|null}}
     */
    tmmPhaseOne(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options = {}) {
        const N = layers.length;
        const omega = options.omega ?? omegaFromLambdaNm(lambda_nm);
        const sinJet = options.sinTheta0Jet || null;
        // arena: layerJets[8N] | thick[N] | n0[8] | ns[8] | sin[8] | out[10]
        const oLay = 0, oThick = 8 * N, oN0 = oThick + N, oNs = oN0 + 8,
              oSin = oNs + 8, oOut = oSin + 8;
        const need = oOut + 10;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            writeJet(buf, oLay + 8 * i, layers[i].nJet);
            buf[oThick + i] = layers[i].d;
        }
        writeJet(buf, oN0, n0Jet);
        writeJet(buf, oNs, nsJet);
        if (sinJet) writeJet(buf, oSin, sinJet);
        const P = (off) => ptr + off * 8;
        this._tmm_phase_one(lambda_nm, omega, theta_deg, polCode | 0,
            P(oN0), P(oNs), P(oLay), P(oThick), N, sinJet ? P(oSin) : 0, P(oOut));
        const out = this._view(ptr, need);
        return { r: readPhase(out, oOut), t: readPhase(out, oOut + 5) };
    }

    /**
     * Batched phase dispersion over a λ grid : the boundary-amortizing path.
     *
     * Unlike `tmmSpectrum` this takes one polarization, because the kernel is an
     * order of magnitude dearer per sample and callers at normal incidence would
     * otherwise pay twice for the same numbers.
     *
     * @param {number[]} lambdas
     * @param {number[][][]} n0Jets  index jet per λ
     * @param {number[][][]} nsJets  index jet per λ
     * @param {number[][][][]} layerJets  [layer][λ] = index jet
     * @param {number[]} thick  layer thicknesses (nm), length N
     * @param {{omegas?:number[], sinJets?:number[][][]}} [options]
     * @returns {{r: object, t: object}} each `{phaseRad, gd, gdd, tod,
     *   magnitudeSquared}` of Float64Array(nLam). Failed samples hold NaN.
     */
    tmmPhaseSpectrum(lambdas, n0Jets, nsJets, layerJets, thick, theta_deg, polCode, options = {}) {
        const nLam = lambdas.length;
        const N = thick.length;
        const omegas = options.omegas
            || lambdas.map(lambda => omegaFromLambdaNm(lambda));
        const sinJets = options.sinJets || null;

        const lamPtr = this._alloc(nLam);
        const omPtr = this._alloc(nLam);
        const n0Ptr = this._alloc(8 * nLam);
        const nsPtr = this._alloc(8 * nLam);
        const matPtr = this._alloc(Math.max(1, 8 * N * nLam));
        const thPtr = this._alloc(Math.max(1, N));
        const sinPtr = sinJets ? this._alloc(8 * nLam) : 0;
        const outPtr = this._alloc(10 * nLam);

        // Views created after all mallocs (the buffer may have grown/detached).
        const lam = this._view(lamPtr, nLam);
        const om = this._view(omPtr, nLam);
        const n0v = this._view(n0Ptr, 8 * nLam);
        const nsv = this._view(nsPtr, 8 * nLam);
        const matv = this._view(matPtr, Math.max(1, 8 * N * nLam));
        const thv = this._view(thPtr, Math.max(1, N));
        const sinv = sinJets ? this._view(sinPtr, 8 * nLam) : null;
        for (let i = 0; i < nLam; i++) {
            lam[i] = lambdas[i];
            om[i] = omegas[i];
            writeJet(n0v, 8 * i, n0Jets[i]);
            writeJet(nsv, 8 * i, nsJets[i]);
            if (sinv) writeJet(sinv, 8 * i, sinJets[i]);
        }
        for (let k = 0; k < N; k++) {
            thv[k] = thick[k];
            const row = layerJets[k];
            const base = k * nLam * 8;
            for (let i = 0; i < nLam; i++) writeJet(matv, base + 8 * i, row[i]);
        }

        this._tmm_phase_spectrum(lamPtr, omPtr, nLam, n0Ptr, nsPtr, matPtr, thPtr, N,
            theta_deg, polCode | 0, sinPtr, outPtr);

        // De-interleave into one array per quantity before freeing.
        const out = this._view(outPtr, 10 * nLam);
        const side = (base) => {
            const q = {
                phaseRad: new Float64Array(nLam), gd: new Float64Array(nLam),
                gdd: new Float64Array(nLam), tod: new Float64Array(nLam),
                magnitudeSquared: new Float64Array(nLam),
            };
            const keys = ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared'];
            for (let i = 0; i < nLam; i++) {
                for (let j = 0; j < 5; j++) q[keys[j]][i] = out[10 * i + base + j];
            }
            return q;
        };
        const result = { r: side(0), t: side(5) };
        for (const p of [lamPtr, omPtr, n0Ptr, nsPtr, matPtr, thPtr, outPtr]) this.free(p);
        if (sinPtr) this.free(sinPtr);
        return result;
    }

    /**
     * Phase dispersion plus exact thickness derivatives : mirrors
     * tmmPhaseThicknessJacobian(). Layers used AS-IS (index parity).
     *
     * @returns {{r, t}} each the phase quantities plus `dPhaseDeg`, `dGd`,
     *   `dGdd`, `dTod` and `dLogMagnitudeSquared` as Float64Array(N).
     */
    tmmPhaseJacobian(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options = {}) {
        const N = layers.length;
        const M = Math.max(1, N);
        const omega = options.omega ?? omegaFromLambdaNm(lambda_nm);
        const sinJet = options.sinTheta0Jet || null;
        // arena: layerJets[8N] | thick[N] | n0[8] | ns[8] | sin[8] | out[10] | deriv[10M]
        const oLay = 0, oThick = 8 * N, oN0 = oThick + N, oNs = oN0 + 8,
              oSin = oNs + 8, oOut = oSin + 8, oDeriv = oOut + 10;
        const need = oDeriv + 10 * M;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            writeJet(buf, oLay + 8 * i, layers[i].nJet);
            buf[oThick + i] = layers[i].d;
        }
        writeJet(buf, oN0, n0Jet);
        writeJet(buf, oNs, nsJet);
        if (sinJet) writeJet(buf, oSin, sinJet);
        const P = (off) => ptr + off * 8;
        this._tmm_phase_jacobian(lambda_nm, omega, theta_deg, polCode | 0,
            P(oN0), P(oNs), P(oLay), P(oThick), N, sinJet ? P(oSin) : 0,
            P(oOut), P(oDeriv));
        const out = this._view(ptr, need);
        const side = (phaseBase, derivBase) => {
            const base = readPhase(out, phaseBase);
            if (!base) return null;
            const take = (q) => out.slice(oDeriv + (derivBase + q) * N, oDeriv + (derivBase + q) * N + N);
            return {
                ...base,
                dPhaseDeg: take(0), dGd: take(1), dGdd: take(2), dTod: take(3),
                dLogMagnitudeSquared: take(4),
            };
        };
        return { r: side(oOut, 0), t: side(oOut + 5, 5) };
    }

    /** True if the loaded module carries the batched phase Jacobian. */
    hasPhaseJacobianSpectrum() { return !!this._tmm_phase_jacobian_spectrum; }

    /**
     * Batched phase Jacobian over a λ grid : tmmPhaseJacobian at every
     * wavelength in one call, for fitting a whole measured spectrum of phase
     * quantities. Arguments as tmmPhaseSpectrum.
     *
     * @returns {{r: object, t: object}} each the five Float64Array(nLam) of
     *   tmmPhaseSpectrum plus `dPhaseDeg`, `dGd`, `dGdd`, `dTod` and
     *   `dLogMagnitudeSquared` as Float64Array(nLam × N), the derivative for
     *   wavelength i and layer k at `[i * N + k]`.
     */
    tmmPhaseJacobianSpectrum(lambdas, n0Jets, nsJets, layerJets, thick, theta_deg, polCode, options = {}) {
        const nLam = lambdas.length;
        const N = thick.length;
        const omegas = options.omegas
            || lambdas.map(lambda => omegaFromLambdaNm(lambda));
        const sinJets = options.sinJets || null;
        const nDeriv = Math.max(1, 10 * N * nLam);

        const lamPtr = this._alloc(nLam);
        const omPtr = this._alloc(nLam);
        const n0Ptr = this._alloc(8 * nLam);
        const nsPtr = this._alloc(8 * nLam);
        const matPtr = this._alloc(Math.max(1, 8 * N * nLam));
        const thPtr = this._alloc(Math.max(1, N));
        const sinPtr = sinJets ? this._alloc(8 * nLam) : 0;
        const outPtr = this._alloc(10 * nLam);
        const derivPtr = this._alloc(nDeriv);

        const lam = this._view(lamPtr, nLam);
        const om = this._view(omPtr, nLam);
        const n0v = this._view(n0Ptr, 8 * nLam);
        const nsv = this._view(nsPtr, 8 * nLam);
        const matv = this._view(matPtr, Math.max(1, 8 * N * nLam));
        const thv = this._view(thPtr, Math.max(1, N));
        const sinv = sinJets ? this._view(sinPtr, 8 * nLam) : null;
        for (let i = 0; i < nLam; i++) {
            lam[i] = lambdas[i];
            om[i] = omegas[i];
            writeJet(n0v, 8 * i, n0Jets[i]);
            writeJet(nsv, 8 * i, nsJets[i]);
            if (sinv) writeJet(sinv, 8 * i, sinJets[i]);
        }
        for (let k = 0; k < N; k++) {
            thv[k] = thick[k];
            const row = layerJets[k];
            const base = k * nLam * 8;
            for (let i = 0; i < nLam; i++) writeJet(matv, base + 8 * i, row[i]);
        }

        this._tmm_phase_jacobian_spectrum(lamPtr, omPtr, nLam, n0Ptr, nsPtr, matPtr, thPtr, N,
            theta_deg, polCode | 0, sinPtr, outPtr, derivPtr);

        // De-interleave: per λ the kernel writes 10 phase values and then a
        // [side][quantity][layer] block of 10 × N derivatives.
        const out = this._view(outPtr, 10 * nLam);
        const deriv = this._view(derivPtr, nDeriv);
        const side = (base, derivBase) => {
            const q = {
                phaseRad: new Float64Array(nLam), gd: new Float64Array(nLam),
                gdd: new Float64Array(nLam), tod: new Float64Array(nLam),
                magnitudeSquared: new Float64Array(nLam),
                dPhaseDeg: new Float64Array(nLam * N), dGd: new Float64Array(nLam * N),
                dGdd: new Float64Array(nLam * N), dTod: new Float64Array(nLam * N),
                dLogMagnitudeSquared: new Float64Array(nLam * N),
            };
            const keys = ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared'];
            const derivKeys = ['dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared'];
            for (let i = 0; i < nLam; i++) {
                for (let j = 0; j < 5; j++) q[keys[j]][i] = out[10 * i + base + j];
                for (let j = 0; j < derivKeys.length; j++) {
                    const from = 10 * N * i + (derivBase + j) * N;
                    for (let k = 0; k < N; k++) q[derivKeys[j]][i * N + k] = deriv[from + k];
                }
            }
            return q;
        };
        const result = { r: side(0, 0), t: side(5, 5) };
        for (const p of [lamPtr, omPtr, n0Ptr, nsPtr, matPtr, thPtr, outPtr, derivPtr]) this.free(p);
        if (sinPtr) this.free(sinPtr);
        return result;
    }

    /**
     * Analytic needle P-function scan : mirrors tmmNeedleScan() in
     * thinFilmMath.js, reshaping the flat WASM output into the SAME nested
     * structure the synthesis scanners consume.
     * @param {{n:[number,number],d:number}[]} layers   used AS-IS (index parity)
     * @param {[number,number][]} candidateNs           candidate ñ
     * @param {number[]} intraFracs                      intra-layer split fractions
     * @returns {{R,T,A,N, gaps:Array, intra:Array}}
     *   gaps[pos][ci] = {dR,dT,dA}  (pos = 0..N)
     *   intra[k][fi]  = {frac, perCand:[{dR,dT,dA}]}
     */
    tmmNeedleScan(lambda_nm, theta_deg, polCode, n0, ns, layers, candidateNs, intraFracs = []) {
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
    }
}

/** Instantiate from raw bytes (ArrayBuffer / Uint8Array). Sets the singleton. */
export async function instantiateTmmWasm(bytes) {
    const { instance } = await WebAssembly.instantiate(
        bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes), wasmImports());
    _instance = new TmmWasmInstance(instance);
    return _instance;
}

/** Renderer/Node helper: fetch the `.wasm` at `url` and instantiate it. */
export function initTmmWasmFromUrl(url) {
    if (_initPromise) return _initPromise;
    _initPromise = (async () => {
        try {
            const resp = await fetch(url);
            if (!resp.ok) throw new Error(`fetch ${url} → ${resp.status}`);
            const buf = await resp.arrayBuffer();
            await instantiateTmmWasm(buf);
            return true;
        } catch (e) {
            // Not built yet / not found : silent fallback to JS.
            _instance = null;
            return false;
        }
    })();
    return _initPromise;
}

export function setTmmWasmEnabled(on) { _enabled = !!on; }
/** TEST-ONLY: inject a TmmWasmInstance (or null) directly, bypassing the .wasm
 *  fetch/instantiate, so the integration seam can be exercised with a mock. */
export function __setTmmWasmInstanceForTest(inst) { _instance = inst; }

// ── Cross-thread plumbing ────────────────────────────────────────────────────
// The renderer (main thread) loads the .wasm bytes once via IPC, instantiates
// its own module, and BROADCASTS the same bytes to each pool/worker (workers
// can't fetch a file:// asset under contextIsolation). Each worker instantiates
// its OWN module (no shared memory) and enables the flag.

let _workerBytes = null;       // main-side: raw bytes to hand to workers
let _workerInitPromise = null; // worker-side: in-flight instantiation

/**
 * MAIN THREAD bootstrap. Store the bytes for worker broadcast and, if the user
 * enabled the feature, instantiate the main-thread module and flip the flag.
 * Safe to call once at startup; failures fall back to JS silently.
 */
export async function initTmmWasmMainThread(bytes, enabled) {
    if (bytes) _workerBytes = bytes;            // remember for workers + later toggles
    if (!enabled) { _enabled = false; return false; }   // toggle off → JS everywhere
    if (!_workerBytes) return false;            // artifact never loaded
    if (!_instance) {                           // instantiate once; reuse on re-toggle
        try { await instantiateTmmWasm(_workerBytes); }
        catch (_) { _instance = null; _enabled = false; return false; }
    }
    _enabled = true;
    return true;
}

/** MAIN THREAD: bytes to ship to a worker : only when the feature is active. */
export function getTmmWasmBytesForWorker() {
    return (_enabled && _workerBytes) ? _workerBytes : null;
}

/**
 * WORKER side: kick off one-time instantiation from broadcast bytes and enable
 * the flag in this worker. Idempotent; no-op without bytes or once instantiated.
 */
export function noteTmmWasmBytes(bytes) {
    if (!bytes || _instance || _workerInitPromise) return;
    _workerInitPromise = instantiateTmmWasm(bytes)
        .then(() => { _enabled = true; return true; })
        .catch(() => { _instance = null; _enabled = false; return false; });
}

/** WORKER side: await any in-flight instantiation before processing a job. */
export function awaitTmmWasmReady() {
    return _workerInitPromise || Promise.resolve(_instance !== null);
}
export function isTmmWasmEnabled() { return _enabled; }
export function isTmmWasmReady() { return _instance !== null; }
/** Active iff the feature flag is on AND a module is instantiated. */
export function tmmWasmActive() { return _enabled && _instance !== null; }
export function getTmmWasm() { return _instance; }
