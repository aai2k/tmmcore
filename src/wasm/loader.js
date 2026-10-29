/**
 * loader.js : loads the WebAssembly kernel and wraps an instance of it. A
 * TmmWasmInstance binds the kernel's exports, carries the bindings of
 * memory.js, rta.js, growing.js, derivatives.js, phase.js and jet.js, and
 * offers the package's calls on them.
 *
 * Acceleration is opt-in, and JavaScript is the fallback: initTmmWasmFromUrl()
 * resolves false when the `.wasm` cannot be fetched or instantiated,
 * tmmWasmActive() tells a caller whether to use the instance, and getTmmWasm()
 * is null until one is loaded. Results agree with the JavaScript calls to
 * float64 round-off.
 *
 * Instances are not shared across threads : there is no shared memory, so each
 * context instantiates its own from the same bytes. Use `instantiateTmmWasm()`
 * where you already hold the bytes (a worker receives them in its init message)
 * and `initTmmWasmFromUrl()` where the artifact is fetchable.
 */

import { stackCalls } from '../js/calls.js';
import { memoryMethods } from './memory.js';
import { wasmJetCalls } from './jet.js';
import { rtaMethods, rtaSpectrumWasm } from './rta.js';
import { derivativeMethods, jacobianSpectrumWasm, hessianSpectrumWasm } from './derivatives.js';
import { growingMethods } from './growing.js';
import { phaseMethods, phaseSpectrumWasm, phaseJacobianSpectrumWasm } from './phase.js';

let _instance = null;     // TmmWasmInstance | null
let _enabled = false;     // feature flag (default OFF)
let _initPromise = null;  // de-dupe concurrent init

// Permissive imports: a STANDALONE_WASM build of pure-math C usually needs no
// imports, but ALLOW_MEMORY_GROWTH may emit `emscripten_notify_memory_growth`,
// and some toolchains emit WASI stubs. Cover them so instantiation never throws.
function wasmImports() {
    return {
        env: { emscripten_notify_memory_growth: () => {} },
        wasi_snapshot_preview1: new Proxy({}, { get: () => () => 0 }),
    };
}

// The instance's kernels in the form calls.js builds on: pol 's' or 'p', the
// stack at one wavelength as an object, a grid as stack.js lays it out.
function wasmKernels(wasm) {
    const code = pol => (pol === 'p' ? 1 : 0);
    return {
        rta: (lambda, theta, pol, media) => wasm._rta(lambda, theta, code(pol), media),
        jacobian: (lambda, theta, pol, media) => wasm._jacobian(lambda, theta, code(pol), media),
        hessian: (lambda, theta, pol, media) => wasm._hessian(lambda, theta, code(pol), media),
        phase: (lambda, theta, pol, media) => wasm._phase(lambda, theta, code(pol), media),
        phaseJacobian: (lambda, theta, pol, media) => wasm._phaseJacobian(lambda, theta, code(pol), media),
        rtaSpectrum: grid => rtaSpectrumWasm(wasm, grid),
        jacobianSpectrum: (grid, pol) => jacobianSpectrumWasm(wasm, grid, code(pol)),
        hessianSpectrum: (grid, pol) => hessianSpectrumWasm(wasm, grid, code(pol)),
        phaseSpectrum: (grid, pol) => phaseSpectrumWasm(wasm, grid, code(pol)),
        phaseJacobianSpectrum: (grid, pol) => phaseJacobianSpectrumWasm(wasm, grid, code(pol)),
    };
}

// Every kernel export of the C kernel. Each is bound as `_<name>` on the
// instance; a module missing any of them is refused.
const KERNEL_EXPORTS = [
    'tmm_rta', 'tmm_rta_spectrum',
    'tmm_jacobian', 'tmm_jacobian_spectrum',
    'tmm_hessian', 'tmm_hessian_spectrum',
    'tmm_needle_scan', 'tmm_needle_scan_spectrum',
    'tmm_phase', 'tmm_phase_spectrum',
    'tmm_phase_jacobian', 'tmm_phase_jacobian_spectrum',
    'tmm_monitor_curve', 'tmm_deposition_spectra',
    'tmm_growing_stack_create', 'tmm_growing_stack_set_top',
    'tmm_growing_stack_sample', 'tmm_growing_stack_free',
    'tmm_jet_constant', 'tmm_jet_from_derivatives', 'tmm_jet_derivatives',
    'tmm_jet_add', 'tmm_jet_subtract', 'tmm_jet_scale', 'tmm_jet_multiply',
    'tmm_jet_reciprocal', 'tmm_jet_divide', 'tmm_jet_sqrt', 'tmm_jet_exp',
    'tmm_jet_log', 'tmm_jet_power', 'tmm_jet_sin_cos', 'tmm_jet_sin', 'tmm_jet_cos',
    'tmm_jet_compose', 'tmm_wavelength_omega_jet', 'tmm_jet_with_imaginary_part',
    'tmm_jet_real_part', 'tmm_jet_clamp_real_minimum', 'tmm_jet_clamp_imaginary',
    'tmm_jet_is_finite', 'tmm_omega_from_lambda_nm',
];

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
        const missing = KERNEL_EXPORTS.filter(name => typeof ex[name] !== 'function');
        if (!this.malloc || !this.free) missing.push('malloc/free');
        if (missing.length) {
            throw new Error(`tmmWasm: exports missing from module: ${missing.join(', ')}`);
        }
        for (const name of KERNEL_EXPORTS) this[`_${name}`] = ex[name];
        this._scratchPtr = 0;   // persistent per-call scratch arena (lazy)
        this._scratchN = 0;
        Object.assign(this, wasmJetCalls(this), stackCalls(wasmKernels(this)));
    }
}

// The bindings of each area, as methods of every instance, not enumerable, as
// the methods of a class are.
for (const methods of [memoryMethods, rtaMethods, growingMethods, derivativeMethods, phaseMethods]) {
    for (const [name, method] of Object.entries(methods)) {
        Object.defineProperty(TmmWasmInstance.prototype, name, { value: method, writable: true, configurable: true });
    }
}

/** Instantiate from raw bytes (ArrayBuffer / Uint8Array). Sets the singleton. */
export async function instantiateTmmWasm(bytes) {
    const { instance } = await WebAssembly.instantiate(
        bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes), wasmImports());
    _instance = new TmmWasmInstance(instance);
    return _instance;
}

/** Fetch the `.wasm` at `url` and instantiate it. */
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
// The main thread holds the .wasm bytes, instantiates its own module, and hands
// the same bytes to each worker, for workers that cannot fetch the file
// themselves. Each worker instantiates its OWN module (no shared memory) and
// enables the flag.

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
