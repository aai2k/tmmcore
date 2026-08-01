/**
 * Benchmark tmmcore's two TMM kernels on the shared cases, cross-check both
 * against Byrnes' results, and print the comparison table.
 *
 * Three timed paths:
 *   js        — the JavaScript `tmm()`, called per (lambda, pol). Same call
 *               granularity as Byrnes' coh_tmm, so this is the like-for-like row.
 *   wasm-one  — kernel `tmm_one()` per (lambda, pol). Same granularity again,
 *               but pays a JS->WASM boundary crossing per call.
 *   wasm-spec — kernel `tmm_spectrum()`, one call for the whole grid and both
 *               polarizations. The path a caller should actually use.
 *
 * Run: node bench_js.mjs   (after bench_py.py has written results_py.json)
 */

import { readFileSync, existsSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));

const wasmPath = join(HERE, '..', 'src', 'tmm_kernel.wasm');
if (!existsSync(wasmPath)) {
    console.error('tmm_kernel.wasm not found — build it with `npm run build:wasm`.');
    process.exit(1);
}
const wasm = await instantiateTmmWasm(readFileSync(wasmPath));

const { cases } = JSON.parse(readFileSync(join(HERE, 'cases.json'), 'utf8'));
const py = existsSync(join(HERE, 'results_py.json'))
    ? JSON.parse(readFileSync(join(HERE, 'results_py.json'), 'utf8'))
    : null;

// ── The three spectrum evaluators ────────────────────────────────────────────
// Each returns { Rs, Ts, Rp, Tp } over the case's lambda grid.

function layersAt(c, i) {
    const out = new Array(c.nLayers);
    for (let k = 0; k < c.nLayers; k++) out[k] = { n: c.layerNK[k][i], d: c.thick[k] };
    return out;
}

function specJS(c) {
    const nLam = c.lambdas.length;
    const Rs = new Float64Array(nLam), Ts = new Float64Array(nLam);
    const Rp = new Float64Array(nLam), Tp = new Float64Array(nLam);
    for (let i = 0; i < nLam; i++) {
        const L = layersAt(c, i);
        const s = tmm(c.lambdas[i], c.theta_deg, 's', c.n0[i], c.ns[i], L);
        const p = tmm(c.lambdas[i], c.theta_deg, 'p', c.n0[i], c.ns[i], L);
        Rs[i] = s.R; Ts[i] = s.T; Rp[i] = p.R; Tp[i] = p.T;
    }
    return { Rs, Ts, Rp, Tp };
}

function specWasmOne(c) {
    const nLam = c.lambdas.length;
    const Rs = new Float64Array(nLam), Ts = new Float64Array(nLam);
    const Rp = new Float64Array(nLam), Tp = new Float64Array(nLam);
    for (let i = 0; i < nLam; i++) {
        const L = layersAt(c, i);
        const s = wasm.tmmOne(c.lambdas[i], c.theta_deg, 0, c.n0[i], c.ns[i], L);
        const p = wasm.tmmOne(c.lambdas[i], c.theta_deg, 1, c.n0[i], c.ns[i], L);
        Rs[i] = s.R; Ts[i] = s.T; Rp[i] = p.R; Tp[i] = p.T;
    }
    return { Rs, Ts, Rp, Tp };
}

function specWasmBatch(c) {
    const r = wasm.tmmSpectrum(c.lambdas, c.n0, c.ns, c.layerNK, c.thick, c.theta_deg);
    return { Rs: r.Rs, Ts: r.Ts, Rp: r.Rp, Tp: r.Tp };
}

// ── Timing: best-of, with rep and wall-time floors (mirrors bench_py.py) ──────

function timeIt(fn, minReps = 3, minSeconds = 1.5) {
    let best = Infinity, reps = 0;
    const t0 = Number(process.hrtime.bigint());
    while (reps < minReps || (Number(process.hrtime.bigint()) - t0) / 1e9 < minSeconds) {
        const a = process.hrtime.bigint();
        fn();
        const dt = Number(process.hrtime.bigint() - a) / 1e9;
        if (dt < best) best = dt;
        reps++;
    }
    return { best, reps };
}

// ── Agreement check ──────────────────────────────────────────────────────────

function maxAbsDiff(a, b) {
    let m = 0;
    for (let i = 0; i < a.length; i++) {
        const d = Math.abs(a[i] - b[i]);
        if (d > m) m = d;
    }
    return m;
}

function compare(x, y) {
    return Math.max(maxAbsDiff(x.Rs, y.Rs), maxAbsDiff(x.Ts, y.Ts),
                    maxAbsDiff(x.Rp, y.Rp), maxAbsDiff(x.Tp, y.Tp));
}

// ── Run ──────────────────────────────────────────────────────────────────────

const rows = [];
console.log('\n=== agreement (max |diff| in R/T over the grid, both pols) ===');
console.log('case            js-vs-wasm1   js-vs-wasmSpec   js-vs-byrnes');
console.log('-'.repeat(64));

for (const c of cases) {
    const j = specJS(c), w1 = specWasmOne(c), wb = specWasmBatch(c);
    const dW1 = compare(j, w1), dWB = compare(j, wb);
    let dPY = null;
    if (py) {
        const p = py.results.find(r => r.name === c.name);
        if (p) dPY = compare(j, { Rs: p.Rs, Ts: p.Ts, Rp: p.Rp, Tp: p.Tp });
    }
    console.log(`${c.name.padEnd(14)}  ${dW1.toExponential(2).padStart(10)}   `
        + `${dWB.toExponential(2).padStart(12)}   `
        + `${dPY === null ? '     n/a' : dPY.toExponential(2).padStart(12)}`);
    rows.push({ c, dW1, dWB, dPY });
}

console.log('\n=== timing: one full spectrum (all lambda, s+p) ===');
console.log('case            N  nLam        js ms   wasm1 ms   wasmSpec ms');
console.log('-'.repeat(64));

const out = [];
for (const { c } of rows) {
    const tj = timeIt(() => specJS(c));
    const t1 = timeIt(() => specWasmOne(c));
    const tb = timeIt(() => specWasmBatch(c));
    console.log(`${c.name.padEnd(14)} ${String(c.nLayers).padStart(2)} `
        + `${String(c.lambdas.length).padStart(5)} `
        + `${(tj.best * 1e3).toFixed(3).padStart(12)} `
        + `${(t1.best * 1e3).toFixed(3).padStart(10)} `
        + `${(tb.best * 1e3).toFixed(3).padStart(13)}`);
    out.push({
        name: c.name, nLayers: c.nLayers, nLam: c.lambdas.length,
        js_s: tj.best, wasm_one_s: t1.best, wasm_spec_s: tb.best,
        js_reps: tj.reps, wasm_one_reps: t1.reps, wasm_spec_reps: tb.reps,
    });
}

writeFileSync(join(HERE, 'results_js.json'),
    JSON.stringify({ node: process.version, results: out }, null, 1));
console.log('\nwrote results_js.json');
