/**
 * Benchmark tmmcore on the shared cases, cross-check it against Byrnes'
 * results, and write results_js.json for summarize.mjs.
 *
 * Both timed paths make the call a caller makes, rtaSpectrum, once per
 * spectrum: the whole wavelength grid, s and p together.
 *   js    the JavaScript module, which runs the one-wavelength kernel at each
 *         point of the grid
 *   wasm  a WebAssembly instance, which runs the C kernel's
 *         tmm_rta_spectrum() once for the whole grid
 *
 * Run: node bench_js.mjs   (after bench_py.py has written results_py.json)
 */

import { readFileSync, existsSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum } from '../src/js/index.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const { version } = JSON.parse(readFileSync(join(HERE, '..', 'package.json'), 'utf8'));

const wasmPath = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');
if (!existsSync(wasmPath)) {
    console.error('tmm_kernel.wasm not found; build it with `npm run build:wasm`.');
    process.exit(1);
}
const wasm = await instantiateTmmWasm(readFileSync(wasmPath));

const { cases } = JSON.parse(readFileSync(join(HERE, 'cases.json'), 'utf8'));
const py = existsSync(join(HERE, 'results_py.json'))
    ? JSON.parse(readFileSync(join(HERE, 'results_py.json'), 'utf8'))
    : null;

// ── The two spectrum evaluators ──────────────────────────────────────────────
// A case as the stack rtaSpectrum takes, each n and k one value per
// wavelength. It is built once per case, outside the timed region, as the
// other drivers build their inputs.

const medium = pairs => ({ n: pairs.map(p => p[0]), k: pairs.map(p => p[1]) });

function stackOf(c) {
    return {
        incident: medium(c.n0),
        substrate: medium(c.ns),
        layers: c.layerNK.map((row, k) => ({ ...medium(row), d: c.thick[k] })),
    };
}

// Each returns { Rs, Ts, Rp, Tp, ... } over the case's wavelength grid.
const specJS = (c, stack) => rtaSpectrum(stack, c.lambdas, { angle: c.theta_deg });
const specWasm = (c, stack) => wasm.rtaSpectrum(stack, c.lambdas, { angle: c.theta_deg });

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
console.log('case              js-vs-wasm   js-vs-byrnes');
console.log('-'.repeat(46));

for (const c of cases) {
    const stack = stackOf(c);
    const j = specJS(c, stack), w = specWasm(c, stack);
    const dW = compare(j, w);
    let dPY = null;
    if (py) {
        const p = py.results.find(r => r.name === c.name);
        if (p) dPY = compare(j, { Rs: p.Rs, Ts: p.Ts, Rp: p.Rp, Tp: p.Tp });
    }
    console.log(`${c.name.padEnd(16)}  ${dW.toExponential(2).padStart(10)}   `
        + `${dPY === null ? '         n/a' : dPY.toExponential(2).padStart(12)}`);
    rows.push({ c, stack });
}

console.log('\n=== timing: one full spectrum (all lambda, s+p) ===');
console.log('case              N  nLam        js ms      wasm ms');
console.log('-'.repeat(52));

const out = [];
for (const { c, stack } of rows) {
    const tj = timeIt(() => specJS(c, stack));
    const tw = timeIt(() => specWasm(c, stack));
    console.log(`${c.name.padEnd(16)} ${String(c.nLayers).padStart(2)} `
        + `${String(c.lambdas.length).padStart(5)} `
        + `${(tj.best * 1e3).toFixed(3).padStart(12)} `
        + `${(tw.best * 1e3).toFixed(3).padStart(12)}`);
    out.push({
        name: c.name, nLayers: c.nLayers, nLam: c.lambdas.length,
        js_s: tj.best, wasm_spec_s: tw.best,
        js_reps: tj.reps, wasm_spec_reps: tw.reps,
    });
}

writeFileSync(join(HERE, 'results_js.json'),
    JSON.stringify({ versions: { node: process.version, tmmcore: version }, results: out }, null, 1));
console.log('\nwrote results_js.json');
