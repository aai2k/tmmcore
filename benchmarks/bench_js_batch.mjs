/**
 * Batch workload for tmmcore: NBATCH stacks of each 71-point case, each over
 * the full wavelength grid, s and p together. Writes results_js_batch.json for
 * summarize.mjs.
 *
 * Every tmmcore call takes one stack, so a batch is NBATCH rtaSpectrum calls,
 * on the JavaScript module and on a WebAssembly instance. Of the Python
 * libraries in bench_fast.py, only tmm_fast takes all the stacks in one call;
 * the driver there loops over the stacks for Byrnes' tmm and for tmm_faster.
 * What the caller pays is the total time for the same work.
 *
 * The thicknesses are jittered by 10% (one standard deviation), the spread
 * bench_fast.py applies, from this file's own seeded generator: the stacks
 * have the same spread as those of bench_fast.py but are not the same draws.
 *
 * Run: node bench_js_batch.mjs
 */

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum } from '../src/js/index.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const { version } = JSON.parse(readFileSync(join(HERE, '..', 'package.json'), 'utf8'));

const wasm = await instantiateTmmWasm(readFileSync(join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm')));
const { cases } = JSON.parse(readFileSync(join(HERE, 'cases.json'), 'utf8'));

const NBATCH = 128;

// Seeded jitter: mulberry32 for the uniform draws, Box-Muller for the normal.
function makeRng(seed) {
    let a = seed >>> 0;
    return () => {
        a |= 0; a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
function gauss(rng) {
    const u = Math.max(rng(), 1e-12), v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

// A case as the stack rtaSpectrum takes, each n and k one value per
// wavelength, with the thicknesses given.
const medium = pairs => ({ n: pairs.map(p => p[0]), k: pairs.map(p => p[1]) });

function stackOf(c, thick) {
    return {
        incident: medium(c.n0),
        substrate: medium(c.ns),
        layers: c.layerNK.map((row, k) => ({ ...medium(row), d: thick[k] })),
    };
}

// The NBATCH stacks of a case, built before any timing.
function makeBatch(c) {
    const rng = makeRng(0);
    const out = [];
    for (let s = 0; s < NBATCH; s++) {
        out.push(stackOf(c, c.thick.map(d => d * (1 + 0.1 * gauss(rng)))));
    }
    return out;
}

function timeIt(fn, minReps = 2, minSeconds = 1.0) {
    let best = Infinity, reps = 0;
    const t0 = Number(process.hrtime.bigint());
    while (reps < minReps || (Number(process.hrtime.bigint()) - t0) / 1e9 < minSeconds) {
        const a = process.hrtime.bigint();
        fn();
        const dt = Number(process.hrtime.bigint() - a) / 1e9;
        if (dt < best) best = dt;
        reps++;
    }
    return best;
}

console.log(`\n=== BATCH of ${NBATCH} stacks, full grid, both pols (single-threaded) ===`);
console.log('case            N  nLam       js ms     wasm ms   wasm per-stack ms');
console.log('-'.repeat(70));

const out = [];
for (const c of cases) {
    if (!c.name.endsWith('g71')) continue;
    const batch = makeBatch(c);
    const options = { angle: c.theta_deg };

    const tj = timeIt(() => {
        for (const stack of batch) rtaSpectrum(stack, c.lambdas, options);
    });

    const tw = timeIt(() => {
        for (const stack of batch) wasm.rtaSpectrum(stack, c.lambdas, options);
    });

    console.log(`${c.name.padEnd(14)} ${String(c.nLayers).padStart(2)} `
        + `${String(c.lambdas.length).padStart(5)} ${(tj * 1e3).toFixed(1).padStart(11)} `
        + `${(tw * 1e3).toFixed(1).padStart(11)} ${(tw / NBATCH * 1e3).toFixed(3).padStart(19)}`);
    out.push({ name: c.name, nLayers: c.nLayers, nLam: c.lambdas.length,
               nbatch: NBATCH, js_s: tj, wasm_s: tw });
}

writeFileSync(join(HERE, 'results_js_batch.json'),
    JSON.stringify({ versions: { node: process.version, tmmcore: version }, nbatch: NBATCH,
                     results: out }, null, 1));
console.log('\nwrote results_js_batch.json');
