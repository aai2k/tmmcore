/**
 * Batch workload for the tmmcore kernels: NBATCH distinct stacks, each over
 * the full wavelength grid, both polarizations.
 *
 * The vectorized Python libraries fold this into one call; tmmcore has no
 * batch entry point, so it simply loops. That is the honest comparison — what
 * matters to the caller is total time for the same work.
 *
 * Stack thicknesses are jittered with the same generator and seed logic as
 * bench_fast.py so both sides evaluate comparably diverse stacks.
 */

import { readFileSync, existsSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));

const wasm = await instantiateTmmWasm(readFileSync(join(HERE, '..', 'src', 'tmm_kernel.wasm')));
const { cases } = JSON.parse(readFileSync(join(HERE, 'cases.json'), 'utf8'));

const NBATCH = 128;

// Deterministic jitter (mulberry32 + Box-Muller) — matches the 10% relative
// spread bench_fast.py applies, without needing the exact same RNG stream.
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

function makeBatch(c) {
    const rng = makeRng(0);
    const out = [];
    for (let s = 0; s < NBATCH; s++) {
        out.push(c.thick.map(d => d * (1 + 0.1 * gauss(rng))));
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

    const tj = timeIt(() => {
        for (const th of batch) {
            for (let i = 0; i < c.lambdas.length; i++) {
                const L = c.layerNK.map((row, k) => ({ n: row[i], d: th[k] }));
                tmm(c.lambdas[i], c.theta_deg, 's', c.n0[i], c.ns[i], L);
                tmm(c.lambdas[i], c.theta_deg, 'p', c.n0[i], c.ns[i], L);
            }
        }
    });

    const tw = timeIt(() => {
        for (const th of batch) {
            wasm.tmmSpectrum(c.lambdas, c.n0, c.ns, c.layerNK, th, c.theta_deg);
        }
    });

    console.log(`${c.name.padEnd(14)} ${String(c.nLayers).padStart(2)} `
        + `${String(c.lambdas.length).padStart(5)} ${(tj * 1e3).toFixed(1).padStart(11)} `
        + `${(tw * 1e3).toFixed(1).padStart(11)} ${(tw / NBATCH * 1e3).toFixed(3).padStart(19)}`);
    out.push({ name: c.name, nLayers: c.nLayers, nLam: c.lambdas.length,
               nbatch: NBATCH, js_s: tj, wasm_s: tw });
}

writeFileSync(join(HERE, 'results_js_batch.json'),
    JSON.stringify({ node: process.version, nbatch: NBATCH, results: out }, null, 1));
console.log('\nwrote results_js_batch.json');
