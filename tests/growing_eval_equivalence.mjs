/**
 * Persistent growing-layer evaluator vs the plain JavaScript TMM.
 *
 * growingEval() folds the completed stack once per (λ, pol) and keeps it in
 * kernel memory; every sample() answers one thickness of the growing layer
 * across the whole wavelength grid. The oracle is the naive path: tmm() over
 * the full stack per (λ, thickness), with the reverse pass evaluated
 * explicitly on the reversed stack at the refracted angle. Also exercised:
 * changing the growing material via setTop() on a live handle, d = 0 (bare
 * completed stack), buffer reuse through sample(d, out), and free().
 *
 * Skips cleanly if tmm_kernel.wasm has not been built or predates the kernel.
 *
 *   node tests/growing_eval_equivalence.mjs
 */

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');

if (!existsSync(WASM)) {
    console.log('SKIP : src/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}

const k = await instantiateTmmWasm(readFileSync(WASM));
if (!k) { console.log('FAIL : instantiation returned null'); process.exit(1); }
if (!k.hasGrowingEval()) {
    console.log('SKIP : tmm_kernel.wasm predates the growing evaluator. Rebuild it.');
    process.exit(0);
}

const ABS = 1e-9;   // R and T live in [0, 1]

let checks = 0;
let failures = 0;
let worstDiff = 0;

function near(js, wa, label) {
    checks++;
    const diff = Math.abs(js - wa);
    if (diff > worstDiff) worstDiff = diff;
    if (diff <= ABS) return;
    failures++;
    if (failures <= 10) console.log(`  FAIL ${label}: js=${js} wasm=${wa} diff=${diff}`);
}

const n0c = [1, 0], nsc = [1.52, 0];
const H = [2.35, 0.0005], L = [1.46, 0], AG = [0.15, 3.2];
const LAMBDAS = [420, 550, 700, 980, 1550];

function thetaSubDeg(theta_deg) {
    const sinSub = (n0c[0] * Math.sin(theta_deg * Math.PI / 180)) / nsc[0];
    return Math.asin(sinSub) * 180 / Math.PI;
}

const BASES = {
    'bare substrate': [],
    '3-layer dielectric': [{ n: H, d: 60 }, { n: L, d: 95 }, { n: H, d: 110 }],
    'stack with silver + zero': [{ n: L, d: 140 }, { n: H, d: 0 }, { n: AG, d: 22 }, { n: H, d: 88 }],
    '40-layer QW': Array.from({ length: 40 }, (_, i) => (i % 2 ? { n: H, d: 58.5 } : { n: L, d: 94.2 })),
};
const D_SAMPLES = [0, 0.5, 12, 58.5, 94.2, 275];

const n0List = LAMBDAS.map(() => n0c);
const nsList = LAMBDAS.map(() => nsc);

// A d > 0 sample before setTop must throw, not return unwritten memory;
// d ≤ 0 (the bare completed stack) needs no growing layer.
{
    const ev = k.growingEval(LAMBDAS, n0List, nsList, [], [], 0);
    let threw = false;
    try { ev.sample(10); } catch (_) { threw = true; }
    checks++;
    if (!threw) { failures++; console.log('  FAIL sample() before setTop did not throw'); }
    const bare = ev.sample(0);
    LAMBDAS.forEach((lam, li) => {
        const s = tmm(lam, 0, 's', n0c, nsc, []);
        near(s.R, bare.Rs[li], `bare stack before setTop λ=${lam} Rs`);
    });
    ev.free();
}

for (const [baseName, base] of Object.entries(BASES)) {
    for (const th of [0, 30, 60]) {
        const layerNK = base.map(l => LAMBDAS.map(() => l.n));
        const thick = base.map(l => l.d);
        const ev = k.growingEval(LAMBDAS, n0List, nsList, layerNK, thick, th);
        const thSub = thetaSubDeg(th);
        let out = null;   // reused across samples, as a scan loop would
        // Two growing materials on the SAME handle: setTop must fully replace
        // the previous top.
        for (const [ngName, ng] of [['H', H], ['Ag', AG]]) {
            ev.setTop(LAMBDAS.map(() => ng));
            for (const d of D_SAMPLES) {
                out = ev.sample(d, out);
                const full = [{ n: ng, d }, ...base];
                LAMBDAS.forEach((lam, li) => {
                    const label = `${baseName} +${ngName} λ=${lam} θ=${th} d=${d}`;
                    const s = tmm(lam, th, 's', n0c, nsc, full);
                    const p = tmm(lam, th, 'p', n0c, nsc, full);
                    near(s.R, out.Rs[li], `${label} Rs`);
                    near(s.T, out.Ts[li], `${label} Ts`);
                    near(p.R, out.Rp[li], `${label} Rp`);
                    near(p.T, out.Tp[li], `${label} Tp`);
                    const rs = tmm(lam, thSub, 's', nsc, n0c, [...full].reverse());
                    const rp = tmm(lam, thSub, 'p', nsc, n0c, [...full].reverse());
                    near(rs.R, out.Rrs[li], `${label} Rrs`);
                    near(rp.R, out.Rrp[li], `${label} Rrp`);
                });
            }
        }
        ev.free();
        let threw = false;
        try { ev.sample(10); } catch (_) { threw = true; }
        if (!threw) { failures++; console.log(`  FAIL ${baseName}: sample() after free() did not throw`); }
        checks++;
    }
}

console.log(`growing-eval equivalence: ${checks} checks, ${failures} failures, worst diff ${worstDiff.toExponential(2)}`);
process.exit(failures ? 1 : 0);
