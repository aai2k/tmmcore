/**
 * Growing-stack kernels vs the plain JavaScript TMM.
 *
 * tmm_monitor_curve and tmm_deposition_spectra never rebuild the full stack:
 * they hold the completed product and fold or top it incrementally, taking the
 * reverse-direction reflectance from the anti-transposed matrix. The oracle
 * here is the naive path: tmm() over the full stack at every sample, and the
 * reverse pass evaluated explicitly on the reversed stack at the refracted
 * angle. Reciprocity makes the two formulations equal analytically; float64
 * grouping leaves ULP-level noise, far under the tolerance.
 *
 * Skips cleanly if tmm_kernel.wasm has not been built or predates the kernels.
 *
 *   node tests/growing_equivalence.mjs
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
if (!k.hasGrowingKernels()) {
    console.log('SKIP : tmm_kernel.wasm predates the growing-stack kernels. Rebuild it.');
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

const n0 = [1, 0], ns = [1.52, 0];
const H = [2.35, 0.0005], L = [1.46, 0], AG = [0.15, 3.2];

/** Refracted-ray angle in the substrate, real-part Snell (cases stay far from grazing). */
function thetaSubDeg(theta_deg) {
    const sinSub = (n0[0] * Math.sin(theta_deg * Math.PI / 180)) / ns[0];
    return Math.asin(sinSub) * 180 / Math.PI;
}

// ── tmm_monitor_curve ────────────────────────────────────────────────────────
// A layer grows on each base stack; every sample is checked against the naive
// full-stack evaluation, forward and reverse.

const BASES = {
    'bare substrate': [],
    '3-layer dielectric': [{ n: H, d: 60 }, { n: L, d: 95 }, { n: H, d: 110 }],
    'stack with silver': [{ n: L, d: 140 }, { n: AG, d: 22 }, { n: H, d: 88 }],
    '40-layer QW': Array.from({ length: 40 }, (_, i) => (i % 2 ? { n: H, d: 58.5 } : { n: L, d: 94.2 })),
};
const GROWING = [['H', H], ['L', L], ['Ag', AG]];
const D_SAMPLES = [0, 0.5, 12, 58.5, 94.2, 130, 275];

for (const [baseName, base] of Object.entries(BASES)) {
    for (const lam of [400, 550, 1064]) {
        for (const th of [0, 30, 60]) {
            for (const [ngName, ng] of GROWING) {
                const label = `${baseName} +${ngName} λ=${lam} θ=${th}`;
                const got = k.monitorCurve(lam, th, n0, ns, base, ng, D_SAMPLES);
                const thSub = thetaSubDeg(th);
                D_SAMPLES.forEach((d, i) => {
                    const full = [{ n: ng, d }, ...base];
                    const s = tmm(lam, th, 's', n0, ns, full);
                    const p = tmm(lam, th, 'p', n0, ns, full);
                    near(s.R, got.Rs[i], `${label} d=${d} Rs`);
                    near(s.T, got.Ts[i], `${label} d=${d} Ts`);
                    near(p.R, got.Rp[i], `${label} d=${d} Rp`);
                    near(p.T, got.Tp[i], `${label} d=${d} Tp`);
                    const rs = tmm(lam, thSub, 's', ns, n0, [...full].reverse());
                    const rp = tmm(lam, thSub, 'p', ns, n0, [...full].reverse());
                    near(rs.R, got.Rrs[i], `${label} d=${d} Rrs`);
                    near(rp.R, got.Rrp[i], `${label} d=${d} Rrp`);
                });
            }
        }
    }
}

// ── tmm_deposition_spectra ───────────────────────────────────────────────────
// A run deposited layer by layer; after each step the spectrum is checked
// against the naive evaluation of the prefix, reversed into incident-first
// order. A zero-thickness step must repeat the previous spectrum.

const RUNS = {
    'QW mirror 24': Array.from({ length: 24 }, (_, i) => ({ n: i % 2 ? L : H, d: i % 2 ? 94.2 : 58.5 })),
    'with zero step': [
        { n: H, d: 60 }, { n: L, d: 0 }, { n: L, d: 95 }, { n: AG, d: 18 }, { n: H, d: 110 },
    ],
};
const LAMBDAS = [420, 550, 700, 980];

for (const [runName, dep] of Object.entries(RUNS)) {
    for (const th of [0, 45]) {
        const layerNK = dep.map(l => LAMBDAS.map(() => l.n));
        const thick = dep.map(l => l.d);
        const n0List = LAMBDAS.map(() => n0);
        const nsList = LAMBDAS.map(() => ns);
        const got = k.depositionSpectra(LAMBDAS, n0List, nsList, layerNK, thick, th);
        const thSub = thetaSubDeg(th);
        for (let step = 1; step <= dep.length; step++) {
            const stack = dep.slice(0, step).filter(l => l.d > 0).reverse();
            LAMBDAS.forEach((lam, li) => {
                const at = (step - 1) * LAMBDAS.length + li;
                const label = `${runName} θ=${th} step=${step} λ=${lam}`;
                const s = tmm(lam, th, 's', n0, ns, stack);
                const p = tmm(lam, th, 'p', n0, ns, stack);
                near(s.R, got.Rs[at], `${label} Rs`);
                near(s.T, got.Ts[at], `${label} Ts`);
                near(p.R, got.Rp[at], `${label} Rp`);
                near(p.T, got.Tp[at], `${label} Tp`);
                const rs = tmm(lam, thSub, 's', ns, n0, [...stack].reverse());
                const rp = tmm(lam, thSub, 'p', ns, n0, [...stack].reverse());
                near(rs.R, got.Rrs[at], `${label} Rrs`);
                near(rp.R, got.Rrp[at], `${label} Rrp`);
            });
        }
    }
}

console.log(`growing-stack equivalence: ${checks} checks, ${failures} failures, worst diff ${worstDiff.toExponential(2)}`);
process.exit(failures ? 1 : 0);
