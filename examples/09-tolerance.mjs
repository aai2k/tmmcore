/**
 * What random thickness errors cost a design, from the exact Hessian.
 *
 * The three-layer antireflection coating of 08-refine.mjs, with its materials and
 * its 420-680 nm grid, at normal incidence. The merit here is the mean
 * reflectance over the grid,
 *
 *     F = (1/M) Σᵢ R(λᵢ),
 *
 * Macleod's figure of merit with every target zero, equal weights and q = 1
 * (Thin-Film Optical Filters 5th ed., §3.4.6, Eq. 3.31). Its gradient and Hessian
 * are the grid means of ∂R/∂dⱼ and ∂²R/∂dⱼ∂dₖ, which hessianSpectrum returns
 * exactly, the first derivatives alongside the second.
 *
 * Newton steps with that Hessian take 08's thicknesses, which minimize the RMS
 * reflectance, to the minimum of F. There the gradient is zero, so an error e in
 * the thicknesses raises the merit by
 *
 *     ΔF = ½ eᵀHe + O(|e|³).
 *
 * If every layer has an independent error of zero mean and standard deviation
 * σ, then E[eⱼeₖ] = σ² when j = k and 0 otherwise, and the mean rise is
 *
 *     E[ΔF] = ½ σ² tr(H),
 *
 * with a correction of order σ⁴, since the odd moments of e average to zero.
 * The eigenvectors of H say which combination of errors costs the most.
 *
 * The check is a Monte Carlo of exact evaluations with errors of that kind,
 * drawn from a normal population: Ritchie's method, as Macleod describes it for
 * tolerancing (§13.7, Tolerances).
 *
 * Writes the plot to docs/img/tolerance.svg.
 *
 *   node examples/09-tolerance.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum, hessianSpectrum } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';
import { solve, symmetricEigen } from './_linalg.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const low = lam => 1.377 + 2000 / (lam * lam);
const high = lam => 2.170 + 8000 / (lam * lam);
const mid = lam => 1.680 + 5000 / (lam * lam);
const glass = lam => 1.504 + 4200 / (lam * lam);

const MATS = [low, high, mid];
const REFINED = [93.62, 117.56, 75.42];      // nm, as 08-refine.mjs prints them
const GRID = [];
for (let lam = 420; lam <= 680; lam += 10) GRID.push(lam);
const M = GRID.length;
const N = MATS.length;

const stackAt = d => ({
    incident: 1,
    substrate: { n: GRID.map(glass) },
    layers: MATS.map((mat, i) => ({ n: GRID.map(mat), d: d[i] })),
});

// s alone, since at normal incidence s and p coincide.
const merit = d => rtaSpectrum(stackAt(d), GRID, { pol: 's' }).R.reduce((a, R) => a + R, 0) / M;

// F, its gradient g and its Hessian H, as fractions per nm and per nm².
function meritWithDerivatives(d) {
    const { R, dRdd, d2Rdd } = hessianSpectrum(stackAt(d), GRID, { pol: 's' });
    const g = new Array(N).fill(0);
    const H = Array.from({ length: N }, () => new Array(N).fill(0));
    for (let i = 0; i < M; i++) {
        for (let j = 0; j < N; j++) {
            g[j] += dRdd[i * N + j] / M;
            for (let k = 0; k < N; k++) H[j][k] += d2Rdd[(i * N + j) * N + k] / M;
        }
    }
    return { F: R.reduce((a, v) => a + v, 0) / M, g, H };
}

const norm = x => Math.hypot(...x);
const percent = v => v * 100;
const columns = (cells, widths) => '  ' + cells.map((c, i) => String(c).padStart(widths[i])).join('');

// ── Newton steps to the minimum of F ─────────────────────────────────────────

console.log('Newton steps on the mean reflectance, from the thicknesses of 08-refine.mjs\n');
const NEWTON = [4, 17, 14, 10, 9, 9];
console.log(columns(['step', 'mean R (%)', '|∇F| (%/nm)', 'd1 (nm)', 'd2', 'd3'], NEWTON));
let d = [...REFINED];
let at = meritWithDerivatives(d);
let previous = Infinity;
for (let step = 0; step <= 10; step++) {
    const gradient = norm(at.g);
    console.log(columns([step, percent(at.F).toFixed(12), percent(gradient).toExponential(2),
        ...d.map(x => x.toFixed(3))], NEWTON));
    // Stop once a step no longer shrinks the gradient: it is at the rounding
    // floor. Newton converges quadratically here, so ten steps is only a guard.
    if (gradient > previous / 10) break;
    previous = gradient;
    const delta = solve(at.H, at.g.map(x => -x));
    d = d.map((x, i) => x + delta[i]);
    at = meritWithDerivatives(d);
}
const F0 = at.F;
const H = at.H;

// ── The shape of the minimum ─────────────────────────────────────────────────

const eigen = symmetricEigen(H);
const trace = H.reduce((sum, row, i) => sum + row[i], 0);
console.log('\nHessian at the minimum: rise in mean R for a 1 nm error along each eigenvector\n');
const EIGEN = [9, 9, 8, 8];
console.log(columns(['rise (%)', 'd1', 'd2', 'd3'], EIGEN));
for (const { value, vector } of eigen) {
    console.log(columns([percent(value / 2).toExponential(2), ...vector.map(x => x.toFixed(3))], EIGEN));
}

// ── Second-order prediction against a Monte Carlo ────────────────────────────

// mulberry32, a small seeded generator, so every run draws the same errors.
function mulberry32(seed) {
    return () => {
        seed = (seed + 0x6D2B79F5) | 0;
        let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}
const uniform = mulberry32(9);
// Box-Muller; 1 − u is in (0, 1], so the logarithm is finite.
const normal = () => Math.sqrt(-2 * Math.log(1 - uniform())) * Math.cos(2 * Math.PI * uniform());

// One set of unit draws, scaled by each σ, so the rows differ only in σ.
const SAMPLES = 10000;
const draws = Array.from({ length: SAMPLES }, () => Array.from({ length: N }, normal));
const SIGMAS = [0.5, 1, 2, 3, 4, 5];   // nm

const predicted = sigma => 0.5 * sigma * sigma * trace;
const rows = SIGMAS.map(sigma => {
    let sum = 0, sumSquares = 0;
    for (const z of draws) {
        const rise = merit(d.map((x, i) => x + sigma * z[i])) - F0;
        sum += rise; sumSquares += rise * rise;
    }
    const mean = sum / SAMPLES;
    const error = Math.sqrt((sumSquares / SAMPLES - mean * mean) / (SAMPLES - 1));
    return { sigma, mean, error };
});

console.log(`\nmean rise in mean R, errors of standard deviation σ in every layer, ${SAMPLES} designs per row\n`);
const RISE = [6, 16, 18, 11, 8];
console.log(columns(['σ (nm)', '½σ² tr(H) (%)', 'Monte Carlo (%)', '', 'ratio'], RISE));
for (const { sigma, mean, error } of rows) {
    console.log(columns([sigma.toFixed(1), percent(predicted(sigma)).toExponential(3),
        percent(mean).toExponential(3), `± ${percent(error).toExponential(1)}`,
        (mean / predicted(sigma)).toFixed(3)], RISE));
}

// The y range leaves the top of the plot clear for the legend.
const fine = Array.from({ length: 101 }, (_, i) => i * 0.05);
mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'tolerance.svg'), linePlot({
    title: 'Mean reflectance rise from random thickness errors',
    xLabel: 'standard deviation of the thickness errors, σ (nm)',
    yLabel: 'rise in mean R (%)',
    yMin: 0, yMax: 0.2,
    series: [
        { label: '½σ² tr(H)', x: fine, y: fine.map(s => percent(predicted(s))) },
        { label: 'Monte Carlo', x: [0, ...SIGMAS], y: [0, ...rows.map(r => percent(r.mean))], dash: '5 4' },
    ],
}));
console.log('\nwrote docs/img/tolerance.svg');
