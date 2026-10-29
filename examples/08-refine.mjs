/**
 * Optimizing a coating with the analytic Jacobian.
 *
 * `jacobianSpectrum` returns ∂R/∂d, ∂T/∂d and ∂A/∂d for every layer at every
 * wavelength, exactly, from the same characteristic-matrix product that
 * produced R, not from finite differences and not from automatic
 * differentiation.
 *
 * Below, a deliberately detuned antireflection stack is refined by damped
 * Gauss-Newton (Levenberg-Marquardt) against a target of zero reflectance
 * across the visible. The residuals are R(λᵢ); the Jacobian rows are ∂R(λᵢ)/∂dⱼ.
 *
 * A finite-difference Jacobian would need one extra spectrum evaluation per
 * layer per iteration. Here it comes back with the spectrum.
 *
 *   node examples/08-refine.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum, jacobianSpectrum } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';
import { solve } from './_linalg.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const low = lam => 1.377 + 2000 / (lam * lam);
const high = lam => 2.170 + 8000 / (lam * lam);
const mid = lam => 1.680 + 5000 / (lam * lam);
const glass = lam => 1.504 + 4200 / (lam * lam);

const MATS = [low, high, mid];
const START = [70, 150, 40];        // nm, a poor starting guess
const GRID = [];
for (let lam = 420; lam <= 680; lam += 10) GRID.push(lam);

// The design with thicknesses d, its indices evaluated on a wavelength grid.
const stackOver = (lams, d) => ({
    incident: 1,
    substrate: { n: lams.map(glass) },
    layers: MATS.map((mat, i) => ({ n: lams.map(mat), d: d[i] })),
});

// Merit: root-mean-square reflectance over the band, unpolarized at normal
// incidence (where s and p coincide, so one evaluation suffices).
function residuals(d) {
    return Array.from(rtaSpectrum(stackOver(GRID, d), GRID, { pol: 's' }).R);
}
const merit = r => Math.sqrt(r.reduce((a, v) => a + v * v, 0) / r.length);

// One row of ∂R/∂dⱼ per wavelength, from the grid's wavelength-major array.
function jacobianRows(d) {
    const { dRdd, N } = jacobianSpectrum(stackOver(GRID, d), GRID, { pol: 's' });
    return GRID.map((_, i) => dRdd.subarray(i * N, (i + 1) * N));
}

// The damped step δ from (JᵀJ + λ diag(JᵀJ)) δ = −Jᵀr; null if singular.
function solveLM(J, r, damping) {
    const n = J[0].length;
    const A = Array.from({ length: n }, () => new Array(n).fill(0));
    const b = new Array(n).fill(0);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            let s = 0;
            for (let k = 0; k < J.length; k++) s += J[k][i] * J[k][j];
            A[i][j] = s + (i === j ? damping * s + 1e-18 : 0);
        }
        let g = 0;
        for (let k = 0; k < J.length; k++) g += J[k][i] * r[k];
        b[i] = -g;
    }
    return solve(A, b);
}

let d = [...START];
let r = residuals(d);
let F = merit(r);
let damping = 1e-3;

console.log(`start        rms R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);

for (let iter = 1; iter <= 40; iter++) {
    const step = solveLM(jacobianRows(d), r, damping);
    if (!step) break;

    const trial = d.map((v, i) => Math.max(1, v + step[i]));
    const rt = residuals(trial);
    const Ft = merit(rt);

    if (Ft < F) {
        d = trial; r = rt;
        const gain = F - Ft;
        F = Ft;
        damping = Math.max(damping * 0.5, 1e-9);
        if (iter <= 3 || gain / F > 1e-3) {
            console.log(`iteration ${String(iter).padStart(2)}  rms R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);
        }
        if (gain < 1e-12) break;
    } else {
        damping *= 4;                 // uphill: shorten the step and retry
        if (damping > 1e8) break;
    }
}

console.log(`\nrefined      rms R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);

const plotGrid = [];
for (let lam = 400; lam <= 700; lam += 2) plotGrid.push(lam);
const curve = dv => Array.from(rtaSpectrum(stackOver(plotGrid, dv), plotGrid, { pol: 's' }).R, R => R * 100);

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'refinement.svg'), linePlot({
    title: 'Refinement driven by the analytic thickness Jacobian',
    xLabel: 'wavelength (nm)',
    yLabel: 'reflectance (%)',
    yMin: 0,
    series: [
        { label: 'starting guess', x: plotGrid, y: curve(START), dash: '5 4' },
        { label: 'refined', x: plotGrid, y: curve(d) },
    ],
}));
console.log('wrote docs/img/refinement.svg');
