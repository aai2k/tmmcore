/**
 * Optimizing a coating with the analytic Jacobian.
 *
 * This is what tmmcore offers that other transfer-matrix packages do not:
 * `tmmThicknessJacobian` returns ∂R/∂d, ∂T/∂d and ∂A/∂d for every layer,
 * exactly — from the same characteristic-matrix product that produced R, not
 * from finite differences and not from automatic differentiation.
 *
 * Below, a deliberately detuned antireflection stack is refined by damped
 * Gauss–Newton (Levenberg–Marquardt) against a target of zero reflectance
 * across the visible. The residuals are R(λᵢ); the Jacobian rows are ∂R(λᵢ)/∂dⱼ.
 *
 * A finite-difference Jacobian would need one extra spectrum evaluation per
 * layer per iteration. Here it comes back with the spectrum.
 *
 *   node examples/04-refine.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm, tmmThicknessJacobian } from '../src/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const low = lam => [1.377 + 2000 / (lam * lam), 0];
const high = lam => [2.170 + 8000 / (lam * lam), 0];
const mid = lam => [1.680 + 5000 / (lam * lam), 0];
const glass = lam => [1.504 + 4200 / (lam * lam), 0];

const MATS = [low, high, mid];
const START = [70, 150, 40];        // nm — a poor starting guess
const GRID = [];
for (let lam = 420; lam <= 680; lam += 10) GRID.push(lam);

const stackAt = (lam, d) => MATS.map((mat, i) => ({ n: mat(lam), d: d[i] }));

// Merit: root-mean-square reflectance over the band, unpolarized at normal
// incidence (where s and p coincide, so one evaluation suffices).
function residuals(d) {
    return GRID.map(lam => tmm(lam, 0, 's', [1, 0], glass(lam), stackAt(lam, d)).R);
}
const merit = r => Math.sqrt(r.reduce((a, v) => a + v * v, 0) / r.length);

function jacobian(d) {
    return GRID.map(lam =>
        tmmThicknessJacobian(lam, 0, 's', [1, 0], glass(lam), stackAt(lam, d)).dRdd);
}

// Solve (JᵀJ + λI) δ = −Jᵀr by Gaussian elimination with partial pivoting.
function solveLM(J, r, damping) {
    const n = J[0].length;
    const A = Array.from({ length: n }, () => new Float64Array(n + 1));
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            let s = 0;
            for (let k = 0; k < J.length; k++) s += J[k][i] * J[k][j];
            A[i][j] = s + (i === j ? damping * s + 1e-18 : 0);
        }
        let g = 0;
        for (let k = 0; k < J.length; k++) g += J[k][i] * r[k];
        A[i][n] = -g;
    }
    for (let c = 0; c < n; c++) {
        let piv = c;
        for (let i = c + 1; i < n; i++) if (Math.abs(A[i][c]) > Math.abs(A[piv][c])) piv = i;
        [A[c], A[piv]] = [A[piv], A[c]];
        if (Math.abs(A[c][c]) < 1e-30) return null;
        for (let i = c + 1; i < n; i++) {
            const f = A[i][c] / A[c][c];
            for (let j = c; j <= n; j++) A[i][j] -= f * A[c][j];
        }
    }
    const x = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let s = A[i][n];
        for (let j = i + 1; j < n; j++) s -= A[i][j] * x[j];
        x[i] = s / A[i][i];
    }
    return x;
}

let d = [...START];
let r = residuals(d);
let F = merit(r);
let damping = 1e-3;

console.log(`start        mean R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);

for (let iter = 1; iter <= 40; iter++) {
    const step = solveLM(jacobian(d), r, damping);
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
            console.log(`iteration ${String(iter).padStart(2)}  mean R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);
        }
        if (gain < 1e-12) break;
    } else {
        damping *= 4;                 // uphill — shorten the step and retry
        if (damping > 1e8) break;
    }
}

console.log(`\nrefined      mean R = ${(F * 100).toFixed(4)} %   d = [${d.map(v => v.toFixed(2))}]`);

const plotGrid = [];
for (let lam = 400; lam <= 700; lam += 2) plotGrid.push(lam);
const curve = dv => plotGrid.map(lam =>
    tmm(lam, 0, 's', [1, 0], glass(lam), stackAt(lam, dv)).R * 100);

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
