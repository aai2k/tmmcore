/**
 * Which layer controls which part of the spectrum.
 *
 * tmmThicknessJacobian returns dR/dd for every layer in the same call that
 * produces R, from the same characteristic-matrix product. Sweeping that over
 * wavelength gives a sensitivity map: one row per layer, one column per
 * wavelength, signed.
 *
 * The stack is the quarter-wave reflector from example 05. Deep inside the
 * high-reflectance zone R is pinned near 1, so no layer can move it and the map
 * goes neutral; all the leverage sits at the band edges and in the sidelobes.
 *
 * The script also checks the analytic derivative against a central difference,
 * because a derivative nobody verified is just an assertion.
 *
 *   node examples/07-jacobian-map.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm, tmmThicknessJacobian } from '../src/index.js';
import { heatMap } from './_heatmap.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const nH = [2.35, 0];
const nL = [1.46, 0];
const AIR = [1, 0];
const GLASS = [1.52, 0];

const LAMBDA0 = 700;
const PAIRS = 8;
const dH = LAMBDA0 / (4 * nH[0]);
const dL = LAMBDA0 / (4 * nL[0]);

const stack = [];
for (let i = 0; i < PAIRS; i++) {
    stack.push({ n: nH, d: dH });
    stack.push({ n: nL, d: dL });
}
stack.push({ n: nH, d: dH });

const lambdas = [];
for (let lam = 400; lam <= 1000; lam += 5) lambdas.push(lam);

// Z[layer][wavelength]. Layer 1 faces the incident medium; the heat map draws
// row 0 at the bottom, so the stack reads top-down as depth increases upward.
const N = stack.length;
const Z = Array.from({ length: N }, () => []);
let peak = 0;
for (const lam of lambdas) {
    const { dRdd } = tmmThicknessJacobian(lam, 0, 's', AIR, GLASS, stack);
    for (let j = 0; j < N; j++) {
        Z[j].push(dRdd[j]);
        if (Math.abs(dRdd[j]) > peak) peak = Math.abs(dRdd[j]);
    }
}

console.log(`quarter-wave stack  (HL)^${PAIRS} H   lambda0 = ${LAMBDA0} nm   ${N} layers`);
console.log(`largest |dR/dd| over the grid: ${peak.toExponential(4)} /nm\n`);

// Where the stack is most and least sensitive, as numbers rather than a colour.
const rowPeak = Z.map(row => Math.max(...row.map(Math.abs)));
const bandIdx = lambdas.map((l, i) => (l >= 620 && l <= 800 ? i : -1)).filter(i => i >= 0);
const inBand = Math.max(...Z.flatMap(row => bandIdx.map(i => Math.abs(row[i]))));
console.log(`peak |dR/dd| inside the 620-800 nm zone : ${inBand.toExponential(4)} /nm`);
console.log(`peak |dR/dd| anywhere on the grid       : ${peak.toExponential(4)} /nm`);
console.log(`ratio                                   : ${(peak / inBand).toFixed(1)}x\n`);

console.log('most sensitive layer, by peak |dR/dd|');
const order = rowPeak.map((v, j) => [j + 1, v]).sort((a, b) => b[1] - a[1]);
for (const [layer, v] of order.slice(0, 3)) {
    console.log(`  layer ${String(layer).padStart(2)}  ${v.toExponential(4)} /nm`);
}

// Analytic against a central difference, at a wavelength on the band edge where
// the derivative is large.
const LAM_CHK = 615;
const { dRdd } = tmmThicknessJacobian(LAM_CHK, 0, 's', AIR, GLASS, stack);
const bump = (j, h) => {
    const s = stack.map((l, k) => ({ n: l.n, d: k === j ? l.d + h : l.d }));
    return tmm(LAM_CHK, 0, 's', AIR, GLASS, s).R;
};
console.log(`\nanalytic dR/dd against a central difference at ${LAM_CHK} nm, h = 1e-4 nm`);
console.log('  layer    analytic          central diff      rel. difference');
for (const j of [0, 4, 8, 16]) {
    const fd = (bump(j, 1e-4) - bump(j, -1e-4)) / 2e-4;
    const rel = Math.abs(fd - dRdd[j]) / Math.abs(dRdd[j]);
    console.log(`  ${String(j + 1).padStart(5)}    ${dRdd[j].toExponential(6)}    ${fd.toExponential(6)}    ${rel.toExponential(1)}`);
}

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'jacobian-map.svg'), heatMap({
    title: `Thickness sensitivity dR/dd of a (HL)⁸H stack, λ₀ = ${LAMBDA0} nm`,
    panels: [{ title: 'normal incidence, s-polarization', z: Z }],
    x: lambdas,
    y: Array.from({ length: N }, (_, j) => j + 1),
    xLabel: 'wavelength (nm)',
    yLabel: 'layer (1 = facing air)',
    zLabel: 'dR/dd',
    zMin: -peak, zMax: peak,
    ramp: 'diverging',
    // Two narrow band edges carry derivatives 150 times larger than anything
    // else, so a linear scale renders 87% of the map as blank neutral. The
    // signed cube root keeps the sign, compresses those peaks and opens up the
    // sidelobe structure. Cube root also puts +-peak/8 exactly a quarter of the
    // way along the bar, so the five ticks stay evenly spaced and honest.
    zTransform: v => 0.5 + 0.5 * Math.sign(v) * Math.cbrt(Math.abs(v) / peak),
    zTickValues: [-peak, -peak / 8, 0, peak / 8, peak],
    zFormat: v => (v === 0 ? '0' : v.toExponential(1)),
    cellW: 3, cellH: 12,
}));
console.log('\nwrote docs/img/jacobian-map.svg');
