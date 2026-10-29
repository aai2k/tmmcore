/**
 * Which layer controls which part of the spectrum.
 *
 * jacobianSpectrum returns dR/dd for every layer at every wavelength, in the
 * same call that produces R, from the same characteristic-matrix product. That
 * is a sensitivity map: one row per layer, one column per wavelength, signed.
 *
 * The stack is the quarter-wave reflector from example 05. Deep inside the
 * high-reflectance zone R is pinned near 1, so no layer can move it and the map
 * goes neutral; all the sensitivity sits at the band edges and in the sidelobes.
 *
 * The script also checks the analytic derivative against a central difference.
 *
 *   node examples/07-jacobian-map.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rta, jacobian, jacobianSpectrum } from '../src/js/index.js';
import { heatMap } from './_heatmap.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const nH = 2.35;
const nL = 1.46;

const LAMBDA0 = 700;
const PAIRS = 8;
const dH = LAMBDA0 / (4 * nH);
const dL = LAMBDA0 / (4 * nL);

const layers = [];
for (let i = 0; i < PAIRS; i++) {
    layers.push({ n: nH, d: dH });
    layers.push({ n: nL, d: dL });
}
layers.push({ n: nH, d: dH });
const stack = { incident: 1, substrate: 1.52, layers };

const lambdas = [];
for (let lam = 400; lam <= 1000; lam += 5) lambdas.push(lam);

// Z[layer][wavelength]. Layer 1 faces the incident medium; the heat map draws
// row 0 at the bottom, so the stack reads top-down as depth increases upward.
// The grid's derivatives are wavelength-major: wavelength i, layer j at [i * N + j].
const N = layers.length;
const Z = Array.from({ length: N }, () => []);
let peak = 0;
const sensitivity = jacobianSpectrum(stack, lambdas, { pol: 's' }).dRdd;
for (let i = 0; i < lambdas.length; i++) {
    for (let j = 0; j < N; j++) {
        const value = sensitivity[i * N + j];
        Z[j].push(value);
        if (Math.abs(value) > peak) peak = Math.abs(value);
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
const { dRdd } = jacobian(stack, LAM_CHK, { pol: 's' });
const bump = (j, h) => {
    const bumped = layers.map((l, k) => ({ n: l.n, d: k === j ? l.d + h : l.d }));
    return rta({ ...stack, layers: bumped }, LAM_CHK, { pol: 's' }).R;
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
    // way along the bar, so the five ticks stay evenly spaced.
    zTransform: v => 0.5 + 0.5 * Math.sign(v) * Math.cbrt(Math.abs(v) / peak),
    zTickValues: [-peak, -peak / 8, 0, peak / 8, peak],
    zFormat: v => (v === 0 ? '0' : v.toExponential(1)),
    cellW: 3, cellH: 12,
}));
console.log('\nwrote docs/img/jacobian-map.svg');
