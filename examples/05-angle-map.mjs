/**
 * Reflectance over wavelength and angle of incidence.
 *
 * A quarter-wave stack tilted away from normal incidence: its high-reflectance
 * zone moves to shorter wavelengths, and s and p separate as the angle grows.
 *
 * The mechanism is the phase thickness at oblique incidence, delta =
 * 2*pi*n*d*cos(theta)/lambda, which reads as an apparent optical thickness
 * n*d*cos(theta), so the layers behave as though they were thinner when tilted
 * (Macleod, Thin-Film Optical Filters 5th ed., §8.2.5.1, "Simple Tilts in
 * Collimated Light").
 *
 * Indices are non-dispersive here so the shift is the only effect in view.
 *
 *   node examples/05-angle-map.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rta, rtaSpectrum } from '../src/js/index.js';
import { heatMap } from './_heatmap.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const nH = 2.35;               // high index, TiO2-like
const nL = 1.46;               // low index, SiO2-like

const LAMBDA0 = 700;           // nm, reference wavelength of the quarter waves
const PAIRS = 8;               // (HL)^8 H  ->  17 layers

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
const angles = [];
for (let a = 0; a <= 80; a += 1) angles.push(a);

// Point-sampled R, on the grid exactly. These are the numbers reported below,
// so anything printed here is directly reproducible in another tool.
// z[angle][wavelength] is the row-major order heatMap expects.
const Zs = [], Zp = [];
for (const theta of angles) {
    const { Rs, Rp } = rtaSpectrum(stack, lambdas, { angle: theta });
    Zs.push(Array.from(Rs)); Zp.push(Array.from(Rp));
}

// For display only. The sidelobe fringes on the short-wavelength side are
// narrower than one cell, so point sampling aliases them into speckle. Each
// cell of the figure instead averages R over its own wavelength bin, which is
// what a spectrophotometer of finite bandwidth measures.
const SUB = 4;
const STEP = lambdas[1] - lambdas[0];
function binned(pol) {
    return angles.map(theta => {
        const acc = new Array(lambdas.length).fill(0);
        for (let k = 0; k < SUB; k++) {
            const shifted = lambdas.map(lam => lam + ((k + 0.5) / SUB - 0.5) * STEP);
            rtaSpectrum(stack, shifted, { angle: theta, pol }).R.forEach((R, i) => { acc[i] += R; });
        }
        return acc.map(sum => sum / SUB);
    });
}

// Band edges: where R crosses 0.9.
function band(row) {
    const hits = row.map((r, i) => (r >= 0.9 ? lambdas[i] : null)).filter(v => v !== null);
    return hits.length ? [hits[0], hits[hits.length - 1]] : null;
}

console.log(`quarter-wave stack  (HL)^${PAIRS} H   lambda0 = ${LAMBDA0} nm`);
console.log(`  H  n = ${nH}  d = ${dH.toFixed(2)} nm`);
console.log(`  L  n = ${nL}  d = ${dL.toFixed(2)} nm\n`);
console.log('R >= 0.9 band, nm');
console.log('  AOI      s-pol           p-pol');
for (const theta of [0, 30, 45, 60]) {
    const i = angles.indexOf(theta);
    const bs = band(Zs[i]), bp = band(Zp[i]);
    const f = b => (b ? `${b[0]}-${b[1]}` : 'none');
    console.log(`  ${String(theta).padStart(2)}°   ${f(bs).padEnd(14)}  ${f(bp)}`);
}

const c0 = band(Zs[0]);
const c60 = band(Zs[angles.indexOf(60)]);
if (c0 && c60) {
    const mid = b => (b[0] + b[1]) / 2;
    console.log(`\ns-pol band centre moves ${(mid(c0) - mid(c60)).toFixed(0)} nm to the blue ` +
                `between 0° and 60° (${mid(c0)} -> ${mid(c60)} nm).`);
}

// Single points, at full precision, for checking this design against another
// implementation without depending on where a band edge lands on the grid.
console.log('\nreference points');
for (const [lam, theta, pol] of [
    [700, 0, 's'], [550, 0, 's'], [450, 0, 's'],
    [700, 45, 's'], [700, 45, 'p'],
    [600, 60, 's'], [600, 60, 'p'],
]) {
    const { R } = rta(stack, lam, { angle: theta, pol });
    console.log(`  ${String(lam).padStart(4)} nm  ${String(theta).padStart(2)}°  ${pol}   R = ${R.toFixed(12)}`);
}

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'angle-map.svg'), heatMap({
    title: `Reflectance of a (HL)⁸H quarter-wave stack, λ₀ = ${LAMBDA0} nm`,
    panels: [
        { title: 's-polarization', z: binned('s') },
        { title: 'p-polarization', z: binned('p') },
    ],
    x: lambdas, y: angles,
    xLabel: 'wavelength (nm)',
    yLabel: 'angle of incidence (degrees)',
    zLabel: 'R',
    cellW: 3.0, cellH: 3.0,
}));
console.log('wrote docs/img/angle-map.svg');
