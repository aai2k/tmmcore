/**
 * Reflectance over wavelength and angle of incidence.
 *
 * A quarter-wave stack tilted away from normal incidence: its high-reflectance
 * zone moves to shorter wavelengths, and s and p separate as the angle grows.
 *
 * The mechanism is the phase thickness at oblique incidence, delta =
 * 2*pi*n*d*cos(theta)/lambda, which reads as an apparent optical thickness
 * n*d*cos(theta), so the layers behave as though they were thinner when tilted
 * (Macleod, Thin-Film Optical Filters 5th ed., ch. 8, "Simple Tilts in
 * Collimated Light"). For the quarter-wave stack specifically: "As the
 * multilayer is tilted to greater angles of incidence, the characteristic moves
 * to a shorter wavelength" (ch. 10, One-Dimensional Photonic Crystals).
 *
 * Indices are non-dispersive here so the shift is the only effect in view.
 *
 *   node examples/05-angle-map.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/index.js';
import { heatMap } from './_heatmap.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const nH = [2.35, 0];          // high index, TiO2-like
const nL = [1.46, 0];          // low index, SiO2-like
const AIR = [1, 0];
const GLASS = [1.52, 0];

const LAMBDA0 = 700;           // nm, reference wavelength of the quarter waves
const PAIRS = 8;               // (HL)^8 H  ->  17 layers

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
const angles = [];
for (let a = 0; a <= 80; a += 1) angles.push(a);

// Point-sampled R, on the grid exactly. These are the numbers reported below,
// so anything printed here is directly reproducible in another tool.
// z[angle][wavelength] is the row-major order heatMap expects.
const Zs = [], Zp = [];
for (const theta of angles) {
    const rs = [], rp = [];
    for (const lam of lambdas) {
        rs.push(tmm(lam, theta, 's', AIR, GLASS, stack).R);
        rp.push(tmm(lam, theta, 'p', AIR, GLASS, stack).R);
    }
    Zs.push(rs); Zp.push(rp);
}

// For display only. The sidelobe fringes on the short-wavelength side are
// narrower than one cell, so point sampling aliases them into speckle. Each
// cell of the figure instead averages R over its own wavelength bin, which is
// what a spectrophotometer of finite bandwidth measures.
const SUB = 4;
const STEP = lambdas[1] - lambdas[0];
function binned(pol) {
    return angles.map(theta => lambdas.map(lam => {
        let acc = 0;
        for (let k = 0; k < SUB; k++) {
            acc += tmm(lam + ((k + 0.5) / SUB - 0.5) * STEP, theta, pol, AIR, GLASS, stack).R;
        }
        return acc / SUB;
    }));
}

// Band edges, defined as where R crosses 0.9, reported rather than asserted:
// "high-reflectance zone" is a judgement about a number.
function band(row) {
    const hits = row.map((r, i) => (r >= 0.9 ? lambdas[i] : null)).filter(v => v !== null);
    return hits.length ? [hits[0], hits[hits.length - 1]] : null;
}

console.log(`quarter-wave stack  (HL)^${PAIRS} H   lambda0 = ${LAMBDA0} nm`);
console.log(`  H  n = ${nH[0]}  d = ${dH.toFixed(2)} nm`);
console.log(`  L  n = ${nL[0]}  d = ${dL.toFixed(2)} nm\n`);
console.log('R >= 0.9 band, nm');
console.log('  AOI      s-pol           p-pol');
for (const theta of [0, 30, 45, 60]) {
    const i = angles.indexOf(theta);
    const bs = band(Zs[i]), bp = band(Zp[i]);
    const f = b => (b ? `${b[0]}-${b[1]}`.padEnd(14) : 'none'.padEnd(14));
    console.log(`  ${String(theta).padStart(2)}°   ${f(bs)}  ${f(bp)}`);
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
    const { R } = tmm(lam, theta, pol, AIR, GLASS, stack);
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
