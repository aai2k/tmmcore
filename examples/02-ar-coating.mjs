/**
 * Broadband antireflection coating: a spectrum, and a figure.
 *
 * The quarter-half-quarter design: a quarter-wave of low index facing
 * air, a half-wave of high index, and a quarter-wave of intermediate index
 * against the substrate. The half-wave layer is an absentee at the reference
 * wavelength and flattens the response either side of it.
 *
 * Macleod, Thin-Film Optical Filters 5th ed., ch. 4, "Multilayer Antireflection
 * Coatings".
 *
 * Writes the plot to docs/img/ar-coating.svg.
 *
 *   node examples/02-ar-coating.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

// Dispersion: n(lambda) with lambda in nm. Cauchy-style, chosen to be
// representative rather than exact; in practice you would supply measured n,k.
const low = lam => 1.377 + 2000 / (lam * lam);      // MgF2-like
const high = lam => 2.170 + 8000 / (lam * lam);     // ZrO2/TiO2-like
const mid = lam => 1.680 + 5000 / (lam * lam);      // Al2O3-like
const glass = lam => 1.504 + 4200 / (lam * lam);    // BK7-like

const LAMBDA0 = 510;   // reference wavelength, nm

// Quarter, half, quarter, outermost layer first.
const design = [
    { mat: low, waves: 0.25 },
    { mat: high, waves: 0.5 },
    { mat: mid, waves: 0.25 },
].map(({ mat, waves }) => ({ mat, d: (waves * LAMBDA0) / mat(LAMBDA0) }));

const lambdas = [];
for (let lam = 400; lam <= 700; lam += 2) lambdas.push(lam);

// Each index is a list with one value per wavelength.
const substrate = { n: lambdas.map(glass) };
const layers = design.map(({ mat, d }) => ({ n: lambdas.map(mat), d }));

// R is the mean of s and p, natural light. At normal incidence the two are
// identical, but the mean keeps the code correct if you change the angle.
const percent = values => Array.from(values, R => R * 100);
const coated = percent(rtaSpectrum({ incident: 1, substrate, layers }, lambdas).R);
const bare = percent(rtaSpectrum({ incident: 1, substrate, layers: [] }, lambdas).R);

const mean = a => a.reduce((x, y) => x + y, 0) / a.length;
console.log(`mean R, bare glass   ${mean(bare).toFixed(3)} %`);
console.log(`mean R, coated       ${mean(coated).toFixed(3)} %`);
console.log(`peak R, coated       ${Math.max(...coated).toFixed(3)} %`);

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'ar-coating.svg'), linePlot({
    title: 'Quarter-half-quarter antireflection coating on BK7, normal incidence',
    xLabel: 'wavelength (nm)',
    yLabel: 'reflectance (%)',
    yMin: 0,
    series: [
        { label: 'uncoated', x: lambdas, y: bare, dash: '5 4' },
        { label: 'coated', x: lambdas, y: coated },
    ],
}));
console.log('\nwrote docs/img/ar-coating.svg');
