/**
 * A coated plate as a spectrophotometer measures it.
 *
 * The antireflection coating of examples/02-ar-coating.mjs on the front face of
 * a glass plate, with the back face bare. The plate is thick: the trip to the
 * back face and back again is far longer than the coherence length of the
 * measuring beam, so the light from the two faces adds in irradiance, not in
 * amplitude, and the plate is not one more layer of the stack.
 *
 * Instead each face is computed on its own and the beams bouncing between them
 * are summed. With a the coated front face and b the bare back face, + for
 * light going into the plate and - for light coming back out (Macleod,
 * Thin-Film Optical Filters 5th ed., §2.14 "Incoherent Reflection at Two or More
 * Surfaces", Figure 2.12):
 *
 *     R = Ra+ + Ta Rb+ Ta / (1 - Ra- Rb+)
 *
 * Ta is the same in both directions. Ra- is the coating seen from inside the
 * glass: the same layers in reverse order, with the glass as incident medium
 * and air as substrate. This form has no loss in the glass between the faces,
 * so the plate's internal transmittance is 1 here, as it is for a glass of k = 0.
 *
 * The check is the uncoated plate. Both faces then reflect R0 = ((1 - ns)/(1 + ns))²
 * (Eq. 2.53), and Macleod's lossless form of the sum in §2.14,
 * R = (Ra + Rb - 2 Ra Rb) / (1 - Ra Rb), becomes R = 2 R0 / (1 + R0).
 *
 * Writes the plot to docs/img/coated-plate.svg.
 *
 *   node examples/03-coated-plate.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rtaSpectrum } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

// The materials and design of examples/02-ar-coating.mjs, n(lambda) with lambda in nm.
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

const air = 1;
const plate = { n: lambdas.map(glass) };
const layers = design.map(({ mat, d }) => ({ n: lambdas.map(mat), d }));

// Normal incidence, so the unpolarized default is the same as s or p.
const front = rtaSpectrum({ incident: air, substrate: plate, layers }, lambdas);
const frontFromGlass = rtaSpectrum({ incident: plate, substrate: air, layers: [...layers].reverse() }, lambdas);
const back = rtaSpectrum({ incident: plate, substrate: air, layers: [] }, lambdas);
const bareFront = rtaSpectrum({ incident: air, substrate: plate, layers: [] }, lambdas);

// The sum over the beams reflected back and forth inside the plate, Macleod §2.14.
const plateR = (RaIn, Ta, RaOut, Rb) => RaIn + (Ta * Rb * Ta) / (1 - RaOut * Rb);

const coating = Array.from(front.R);
const coated = lambdas.map((_, i) => plateR(front.R[i], front.T[i], frontFromGlass.R[i], back.R[i]));
const uncoated = lambdas.map((_, i) => plateR(bareFront.R[i], bareFront.T[i], back.R[i], back.R[i]));

// Without loss in the layers the coating reflects the same from either side.
const sides = Math.max(...lambdas.map((_, i) => Math.abs(frontFromGlass.R[i] - front.R[i])));

console.log('coating of example 02 on the front of a thick plate, back face bare, normal incidence');
console.log(`coating from the glass side against the air side: largest |Ra- - Ra+| = ${sides.toExponential(1)}\n`);

const at = lam => lambdas.indexOf(lam);
const col = (s, width) => String(s).padStart(width);
console.log('R, %    ' + col('coating alone', 16) + col('coated plate', 16) + col('uncoated plate', 16));
for (const lam of [400, 450, 510, 600, 700]) {
    const i = at(lam);
    const pct = v => col((v * 100).toFixed(3), 16);
    console.log(`  ${lam} nm` + pct(coating[i]) + pct(coated[i]) + pct(uncoated[i]));
}

console.log('\nuncoated plate against the closed form 2 R0 / (1 + R0)');
console.log('        ' + col('incoherent sum', 23) + col('closed form', 23) + col('difference', 12));
for (const lam of [400, 550, 700]) {
    const ns = glass(lam);
    const R0 = ((1 - ns) / (1 + ns)) ** 2;
    const exact = (2 * R0) / (1 + R0);
    const sum = uncoated[at(lam)];
    console.log(`  ${lam} nm` + col(sum.toPrecision(17), 23) + col(exact.toPrecision(17), 23) +
        col(Math.abs(sum - exact).toExponential(1), 12));
}

const percent = values => values.map(R => R * 100);
mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'coated-plate.svg'), linePlot({
    title: 'Coated plate with a bare back face, normal incidence',
    xLabel: 'wavelength (nm)',
    yLabel: 'reflectance (%)',
    // Headroom above the uncoated plate's 8 % keeps its curve clear of the legend.
    yMin: 0,
    yMax: 12,
    series: [
        { label: 'uncoated plate', x: lambdas, y: percent(uncoated), dash: '5 4' },
        { label: 'coating alone', x: lambdas, y: percent(coating) },
        { label: 'coated plate', x: lambdas, y: percent(coated) },
    ],
}));
console.log('\nwrote docs/img/coated-plate.svg');
