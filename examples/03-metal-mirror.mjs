/**
 * A silver mirror at 45° — absorbing media and the s/p split.
 *
 * Metals are where transfer-matrix implementations are most likely to diverge:
 * k is large, the phase thickness is strongly complex, and the field decays
 * within a few tens of nanometres. Handling them correctly is the interesting
 * part of a TMM, not an afterthought.
 *
 * At oblique incidence s and p behave differently, and the absorptance is where
 * that shows most clearly.
 *
 * Silver indices here are a linear approximation over the visible, adequate for
 * a demonstration; use measured n,k for real work.
 *
 *   node examples/03-metal-mirror.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const silver = lam => [0.15 + 0.0006 * lam, 3.2 + 0.004 * lam];
const glass = lam => [1.504 + 4200 / (lam * lam), 0];

const THETA = 45;
const THICKNESS = 120;   // nm — optically opaque at these wavelengths

const lambdas = [];
for (let lam = 400; lam <= 800; lam += 2) lambdas.push(lam);

const Rs = [], Rp = [], As = [], Ap = [];
for (const lam of lambdas) {
    const layers = [{ n: silver(lam), d: THICKNESS }];
    const s = tmm(lam, THETA, 's', [1, 0], glass(lam), layers);
    const p = tmm(lam, THETA, 'p', [1, 0], glass(lam), layers);
    Rs.push(s.R * 100); Rp.push(p.R * 100);
    As.push(s.A * 100); Ap.push(p.A * 100);
}

const at = lam => lambdas.indexOf(lam);
console.log('           Rs %     Rp %     As %     Ap %');
for (const lam of [400, 550, 700]) {
    const i = at(lam);
    console.log(`${lam} nm   ${Rs[i].toFixed(3)}   ${Rp[i].toFixed(3)}   ` +
                `${As[i].toFixed(3)}   ${Ap[i].toFixed(3)}`);
}

// At this thickness silver is near-opaque: whatever is not reflected is
// absorbed, and the leftover 1 − R − A is the light still reaching the
// substrate. Worth printing rather than assuming — "opaque" is a judgement
// about a number, and how small that number is depends on the design.
const maxT = Math.max(...lambdas.map((_, i) => 1 - (Rs[i] + As[i]) / 100));
console.log(`\nmax transmittance through ${THICKNESS} nm of silver: ${maxT.toExponential(2)}`);

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'metal-mirror.svg'), linePlot({
    title: `Opaque silver mirror at ${THETA}°`,
    xLabel: 'wavelength (nm)',
    yLabel: 'percent',
    series: [
        { label: 'R (s-pol)', x: lambdas, y: Rs },
        { label: 'R (p-pol)', x: lambdas, y: Rp },
        { label: 'A (s-pol)', x: lambdas, y: As, dash: '5 4' },
        { label: 'A (p-pol)', x: lambdas, y: Ap, dash: '5 4' },
    ],
}));
console.log('wrote docs/img/metal-mirror.svg');
