/**
 * The needle P-function: where should the next layer go?
 *
 * Given a design and a merit function, tmmNeedleScan returns the derivative of
 * each spectral quantity with respect to inserting an infinitesimally thin
 * layer of a candidate material, at every position in the stack at once. Chain
 * that through the merit function and you get P(z), whose most negative point
 * is the best place to insert.
 *
 * This is the d -> 0 limit of Sullivan and Dobrowolski's numerical pre/post
 * method (Appl. Opt. 35, 5484, 1996), the analytic P-function of Tikhonravov
 * et al. (Appl. Opt. 35, 5493, 1996). Unlike the numerical form it needs no
 * trial thickness and no second spectrum evaluation.
 *
 * The script also checks P against a finite difference, which is the point:
 * the analytic value has to be the limit the numerical one approaches.
 *
 *   node examples/06-needle.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm, tmmNeedleScan } from '../src/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const H = [2.35, 0];
const L = [1.46, 0];
const AIR = [1, 0];
const GLASS = [1.52, 0];

// A deliberately mediocre antireflection start, so there is something to find.
const design = [
    { n: L, d: 120 },
    { n: H, d: 60 },
    { n: L, d: 140 },
    { n: H, d: 40 },
];

// Merit function: mean square reflectance across the visible. Targeting zero
// makes dMF/dx = (2/K) * sum R * dR/dx, so the needle metric chains directly.
const band = [];
for (let lam = 450; lam <= 650; lam += 20) band.push(lam);

const meritOf = layers =>
    band.reduce((acc, lam) => acc + tmm(lam, 0, 's', AIR, GLASS, layers).R ** 2, 0) / band.length;

const FRACS = Array.from({ length: 49 }, (_, i) => (i + 1) / 50);
const CANDS = [H, L];

// P(z) for each candidate, sampled inside every layer. Depth is measured from
// the incident medium, which is the order the layers are stored in.
const depthOf = (k, frac) =>
    design.slice(0, k).reduce((a, l) => a + l.d, 0) + frac * design[k].d;

const P = CANDS.map(() => []);
const Z = [];
for (let k = 0; k < design.length; k++) {
    for (const frac of FRACS) Z.push(depthOf(k, frac));
}

for (let c = 0; c < CANDS.length; c++) {
    const acc = new Array(Z.length).fill(0);
    for (const lam of band) {
        const scan = tmmNeedleScan(lam, 0, 's', AIR, GLASS, design, CANDS, FRACS);
        let i = 0;
        for (let k = 0; k < design.length; k++) {
            for (let f = 0; f < FRACS.length; f++) {
                acc[i++] += (2 / band.length) * scan.R * scan.intra[k][f].perCand[c].dR;
            }
        }
    }
    P[c] = acc;
}

const total = design.reduce((a, l) => a + l.d, 0);
const bounds = design.map((_, k) => design.slice(0, k + 1).reduce((a, l) => a + l.d, 0));

console.log(`start design, ${design.length} layers, ${total} nm total`);
console.log(`merit (mean R^2 over ${band[0]}-${band[band.length - 1]} nm) = ${meritOf(design).toExponential(4)}\n`);

// Best insertion point overall.
let best = { p: Infinity };
for (let c = 0; c < CANDS.length; c++) {
    for (let i = 0; i < Z.length; i++) {
        if (P[c][i] < best.p) best = { p: P[c][i], c, i, z: Z[i] };
    }
}
const nameOf = c => (c === 0 ? 'H (n=2.35)' : 'L (n=1.46)');
console.log(`most negative P: ${best.p.toExponential(4)} /nm`);
console.log(`  candidate ${nameOf(best.c)} at depth ${best.z.toFixed(1)} nm from the incident medium\n`);

// The analytic P must be the limit of the numerical slope as the trial
// thickness goes to zero. Insert a needle of thickness eps at the best spot and
// watch the finite difference converge.
function withNeedle(z, nA, eps) {
    const out = [];
    let acc = 0;
    for (const layer of design) {
        if (z > acc && z < acc + layer.d) {
            out.push({ n: layer.n, d: z - acc });
            out.push({ n: nA, d: eps });
            out.push({ n: layer.n, d: acc + layer.d - z });
        } else {
            out.push({ ...layer });
        }
        acc += layer.d;
    }
    return out;
}

const mf0 = meritOf(design);
console.log('analytic P against a finite difference at the same point');
console.log('   eps (nm)   (MF(eps) - MF(0)) / eps        ratio to analytic');
for (const eps of [1, 0.1, 0.01, 0.001]) {
    const fd = (meritOf(withNeedle(best.z, CANDS[best.c], eps)) - mf0) / eps;
    console.log(`   ${String(eps).padEnd(9)}  ${fd.toExponential(6).padEnd(24)}  ${(fd / best.p).toFixed(6)}`);
}

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'needle.svg'), linePlot({
    title: 'Needle P-function through a 4-layer antireflection start',
    xLabel: 'depth from incident medium (nm)',
    yLabel: 'dMF / dd  (1/nm)',
    zeroLine: true,
    vlines: bounds.slice(0, -1).map((x, i) => ({ x, label: `${i + 1}|${i + 2}` })),
    series: [
        { label: 'insert H (n=2.35)', x: Z, y: P[0] },
        { label: 'insert L (n=1.46)', x: Z, y: P[1] },
    ],
}));
console.log('\nwrote docs/img/needle.svg');
