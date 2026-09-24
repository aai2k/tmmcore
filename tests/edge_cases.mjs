/**
 * Inputs at the edge of the physics, where a formula that is right everywhere
 * else divides zero by zero or cancels to nothing.
 *
 *  1. Grazing incidence. R and T of a bare interface against Fresnel's
 *     equations with an exact cosine, up to and at 90°; a layer of the incident
 *     index changes nothing there; tmm() and the phase kernel agree.
 *  2. The exact critical angle, inside a layer and at the substrate. The layer
 *     matrix is 0/0 as written and finite in the limit, and R and T match that
 *     limit in closed form and the neighbouring angles on both sides.
 *  3. Invalid thickness. A layer whose thickness is negative or NaN is absent
 *     from every evaluator and has zero derivative in every derivative kernel,
 *     so the thickness Jacobian and the phase Jacobian agree about it.
 *  4. The p sign convention: r_p = r_s at normal incidence.
 *  5. The WebAssembly kernel on all of the above.
 *
 *   node tests/edge_cases.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import {
    tmm, tmmThicknessJacobian, tmmThicknessHessian, tmmNeedleScan, snellCosTheta,
} from '../src/tmm.js';
import { tmmPhaseDispersion, tmmPhaseThicknessJacobian } from '../src/phase.js';
import { jetConstant } from '../src/taylorJet.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));

const near = (actual, expected, tolerance, message) => assert.ok(
    Math.abs(actual - expected) <= tolerance,
    `${message}: got ${actual}, expected ${expected}, tolerance ${tolerance}`);
const relNear = (actual, expected, tolerance, message) => assert.ok(
    Math.abs(actual - expected) <= tolerance * Math.abs(expected),
    `${message}: got ${actual}, expected ${expected}, relative tolerance ${tolerance}`);

// Every combination of one entry from each list.
const grid = (...lists) => lists.reduce(
    (combos, list) => combos.flatMap(combo => list.map(item => [...combo, item])), [[]]);

const AIR = [1, 0], GLASS = [1.52, 0];
const POLS = ['s', 'p'];
const asJets = layers => layers.map(({ n, d }) => ({ nJet: jetConstant(n[0], n[1]), d }));

// ── 1. Grazing incidence ─────────────────────────────────────────────────────
// Fresnel for a bare interface with real indices: T = 4 η0 ηs / (η0 + ηs)²,
// η = n cosθ (s) or n / cosθ (p), cosθ0 taken directly and
// ns cosθs = sqrt(ns² − n0² sin²θ0), which has no cancellation here.

const ADMITTANCE = { s: (n, c) => n * c, p: (n, c) => n / c };

function fresnelT(thetaDeg, pol, n0, ns) {
    const rad = thetaDeg * Math.PI / 180;
    const cs = Math.sqrt(ns * ns - (n0 * Math.sin(rad)) ** 2) / ns;
    const eta0 = ADMITTANCE[pol](n0, Math.cos(rad));
    const etaS = ADMITTANCE[pol](ns, cs);
    return 4 * eta0 * etaS / (eta0 + etaS) ** 2;
}

const GRAZING = [89.9, 89.99999, 89.9999999, 90];
const AIR_LAYER = [{ n: AIR, d: 50 }];
for (const [theta, pol] of grid(GRAZING, POLS)) {
    const at = `${theta}° ${pol}`;
    const bare = tmm(550, theta, pol, AIR, GLASS, []);
    const exact = fresnelT(theta, pol, 1, 1.52);
    relNear(bare.T, exact, 1e-12, `${at}: T of a bare interface`);
    near(bare.R, 1 - exact, 1e-15, `${at}: R of a bare interface`);
    // An air layer under air is a phase delay and nothing else.
    const withAir = tmm(550, theta, pol, AIR, GLASS, AIR_LAYER);
    relNear(withAir.T, bare.T, 1e-12, `${at}: T with a layer of the incident index`);
    near(withAir.R, bare.R, 1e-15, `${at}: R with a layer of the incident index`);
    // The phase kernel carries the same geometry.
    const phase = tmmPhaseDispersion(550, theta, pol, jetConstant(1), jetConstant(1.52), asJets(AIR_LAYER));
    near(phase.r.magnitudeSquared, withAir.R, 1e-15, `${at}: |r|² of the phase kernel`);
}

// ── 2. The exact critical angle ──────────────────────────────────────────────
// From glass of index 2 into a medium of index 1 the critical angle is 30°.
// The doubles near 30° whose sine is exactly 1/2 put cosθ = 0 exactly in the
// medium of index 1; they are searched for rather than assumed, and the
// doubles either side of them give cosθ tiny but nonzero.

const DENSE = [2, 0], RARE = [1, 0];
const cosInRare = theta => {
    const rad = theta * Math.PI / 180;
    return snellCosTheta(DENSE, [Math.sin(rad), 0], RARE, [Math.cos(rad), 0]);
};
const ULP = 2 ** -48;   // one unit in the last place of doubles near 30
const exactZeros = Array.from({ length: 2001 }, (_, i) => 30 + (i - 1000) * ULP)
    .filter(theta => cosInRare(theta).every(part => part === 0));
assert.ok(exactZeros.length > 0, 'a double near 30° puts cosθ exactly at zero in the rarer medium');
const critical = exactZeros[0];
const below = critical - ULP, above = exactZeros[exactZeros.length - 1] + ULP;
assert.ok(cosInRare(below)[0] > 0, 'just below the critical angle the wave propagates');
assert.ok(cosInRare(above)[1] > 0, 'just above it the wave is evanescent');

// In a layer of the rarer medium between two dense ones the matrix at cosθ = 0
// is I + d·G, and with x = (2π/λ) d and η the admittance of the dense media,
//   s: r = −i x η / (2 − i x η),        η = 2 cosθ0
//   p: r =  i x n² / (2η − i x n²),     η = 2 / cosθ0,   n = 1
// so R = u² / (v² + u²) with (u, v) below, and T = 1 − R, the layer being
// lossless.
const LAMBDA = 550, GAP = 120;
const gapLayer = [{ n: RARE, d: GAP }];
const LIMIT_TERMS = {
    s: (x, c0) => [x * 2 * c0, 2],
    p: (x, c0) => [x, 2 * (2 / c0)],
};
for (const pol of POLS) {
    const [u, v] = LIMIT_TERMS[pol](2 * Math.PI / LAMBDA * GAP, Math.cos(critical * Math.PI / 180));
    const R = u * u / (v * v + u * u);
    const at = `critical angle in a ${GAP} nm layer, ${pol}`;
    const k = tmm(LAMBDA, critical, pol, DENSE, DENSE, gapLayer);
    near(k.R, R, 1e-15, `${at}: R against the limit in closed form`);
    near(k.T, 1 - R, 1e-15, `${at}: T against the limit in closed form`);
}
for (const [pol, theta] of grid(POLS, [below, above])) {
    const at = `critical angle in a ${GAP} nm layer, ${pol}, one ulp away at ${theta}`;
    const k = tmm(LAMBDA, critical, pol, DENSE, DENSE, gapLayer);
    const side = tmm(LAMBDA, theta, pol, DENSE, DENSE, gapLayer);
    near(side.R, k.R, 1e-12, `${at}: R`);
    near(side.T, k.T, 1e-12, `${at}: T`);
}

// The substrate at its critical angle reflects everything, with a layer on
// top or without, and transmits nothing.
for (const [pol, layers] of grid(POLS, [[], [{ n: [1.5, 0], d: 90 }]])) {
    const at = `critical angle at the substrate, ${layers.length} layers, ${pol}`;
    const k = tmm(LAMBDA, critical, pol, DENSE, RARE, layers);
    near(k.R, 1, 1e-15, `${at}: R`);
    assert.equal(k.T, 0, `${at}: T`);
    near(k.A, 0, 1e-15, `${at}: A`);
}

// The derivative kernels are finite there too. Their agreement with finite
// differences is in tests/derivatives_fd.mjs.
const CRITICAL_STACK = [{ n: [1.5, 0], d: 90 }, ...gapLayer, { n: [1.5, 0], d: 70 }];
for (const pol of POLS) {
    const args = [LAMBDA, critical, pol, DENSE, DENSE, CRITICAL_STACK];
    const jac = tmmThicknessJacobian(...args);
    const hess = tmmThicknessHessian(...args);
    const scan = tmmNeedleScan(...args, [RARE, [1.5, 0]], [0.5]);
    const needles = [...scan.gaps.flat(), ...scan.intra.flat().flatMap(({ perCand }) => perCand)];
    const values = [...jac.dRdd, ...jac.dTdd, ...hess.d2Rdd.flat(), ...hess.d2Tdd.flat(),
        ...needles.flatMap(g => [g.dR, g.dT])];
    assert.ok(values.every(Number.isFinite), `critical angle, ${pol}: every derivative is finite`);
}

// ── 3. Invalid thickness ─────────────────────────────────────────────────────

const H = [2.35, 0], L = [1.46, 0];
const VALID = [{ n: H, d: 110 }, { n: H, d: 88 }];
const withThickness = d => [VALID[0], { n: L, d }, VALID[1]];
for (const [bad, theta, pol] of grid([NaN, -25], [0, 40], POLS)) {
    const at = `thickness ${bad}, ${theta}° ${pol}`;
    const layers = withThickness(bad);
    const plain = tmm(550, theta, pol, AIR, GLASS, VALID);
    assert.deepEqual(tmm(550, theta, pol, AIR, GLASS, layers), plain, `${at}: R, T, A as without the layer`);

    const jac = tmmThicknessJacobian(550, theta, pol, AIR, GLASS, layers);
    const phase = tmmPhaseThicknessJacobian(550, theta, pol, jetConstant(1), jetConstant(1.52),
        asJets(layers)).r;
    near(jac.R, plain.R, 1e-15, `${at}: R of the Jacobian`);
    const atLayer = [jac.dRdd[1], jac.dTdd[1], jac.dAdd[1], phase.dLogMagnitudeSquared[1]];
    assert.ok(atLayer.every(value => value === 0), `${at}: derivatives at the layer ${atLayer}`);
    // dR/dd = R d(ln|r|²)/dd: the two Jacobians agree layer by layer.
    relNear(jac.dRdd[0], jac.R * phase.dLogMagnitudeSquared[0], 1e-9, `${at}: dR/dd[0], both Jacobians`);
    relNear(jac.dRdd[2], jac.R * phase.dLogMagnitudeSquared[2], 1e-9, `${at}: dR/dd[2], both Jacobians`);
    const scan = tmmNeedleScan(550, theta, pol, AIR, GLASS, layers, [H, L], [0.5]);
    assert.ok(scan.gaps.flat().every(g => Number.isFinite(g.dR)), `${at}: needle gaps are finite`);
}

// ── 4. The p sign convention ─────────────────────────────────────────────────
// r_p follows the tilted admittances n/cosθ (Macleod 5th ed., §2.4), so at
// normal incidence it is r_s itself, the opposite sign to the Fresnel r_p of
// Born and Wolf. benchmarks/compare.mjs checks the sign against Byrnes' tmm at
// oblique incidence.

{
    const layers = asJets([{ n: H, d: 110 }, { n: [0.15, 3.2], d: 20 }, { n: L, d: 95 }]);
    const [s, p] = POLS.map(pol =>
        tmmPhaseDispersion(550, 0, pol, jetConstant(1), jetConstant(1.52), layers).r);
    assert.equal(p.phaseRad, s.phaseRad, 'phase of r_p equals phase of r_s at normal incidence');
    assert.equal(p.magnitudeSquared, s.magnitudeSquared, '|r_p| equals |r_s| at normal incidence');
}

// ── 5. The WebAssembly kernel ────────────────────────────────────────────────

const CASES = [
    ...GRAZING.map(theta => [`${theta}°`, theta, AIR, GLASS, AIR_LAYER]),
    ['critical layer', critical, DENSE, DENSE, CRITICAL_STACK],
    ['critical substrate', critical, DENSE, RARE, [{ n: [1.5, 0], d: 90 }]],
    ['NaN layer', 40, AIR, GLASS, withThickness(NaN)],
    ['negative layer', 40, AIR, GLASS, withThickness(-25)],
];
const CANDIDATES = [RARE, [1.5, 0]];

// [label, JavaScript value, WebAssembly value] for every output of one case.
function portPairs(kernel, [name, theta, n0, ns, layers], pol) {
    const code = POLS.indexOf(pol);
    const at = `${name} ${pol}`;
    const js = [LAMBDA, theta, pol, n0, ns, layers];
    const wa = [LAMBDA, theta, code, n0, ns, layers];
    const a = tmm(...js), b = kernel.tmmOne(...wa);
    const ja = tmmThicknessJacobian(...js), jb = kernel.tmmJacobian(...wa);
    const ha = tmmThicknessHessian(...js), hb = kernel.tmmHessian(...wa);
    const na = tmmNeedleScan(...js, CANDIDATES, [0.5]), nb = kernel.tmmNeedleScan(...wa, CANDIDATES, [0.5]);
    const N = layers.length;
    const needle = (scan, pos, c, q) => scan.gaps[pos][c][q];
    const inside = (scan, k, c, q) => scan.intra[k][0].perCand[c][q];
    return [
        ...['R', 'T', 'A'].map(q => [`${at} ${q}`, a[q], b[q]]),
        ...grid(['dRdd', 'dTdd', 'dAdd'], [...Array(N).keys()]).map(([q, i]) =>
            [`${at} ${q}[${i}]`, ja[q][i], jb[q][i]]),
        ...grid(['d2Rdd', 'd2Tdd', 'd2Add'], [...Array(N).keys()], [...Array(N).keys()]).map(([q, i, j]) =>
            [`${at} ${q}[${i}][${j}]`, ha[q][i][j], hb[q][i][j]]),
        ...grid([...Array(N + 1).keys()], [0, 1], ['dR', 'dT', 'dA']).map(([pos, c, q]) =>
            [`${at} gap[${pos}][${c}] ${q}`, needle(na, pos, c, q), needle(nb, pos, c, q)]),
        ...grid([...Array(N).keys()], [0, 1], ['dR', 'dT', 'dA']).map(([k, c, q]) =>
            [`${at} intra[${k}][${c}] ${q}`, inside(na, k, c, q), inside(nb, k, c, q)]),
    ];
}

function phasePairs(kernel, theta, pol) {
    const args = [LAMBDA, theta, pol, jetConstant(1), jetConstant(1.52), asJets(AIR_LAYER)];
    const a = tmmPhaseDispersion(...args).r;
    const b = kernel.tmmPhaseOne(...args.slice(0, 2), POLS.indexOf(pol), ...args.slice(3)).r;
    return ['magnitudeSquared', 'phaseRad'].map(q => [`phase ${theta}° ${pol} ${q}`, a[q], b[q]]);
}

const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');
const kernel = existsSync(WASM) ? await instantiateTmmWasm(readFileSync(WASM)) : null;
const pairs = kernel ? [
    ...grid(CASES, POLS).flatMap(([c, pol]) => portPairs(kernel, c, pol)),
    ...grid(GRAZING, POLS).flatMap(([theta, pol]) => phasePairs(kernel, theta, pol)),
] : [];
for (const [label, js, wa] of pairs) {
    const diff = Math.abs(js - wa);
    assert.ok(diff <= 1e-12 || diff <= 1e-7 * Math.abs(js), `wasm ${label}: js=${js} wasm=${wa}`);
}
if (!kernel) console.log('NOTE : src/tmm_kernel.wasm not built; the WebAssembly checks were skipped.');

console.log(`exact critical angle found at θ0 = ${critical}°`);
console.log(`WebAssembly comparisons: ${pairs.length}`);
console.log('PASS : edge cases');
