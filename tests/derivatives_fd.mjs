/**
 * Analytic derivatives against finite differences of the plain evaluators.
 *
 * tests/equivalence.mjs holds the C port to the JavaScript, so it cannot catch
 * a formula the two share. This file holds every derivative kernel to finite
 * differences of tmm() and tmmPhaseDispersion(), which share nothing with the
 * kernels but the layer matrix:
 *
 *   • the thickness Jacobian, dR/dd, dT/dd and dA/dd of every layer
 *   • the thickness Hessian, every d²R/dd_i dd_j and d²T/dd_i dd_j
 *   • the needle P-function at every gap and at positions inside each layer
 *   • the phase thickness Jacobian: phase, GD, GDD and ln|r|²
 *
 * The stacks are chosen to reach the two overflow guards. An opaque metal, and
 * a lossless air gap past the critical angle, have a phase thickness whose
 * imaginary part passes the bound layerMatrix holds it to; a deep Ge/ZnS stack
 * seen where Ge absorbs has prefix and suffix products past 1e77, where a
 * complex division by den² forms |den|⁴ and overflows. An air gap exactly at
 * its critical angle has a layer matrix that is 0/0 as written, and an
 * absorbing incident medium gives A a term that dR and dT do not carry. An
 * ordinary dielectric stack runs alongside.
 *
 * With the WebAssembly kernel built, its derivatives on the same stacks,
 * needles inside layers included, are held to the JavaScript as well.
 *
 *   node tests/derivatives_fd.mjs
 */

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

// Thickness step for the differences, nm. The stencils below leave truncation
// errors of order h⁴ times a fifth or sixth thickness derivative, below 1e-12
// for these stacks, and round-off of order ε/h² at worst, about 1e-12.
const STEP_NM = 0.01;
const FIRST = { abs: 1e-9, rel: 1e-6 };
const SECOND = { abs: 1e-9, rel: 1e-6 };

// GDD is a difference of nearly equal terms and its value carries noise near
// 3e-10 fs² on the phase stacks below, so their differences take a wider step:
// 0.25 nm leaves about 2e-9 of that noise in the derivative, and a truncation
// error near 1e-11 against thickness oscillations of period λ/2n.
const PHASE_STEP_NM = 0.25;
const PHASE = { abs: 1e-8, rel: 1e-6 };

// The port reproduces the same arithmetic up to libm round-off.
const PORT = { abs: 1e-12, rel: 1e-7 };

let checks = 0;
let failures = 0;
let worstFraction = 0;
let worstLabel = '';

function agree(analytic, reference, tolerance, label) {
    checks++;
    const diff = Math.abs(analytic - reference);
    const allowed = tolerance.abs + tolerance.rel * Math.abs(reference);
    const fraction = diff / allowed;
    if (fraction > worstFraction) [worstFraction, worstLabel] = [fraction, label];
    if (diff <= allowed) return;
    failures++;
    if (failures <= 25) console.log(`  FAIL ${label}: analytic=${analytic} reference=${reference}`);
}

// Every combination of indices below the given bounds.
const grid = (...bounds) => bounds.reduce(
    (combos, bound) => combos.flatMap(combo => Array.from({ length: bound }, (_, i) => [...combo, i])),
    [[]]);

// ── Finite differences ───────────────────────────────────────────────────────

const central = (f, h = STEP_NM) =>
    (-f(2 * h) + 8 * f(h) - 8 * f(-h) + f(-2 * h)) / (12 * h);

// One-sided, third order: a needle cannot have negative thickness.
const forward = f =>
    (-11 * f(0) + 18 * f(STEP_NM) - 9 * f(2 * STEP_NM) + 2 * f(3 * STEP_NM)) / (6 * STEP_NM);

// The stack with h_i added to layer i and h_j to layer j.
const shifted = (layers, i, j, hi, hj) => layers.map((layer, k) => ({
    n: layer.n,
    d: layer.d + (k === i ? hi : 0) + (k === j ? hj : 0),
}));

// Second partial by the four-point stencil, Richardson-extrapolated from steps
// h and 2h to fourth order. With i = j it is the second difference at 2h.
function secondPartial(evaluate, layers, i, j) {
    const at = h => (evaluate(shifted(layers, i, j, h, h)) - evaluate(shifted(layers, i, j, h, -h))
        - evaluate(shifted(layers, i, j, -h, h)) + evaluate(shifted(layers, i, j, -h, -h))) / (4 * h * h);
    return (4 * at(STEP_NM) - at(2 * STEP_NM)) / 3;
}

// The stack with a needle of index n and thickness h at gap `pos`.
const withNeedle = (layers, pos, n) => h => [...layers.slice(0, pos), { n, d: h }, ...layers.slice(pos)];

// The stack with layer k split at `frac` of its thickness from the incident side.
const splitAt = (layers, k, frac) => [...layers.slice(0, k),
    { n: layers[k].n, d: frac * layers[k].d },
    { n: layers[k].n, d: (1 - frac) * layers[k].d },
    ...layers.slice(k + 1)];

// ── Stacks ───────────────────────────────────────────────────────────────────

const AIR = [1, 0];
const GLASS = [1.52, 0];
const H = [2.35, 0.0005];
const L = [1.46, 0];
const SIO2 = [1.46, 0];
// Al at 550 nm (Rakić): Im δ grows by 0.076 per nm at normal incidence, so
// 700 nm is just past the bound of 50 and 1500 nm far past it.
const AL = [0.96, 6.69];
// Ge and ZnS at 579 nm (Nunley; Querry), where Ge sits on its E1 absorption.
const GE = [5.564, 1.88];
const ZNS = [2.363, 0.001];

const mk = rows => rows.map(([n, d]) => ({ n, d }));

// From glass of index 2 the critical angle into air is 30°. The double nearest
// 30° whose sine is exactly 1/2 puts cosθ exactly at zero in air, where the
// layer matrix is 0/0 as written and every kernel takes its limit.
const DENSE = [2, 0];
const CRITICAL = Array.from({ length: 2001 }, (_, i) => 30 + (i - 1000) * 2 ** -48).find(theta => {
    const rad = theta * Math.PI / 180;
    const c = snellCosTheta(DENSE, [Math.sin(rad), 0], AIR, [Math.cos(rad), 0]);
    return c[0] === 0 && c[1] === 0;
});

const CASES = [
    {
        name: 'dielectric AR', lambda: 550, angles: [0, 45], n0: AIR, ns: GLASS,
        layers: mk([[H, 110], [L, 95], [H, 60], [L, 130]]),
        candidates: [H, L], fracs: [0, 0.3, 1],
    },
    {
        name: 'Al 700 nm between SiO2', lambda: 550, angles: [0, 30], n0: AIR, ns: GLASS,
        layers: mk([[SIO2, 100], [AL, 700], [SIO2, 100]]),
        candidates: [SIO2, H], fracs: [0, 0.05, 0.3, 0.9, 1],
    },
    {
        name: 'Al 1500 nm between SiO2', lambda: 550, angles: [0], n0: AIR, ns: GLASS,
        layers: mk([[SIO2, 100], [AL, 1500], [SIO2, 100]]),
        candidates: [SIO2], fracs: [0.05, 0.3],
    },
    {
        name: 'k = 5 slab, 3 µm and 12 µm', lambda: 532, angles: [0], n0: AIR, ns: GLASS,
        layers: mk([[[1, 5], 3000], [L, 80], [[1, 5], 12000]]),
        candidates: [L], fracs: [0.5],
    },
    {
        // n0 sinθ0 = 1.43 at 70°, so the gap carries an evanescent wave whose
        // phase thickness is purely imaginary and about 116 over 10 µm.
        name: 'air gap 10 µm under a prism at 70°', lambda: 550, angles: [70], n0: GLASS, ns: GLASS,
        layers: mk([[H, 80], [AIR, 10000], [L, 100]]),
        candidates: [H, L], fracs: [0.02, 0.5],
    },
    {
        // A quarter-wave long-wave stack for 10 µm, ZnS outermost, 16 pairs.
        name: 'Ge/ZnS 16 pairs at 579 nm', lambda: 579.4, angles: [0], n0: AIR, ns: GE,
        layers: Array.from({ length: 32 }, (_, i) => (i % 2 ? { n: GE, d: 625 } : { n: ZNS, d: 1136 })),
        candidates: [ZNS, GE], fracs: [0.5],
    },
    {
        // The middle layer and the first candidate sit at the critical angle.
        name: 'air gap at its exact critical angle', lambda: 550, angles: [CRITICAL], n0: DENSE, ns: DENSE,
        layers: mk([[[1.5, 0], 90], [AIR, 120], [[1.5, 0], 70]]),
        candidates: [AIR, [1.5, 0]], fracs: [0, 0.3, 1],
    },
    {
        // A carries the mixed Poynting term here, so dA is not −(dR + dT).
        name: 'absorbing incident medium', lambda: 550, angles: [0, 50], n0: [1.5, 0.01], ns: GLASS,
        layers: mk([[H, 110], [AL, 18], [L, 95]]),
        candidates: [H, L], fracs: [0.5],
    },
];

// One entry per stack, angle and polarization.
const RUNS = CASES.flatMap(c => c.angles.flatMap(theta => [['s', 0], ['p', 1]].map(([pol, code]) => ({
    ...c, theta, pol, code,
    at: `${c.name} θ=${theta} ${pol}`,
    args: [c.lambda, theta, pol, c.n0, c.ns, c.layers],
    kernelArgs: [c.lambda, theta, code, c.n0, c.ns, c.layers],
}))));

const plainOf = run => (layers, q) => tmm(run.lambda, run.theta, run.pol, run.n0, run.ns, layers)[q];

// T through a layer past the imaginary-phase bound is below 1e-43, where the
// absolute tolerance sees nothing at all, so every derivative of T is held
// relative to T as well.
function agreeWithT(q, analytic, reference, T, tolerance, label) {
    agree(analytic, reference, tolerance, label);
    if (q === 'T' && T > 0) agree(analytic / T, reference / T, tolerance, `${label} relative to T`);
}

// ── Jacobian, Hessian and needles against tmm() ──────────────────────────────

const FIRST_KEYS = [['R', 'dRdd'], ['T', 'dTdd'], ['A', 'dAdd']];
const SECOND_KEYS = [['R', 'd2Rdd'], ['T', 'd2Tdd'], ['A', 'd2Add']];

function checkJacobian(run) {
    const plain = plainOf(run);
    const jac = tmmThicknessJacobian(...run.args);
    for (const [k, which] of grid(run.layers.length, FIRST_KEYS.length)) {
        const [q, key] = FIRST_KEYS[which];
        const reference = central(h => plain(shifted(run.layers, k, -1, h, 0), q));
        agreeWithT(q, jac[key][k], reference, jac.T, FIRST, `${run.at} ${key}[${k}]`);
    }
}

function checkHessian(run) {
    const plain = plainOf(run);
    const hess = tmmThicknessHessian(...run.args);
    const N = run.layers.length;
    for (const [i, j, which] of grid(N, N, SECOND_KEYS.length).filter(([i, j]) => j >= i)) {
        const [q, key] = SECOND_KEYS[which];
        const reference = secondPartial(layers => plain(layers, q), run.layers, i, j);
        agreeWithT(q, hess[key][i][j], reference, hess.T, SECOND, `${run.at} ${key}[${i}][${j}]`);
    }
}

function checkNeedles(run) {
    const plain = plainOf(run);
    const scan = tmmNeedleScan(...run.args, run.candidates, run.fracs);
    const N = run.layers.length;
    const nc = run.candidates.length;
    for (const [pos, ci, which] of grid(N + 1, nc, 3)) {
        const q = ['R', 'T', 'A'][which];
        const reference = forward(h => plain(withNeedle(run.layers, pos, run.candidates[ci])(h), q));
        agreeWithT(q, scan.gaps[pos][ci][`d${q}`], reference, scan.T, FIRST,
            `${run.at} gap[${pos}][${ci}] d${q}`);
    }
    for (const [k, fi, ci, which] of grid(N, run.fracs.length, nc, 3)) {
        const q = ['R', 'T', 'A'][which];
        const split = splitAt(run.layers, k, run.fracs[fi]);
        const reference = forward(h => plain(withNeedle(split, k + 1, run.candidates[ci])(h), q));
        agreeWithT(q, scan.intra[k][fi].perCand[ci][`d${q}`], reference, scan.T, FIRST,
            `${run.at} intra[${k}] at ${run.fracs[fi]} [${ci}] d${q}`);
    }
}

for (const run of RUNS) {
    checkJacobian(run);
    checkHessian(run);
    checkNeedles(run);
}

// ── Phase thickness Jacobian against tmmPhaseDispersion() ────────────────────
// One stack whose middle layer passes the imaginary-phase bound, and one whose
// six such layers take the matrix product past the rescale threshold.

const K = jetConstant(2.0, 0.6);
const PHASE_CASES = {
    'one clamped layer': [[jetConstant(2.35, 0.0005), 110], [K, 20000], [jetConstant(1.46), 90]],
    'six clamped layers': [[jetConstant(2.35, 0.0005), 110],
        ...Array.from({ length: 6 }, () => [K, 20000]), [jetConstant(1.46), 90]],
};
const PHASE_RUNS = Object.entries(PHASE_CASES).flatMap(([name, rows]) => [550, 1064].flatMap(lambda =>
    ['s', 'p'].map(pol => ({ name, lambda, pol, layers: rows.map(([nJet, d]) => ({ nJet, d })) }))));

const wrap = angle => angle - 2 * Math.PI * Math.round(angle / (2 * Math.PI));

// Both coefficients. t is checked where |t|² is still inside the double
// range: through the six clamped layers it is 1e-736 at 550 nm and 1e-381 at
// 1064 nm.
function checkPhase({ name, lambda, pol, layers }) {
    const n0Jet = jetConstant(1);
    const nsJet = jetConstant(1.52);
    const jacobian = tmmPhaseThicknessJacobian(lambda, 30, pol, n0Jet, nsJet, layers);
    const plain = tmmPhaseDispersion(lambda, 30, pol, n0Jet, nsJet, layers);
    for (const side of ['r', 't'].filter(s => plain[s].magnitudeSquared > 0)) {
        const jac = jacobian[side];
        const base = plain[side];
        const coefficientAt = (k, h) => tmmPhaseDispersion(lambda, 30, pol, n0Jet, nsJet,
            layers.map((layer, j) => ({ nJet: layer.nJet, d: layer.d + (j === k ? h : 0) })))[side];
        const quantities = [
            ['dPhaseDeg', c => wrap(c.phaseRad - base.phaseRad) * 180 / Math.PI],
            ['dGd', c => c.gd],
            ['dGdd', c => c.gdd],
            ['dLogMagnitudeSquared', c => Math.log(c.magnitudeSquared)],
        ];
        for (const [k, which] of grid(layers.length, quantities.length)) {
            const [key, read] = quantities[which];
            const reference = central(h => read(coefficientAt(k, h)), PHASE_STEP_NM);
            const scale = Math.max(1, Math.abs(read(base)));
            agree((jac[key] || [])[k], reference, { abs: PHASE.abs * scale, rel: PHASE.rel },
                `${name} λ=${lambda} ${pol} ${side} ${key}[${k}]`);
        }
    }
}

PHASE_RUNS.forEach(checkPhase);

// ── The WebAssembly kernel on the same stacks ────────────────────────────────

// [label, JavaScript value, WebAssembly value, quantity] for every derivative
// of one run, the quantity being 'R', 'T' or 'A'.
function portPairs(kernel, run) {
    const N = run.layers.length;
    const ja = tmmThicknessJacobian(...run.args);
    const jb = kernel.tmmJacobian(...run.kernelArgs);
    const ha = tmmThicknessHessian(...run.args);
    const hb = kernel.tmmHessian(...run.kernelArgs);
    const na = tmmNeedleScan(...run.args, run.candidates, run.fracs);
    const nb = kernel.tmmNeedleScan(...run.kernelArgs, run.candidates, run.fracs);
    const Q = ['dR', 'dT', 'dA'];
    const RTA = ['R', 'T', 'A'];
    return [
        ...grid(N, 3).map(([i, w]) => {
            const key = FIRST_KEYS[w][1];
            return [`${key}[${i}]`, ja[key][i], jb[key][i], RTA[w]];
        }),
        ...grid(N, N, 3).map(([i, j, w]) => {
            const key = ['d2Rdd', 'd2Tdd', 'd2Add'][w];
            return [`${key}[${i}][${j}]`, ha[key][i][j], hb[key][i][j], RTA[w]];
        }),
        ...grid(N + 1, run.candidates.length, 3).map(([pos, ci, w]) =>
            [`gap[${pos}][${ci}] ${Q[w]}`, na.gaps[pos][ci][Q[w]], nb.gaps[pos][ci][Q[w]], RTA[w]]),
        ...grid(N, run.fracs.length, run.candidates.length, 3).map(([k, fi, ci, w]) =>
            [`intra[${k}][${fi}][${ci}] ${Q[w]}`,
                na.intra[k][fi].perCand[ci][Q[w]], nb.intra[k][fi].perCand[ci][Q[w]], RTA[w]]),
    ].map(pair => [...pair, ja.T]);
}

// Derivatives of T are held relative to T here too, or the port of every
// derivative through a clamped layer, where T is below 1e-43, goes unchecked.
const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');
const kernel = existsSync(WASM) ? await instantiateTmmWasm(readFileSync(WASM)) : null;
const pairs = kernel ? RUNS.flatMap(run => portPairs(kernel, run).map(([label, ...rest]) =>
    [`${run.at} ${label}`, ...rest])) : [];
for (const [label, js, wa, q, T] of pairs) agreeWithT(q, wa, js, T, PORT, `wasm ${label}`);
if (!kernel) console.log('NOTE : src/tmm_kernel.wasm not built; the WebAssembly checks were skipped.');

console.log(`${checks} comparisons against finite differences and the port, ${pairs.length} of them WebAssembly.`);
console.log(`worst difference as a fraction of its tolerance: ${worstFraction.toExponential(2)}, ${worstLabel}`);
console.log(failures === 0 ? 'PASS : analytic derivatives match finite differences.'
                           : `FAIL : ${failures} comparisons outside tolerance.`);
process.exit(failures === 0 ? 0 : 1);
