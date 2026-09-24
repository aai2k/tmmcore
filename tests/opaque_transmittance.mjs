/**
 * Transmittance through layers past the opaque-layer bound.
 *
 * layerMatrix holds |Im δ| of a layer at 50 so that cosh δ cannot overflow,
 * and the factor it holds back is carried as a log scale for t. Up to 0.4.1
 * that factor was dropped: R stayed right, and T stopped near 1e-43 for each
 * such layer however thick it was, about 1e23 too large on 1 µm of aluminium
 * and up to 1e91 on the 10 µm air gap below.
 *
 * Every evaluator that reports T or |t| is held here to the Airy recursion,
 * which builds r and t layer by layer from interface coefficients and
 * propagation factors e^{iδ}, never forms cosh, and shares no code with the
 * characteristic matrix:
 *
 *   • tmm(), and T as the Jacobian, Hessian and needle scan report it
 *   • |t|² and the phase of t from the phase kernel
 *   • the WebAssembly kernel: tmm_one, tmm_spectrum, the phase kernels and
 *     the three growing-stack kernels
 *
 * An ordinary dielectric stack runs alongside, where nothing is held.
 *
 *   node tests/opaque_transmittance.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import {
    tmm, tmmThicknessJacobian, tmmThicknessHessian, tmmNeedleScan, layerLogScale,
} from '../src/tmm.js';
import { tmmPhaseDispersion, tmmPhaseThicknessJacobian } from '../src/phase.js';
import { jetConstant } from '../src/taylorJet.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));

let checks = 0;
const relNear = (actual, expected, tolerance, message) => {
    checks++;
    assert.ok(Math.abs(actual - expected) <= tolerance * Math.abs(expected),
        `${message}: got ${actual}, expected ${expected}, relative tolerance ${tolerance}`);
};
const near = (actual, expected, tolerance, message) => {
    checks++;
    assert.ok(Math.abs(actual - expected) <= tolerance,
        `${message}: got ${actual}, expected ${expected}, tolerance ${tolerance}`);
};

// ── The Airy recursion ───────────────────────────────────────────────────────
// Tilted admittances in free-space units, η = n cosθ (s) or n / cosθ (p), with
// n cosθ = sqrt(n² − a²), a = n0 sinθ0 real, on the root whose imaginary part
// is not negative, so that e^{iδ} decays. With ρ = (ηa − ηb)/(ηa + ηb) and
// τ = 2ηa/(ηa + ηb) at each interface, from the substrate outward:
//   r ← (ρ + r e^{2iδ}) / (1 + ρ r e^{2iδ}),   t ← τ t e^{iδ} / (1 + ρ r e^{2iδ})
// and T = Re(ηs)/Re(η0) |t|².

const add = ([a, b], [c, d]) => [a + c, b + d];
const sub = ([a, b], [c, d]) => [a - c, b - d];
const mul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const div = ([a, b], [c, d]) => {
    const q = c * c + d * d;
    return [(a * c + b * d) / q, (b * c - a * d) / q];
};
const abs2 = ([a, b]) => a * a + b * b;
const expI = ([a, b]) => [Math.exp(-b) * Math.cos(a), Math.exp(-b) * Math.sin(a)];

// The larger part of the root comes from the modulus and the smaller one from
// y divided by it, so a small imaginary part keeps its digits.
function normalComponent(n, a) {
    const [x, y] = sub(mul(n, n), [a * a, 0]);
    const m = Math.hypot(x, y);
    let root;
    if (x >= 0) {
        const re = Math.sqrt((m + x) / 2);
        root = [re, y / (2 * re)];
    } else {
        const im = Math.sqrt((m - x) / 2);
        root = [y / (2 * im), im];
    }
    if (root[1] < 0 || (root[1] === 0 && root[0] < 0)) root = [-root[0], -root[1]];
    return root;
}

function airy(lambda, thetaDeg, pol, n0, ns, layers) {
    const rad = thetaDeg * Math.PI / 180;
    const a = n0[0] * Math.sin(rad);
    // The incident medium is transparent here and keeps its own cosine.
    const eta = (n, i) => {
        const q = i === 0 ? [n0[0] * Math.cos(rad), 0] : normalComponent(n, a);
        return pol === 's' ? q : div(mul(n, n), q);
    };
    const media = [n0, ...layers.map(l => l.n), ns];
    const etas = media.map(eta);
    const k0 = 2 * Math.PI / lambda;
    const last = media.length - 2;
    let r = div(sub(etas[last], etas[last + 1]), add(etas[last], etas[last + 1]));
    let t = div(mul([2, 0], etas[last]), add(etas[last], etas[last + 1]));
    for (let j = last - 1; j >= 0; j--) {
        const layer = layers[j];
        const phase = expI(mul(normalComponent(layer.n, a), [k0 * layer.d, 0]));
        const rho = div(sub(etas[j], etas[j + 1]), add(etas[j], etas[j + 1]));
        const tau = div(mul([2, 0], etas[j]), add(etas[j], etas[j + 1]));
        const round = mul(r, mul(phase, phase));
        const den = add([1, 0], mul(rho, round));
        r = div(add(rho, round), den);
        t = div(mul(tau, mul(t, phase)), den);
    }
    return { R: abs2(r), T: etas[media.length - 1][0] / etas[0][0] * abs2(t), t, eta0: etas[0], etaS: etas[media.length - 1] };
}

// ── Stacks ───────────────────────────────────────────────────────────────────

const AIR = [1, 0], GLASS = [1.52, 0], SIO2 = [1.46, 0], MGF2 = [1.38, 0], TIO2 = [2.35, 0.0005];
// Al at 550 nm (Rakić): Im δ grows by 0.0764 per nm at normal incidence.
const AL = [0.96, 6.69];
const mk = rows => rows.map(([n, d]) => ({ n, d }));

const CASES = [
    // Im δ of the Al is 76 at normal incidence: 26 past the bound.
    { name: 'SiO2 90 nm over Al 1000 nm', lambdas: [550], angles: [0, 30, 60],
        n0: AIR, ns: GLASS, layers: mk([[SIO2, 90], [AL, 1000]]), held: true },
    // Im δ 191: T near 1e-166.
    { name: 'Al 2500 nm', lambdas: [550], angles: [0, 45],
        n0: AIR, ns: GLASS, layers: mk([[AL, 2500]]), held: true },
    // Evanescent past 41.1°: Im δ 55 at 45° and 450 nm, 127 at 80° and 550 nm.
    { name: 'air gap 10 µm in glass', lambdas: [450, 550], angles: [45, 80],
        n0: GLASS, ns: GLASS, layers: mk([[TIO2, 80], [AIR, 10000], [MGF2, 100]]), held: true },
    { name: 'dielectric AR', lambdas: [450, 550, 700], angles: [0, 45],
        n0: AIR, ns: GLASS, layers: mk([[TIO2, 110], [SIO2, 95], [TIO2, 60], [MGF2, 130]]), held: false },
];

const RUNS = CASES.flatMap(c => c.lambdas.flatMap(lambda => c.angles.flatMap(theta =>
    ['s', 'p'].map(pol => ({ ...c, lambda, theta, pol, code: pol === 'p' ? 1 : 0,
        at: `${c.name}, ${lambda} nm, ${theta}° ${pol}` })))));

// R, which the held factor cannot reach, agrees to rounding amplified by the
// stack: up to 2e-12 where the absorbing TiO2 sits over the evanescent gap at
// 80°. T past the bound is what is under test.
const R_TOLERANCE = 5e-12;
const T_TOLERANCE = 1e-12;
const asJets = layers => layers.map(({ n, d }) => ({ nJet: jetConstant(n[0], n[1]), d }));
const wrap = angle => angle - 2 * Math.PI * Math.round(angle / (2 * Math.PI));

// ── The JavaScript evaluators ────────────────────────────────────────────────

let worstBefore = 1;
const heldBy = new Map();
for (const run of RUNS) {
    const args = [run.lambda, run.theta, run.pol, run.n0, run.ns, run.layers];
    const exact = airy(...args);
    const plain = tmm(...args);
    near(plain.R, exact.R, R_TOLERANCE, `${run.at}: R`);
    relNear(plain.T, exact.T, T_TOLERANCE, `${run.at}: T`);

    // What the held matrix alone gives, as up to 0.4.1: every held layer
    // puts T too high by e^{2·excess}.
    const excess = run.layers.reduce((sum, { n, d }) => {
        const rad = run.theta * Math.PI / 180;
        const cosTheta = div(normalComponent(n, run.n0[0] * Math.sin(rad)), n);
        return sum + layerLogScale(n, d, run.lambda, cosTheta);
    }, 0);
    if (!run.held) assert.equal(excess, 0, `${run.at}: a layer past the bound`);
    heldBy.set(run.name, Math.max(heldBy.get(run.name) ?? 0, excess));
    worstBefore = Math.max(worstBefore, Math.exp(2 * excess));

    for (const [name, result] of [
        ['Jacobian', tmmThicknessJacobian(...args)],
        ['Hessian', tmmThicknessHessian(...args)],
        ['needle scan', tmmNeedleScan(...args, [SIO2], [0.5])],
    ]) {
        near(result.R, exact.R, R_TOLERANCE, `${run.at}: R of the ${name}`);
        relNear(result.T, exact.T, T_TOLERANCE, `${run.at}: T of the ${name}`);
    }

    // The phase kernel: |t|² times the admittance ratio is T, and the phase
    // of t is the one the recursion builds (reported with Macleod's sign).
    const Tfac = exact.etaS[0] / exact.eta0[0];
    const n0Jet = jetConstant(run.n0[0], run.n0[1]);
    const nsJet = jetConstant(run.ns[0], run.ns[1]);
    for (const [name, phase] of [
        ['phase kernel', tmmPhaseDispersion(run.lambda, run.theta, run.pol, n0Jet, nsJet, asJets(run.layers))],
        ['phase Jacobian', tmmPhaseThicknessJacobian(run.lambda, run.theta, run.pol, n0Jet, nsJet, asJets(run.layers))],
    ]) {
        relNear(Tfac * phase.t.magnitudeSquared, exact.T, T_TOLERANCE, `${run.at}: |t|² of the ${name}`);
        near(wrap(phase.t.phaseRad + Math.atan2(exact.t[1], exact.t[0])), 0, 1e-9,
            `${run.at}: phase of t from the ${name}`);
    }
}

for (const c of CASES.filter(c => c.held)) {
    assert.ok(heldBy.get(c.name) > 0, `${c.name}: no layer past the bound`);
}

// ── The WebAssembly kernel ───────────────────────────────────────────────────

const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');
const kernel = existsSync(WASM) ? await instantiateTmmWasm(readFileSync(WASM)) : null;
if (kernel) {
    for (const run of RUNS) {
        const exact = airy(run.lambda, run.theta, run.pol, run.n0, run.ns, run.layers);
        const kargs = [run.lambda, run.theta, run.code, run.n0, run.ns, run.layers];
        const one = kernel.tmmOne(...kargs);
        near(one.R, exact.R, R_TOLERANCE, `wasm ${run.at}: R`);
        relNear(one.T, exact.T, T_TOLERANCE, `wasm ${run.at}: T`);
        relNear(kernel.tmmJacobian(...kargs).T, exact.T, T_TOLERANCE, `wasm ${run.at}: T of the Jacobian`);
        relNear(kernel.tmmHessian(...kargs).T, exact.T, T_TOLERANCE, `wasm ${run.at}: T of the Hessian`);
        relNear(kernel.tmmNeedleScan(...kargs, [SIO2], [0.5]).T, exact.T, T_TOLERANCE,
            `wasm ${run.at}: T of the needle scan`);

        const Tfac = exact.etaS[0] / exact.eta0[0];
        const jets = [jetConstant(run.n0[0], run.n0[1]), jetConstant(run.ns[0], run.ns[1]), asJets(run.layers)];
        for (const [name, phase] of [
            ['phase kernel', kernel.tmmPhaseOne(run.lambda, run.theta, run.code, ...jets)],
            ['phase Jacobian', kernel.tmmPhaseJacobian(run.lambda, run.theta, run.code, ...jets)],
        ]) {
            relNear(Tfac * phase.t.magnitudeSquared, exact.T, T_TOLERANCE, `wasm ${run.at}: |t|² of the ${name}`);
            near(wrap(phase.t.phaseRad + Math.atan2(exact.t[1], exact.t[0])), 0, 1e-9,
                `wasm ${run.at}: phase of t from the ${name}`);
        }
    }

    // The batched spectrum and the three growing-stack kernels, over each
    // stack's wavelengths. The growing kernels take the outermost layer as the
    // one being deposited on the rest.
    for (const c of CASES) {
        for (const theta of c.angles) {
            const exact = pol => c.lambdas.map(lambda => airy(lambda, theta, pol, c.n0, c.ns, c.layers));
            const [s, p] = [exact('s'), exact('p')];
            const at = `${c.name}, ${theta}°`;
            const nkRows = c.layers.map(({ n }) => c.lambdas.map(() => n));
            const n0List = c.lambdas.map(() => c.n0), nsList = c.lambdas.map(() => c.ns);
            const spectrum = kernel.tmmSpectrum(c.lambdas, n0List, nsList, nkRows,
                c.layers.map(l => l.d), theta);
            const [top, ...below] = c.layers;
            const growing = kernel.growingEval(c.lambdas, n0List, nsList, below.map(({ n }) => c.lambdas.map(() => n)),
                below.map(l => l.d), theta);
            growing.setTop(c.lambdas.map(() => top.n));
            const sampled = growing.sample(top.d);
            growing.free();
            // Deposition order: substrate side first, the finished stack last.
            const reversed = [...c.layers].reverse();
            const deposited = kernel.depositionSpectra(c.lambdas, n0List, nsList,
                reversed.map(({ n }) => c.lambdas.map(() => n)), reversed.map(l => l.d), theta);
            const lastStep = (c.layers.length - 1) * c.lambdas.length;
            c.lambdas.forEach((lambda, i) => {
                const curve = kernel.monitorCurve(lambda, theta, c.n0, c.ns, below, top.n, [top.d]);
                for (const [name, Ts, Tp] of [
                    ['tmm_spectrum', spectrum.Ts[i], spectrum.Tp[i]],
                    ['growing evaluator', sampled.Ts[i], sampled.Tp[i]],
                    ['deposition spectra', deposited.Ts[lastStep + i], deposited.Tp[lastStep + i]],
                    ['monitor curve', curve.Ts[0], curve.Tp[0]],
                ]) {
                    relNear(Ts, s[i].T, T_TOLERANCE, `wasm ${at}, ${lambda} nm: T s of the ${name}`);
                    relNear(Tp, p[i].T, T_TOLERANCE, `wasm ${at}, ${lambda} nm: T p of the ${name}`);
                }
            });
        }
    }
} else {
    console.log('NOTE : src/tmm_kernel.wasm not built; the WebAssembly checks were skipped.');
}

console.log(`${checks} checks against the Airy recursion${kernel ? ', WebAssembly included' : ''}.`);
console.log(`without the held factor, T would be off by up to ${worstBefore.toExponential(1)}`);
console.log('PASS : transmittance past the opaque-layer bound.');
