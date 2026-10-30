/**
 * The calls of 0.5 on the package index and on a WebAssembly instance: rta,
 * rtaSpectrum, jacobian, jacobianSpectrum, hessian, hessianSpectrum,
 * needleScan, needleScanSpectrum, phase, phaseSpectrum, phaseJacobian,
 * phaseJacobianSpectrum, monitorCurve, depositionSpectra and growingStack.
 *
 *   • module and instance agree on every call, at one wavelength and over a
 *     grid, with one angle and one per wavelength, in s, p and avg, to the
 *     tolerances of equivalence.mjs, and return the same shapes with their
 *     keys in the same order
 *   • within each, a grid call gives at every wavelength exactly what the
 *     one-wavelength call gives there, and a growing stack what monitorCurve
 *     gives there
 *   • depositionSpectra deposits the layers from the substrate outward: each
 *     step is the spectrum of the stack's last layers, as many as deposited
 *   • a growing stack writes into the result it is given, and refuses a
 *     sample after free() with the same message on both
 *   • the short forms of a stack read as the long ones, and lists of numbers
 *     and the terms of a jet may be typed arrays
 *   • the phase calls may leave pol out at normal incidence, where s and p
 *     are the same, and take the incident sine as jets in place of the angle:
 *     a constant sine is the angle, and a sine that follows a dispersive
 *     incident medium gives the GD and GDD of differences of the phase in ω
 *   • a misspelled option, a malformed stack, a key the stack, a layer or a
 *     medium does not take, a wavelength list holding
 *     anything but numbers, missing or malformed needle candidates, malformed
 *     fractions, a missing or malformed growing layer or list of thicknesses,
 *     and an index per wavelength with no derivatives in a phase call are
 *     refused; two index values read as two wavelengths in a grid call and are
 *     refused at one wavelength
 *
 * The names of 0.4 are in old_names.mjs. Skips cleanly if tmm_kernel.wasm has
 * not been built.
 *
 *   node tests/calls.mjs
 */

import assert from 'node:assert/strict';
import * as tmmcore from '../src/js/index.js';
import * as J from '../src/js/jet.js';
import {
    loadKernel, LAMBDAS, ANGLES, N, GRID_STACK, GRID_JET_STACK, pointStack, pointJetStack, sameBits, counts,
} from './design.mjs';

const wasm = await loadKernel();
const warnings = [];
console.warn = message => warnings.push(String(message));

const ABS_RTA = 1e-9;      // R, T, A live in [0,1]
const REL_DERIV = 1e-7;    // derivatives and phase quantities span many orders of magnitude
const ABS_DERIV = 1e-12;
const BACKENDS = [['module', tmmcore], ['instance', wasm]];

// ── Module against instance ──────────────────────────────────────────────────

let compared = 0, differing = 0;
const worst = { rta: 0, rel: 0 };
const RTA_KEYS = new Set(['R', 'T', 'A', 'Rs', 'Ts', 'As', 'Rp', 'Tp', 'Ap', 'Rrs', 'Rrp']);

function closeNumber(js, wa, what, key) {
    compared++;
    assert.equal(typeof wa, 'number', `${what}: ${wa} for a number`);
    if (Object.is(js, wa)) return;
    differing++;
    const diff = Math.abs(js - wa);
    if (RTA_KEYS.has(key)) {
        worst.rta = Math.max(worst.rta, diff);
        assert.ok(diff <= ABS_RTA, `${what}: ${js} against ${wa}`);
        return;
    }
    const rel = diff / Math.max(Math.abs(js), Math.abs(wa), Number.MIN_VALUE);
    if (diff > ABS_DERIV) worst.rel = Math.max(worst.rel, rel);
    assert.ok(diff <= ABS_DERIV || rel <= REL_DERIV, `${what}: ${js} against ${wa}`);
}

// The same keys in the same order and the same containers, the numbers within
// round-off: R, T and A to an absolute tolerance, the rest relative.
function close(js, wa, what, key = '') {
    if (typeof js === 'number') return closeNumber(js, wa, what, key);
    if (js === null) return assert.equal(wa, null, what);
    const kind = Object.prototype.toString;
    assert.equal(kind.call(wa), kind.call(js), `${what}: container`);
    if (ArrayBuffer.isView(js) || Array.isArray(js)) {
        assert.equal(wa.length, js.length, `${what}: length`);
        return js.forEach((value, i) => close(value, wa[i], `${what}[${i}]`, key));
    }
    assert.deepEqual(Object.keys(wa), Object.keys(js), `${what}: keys`);
    return Object.keys(js).forEach(name => close(js[name], wa[name], `${what}.${name}`, name));
}

// ── Grids against points ─────────────────────────────────────────────────────
// The needle scans try two candidates, a number and the first layer's material,
// one value per wavelength over the grid, and split each layer at three places,
// the first and last of them its faces.
const FRACTIONS = [0, 0.37, 1];
const needles = i => {
    const { n, k } = i === null ? GRID_STACK.layers[0] : pointStack(i).layers[0];
    return { candidates: [1.38, { n, k }], fractions: FRACTIONS };
};
const N_CAND = 2;

// Wavelength i of a grid result, in the shape of the one-wavelength result.
const SLICERS = { r: (value, i) => sliceAt(value, i), t: (value, i) => sliceAt(value, i), N: value => value };
for (const key of ['R', 'T', 'A', 'phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared']) {
    SLICERS[key] = (value, i) => value[i];
}
for (const key of ['dRdd', 'dTdd', 'dAdd', 'dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared']) {
    SLICERS[key] = (value, i) => value.subarray(i * N, (i + 1) * N);
}
for (const key of ['d2Rdd', 'd2Tdd', 'd2Add']) {
    SLICERS[key] = (value, i) =>
        Array.from({ length: N }, (_, j) => value.subarray((i * N + j) * N, (i * N + j + 1) * N));
}
// A needle scan's gaps as [position][candidate] and intra as [layer][fraction],
// each entry from [(i * (N + 1) + position) * nCand + candidate] and
// [((i * N + layer) * nFrac + fraction) * nCand + candidate] of the flat lists.
const triple = ({ dR, dT, dA }, at) => ({ dR: dR[at], dT: dT[at], dA: dA[at] });
const candidatesAt = (values, at) => Array.from({ length: N_CAND }, (_, c) => triple(values, at * N_CAND + c));
SLICERS.gaps = (values, i) => Array.from({ length: N + 1 }, (_, pos) => candidatesAt(values, i * (N + 1) + pos));
SLICERS.intra = (values, i) => Array.from({ length: N }, (_, k) => FRACTIONS.map((frac, f) =>
    ({ frac, perCand: candidatesAt(values, (i * N + k) * FRACTIONS.length + f) })));
function sliceAt(grid, i) {
    const keys = Object.keys(grid).filter(key => key in SLICERS);
    return Object.fromEntries(keys.map(key => [key, SLICERS[key](grid[key], i)]));
}

const CALLS = [
    ['rta', 'rtaSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['jacobian', 'jacobianSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['hessian', 'hessianSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['needleScan', 'needleScanSpectrum', [undefined, 's', 'p', 'avg'], false, needles],
    ['phase', 'phaseSpectrum', ['s', 'p'], true],
    ['phaseJacobian', 'phaseJacobianSpectrum', ['s', 'p'], true],
];
const optionsOf = (angle, pol) => (pol === undefined ? { angle } : { angle, pol });
const angleAt = (angle, i) => (typeof angle === 'number' ? angle : angle[i]);

// One call over the grid and at each of its wavelengths, on both backends;
// more(i) gives the call's other options, at wavelength i or, for null, over the grid.
function checkCall([point, spectrum, , jets, more = () => ({})], angle, pol) {
    const tag = `${spectrum} θ ${typeof angle === 'number' ? angle : 'per λ'} pol ${pol}`;
    const grids = BACKENDS.map(([, b]) => b[spectrum](jets ? GRID_JET_STACK : GRID_STACK, LAMBDAS,
        { ...more(null), ...optionsOf(angle, pol) }));
    close(grids[0], grids[1], tag);
    LAMBDAS.forEach((lambda, i) => {
        const stack = jets ? pointJetStack(i) : pointStack(i);
        const points = BACKENDS.map(([, b]) => b[point](stack, lambda, { ...more(i), ...optionsOf(angleAt(angle, i), pol) }));
        close(points[0], points[1], `${point} λ ${lambda} ${tag}`);
        BACKENDS.forEach(([name], b) =>
            sameBits(points[b], sliceAt(grids[b], i), `${name} ${spectrum} λ ${lambda} against ${point} ${tag}`));
    });
}

for (const call of CALLS) {
    for (const angle of [0, 45, ANGLES]) call[2].forEach(pol => checkCall(call, angle, pol));
}

// rtaSpectrum also returns s and p, whatever pol asked for; the pol asked for
// comes as its own copy.
for (const [name, b] of BACKENDS) {
    const grid = b.rtaSpectrum(GRID_STACK, LAMBDAS, { angle: ANGLES, pol: 's' });
    assert.notEqual(grid.R, grid.Rs, `${name} rtaSpectrum R is a copy`);
    sameBits(grid.R, grid.Rs, `${name} rtaSpectrum R in s`);
    LAMBDAS.forEach((lambda, i) => ['s', 'p'].forEach(pol => {
        const one = b.rta(pointStack(i), lambda, { angle: ANGLES[i], pol });
        sameBits([grid[`R${pol}`][i], grid[`T${pol}`][i], grid[`A${pol}`][i]], [one.R, one.T, one.A],
            `${name} rtaSpectrum ${pol} λ ${lambda}`);
    }));
}

// ── Short forms ──────────────────────────────────────────────────────────────
// A medium as a number reads as { n, k: 0 }, a missing k as 0, a number in a
// phase call as a constant jet, and lists of wavelengths, angles and indices,
// and the terms of a jet, may be typed arrays.
const SHORT = { incident: 1, substrate: 1.52, layers: [{ n: 2.35, d: 58.5 }, { n: 1.46, k: 0.001, d: 94.2 }] };
const LONG = {
    incident: { n: 1, k: 0 }, substrate: { n: 1.52, k: 0 },
    layers: [{ n: 2.35, k: 0, d: 58.5 }, { n: 1.46, k: 0.001, d: 94.2 }],
};
const LONG_JETS = {
    incident: { n: J.jetConstant(1) }, substrate: { n: J.jetConstant(1.52) },
    layers: [{ n: J.jetConstant(2.35), d: 58.5 }, { n: J.jetConstant(1.46, 0.001), d: 94.2 }],
};
const TYPED_STACK = {
    incident: 1, substrate: { n: Float64Array.from(GRID_STACK.substrate.n) },
    layers: GRID_STACK.layers.map(({ n, k, d }) => ({ n: Float64Array.from(n), k: Float64Array.from(k), d })),
};
const typedTerms = jet => jet.map(term => Float64Array.from(term));
const withTypedTerms = (stack, typed) => ({
    incident: { n: typed(stack.incident.n) }, substrate: { n: typed(stack.substrate.n) },
    layers: stack.layers.map(({ n, d }) => ({ n: typed(n), d })),
});
const TYPED_TERMS = withTypedTerms(pointJetStack(2), typedTerms);
const TYPED_TERMS_GRID = withTypedTerms(GRID_JET_STACK, jets => jets.map(typedTerms));
for (const [name, b] of BACKENDS) {
    sameBits(b.rta(SHORT, 550, { angle: 30 }), b.rta(LONG, 550, { angle: 30 }), `${name} short form, rta`);
    sameBits(b.rtaSpectrum(SHORT, LAMBDAS), b.rtaSpectrum(LONG, LAMBDAS), `${name} short form, rtaSpectrum`);
    sameBits(b.phase(SHORT, 550, { pol: 'p', angle: 30 }), b.phase(LONG_JETS, 550, { pol: 'p', angle: 30 }),
        `${name} short form, phase`);
    sameBits(b.rtaSpectrum(TYPED_STACK, Float64Array.from(LAMBDAS), { angle: Float64Array.from(ANGLES) }),
        b.rtaSpectrum(GRID_STACK, LAMBDAS, { angle: ANGLES }), `${name} typed arrays`);
    sameBits(b.phaseJacobian(TYPED_TERMS, LAMBDAS[2], { pol: 'p', angle: 30 }),
        b.phaseJacobian(pointJetStack(2), LAMBDAS[2], { pol: 'p', angle: 30 }), `${name} typed jet terms`);
    sameBits(b.phaseJacobianSpectrum(TYPED_TERMS_GRID, LAMBDAS, { pol: 's', angle: ANGLES }),
        b.phaseJacobianSpectrum(GRID_JET_STACK, LAMBDAS, { pol: 's', angle: ANGLES }),
        `${name} typed jet terms over the grid`);
    sameBits(b.needleScan(SHORT, 550, { candidates: Float64Array.of(1.38, 2.35), fractions: Float64Array.of(0.5) }),
        b.needleScan(LONG, 550, { candidates: [{ n: 1.38, k: 0 }, { n: 2.35, k: 0 }], fractions: [0.5] }),
        `${name} short forms, needleScan`);
    sameBits(b.needleScan(SHORT, 550, { candidates: [1.38] }), b.needleScan(SHORT, 550, { candidates: [1.38], fractions: [] }),
        `${name} needleScan without fractions`);
}

// Two values over a grid of two wavelengths are two indices, not one complex
// number; at one wavelength they are refused.
const TWO = { incident: 1, substrate: 1.52, layers: [{ n: [2.35, 2.30], d: 58.5 }] };
for (const [name, b] of BACKENDS) {
    const grid = b.rtaSpectrum(TWO, [500, 550]);
    const each = [[2.35, 500], [2.30, 550]].map(([n, lambda]) =>
        b.rta({ incident: 1, substrate: 1.52, layers: [{ n, d: 58.5 }] }, lambda).R);
    sameBits(grid.R, each, `${name} two values over two wavelengths`);
    assert.throws(() => b.rta(TWO, 550), /layer 0 n must be a number at one wavelength/, `${name} two values at one λ`);
}

// ── Polarization of the phase calls at normal incidence ─────────────────────
// s and p are the same there, and pol may be left out.
const PHASE_CALLS = CALLS.filter(([, , , jets]) => jets);
for (const [name, b] of BACKENDS) {
    for (const [point, spectrum] of PHASE_CALLS) {
        const s = b[point](pointJetStack(2), LAMBDAS[2], { pol: 's' });
        sameBits(b[point](pointJetStack(2), LAMBDAS[2], { pol: 'p' }), s, `${name} ${point} p at normal incidence`);
        sameBits(b[point](pointJetStack(2), LAMBDAS[2]), s, `${name} ${point} without pol`);
        sameBits(b[point](pointJetStack(2), LAMBDAS[2], { sine: J.jetConstant(0) }), s,
            `${name} ${point} without pol, a sine of 0`);
        const grid = b[spectrum](GRID_JET_STACK, LAMBDAS, { pol: 'p' });
        sameBits(b[spectrum](GRID_JET_STACK, LAMBDAS, { angle: LAMBDAS.map(() => 0) }), grid,
            `${name} ${spectrum} without pol`);
        sameBits(b[spectrum](GRID_JET_STACK, LAMBDAS, { sine: LAMBDAS.map(() => J.jetConstant(0)) }), grid,
            `${name} ${spectrum} without pol, sines of 0`);
    }
}

// ── The incident sine ────────────────────────────────────────────────────────
// A stack seen from inside Cauchy glass, air beyond it, lit at 45° in air: the
// angle inside the glass follows the glass's index, and `sine` carries it as
// a jet.
const glassAt = lambda => 1.5 + 4000 / (lambda * lambda);
const glassJet = lambda => {
    const lambdaJet = J.wavelengthOmegaJet(lambda, tmmcore.omegaFromLambdaNm(lambda));
    return J.jetAdd(J.jetConstant(1.5), J.jetScale(J.jetReciprocal(J.jetMultiply(lambdaJet, lambdaJet)), 4000));
};
const AIR_SINE = Math.sin(Math.PI / 4);
const sineAt = lambda => J.jetDivide(J.jetConstant(AIR_SINE), glassJet(lambda));
const fromGlass = (glass, layers) => ({ incident: { n: glass }, substrate: 1, layers: [...layers].reverse() });
const SINES = LAMBDAS.map(sineAt);

for (const [point, spectrum] of PHASE_CALLS) {
    for (const pol of ['s', 'p']) {
        const tag = `${spectrum} with the sine, pol ${pol}`;
        const grids = BACKENDS.map(([, b]) =>
            b[spectrum](fromGlass(LAMBDAS.map(glassJet), GRID_JET_STACK.layers), LAMBDAS, { pol, sine: SINES }));
        close(grids[0], grids[1], tag);
        LAMBDAS.forEach((lambda, i) => {
            const stack = fromGlass(glassJet(lambda), pointJetStack(i).layers);
            const points = BACKENDS.map(([, b]) => b[point](stack, lambda, { pol, sine: SINES[i] }));
            close(points[0], points[1], `${point} λ ${lambda} ${tag}`);
            BACKENDS.forEach(([name], b) =>
                sameBits(points[b], sliceAt(grids[b], i), `${name} ${spectrum} λ ${lambda} against ${point} ${tag}`));
        });
    }
}

// GD and GDD with the sine are the frequency derivatives of the phase with the
// angle inside set by Snell's law at each frequency. With the inside angle held
// at its value at 800 nm, which is what a sine without its derivative terms
// gives, they are not, by about 1e-3 in GD and 1e-2 in GDD. The two-layer
// antireflection coating turned round, central differences in ω with step H
// rad/fs: they agree to about 2e-9 in GD, the truncation error, and 5e-8 in
// GDD, the round-off, about 1e-15/H² in fs². The checks allow 1e-6 and 1e-5,
// far above the agreement and far below the gap.
const AR = [{ n: 2.35, d: 58.5 }, { n: 1.46, d: 94.2 }];
const H = 3e-4;
const OMEGA = tmmcore.omegaFromLambdaNm(800);
const insideDeg = lambda => Math.asin(AIR_SINE / glassAt(lambda)) * 180 / Math.PI;
const atOmega = omega => 2 * Math.PI * tmmcore.C_NM_PER_FS / omega;
const wrap = angle => angle - 2 * Math.PI * Math.round(angle / (2 * Math.PI));
for (const [name, b] of BACKENDS) {
    for (const pol of ['s', 'p']) {
        const phi = omega => b.phase(fromGlass(glassAt(atOmega(omega)), AR), atOmega(omega),
            { pol, angle: insideDeg(atOmega(omega)) }).r.phaseRad;
        const [down, mid, up] = [OMEGA - H, OMEGA, OMEGA + H].map(phi);
        const gd = -wrap(up - down) / (2 * H);
        const gdd = -(wrap(up - mid) - wrap(mid - down)) / (H * H);
        const { r } = b.phase(fromGlass(glassJet(800), AR), 800, { pol, sine: sineAt(800) });
        const fixed = b.phase(fromGlass(glassJet(800), AR), 800, { pol, angle: insideDeg(800) }).r;
        const tag = `${name} 800 nm ${pol}`;
        assert.ok(Math.abs(r.gd - gd) < 1e-6 * Math.abs(gd), `${tag}: GD ${r.gd} against ${gd}`);
        assert.ok(Math.abs(r.gdd - gdd) < 1e-5 * Math.abs(gdd), `${tag}: GDD ${r.gdd} against ${gdd}`);
        assert.ok(Math.abs(fixed.gd - gd) > 1e-4 * Math.abs(gd), `${tag}: GD at a fixed inside angle`);
    }
}

// A constant sine is the angle whose sine it is. The kernels take cosθ0 from
// the angle but √(1 − sin²θ0) from a sine, which differ in the last bit, and
// the results by about 1e-15 relative. A real difference between the two paths
// is far larger than the 1e-12 allowed.
for (const [name, b] of BACKENDS) {
    const bySine = b.phase(pointJetStack(2), LAMBDAS[2], { pol: 'p', sine: J.jetConstant(Math.sin(30 * Math.PI / 180)) });
    const byAngle = b.phase(pointJetStack(2), LAMBDAS[2], { pol: 'p', angle: 30 });
    for (const side of ['r', 't']) {
        for (const [key, value] of Object.entries(byAngle[side])) {
            const tag = `${name} constant sine ${side}.${key}: ${bySine[side][key]} against ${value}`;
            assert.ok(Math.abs(bySine[side][key] - value) <= 1e-12 * Math.abs(value), tag);
        }
    }
}

// ── Monitoring calls ─────────────────────────────────────────────────────────
// The first layer of the design grows on the other three: a monitoring curve
// at each wavelength, and a growing stack over the grid.
const GROWN = [0, 12.5, 41];
const finished = stack => ({ ...stack, layers: stack.layers.slice(1) });
const topOf = ({ layers: [{ n, k }] }) => ({ n, k });

for (const angle of [0, 45, ANGLES]) {
    const tag = `θ ${typeof angle === 'number' ? angle : 'per λ'}`;
    const steps = BACKENDS.map(([, b]) => b.depositionSpectra(GRID_STACK, LAMBDAS, { angle }));
    close(steps[0], steps[1], `depositionSpectra ${tag}`);
    const grown = BACKENDS.map(([, b]) => {
        const growing = b.growingStack(finished(GRID_STACK), LAMBDAS, { angle });
        growing.setTop(topOf(GRID_STACK));
        const samples = GROWN.map(d => growing.sample(d));
        growing.free();
        return samples;
    });
    close(grown[0], grown[1], `growingStack ${tag}`);
    LAMBDAS.forEach((lambda, i) => {
        const options = { angle: angleAt(angle, i), top: topOf(pointStack(i)), thicknesses: GROWN };
        const curves = BACKENDS.map(([, b]) => b.monitorCurve(finished(pointStack(i)), lambda, options));
        close(curves[0], curves[1], `monitorCurve λ ${lambda} ${tag}`);
        BACKENDS.forEach(([name], b) => GROWN.forEach((d, s) => {
            const sampled = Object.fromEntries(Object.entries(grown[b][s]).map(([key, values]) => [key, values[i]]));
            const curve = Object.fromEntries(Object.entries(curves[b]).map(([key, values]) => [key, values[s]]));
            sameBits(curve, sampled, `${name} monitorCurve λ ${lambda} d ${d} against growingStack ${tag}`);
        }));
    });
}

// Step k of the deposition, k + 1 layers deposited, is the spectrum of the
// stack's last k + 1 layers.
const lastLayers = count => ({ ...GRID_STACK, layers: GRID_STACK.layers.slice(N - count) });
for (const [name, b] of BACKENDS) {
    const steps = b.depositionSpectra(GRID_STACK, LAMBDAS, { angle: ANGLES });
    const spectra = Array.from({ length: N }, (_, k) => b.rtaSpectrum(lastLayers(k + 1), LAMBDAS, { angle: ANGLES }));
    for (const key of ['Rs', 'Ts', 'Rp', 'Tp']) {
        const expected = spectra.flatMap(spectrum => [...spectrum[key]]);
        steps[key].forEach((step, at) => assert.ok(Math.abs(step - expected[at]) <= ABS_RTA,
            `${name} depositionSpectra ${key} step ${Math.floor(at / LAMBDAS.length)}: ${step} against ${expected[at]}`));
    }
}

// A growing stack's handle, on both with the same messages.
const ONE = { incident: 1, substrate: 1.52, layers: [{ n: 2.35, d: 58.5 }] };
for (const [name, b] of BACKENDS) {
    const growing = b.growingStack(ONE, [500, 550]);
    assert.throws(() => growing.sample(10), /^Error: tmmcore: growing stack sampled before setTop$/,
        `${name} sample before setTop`);
    assert.throws(() => growing.setTop({ n: [1.38, 1.37, 1.36] }),
        /top n must be a number or 2 numbers, one per wavelength/, `${name} setTop`);
    assert.throws(() => growing.sample('10'), /d must be a number \(nm\)/, `${name} sample of text`);
    assert.throws(() => growing.setTop({ n: 1.38, kappa: 0 }), /top has 'kappa'; a medium takes n and k/,
        `${name} setTop with a key a medium does not take`);
    growing.setTop(1.38);
    const out = growing.sample(10);
    assert.equal(growing.sample(20, out), out, `${name} writes into out`);
    sameBits(out, growing.sample(20), `${name} written into out`);
    for (const bad of [{ Rs: new Float64Array(2) }, { ...out, Rrp: new Float64Array(3) }, { ...out, Ts: [0, 0] }]) {
        assert.throws(() => growing.sample(10, bad), /out must be \{ Rs, Ts, Rp, Tp, Rrs, Rrp \}, each a Float64Array of 2/,
            `${name} sample into a wrong out`);
    }
    growing.free();
    assert.throws(() => growing.sample(10), /^Error: tmmcore: growing stack already freed$/, `${name} sample after free`);
    assert.throws(() => growing.setTop(1.38), /^Error: tmmcore: growing stack already freed$/, `${name} setTop after free`);
    growing.free();
}

// ── Refused ──────────────────────────────────────────────────────────────────
const withLayer = layer => ({ ...ONE, layers: [layer] });
const REFUSED = [
    [b => b.rta(ONE, 550, { angel: 30 }), /unknown option 'angel'; the options are angle and pol/],
    [b => b.rtaSpectrum(ONE, [550], { polarization: 's' }), /unknown option 'polarization'/],
    [b => b.rta(ONE, 550, null), /options must be \{ angle, pol \}/],
    [b => b.rta(ONE, 550, { pol: 'x' }), /pol must be 's', 'p' or 'avg'/],
    [b => b.phase(ONE, 550, { angle: 30 }), /pol must be 's' or 'p'; it may be left out at normal incidence/],
    [b => b.phaseSpectrum(ONE, [500, 550], { angle: [0, 10] }), /pol must be 's' or 'p'/],
    [b => b.phase(ONE, 550, { sine: J.jetConstant(0.5) }), /pol must be 's' or 'p'/],
    [b => b.phaseSpectrum(ONE, [500, 550], { sine: [J.jetConstant(0), J.jetConstant(0.5)] }), /pol must be 's' or 'p'/],
    [b => b.phaseSpectrum(ONE, [550], { pol: 'avg' }), /pol must be 's' or 'p'/],
    [b => b.phase(ONE, 550, { angel: 30 }), /unknown option 'angel'; the options are angle, sine and pol/],
    [b => b.rta(ONE, 550, { sine: J.jetConstant(0.5) }), /unknown option 'sine'; the options are angle and pol/],
    [b => b.phase(ONE, 550, { pol: 's', angle: 0, sine: J.jetConstant(0) }), /a phase call takes angle or sine, not both/],
    [b => b.phase(ONE, 550, { pol: 's', sine: 0.5 }), /sine must be a jet/],
    [b => b.phaseSpectrum(ONE, [500, 550], { pol: 's', sine: J.jetConstant(0.5) }),
        /sine must be 2 jets, one per wavelength/],
    [b => b.phaseSpectrum(ONE, [450, 550, 650], { pol: 's', sine: Object.assign(Array(3), { 0: J.jetConstant(0.5),
        2: J.jetConstant(0.5) }) }), /sine must be 3 jets, one per wavelength/],
    [b => b.needleScan(ONE, 550), /candidates must be a list of indices, each a number or \{ n, k \}/],
    [b => b.needleScan(ONE, 550, { candidates: 1.38 }), /candidates must be a list of indices/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38, [1.46, 0]] }), /candidate 1 must be a number or \{ n, k \}/],
    [b => b.needleScanSpectrum(ONE, [500, 550], { candidates: [{ n: [1.38, 1.37, 1.36] }] }),
        /candidate 0 n must be a number or 2 numbers, one per wavelength/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38], fractions: 0.5 }), /fractions must be a list of numbers/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38], fractions: ['0.5'] }), /fractions must be a list of numbers/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38], pol: 'x' }), /pol must be 's', 'p' or 'avg'/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38], sine: J.jetConstant(0) }),
        /unknown option 'sine'; the options are angle, pol, candidates and fractions/],
    [b => b.monitorCurve(ONE, 550), /top must be a number or \{ n, k \}/],
    [b => b.monitorCurve(ONE, 550, { top: [1.38, 0], thicknesses: [0] }), /top must be a number or \{ n, k \}/],
    [b => b.monitorCurve(ONE, 550, { top: { n: [1.38, 1.37] }, thicknesses: [0] }),
        /top n must be a number at one wavelength/],
    [b => b.monitorCurve(ONE, 550, { top: 1.38 }), /thicknesses must be a list of numbers \(nm\)/],
    [b => b.monitorCurve(ONE, 550, { top: 1.38, thicknesses: 50 }), /thicknesses must be a list of numbers \(nm\)/],
    [b => b.monitorCurve(ONE, [550], { top: 1.38, thicknesses: [0] }), /wavelength must be a number \(nm\)/],
    [b => b.monitorCurve(ONE, 550, { top: 1.38, thicknesses: [0], pol: 's' }),
        /unknown option 'pol'; the options are angle, top and thicknesses/],
    [b => b.depositionSpectra(ONE, [500, 550], { pol: 's' }), /unknown option 'pol'; the only option is angle/],
    [b => b.depositionSpectra(ONE, [500, 550], { angle: [0, 10, 20] }), /angle must be a number or 2 numbers/],
    [b => b.depositionSpectra(ONE.layers, [500, 550]), /stack must be \{ incident, substrate, layers/],
    [b => b.growingStack(ONE, [500, 550], { pol: 's' }), /unknown option 'pol'; the only option is angle/],
    [b => b.growingStack(withLayer({ n: 2.35, d: '58.5' }), [500, 550]), /layer 0 d must be a number \(nm\)/],
    [b => b.rta(ONE, [550]), /wavelength must be a number \(nm\)/],
    [b => b.rtaSpectrum(ONE, 550), /wavelengths must be a list of numbers \(nm\)/],
    [b => b.rtaSpectrum(ONE, [{}, 550]), /wavelengths must be a list of numbers \(nm\)/],
    [b => b.hessianSpectrum(ONE, ['550', 600]), /wavelengths must be a list of numbers \(nm\)/],
    [b => b.phaseSpectrum(ONE, Object.assign(Array(3), { 0: 450, 2: 650 }), { pol: 's' }),
        /wavelengths must be a list of numbers \(nm\)/],
    [b => b.rta(ONE, 550, { angle: [0] }), /angle must be a number \(degrees\)/],
    [b => b.rtaSpectrum(ONE, [500, 550], { angle: [0, 1, 2] }), /angle must be a number or 2 numbers/],
    [b => b.rtaSpectrum(ONE, [500, 550], { angle: [0, '30'] }), /angle must be a number or 2 numbers/],
    [b => b.rta({ incident: 1, substrate: 1.52 }, 550), /stack must be \{ incident, substrate, layers/],
    [b => b.rta({ ...ONE, incident: [1, 0] }, 550), /incident must be a number or \{ n, k \}/],
    [b => b.rta({ ...ONE, substrate: undefined }, 550), /substrate must be a number or \{ n, k \}/],
    [b => b.rta(withLayer({ n: 2.35, thickness: 58.5 }), 550), /layer 0 needs \{n, d\}/],
    [b => b.rta(withLayer({ index: 2.35, d: 58.5 }), 550), /layer 0 needs \{n, d\}/],
    [b => b.rta(withLayer({ n: 2.35, d: '58.5' }), 550), /layer 0 d must be a number \(nm\)/],
    [b => b.rta(withLayer({ n: 2.35, kappa: 0.1, d: 58.5 }), 550), /layer 0 has 'kappa'; a layer takes n, k and d/],
    [b => b.phaseSpectrum(withLayer({ n: 2.35, d: 58.5, dn: 0 }), [500, 550], { pol: 's' }),
        /layer 0 has 'dn'; a layer takes n, k and d/],
    [b => b.rta({ ...ONE, angle: 30 }, 550), /stack has 'angle'; a stack takes incident, substrate and layers/],
    [b => b.rtaSpectrum({ ...ONE, substrate: { n: 1.52, kappa: 0 } }, [500, 550]),
        /substrate has 'kappa'; a medium takes n and k/],
    [b => b.needleScan(ONE, 550, { candidates: [1.38, { n: 1.46, kappa: 0 }] }),
        /candidate 1 has 'kappa'; a medium takes n and k/],
    [b => b.monitorCurve(ONE, 550, { top: { n: 1.38, kappa: 0 }, thicknesses: [0] }),
        /top has 'kappa'; a medium takes n and k/],
    [b => b.rta(withLayer({ n: 2.35, k: [0, 0.1], d: 58.5 }), 550), /layer 0 k must be a number at one wavelength/],
    [b => b.rtaSpectrum(withLayer({ n: [2.35, 2.3, 2.2], d: 58.5 }), [500, 550]),
        /layer 0 n must be a number or 2 numbers, one per wavelength/],
    [b => b.rtaSpectrum(withLayer({ n: [[2.35, 0.01], [2.3, 0.01]], d: 58.5 }), [500, 550]),
        /layer 0 n must be a number or 2 numbers, one per wavelength; n and k are two lists/],
    [b => b.jacobianSpectrum({ ...ONE, substrate: { n: ['1.52', '1.51'] } }, [500, 550]),
        /substrate n must be a number or 2 numbers/],
    [b => b.phaseSpectrum(withLayer({ n: [2.35, 2.30], d: 58.5 }), [500, 550], { pol: 's' }),
        /layer 0 n per wavelength has no derivatives, so GD would miss the material's dispersion/],
    [b => b.phaseSpectrum({ ...ONE, substrate: { n: Float64Array.of(1.52, 1.51) } }, [500, 550], { pol: 's' }),
        /substrate n per wavelength has no derivatives/],
    [b => b.phaseSpectrum(withLayer({ n: 2.35, k: [0, 0.01], d: 58.5 }), [500, 550], { pol: 's' }),
        /layer 0 k per wavelength has no derivatives/],
    [b => b.phase(withLayer({ n: J.jetConstant(2.35), k: 0.01, d: 58.5 }), 550, { pol: 's' }),
        /layer 0 k goes inside the jet when n is a jet/],
    [b => b.phase(withLayer({ n: [2.35, 0], d: 58.5 }), 550, { pol: 's' }), /layer 0 n must be a number or a jet/],
    [b => b.phaseSpectrum(withLayer({ n: J.jetConstant(2.35), d: 58.5 }), [500, 550], { pol: 's' }),
        /layer 0 n must be a number or 2 jets, one per wavelength/],
];
for (const [name, b] of BACKENDS) {
    REFUSED.forEach(([run, message], i) => assert.throws(() => run(b), message, `${name} refusal ${i}`));
}
assert.deepEqual(warnings, [], 'the calls of 0.5 print nothing, refused or not');

console.log(`${counts.values} values the same bit for bit. Module against instance: ${differing} of ${compared} `
    + `differ in the last bits; worst R/T/A ${worst.rta.toExponential(2)}, worst relative ${worst.rel.toExponential(2)}`);
console.log('PASS : calls');
