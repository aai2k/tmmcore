/**
 * The calls of 0.5 on the package index and on a WebAssembly instance: rta,
 * rtaSpectrum, jacobian, jacobianSpectrum, hessian, hessianSpectrum, phase,
 * phaseSpectrum, phaseJacobian and phaseJacobianSpectrum.
 *
 *   • module and instance agree on every call, at one wavelength and over a
 *     grid, with one angle and one per wavelength, in s, p and avg, to the
 *     tolerances of equivalence.mjs, and return the same shapes with their
 *     keys in the same order
 *   • within each, a grid call gives at every wavelength exactly what the
 *     one-wavelength call gives there
 *   • the short forms of a stack read as the long ones, and lists of numbers
 *     and the terms of a jet may be typed arrays
 *   • a misspelled option, a malformed stack, a wavelength list holding
 *     anything but numbers, and an index per wavelength with no derivatives
 *     in a phase call are refused; two index values read as two wavelengths
 *     in a grid call and are refused at one wavelength
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
const RTA_KEYS = new Set(['R', 'T', 'A', 'Rs', 'Ts', 'As', 'Rp', 'Tp', 'Ap']);

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
function sliceAt(grid, i) {
    const keys = Object.keys(grid).filter(key => key in SLICERS);
    return Object.fromEntries(keys.map(key => [key, SLICERS[key](grid[key], i)]));
}

const CALLS = [
    ['rta', 'rtaSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['jacobian', 'jacobianSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['hessian', 'hessianSpectrum', [undefined, 's', 'p', 'avg'], false],
    ['phase', 'phaseSpectrum', ['s', 'p'], true],
    ['phaseJacobian', 'phaseJacobianSpectrum', ['s', 'p'], true],
];
const optionsOf = (angle, pol) => (pol === undefined ? { angle } : { angle, pol });
const angleAt = (angle, i) => (typeof angle === 'number' ? angle : angle[i]);

// One call over the grid and at each of its wavelengths, on both backends.
function checkCall([point, spectrum, , jets], angle, pol) {
    const tag = `${spectrum} θ ${typeof angle === 'number' ? angle : 'per λ'} pol ${pol}`;
    const grids = BACKENDS.map(([, b]) => b[spectrum](jets ? GRID_JET_STACK : GRID_STACK, LAMBDAS, optionsOf(angle, pol)));
    close(grids[0], grids[1], tag);
    LAMBDAS.forEach((lambda, i) => {
        const stack = jets ? pointJetStack(i) : pointStack(i);
        const points = BACKENDS.map(([, b]) => b[point](stack, lambda, optionsOf(angleAt(angle, i), pol)));
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

assert.deepEqual(warnings, [], 'the calls of 0.5 print nothing');

// ── Refused ──────────────────────────────────────────────────────────────────
const ONE = { incident: 1, substrate: 1.52, layers: [{ n: 2.35, d: 58.5 }] };
const withLayer = layer => ({ ...ONE, layers: [layer] });
const REFUSED = [
    [b => b.rta(ONE, 550, { angel: 30 }), /unknown option 'angel'; the options are angle and pol/],
    [b => b.rtaSpectrum(ONE, [550], { polarization: 's' }), /unknown option 'polarization'/],
    [b => b.rta(ONE, 550, null), /options must be \{ angle, pol \}/],
    [b => b.rta(ONE, 550, { pol: 'x' }), /pol must be 's', 'p' or 'avg'/],
    [b => b.phase(ONE, 550), /pol must be 's' or 'p'/],
    [b => b.phaseSpectrum(ONE, [550], { pol: 'avg' }), /pol must be 's' or 'p'/],
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

console.log(`${counts.values} values the same bit for bit. Module against instance: ${differing} of ${compared} `
    + `differ in the last bits; worst R/T/A ${worst.rta.toExponential(2)}, worst relative ${worst.rel.toExponential(2)}`);
console.log('PASS : calls');
