/**
 * Writes reference.json beside this file: the WebAssembly kernel's results on
 * the cases of the equivalence suite (tests/cases.mjs), with the calls that
 * produced them in the Python package's names and arguments, for
 * test_reference.py to run through the package and compare.
 *
 *   node python/tests/reference.mjs
 *
 * In the file a complex number, or an array of them, is {"complex": pairs},
 * the pairs [re, im] at the innermost level; a jet is four of them. NaN is
 * written as null.
 */

import { readFileSync, writeFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { instantiateTmmWasm } from '../../src/wasm/index.js';
import { omegaFromLambdaNm } from '../../src/js/phase.js';
import { jetConstant, jetDivide } from '../../src/js/jet.js';
import {
    STACKS, LAMBDAS, ANGLES, POLS, CANDIDATES, INTRA_FRACS, DISPERSIVE, PHASE_STACKS,
    PHASE_SUBSTRATE, PHASE_GRID, n0, ns, cauchyJet,
} from '../../tests/cases.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const kernel = await instantiateTmmWasm(readFileSync(join(HERE, '..', '..', 'src', 'wasm', 'tmm_kernel.wasm')));

// The Python names of the result keys that are words; R, T, dRdd and the
// other symbols keep theirs.
const PYTHON_KEYS = {
    phaseRad: 'phase_rad', phaseDeg: 'phase_deg', magnitudeSquared: 'magnitude_squared',
    dPhaseDeg: 'd_phase_deg', dGd: 'd_gd', dGdd: 'd_gdd', dTod: 'd_tod',
    dLogMagnitudeSquared: 'd_log_magnitude_squared',
};

// A result as plain JSON values in the Python names.
function plain(value) {
    if (ArrayBuffer.isView(value)) return Array.from(value);
    if (Array.isArray(value)) return value.map(plain);
    if (value !== null && typeof value === 'object') {
        return Object.fromEntries(Object.entries(value).map(([key, entry]) => [PYTHON_KEYS[key] ?? key, plain(entry)]));
    }
    return value;
}

const complex = pairs => ({ complex: plain(pairs) });

// Every combination of one entry from each list, the last list varying fastest.
const product = (...lists) => lists.reduce(
    (combos, list) => combos.flatMap(combo => list.map(entry => [...combo, entry])), [[]]);

const calls = [];
const add = (call, args, kwargs, expected) => calls.push({ call, args, kwargs, expected: plain(expected) });

// ── R/T and derivative calls ─────────────────────────────────────────────────

const stackOf = rows => ({
    incident: n0[0], substrate: ns[0],
    layers: rows.map(([[n, k], d]) => ({ n, k, d })),
});
const RTA_POLS = ['s', 'p', 'avg'];
const GRID_ANGLES = [...ANGLES, LAMBDAS.map((_, i) => ANGLES[i])];

for (const [rows, lambda, angle, pol] of product(Object.values(STACKS), LAMBDAS, ANGLES, RTA_POLS)) {
    const stack = stackOf(rows);
    const kwargs = { angle, pol };
    add('rta', [stack, lambda], kwargs, kernel.rta(stack, lambda, kwargs));
    add('jacobian', [stack, lambda], kwargs, kernel.jacobian(stack, lambda, kwargs));
    add('hessian', [stack, lambda], kwargs, kernel.hessian(stack, lambda, kwargs));
}
for (const [rows, angle, pol] of product(Object.values(STACKS), GRID_ANGLES, RTA_POLS)) {
    const stack = stackOf(rows);
    const kwargs = { angle, pol };
    add('rta_spectrum', [stack, LAMBDAS], kwargs, kernel.rtaSpectrum(stack, LAMBDAS, kwargs));
    add('jacobian_spectrum', [stack, LAMBDAS], kwargs, kernel.jacobianSpectrum(stack, LAMBDAS, kwargs));
    add('hessian_spectrum', [stack, LAMBDAS], kwargs, kernel.hessianSpectrum(stack, LAMBDAS, kwargs));
}

// Indices per wavelength: the phase stacks' Cauchy materials at their values.
const valueAt = (lambda, coefficients) => cauchyJet(lambda, omegaFromLambdaNm(lambda), coefficients)[0];
const perWavelength = coefficients => {
    const values = PHASE_GRID.map(lambda => valueAt(lambda, coefficients));
    return { n: values.map(([n]) => n), k: values.map(([, k]) => k) };
};
for (const [rows, pol] of product(Object.values(PHASE_STACKS), RTA_POLS)) {
    const stack = {
        incident: 1, substrate: perWavelength(PHASE_SUBSTRATE),
        layers: rows.map(([name, d]) => ({ ...perWavelength(DISPERSIVE[name]), d })),
    };
    const kwargs = { angle: 30, pol };
    add('rta_spectrum', [stack, PHASE_GRID], kwargs, kernel.rtaSpectrum(stack, PHASE_GRID, kwargs));
    add('jacobian_spectrum', [stack, PHASE_GRID], kwargs, kernel.jacobianSpectrum(stack, PHASE_GRID, kwargs));
}

// ── Needle scan ──────────────────────────────────────────────────────────────
// At one wavelength the package returns, for each of dR, dT and dA, an array
// over [interface][candidate] and one over [layer][fraction][candidate], where
// the JavaScript nests lists of { dR, dT, dA }.

const byKey = pick => Object.fromEntries(['dR', 'dT', 'dA'].map(key => [key, pick(key)]));
const needleArrays = ({ R, T, A, gaps, intra, N }) => ({
    R, T, A,
    gaps: byKey(key => gaps.map(row => row.map(entry => entry[key]))),
    intra: byKey(key => intra.map(layer => layer.map(({ perCand }) => perCand.map(entry => entry[key])))),
    N,
});
const NEEDLE_CANDIDATES = CANDIDATES.map(([n, k]) => ({ n, k }));

for (const [rows, lambda, angle, pol] of product(Object.values(STACKS), LAMBDAS, [0, 45], RTA_POLS)) {
    const stack = stackOf(rows);
    const kwargs = { angle, pol, candidates: NEEDLE_CANDIDATES, fractions: INTRA_FRACS };
    add('needle_scan', [stack, lambda], kwargs, needleArrays(kernel.needleScan(stack, lambda, kwargs)));
}
// Over a grid the candidates after the first have an index per wavelength,
// rising by 0.01 from each wavelength to the next.
const GRID_CANDIDATES = NEEDLE_CANDIDATES.map(({ n, k }, c) => (c === 0 ? { n, k }
    : { n: LAMBDAS.map((_, i) => n + 0.01 * i), k: LAMBDAS.map(() => k) }));
for (const [rows, angle, pol] of product(Object.values(STACKS), GRID_ANGLES, RTA_POLS)) {
    const stack = stackOf(rows);
    const kwargs = { angle, pol, candidates: GRID_CANDIDATES, fractions: INTRA_FRACS };
    add('needle_scan_spectrum', [stack, LAMBDAS], kwargs, kernel.needleScanSpectrum(stack, LAMBDAS, kwargs));
}

// ── Phase calls ──────────────────────────────────────────────────────────────

const jetStack = (rows, jetOf) => ({
    incident: 1,
    substrate: { n: jetOf(PHASE_SUBSTRATE) },
    layers: rows.map(([name, d]) => ({ n: jetOf(DISPERSIVE[name]), d })),
});
const encodeMedium = medium => (typeof medium === 'number' ? medium : { n: complex(medium.n) });
const encodeJets = stack => ({
    incident: encodeMedium(stack.incident),
    substrate: encodeMedium(stack.substrate),
    layers: stack.layers.map(({ n, d }) => ({ n: complex(n), d })),
});
const jetsOnGrid = coefficients => PHASE_GRID.map(lambda => cauchyJet(lambda, omegaFromLambdaNm(lambda), coefficients));

for (const [rows, lambda, angle, [pol]] of product(Object.values(PHASE_STACKS), LAMBDAS, ANGLES, POLS)) {
    const stack = jetStack(rows, coefficients => cauchyJet(lambda, omegaFromLambdaNm(lambda), coefficients));
    const kwargs = { angle, pol };
    add('phase', [encodeJets(stack), lambda], kwargs, kernel.phase(stack, lambda, kwargs));
    add('phase_jacobian', [encodeJets(stack), lambda], kwargs, kernel.phaseJacobian(stack, lambda, kwargs));
}
const PHASE_GRID_ANGLES = [22.5, PHASE_GRID.map((_, i) => i)];
for (const [name, angle, [pol]] of product(['chirped 12-layer', 'with a zero layer'], PHASE_GRID_ANGLES, POLS)) {
    const stack = jetStack(PHASE_STACKS[name], jetsOnGrid);
    const kwargs = { angle, pol };
    add('phase_spectrum', [encodeJets(stack), PHASE_GRID], kwargs, kernel.phaseSpectrum(stack, PHASE_GRID, kwargs));
    add('phase_jacobian_spectrum', [encodeJets(stack), PHASE_GRID], kwargs,
        kernel.phaseJacobianSpectrum(stack, PHASE_GRID, kwargs));
}
// Constant indices, a number for each.
for (const [rows, [pol]] of product(Object.values(STACKS), POLS)) {
    const stack = stackOf(rows);
    const kwargs = { angle: 15, pol };
    add('phase', [stack, 550], kwargs, kernel.phase(stack, 550, kwargs));
    add('phase_spectrum', [stack, LAMBDAS], kwargs, kernel.phaseSpectrum(stack, LAMBDAS, kwargs));
}
// The phase stacks seen from inside the substrate glass, air beyond, lit at
// 45° in air, with the incident sine as a jet; and at normal incidence with
// pol left out.
const fromGlass = (rows, jetOf) => ({
    incident: { n: jetOf(PHASE_SUBSTRATE) }, substrate: 1,
    layers: [...rows].reverse().map(([name, d]) => ({ n: jetOf(DISPERSIVE[name]), d })),
});
const sineIn = glass => jetDivide(jetConstant(Math.SQRT1_2), glass);
for (const [rows, lambda, [pol]] of product(Object.values(PHASE_STACKS), LAMBDAS, POLS)) {
    const stack = fromGlass(rows, coefficients => cauchyJet(lambda, omegaFromLambdaNm(lambda), coefficients));
    const sine = sineIn(stack.incident.n);
    add('phase', [encodeJets(stack), lambda], { pol, sine: complex(sine) }, kernel.phase(stack, lambda, { pol, sine }));
    add('phase_jacobian', [encodeJets(stack), lambda], { pol, sine: complex(sine) },
        kernel.phaseJacobian(stack, lambda, { pol, sine }));
}
for (const [pol] of POLS) {
    const stack = fromGlass(PHASE_STACKS['chirped 12-layer'], jetsOnGrid);
    const sine = stack.incident.n.map(sineIn);
    add('phase_spectrum', [encodeJets(stack), PHASE_GRID], { pol, sine: complex(sine) },
        kernel.phaseSpectrum(stack, PHASE_GRID, { pol, sine }));
    add('phase_jacobian_spectrum', [encodeJets(stack), PHASE_GRID], { pol, sine: complex(sine) },
        kernel.phaseJacobianSpectrum(stack, PHASE_GRID, { pol, sine }));
}
for (const rows of Object.values(PHASE_STACKS)) {
    const stack = jetStack(rows, coefficients => cauchyJet(550, omegaFromLambdaNm(550), coefficients));
    add('phase', [encodeJets(stack), 550], {}, kernel.phase(stack, 550));
    add('phase_jacobian', [encodeJets(stack), 550], {}, kernel.phaseJacobian(stack, 550));
}
const normalGrid = jetStack(PHASE_STACKS['chirped 12-layer'], jetsOnGrid);
add('phase_spectrum', [encodeJets(normalGrid), PHASE_GRID], {}, kernel.phaseSpectrum(normalGrid, PHASE_GRID));

add('omega_from_lambda_nm', [LAMBDAS], {}, LAMBDAS.map(omegaFromLambdaNm));

// ── Jets ─────────────────────────────────────────────────────────────────────
// Each call on a grid of jets in Python against the same call jet by jet here.

const JETS = PHASE_GRID.slice(0, 6).map(lambda => cauchyJet(lambda, omegaFromLambdaNm(lambda), DISPERSIVE.M));
const OTHER = PHASE_GRID.slice(6, 12).map(lambda => cauchyJet(lambda, omegaFromLambdaNm(lambda), DISPERSIVE.H));
const REALS = [0.5, -1.25, 2, 3.5, -0.75, 1e-3];
const jetCall = (call, pythonArgs, run) => add(call, pythonArgs, {}, {
    complex: JETS.map((_, i) => run(i)),
});

jetCall('jet_constant', [REALS.map(value => value * 3), REALS], i => kernel.jetConstant(REALS[i] * 3, REALS[i]));
jetCall('jet_from_derivatives', [complex(JETS.map(jet => jet[0])), REALS, complex(JETS.map(jet => jet[2])), 4],
    i => kernel.jetFromDerivatives(JETS[i][0], REALS[i], JETS[i][2], 4));
jetCall('jet_derivatives', [complex(JETS)], i => kernel.jetDerivatives(JETS[i]));
for (const [call, name] of [['jet_add', 'jetAdd'], ['jet_subtract', 'jetSubtract'], ['jet_multiply', 'jetMultiply'],
    ['jet_divide', 'jetDivide'], ['jet_with_imaginary_part', 'jetWithImaginaryPart']]) {
    jetCall(call, [complex(JETS), complex(OTHER)], i => kernel[name](JETS[i], OTHER[i]));
}
for (const [call, name] of [['jet_scale', 'jetScale'], ['jet_power', 'jetPower'],
    ['jet_clamp_real_minimum', 'jetClampRealMinimum'], ['jet_clamp_imaginary', 'jetClampImaginary']]) {
    jetCall(call, [complex(JETS), REALS], i => kernel[name](JETS[i], REALS[i]));
}
for (const [call, name] of [['jet_reciprocal', 'jetReciprocal'], ['jet_sqrt', 'jetSqrt'], ['jet_exp', 'jetExp'],
    ['jet_log', 'jetLog'], ['jet_sin', 'jetSin'], ['jet_cos', 'jetCos'], ['jet_real_part', 'jetRealPart']]) {
    jetCall(call, [complex(JETS)], i => kernel[name](JETS[i]));
}
add('jet_sin_cos', [complex(JETS)], {}, {
    sine: complex(JETS.map(jet => kernel.jetSinCos(jet).sine)),
    cosine: complex(JETS.map(jet => kernel.jetSinCos(jet).cosine)),
});
const DERIVATIVES = REALS.map((value, i) => [value, REALS[(i + 1) % 6], REALS[(i + 2) % 6]]);
jetCall('jet_compose', [REALS, DERIVATIVES, complex(JETS)], i => kernel.jetCompose(REALS[i], DERIVATIVES[i], JETS[i]));
jetCall('wavelength_omega_jet', [PHASE_GRID.slice(0, 6), PHASE_GRID.slice(0, 6).map(omegaFromLambdaNm)],
    i => kernel.wavelengthOmegaJet(PHASE_GRID[i], omegaFromLambdaNm(PHASE_GRID[i])));
add('jet_is_finite', [complex([JETS[0], [[NaN, 0], ...JETS[1].slice(1)], JETS[2]])], {},
    [true, false, true]);

// ── Monitoring calls ─────────────────────────────────────────────────────────
// A layer of the stack's first material grows on the stack; the growing stack
// takes the first candidate.

const THICKNESSES = [0, 12.5, 40, 133.7];
for (const rows of Object.values(STACKS)) {
    const stack = stackOf(rows);
    const [[n, k]] = rows[0];
    const kwargs = { angle: 30, top: { n, k }, thicknesses: THICKNESSES };
    add('monitor_curve', [stack, 632.8], kwargs, kernel.monitorCurve(stack, 632.8, kwargs));
    add('deposition_spectra', [stack, LAMBDAS], { angle: 30 }, kernel.depositionSpectra(stack, LAMBDAS, { angle: 30 }));

    const growing = kernel.growingStack(stack, LAMBDAS, { angle: 30 });
    const top = NEEDLE_CANDIDATES[0];
    growing.setTop(top);
    calls.push({
        call: 'growing_stack', args: [stack, LAMBDAS], kwargs: { angle: 30 }, top, samples: THICKNESSES,
        expected: plain(THICKNESSES.map(d => growing.sample(d))),
    });
    growing.free();
}

writeFileSync(join(HERE, 'reference.json'), JSON.stringify({ calls }));
console.log(`reference.json: ${calls.length} calls`);
