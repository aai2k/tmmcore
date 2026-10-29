/**
 * The names of 0.4, kept through 0.5.x with their old arguments.
 *
 *   • each gives the same bits as the call of 0.5 that replaces it, on the
 *     module and on a WebAssembly instance
 *   • the coefficient-jet calls, which have no replacement, give what the
 *     functions of src/js/phase.js give
 *   • on the module each keeps the name and parameter count of the function
 *     it calls; on an instance the methods are not enumerable, as the methods
 *     of a class are
 *   • each prints one warning, on first use, naming its replacement;
 *     tmmNeedleScan and the growing-layer calls of an instance print none
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/old_names.mjs
 */

import assert from 'node:assert/strict';
import * as tmmcore from '../src/js/index.js';
import * as rtaKernels from '../src/js/rta.js';
import * as derivativeKernels from '../src/js/derivatives.js';
import * as phaseKernels from '../src/js/phase.js';
import * as J from '../src/js/jet.js';
import {
    loadKernel, LAMBDAS, DESIGN, GRID_STACK, GRID_JET_STACK, pointStack, pointJetStack, oldPoint, OLD_GRID,
    sameBits, counts,
} from './design.mjs';

const wasm = await loadKernel();
const warnings = [];
console.warn = message => warnings.push(String(message));

// Each pair: the call of 0.5 and its name of 0.4 on the same design, at one
// wavelength, 30°, in one polarization (code 0 for s and 1 for p).
const i = 2, lambda = LAMBDAS[i], angle = 30;
const old = oldPoint(lambda);
const G = OLD_GRID;
const PAIRS = [
    ['tmm', pol => tmmcore.rta(pointStack(i), lambda, { angle, pol }),
        pol => tmmcore.tmm(lambda, angle, pol, old.n0, old.ns, old.layers)],
    ['tmmThicknessJacobian', pol => tmmcore.jacobian(pointStack(i), lambda, { angle, pol }),
        pol => tmmcore.tmmThicknessJacobian(lambda, angle, pol, old.n0, old.ns, old.layers)],
    ['tmmThicknessHessian', pol => tmmcore.hessian(pointStack(i), lambda, { angle, pol }),
        pol => tmmcore.tmmThicknessHessian(lambda, angle, pol, old.n0, old.ns, old.layers)],
    ['tmmPhaseDispersion', pol => tmmcore.phase(pointJetStack(i), lambda, { angle, pol }),
        pol => tmmcore.tmmPhaseDispersion(lambda, angle, pol, old.n0Jet, old.nsJet, old.jetLayers)],
    ['tmmPhaseThicknessJacobian', pol => tmmcore.phaseJacobian(pointJetStack(i), lambda, { angle, pol }),
        pol => tmmcore.tmmPhaseThicknessJacobian(lambda, angle, pol, old.n0Jet, old.nsJet, old.jetLayers)],
    ['tmmOne', pol => wasm.rta(pointStack(i), lambda, { angle, pol }),
        (pol, code) => wasm.tmmOne(lambda, angle, code, old.n0, old.ns, old.layers)],
    ['tmmSpectrum', () => wasm.rtaSpectrum(GRID_STACK, LAMBDAS, { angle }),
        () => wasm.tmmSpectrum(LAMBDAS, G.n0List, G.nsList, G.layerNK, G.thick, angle)],
    ['tmmJacobian', pol => wasm.jacobian(pointStack(i), lambda, { angle, pol }),
        (pol, code) => wasm.tmmJacobian(lambda, angle, code, old.n0, old.ns, old.layers)],
    ['tmmHessian', pol => wasm.hessian(pointStack(i), lambda, { angle, pol }),
        (pol, code) => wasm.tmmHessian(lambda, angle, code, old.n0, old.ns, old.layers)],
    ['tmmPhaseOne', pol => wasm.phase(pointJetStack(i), lambda, { angle, pol }),
        (pol, code) => wasm.tmmPhaseOne(lambda, angle, code, old.n0Jet, old.nsJet, old.jetLayers)],
    ['tmmPhaseJacobian', pol => wasm.phaseJacobian(pointJetStack(i), lambda, { angle, pol }),
        (pol, code) => wasm.tmmPhaseJacobian(lambda, angle, code, old.n0Jet, old.nsJet, old.jetLayers)],
    ['tmmPhaseSpectrum', pol => wasm.phaseSpectrum(GRID_JET_STACK, LAMBDAS, { angle, pol }),
        (pol, code) => wasm.tmmPhaseSpectrum(LAMBDAS, G.n0Jets, G.nsJets, G.layerJets, G.thick, angle, code)],
    ['tmmPhaseJacobianSpectrum', pol => wasm.phaseJacobianSpectrum(GRID_JET_STACK, LAMBDAS, { angle, pol }),
        (pol, code) => wasm.tmmPhaseJacobianSpectrum(LAMBDAS, G.n0Jets, G.nsJets, G.layerJets, G.thick, angle, code)],
];
for (const [name, current, previous] of PAIRS) {
    for (const [pol, code] of [['s', 0], ['p', 1]]) sameBits(current(pol), previous(pol, code), `${name} ${pol}`);
}

// The grid names of 0.4 read their lists entry by entry, as 0.4 did: typed
// arrays serve as lists, an angle given as text reads as its number, entries
// past the last wavelength go unread, and no layers need no layer list.
const typed = values => Float64Array.from(values);
const SPECTRA = [
    ['tmmSpectrum', (...args) => wasm.tmmSpectrum(...args)],
    ['depositionSpectra', (...args) => wasm.depositionSpectra(...args)],
    ['growingEval', (...args) => {
        const growing = wasm.growingEval(...args);
        growing.setTop(G.layerNK[0]);
        const sample = growing.sample(37.5);
        growing.free();
        return sample;
    }],
];
const SPECTRUM_FORMS = [
    ['typed lists', [typed(LAMBDAS), G.n0List.map(typed), G.nsList.map(typed),
        G.layerNK.map(row => row.map(typed)), typed(G.thick), angle]],
    ['an angle as text', [LAMBDAS, G.n0List, G.nsList, G.layerNK, G.thick, String(angle)]],
    ['rows past the grid', [LAMBDAS, [...G.n0List, [9, 9]], G.nsList,
        G.layerNK.map(row => [...row, [9, 9]]), G.thick, angle]],
];
for (const [name, call] of SPECTRA) {
    const expected = call(LAMBDAS, G.n0List, G.nsList, G.layerNK, G.thick, angle);
    for (const [form, args] of SPECTRUM_FORMS) sameBits(call(...args), expected, `${name} ${form}`);
}
sameBits(wasm.tmmSpectrum(LAMBDAS, G.n0List, G.nsList, undefined, [], angle),
    wasm.tmmSpectrum(LAMBDAS, G.n0List, G.nsList, [], [], angle), 'tmmSpectrum with no layer list');

// An incident sine given as jets stands in for the angle, which may then be
// left null or out.
const SINES = LAMBDAS.map(() => J.jetDivide(J.jetConstant(0.5), J.jetConstant(1.02)));
const typedJets = jets => jets.map(jet => jet.map(typed));
for (const name of ['tmmPhaseSpectrum', 'tmmPhaseJacobianSpectrum']) {
    const call = (lambdas, n0Jets, theta, options) =>
        wasm[name](lambdas, n0Jets, G.nsJets, G.layerJets, G.thick, theta, 1, options);
    sameBits(call(typed(LAMBDAS), typedJets(G.n0Jets), angle), call(LAMBDAS, G.n0Jets, angle), `${name} typed lists`);
    const withSines = call(LAMBDAS, G.n0Jets, angle, { sinJets: SINES });
    for (const theta of [null, undefined]) {
        sameBits(call(LAMBDAS, G.n0Jets, theta, { sinJets: SINES }), withSines, `${name} sines, angle ${theta}`);
    }
}

// The coefficient jets, called through the package and straight from phase.js.
const options = {
    wavelengthJet: J.wavelengthOmegaJet(lambda, phaseKernels.omegaFromLambdaNm(lambda)),
    thetaDeg: angle,
    polarization: 's',
    incidentIndexJet: old.n0Jet,
    substrateIndexJet: old.nsJet,
    layers: DESIGN.layers.map(([m, d]) => ({ indexJet: m.jet(lambda), thicknessNm: d })),
};
const jets = phaseKernels.tmmCoefficientJets(options);
const thickness = phaseKernels.tmmCoefficientThicknessJets(options);
const COEFFICIENT_CALLS = [
    ['tmmCoefficientJets', [options]],
    ['tmmCoefficientThicknessJets', [options]],
    ['coefficientPhaseDispersion', [jets.transmission, jets.transmissionLogScale]],
    ['coefficientPhaseThicknessDerivatives', [thickness.reflection, thickness.reflectionThickness]],
];
for (const [name, args] of COEFFICIENT_CALLS) {
    for (let twice = 0; twice < 2; twice++) sameBits(tmmcore[name](...args), phaseKernels[name](...args), name);
}

// The needle scan and the growing-layer calls have no call of 0.5 in their
// place, and print nothing.
tmmcore.tmmNeedleScan(lambda, angle, 's', old.n0, old.ns, old.layers, [[1.38, 0]], [0.5]);
wasm.tmmNeedleScan(lambda, angle, 0, old.n0, old.ns, old.layers, [[1.38, 0]], [0.5]);
wasm.monitorCurve(lambda, angle, old.n0, old.ns, old.layers, [1.38, 0], [0, 50]);
wasm.depositionSpectra(LAMBDAS, G.n0List, G.nsList, G.layerNK, G.thick, angle);
wasm.growingEval(LAMBDAS, G.n0List, G.nsList, G.layerNK, G.thick, angle).free();

const REPLACED = {
    tmm: 'rta', tmmThicknessJacobian: 'jacobian', tmmThicknessHessian: 'hessian',
    tmmPhaseDispersion: 'phase', tmmPhaseThicknessJacobian: 'phaseJacobian',
    tmmOne: 'rta', tmmSpectrum: 'rtaSpectrum', tmmJacobian: 'jacobian', tmmHessian: 'hessian',
    tmmPhaseOne: 'phase', tmmPhaseJacobian: 'phaseJacobian',
    tmmPhaseSpectrum: 'phaseSpectrum', tmmPhaseJacobianSpectrum: 'phaseJacobianSpectrum',
    tmmCoefficientJets: null, tmmCoefficientThicknessJets: null,
    coefficientPhaseDispersion: null, coefficientPhaseThicknessDerivatives: null,
};
const expected = Object.entries(REPLACED).map(([name, replacement]) =>
    `tmmcore: ${name} is deprecated and will be removed in 0.6.0${replacement ? `; use ${replacement}` : ''}`);
assert.deepEqual([...warnings].sort(), expected.sort(), 'one warning per old name');

const KERNELS = { ...rtaKernels, ...derivativeKernels, ...phaseKernels };
const moduleNames = Object.keys(REPLACED).filter(name => name in tmmcore);
for (const name of moduleNames) {
    assert.equal(tmmcore[name].name, name, `${name}: name`);
    assert.equal(tmmcore[name].length, KERNELS[name].length, `${name}: parameter count`);
}
assert.deepEqual(Object.keys(Object.getPrototypeOf(wasm)), [], 'no enumerable methods on an instance');

console.log(`${counts.values} values the same bit for bit; ${warnings.length} names warned once each`);
console.log('PASS : old names');
