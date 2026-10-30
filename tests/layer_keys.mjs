/**
 * Layers whose keys do not match the call.
 *
 * Each call reads one index key and one thickness key from every layer:
 * { n, d } for tmm() and its derivative and needle calls, { nJet, d } for the
 * phase calls, { indexJet, thicknessNm } for the coefficient jets. A layer
 * without either key must throw a TypeError naming the keys and the layer's
 * position. Read as it stood, a layer given as { n, thickness } had an
 * undefined thickness, failed the d > 0 test and was dropped, so the call
 * returned the stack without it: bare glass for a whole design in the wrong
 * shape.
 *
 * A thickness that is present but zero, negative or NaN keeps its meaning:
 * an absent layer, and in the derivative calls a zero column, where zero
 * thickness keeps its derivative.
 *
 * These are the names of 0.4, which ignore a key a layer does not take, as
 * 0.4 did. The calls of 0.5 refuse it (calls.mjs).
 *
 * The WebAssembly methods are checked too when tmm_kernel.wasm has been built.
 *
 *   node tests/layer_keys.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/js/rta.js';
import { tmmThicknessJacobian, tmmThicknessHessian, tmmNeedleScan } from '../src/js/derivatives.js';
import {
    tmmPhaseDispersion, tmmPhaseThicknessJacobian,
    tmmCoefficientJets, tmmCoefficientThicknessJets, omegaFromLambdaNm,
} from '../src/js/phase.js';
import { jetConstant, wavelengthOmegaJet } from '../src/js/jet.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');

// (HL)^2 quarter waves at 550 nm on glass.
const LAMBDA = 550;
const DESIGN = [[2.35, 550 / 4 / 2.35], [1.46, 550 / 4 / 1.46], [2.35, 550 / 4 / 2.35], [1.46, 550 / 4 / 1.46]];
const N0 = [1, 0], NS = [1.52, 0];
const N0_JET = jetConstant(1), NS_JET = jetConstant(1.52);
const WAVELENGTH_JET = wavelengthOmegaJet(LAMBDA, omegaFromLambdaNm(LAMBDA));

const complexLayers = design => design.map(([n, d]) => ({ n: [n, 0], d }));
const jetLayers = design => design.map(([n, d]) => ({ nJet: jetConstant(n), d }));
const coefficientLayers = design => design.map(([n, d]) => ({ indexJet: jetConstant(n), thicknessNm: d }));

// Every public call that takes a layer list, with the keys it reads.
const CALLS = [
    { name: 'tmm', keys: ['n', 'd'], layers: complexLayers,
      run: layers => tmm(LAMBDA, 0, 's', N0, NS, layers) },
    { name: 'tmmThicknessJacobian', keys: ['n', 'd'], layers: complexLayers,
      run: layers => tmmThicknessJacobian(LAMBDA, 0, 's', N0, NS, layers) },
    { name: 'tmmThicknessHessian', keys: ['n', 'd'], layers: complexLayers,
      run: layers => tmmThicknessHessian(LAMBDA, 0, 's', N0, NS, layers) },
    { name: 'tmmNeedleScan', keys: ['n', 'd'], layers: complexLayers,
      run: layers => tmmNeedleScan(LAMBDA, 0, 's', N0, NS, layers, [[2.35, 0]], [0.5]) },
    { name: 'tmmPhaseDispersion', keys: ['nJet', 'd'], layers: jetLayers,
      run: layers => tmmPhaseDispersion(LAMBDA, 0, 's', N0_JET, NS_JET, layers) },
    { name: 'tmmPhaseThicknessJacobian', keys: ['nJet', 'd'], layers: jetLayers,
      run: layers => tmmPhaseThicknessJacobian(LAMBDA, 0, 's', N0_JET, NS_JET, layers) },
    { name: 'tmmCoefficientJets', keys: ['indexJet', 'thicknessNm'], layers: coefficientLayers,
      run: layers => tmmCoefficientJets({ wavelengthJet: WAVELENGTH_JET, thetaDeg: 0, polarization: 's',
          incidentIndexJet: N0_JET, substrateIndexJet: NS_JET, layers }) },
    { name: 'tmmCoefficientThicknessJets', keys: ['indexJet', 'thicknessNm'], layers: coefficientLayers,
      run: layers => tmmCoefficientThicknessJets({ wavelengthJet: WAVELENGTH_JET, thetaDeg: 0, polarization: 's',
          incidentIndexJet: N0_JET, substrateIndexJet: NS_JET, layers }) },
];

const kernel = existsSync(WASM) ? await instantiateTmmWasm(readFileSync(WASM)) : null;
if (kernel) {
    CALLS.push(
        { name: 'wasm tmmOne', keys: ['n', 'd'], layers: complexLayers,
          run: layers => kernel.tmmOne(LAMBDA, 0, 0, N0, NS, layers) },
        { name: 'wasm tmmJacobian', keys: ['n', 'd'], layers: complexLayers,
          run: layers => kernel.tmmJacobian(LAMBDA, 0, 0, N0, NS, layers) },
        { name: 'wasm tmmNeedleScan', keys: ['n', 'd'], layers: complexLayers,
          run: layers => kernel.tmmNeedleScan(LAMBDA, 0, 0, N0, NS, layers, [[2.35, 0]], [0.5]) },
    );
    if (kernel.hasHessian()) CALLS.push({ name: 'wasm tmmHessian', keys: ['n', 'd'], layers: complexLayers,
        run: layers => kernel.tmmHessian(LAMBDA, 0, 0, N0, NS, layers) });
    if (kernel.hasPhase()) CALLS.push(
        { name: 'wasm tmmPhaseOne', keys: ['nJet', 'd'], layers: jetLayers,
          run: layers => kernel.tmmPhaseOne(LAMBDA, 0, 0, N0_JET, NS_JET, layers) },
        { name: 'wasm tmmPhaseJacobian', keys: ['nJet', 'd'], layers: jetLayers,
          run: layers => kernel.tmmPhaseJacobian(LAMBDA, 0, 0, N0_JET, NS_JET, layers) },
    );
    if (kernel.hasGrowingKernels()) CALLS.push({ name: 'wasm monitorCurve', keys: ['n', 'd'], layers: complexLayers,
        run: layers => kernel.monitorCurve(LAMBDA, 0, N0, NS, layers, [1.46, 0], [0, 50, 100]) });
} else {
    console.log('SKIP : WebAssembly methods, src/wasm/tmm_kernel.wasm not built (npm run build:wasm)');
}

// Layer 1 with one key renamed, the way a caller passing another call's shape
// or a misspelling gives it.
function renameKey(layers, from, to) {
    const copy = layers.map(layer => ({ ...layer }));
    copy[1][to] = copy[1][from];
    delete copy[1][from];
    return copy;
}

for (const call of CALLS) {
    const [indexKey, thicknessKey] = call.keys;
    const layers = call.layers(DESIGN);
    assert.doesNotThrow(() => call.run(layers), `${call.name} with {${indexKey}, ${thicknessKey}}`);
    const expected = { name: 'TypeError', message: `layer 1 needs {${indexKey}, ${thicknessKey}}` };
    assert.throws(() => call.run(renameKey(layers, thicknessKey, 'thickness')), expected,
        `${call.name} with the thickness key missing`);
    assert.throws(() => call.run(renameKey(layers, indexKey, 'index')), expected,
        `${call.name} with the index key missing`);
    assert.deepEqual(call.run(layers.map(layer => ({ ...layer, kappa: 0.1 }))), call.run(layers),
        `${call.name} ignores a key a layer does not take`);
}
console.log(`${CALLS.length} calls throw on a layer missing either key, and ignore one they do not take`);

// ── Zero, negative and NaN thickness keep their meaning ─────────────────────

// The design with a fifth layer of thickness d at position 2.
const withExtra = d => [...DESIGN.slice(0, 2), [2.35, d], ...DESIGN.slice(2)];
// -0 counts as zero.
const isZero = (value, label) => assert.ok(value === 0, `${label}: ${value}`);
const close = (a, b, label) => assert.ok(Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(b)), `${label}: ${a} against ${b}`);

const base = tmm(LAMBDA, 0, 's', N0, NS, complexLayers(DESIGN));
const baseJacobian = tmmThicknessJacobian(LAMBDA, 0, 's', N0, NS, complexLayers(DESIGN));
const basePhase = tmmPhaseDispersion(LAMBDA, 0, 's', N0_JET, NS_JET, jetLayers(DESIGN));
const basePhaseJacobian = tmmPhaseThicknessJacobian(LAMBDA, 0, 's', N0_JET, NS_JET, jetLayers(DESIGN));
const baseCoefficient = tmmCoefficientJets({ wavelengthJet: WAVELENGTH_JET, thetaDeg: 0, polarization: 's',
    incidentIndexJet: N0_JET, substrateIndexJet: NS_JET, layers: coefficientLayers(DESIGN) });

for (const d of [0, -5, NaN]) {
    const label = `extra layer of ${d} nm`;
    // Point evaluators: the layer is absent.
    assert.deepEqual(tmm(LAMBDA, 0, 's', N0, NS, complexLayers(withExtra(d))), base, `tmm, ${label}`);
    assert.deepEqual(tmmPhaseDispersion(LAMBDA, 0, 's', N0_JET, NS_JET, jetLayers(withExtra(d))), basePhase,
        `tmmPhaseDispersion, ${label}`);
    const coefficient = tmmCoefficientJets({ wavelengthJet: WAVELENGTH_JET, thetaDeg: 0, polarization: 's',
        incidentIndexJet: N0_JET, substrateIndexJet: NS_JET, layers: coefficientLayers(withExtra(d)) });
    assert.deepEqual(coefficient.reflection, baseCoefficient.reflection, `tmmCoefficientJets, ${label}`);

    // Derivative calls: the same result, a column for every layer as given, and
    // that column zero unless the thickness is zero.
    const jacobian = tmmThicknessJacobian(LAMBDA, 0, 's', N0, NS, complexLayers(withExtra(d)));
    close(jacobian.R, baseJacobian.R, `tmmThicknessJacobian R, ${label}`);
    assert.equal(jacobian.dRdd.length, DESIGN.length + 1, `tmmThicknessJacobian columns, ${label}`);
    [0, 1, 3, 4].forEach((j, k) => close(jacobian.dRdd[j], baseJacobian.dRdd[k], `dR/dd${j}, ${label}`));
    const phaseJacobian = tmmPhaseThicknessJacobian(LAMBDA, 0, 's', N0_JET, NS_JET, jetLayers(withExtra(d)));
    close(phaseJacobian.t.gd, basePhaseJacobian.t.gd, `tmmPhaseThicknessJacobian GD, ${label}`);
    assert.equal(phaseJacobian.t.dGd.length, DESIGN.length + 1, `tmmPhaseThicknessJacobian columns, ${label}`);
    if (d === 0) {
        assert.ok(jacobian.dRdd[2] !== 0, `dR/dd of a zero-thickness layer is kept`);
        assert.ok(phaseJacobian.t.dPhaseDeg[2] !== 0, `phase derivative of a zero-thickness layer is kept`);
    } else {
        isZero(jacobian.dRdd[2], `dR/dd, ${label}`);
        isZero(jacobian.dTdd[2], `dT/dd, ${label}`);
        isZero(phaseJacobian.t.dPhaseDeg[2], `phase derivative, ${label}`);
        isZero(phaseJacobian.t.dGd[2], `GD derivative, ${label}`);
    }

    if (kernel) {
        const one = kernel.tmmOne(LAMBDA, 0, 0, N0, NS, complexLayers(withExtra(d)));
        close(one.R, base.R, `wasm tmmOne R, ${label}`);
        close(one.T, base.T, `wasm tmmOne T, ${label}`);
        const wasmJacobian = kernel.tmmJacobian(LAMBDA, 0, 0, N0, NS, complexLayers(withExtra(d)));
        if (d === 0) assert.ok(wasmJacobian.dRdd[2] !== 0, `wasm dR/dd of a zero-thickness layer is kept`);
        else isZero(wasmJacobian.dRdd[2], `wasm dR/dd, ${label}`);
    }
}
console.log('zero, negative and NaN thicknesses give the documented results');
console.log('PASS : layer keys');
