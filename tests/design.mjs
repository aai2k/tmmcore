/**
 * The design calls.mjs and old_names.mjs run: four layers of Cauchy materials,
 * one of them a metal, on a dispersive glass, written as a stack of 0.5 and as
 * the arguments of 0.4, with the bit-for-bit comparison both tests use.
 *
 * Cauchy indices n = A + B/λ² (B in nm²) with a constant k, as numbers and as
 * a jet in ω for the phase calls.
 */

import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { omegaFromLambdaNm } from '../src/js/phase.js';
import * as J from '../src/js/jet.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

/** An instance of the built kernel; exits the test with SKIP when it is not built. */
export async function loadKernel() {
    const path = join(dirname(fileURLToPath(import.meta.url)), '..', 'src', 'wasm', 'tmm_kernel.wasm');
    if (!existsSync(path)) {
        console.log('SKIP : src/wasm/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
        process.exit(0);
    }
    return instantiateTmmWasm(readFileSync(path));
}

export const LAMBDAS = [420, 500, 550, 633, 800];
export const ANGLES = [0, 20, 40, 61, 79];

const material = (A, B, k) => ({
    n: lambda => A + B / (lambda * lambda),
    k: () => k,
    jet: lambda => {
        const lambdaJet = J.wavelengthOmegaJet(lambda, omegaFromLambdaNm(lambda));
        const inverseSquare = J.jetReciprocal(J.jetMultiply(lambdaJet, lambdaJet));
        return J.jetAdd(J.jetConstant(A, k), J.jetScale(inverseSquare, B));
    },
});

export const DESIGN = {
    incident: material(1, 0, 0),
    substrate: material(1.5, 4200, 0),
    layers: [
        [material(2.2, 4.0e4, 2e-4), 62],
        [material(1.44, 4.0e3, 0), 105],
        [material(0.15, 0, 3.2), 22],
        [material(2.2, 4.0e4, 2e-4), 41],
    ],
};
export const N = DESIGN.layers.length;

// The design as a stack of 0.5, each value picked by `pick` from a function
// of λ: one value per wavelength for a grid call, one number for a point.
function stackOf(pick, jets) {
    const medium = m => (jets ? { n: pick(m.jet) } : { n: pick(m.n), k: pick(m.k) });
    return {
        incident: medium(DESIGN.incident),
        substrate: medium(DESIGN.substrate),
        layers: DESIGN.layers.map(([m, d]) => ({ ...medium(m), d })),
    };
}
const overGrid = f => LAMBDAS.map(f);
const atPoint = i => f => f(LAMBDAS[i]);
export const GRID_STACK = stackOf(overGrid, false);
export const GRID_JET_STACK = stackOf(overGrid, true);
export const pointStack = i => stackOf(atPoint(i), false);
export const pointJetStack = i => stackOf(atPoint(i), true);

// The same design in the arguments of 0.4.
const nk = (m, lambda) => [m.n(lambda), m.k(lambda)];
export const oldPoint = lambda => ({
    n0: nk(DESIGN.incident, lambda),
    ns: nk(DESIGN.substrate, lambda),
    layers: DESIGN.layers.map(([m, d]) => ({ n: nk(m, lambda), d })),
    n0Jet: DESIGN.incident.jet(lambda),
    nsJet: DESIGN.substrate.jet(lambda),
    jetLayers: DESIGN.layers.map(([m, d]) => ({ nJet: m.jet(lambda), d })),
});
export const OLD_GRID = {
    n0List: LAMBDAS.map(lambda => nk(DESIGN.incident, lambda)),
    nsList: LAMBDAS.map(lambda => nk(DESIGN.substrate, lambda)),
    layerNK: DESIGN.layers.map(([m]) => LAMBDAS.map(lambda => nk(m, lambda))),
    n0Jets: LAMBDAS.map(DESIGN.incident.jet),
    nsJets: LAMBDAS.map(DESIGN.substrate.jet),
    layerJets: DESIGN.layers.map(([m]) => LAMBDAS.map(m.jet)),
    thick: DESIGN.layers.map(([, d]) => d),
};

export const counts = { values: 0 };

/** Every number in `expected` against the same place in `actual`, bit for bit. */
export function sameBits(actual, expected, what) {
    if (typeof expected === 'number') {
        counts.values++;
        if (!Object.is(actual, expected)) assert.fail(`${what}: ${actual} against ${expected}`);
        return;
    }
    if (expected === null) return assert.equal(actual, null, what);
    assert.ok(actual !== null && typeof actual === 'object', `${what}: missing`);
    if (ArrayBuffer.isView(expected) || Array.isArray(expected)) {
        assert.equal(actual.length, expected.length, `${what}: length`);
    }
    for (const name of Object.keys(expected)) sameBits(actual[name], expected[name], `${what}.${name}`);
}
