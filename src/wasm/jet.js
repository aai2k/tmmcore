/**
 * jet.js : the jet arithmetic of src/js/jet.js on a WebAssembly instance,
 * under the same names and with the same arguments and results. Each call
 * runs the matching tmm_jet_* export of src/c/jet.c on a grid of one. A jet
 * that is not 4 terms of [re, im] throws.
 */

import { isJet } from '../js/stack.js';

/** Writes a jet, 4 × [re, im], as 8 doubles from `at`. */
export function writeJet(view, at, jet) {
    for (let order = 0; order < 4; order++) {
        view[at + 2 * order] = jet[order][0];
        view[at + 2 * order + 1] = jet[order][1];
    }
}

const readJet = (view, at) => [0, 1, 2, 3].map(order => [view[at + 2 * order], view[at + 2 * order + 1]]);

// The 8 doubles of a jet, read term by term, so a term may be any list of two.
function jetValues(jet) {
    if (!isJet(jet)) throw new TypeError("a jet must be 4 terms [f, f', f''/2!, f'''/3!], each [re, im]");
    const values = new Float64Array(8);
    writeJet(values, 0, jet);
    return values;
}

const asComplex = value => (Array.isArray(value) ? [value[0], value[1]] : [value, 0]);

// Runs export `name` on one item in the instance's scratch memory: the inputs,
// each a list of doubles, then the count, then one output block per entry of
// `outputs`, that many doubles each. Returns the address of the first output.
function runOnOne(wasm, name, inputs, outputs) {
    const total = [...inputs.map(values => values.length), ...outputs].reduce((sum, size) => sum + size, 0);
    const base = wasm._scratch(total);
    const view = wasm._view(base, total);
    let at = 0;
    const inputAddresses = inputs.map(values => {
        view.set(values, at);
        at += values.length;
        return base + 8 * (at - values.length);
    });
    const outputAddresses = outputs.map(size => {
        at += size;
        return base + 8 * (at - size);
    });
    wasm[`_${name}`](...inputAddresses, 1, ...outputAddresses);
    return outputAddresses[0];
}

/** The jet calls of src/js/jet.js and omegaFromLambdaNm, run on `wasm`. */
export function wasmJetCalls(wasm) {
    // A fresh view after the export has run, since memory growth detaches older ones.
    const run = (name, ...inputs) => readJet(wasm._view(runOnOne(wasm, name, inputs, [8]), 8), 0);
    return {
        jetConstant: (value, imaginary = 0) => run('tmm_jet_constant', [value], [imaginary]),
        jetFromDerivatives: (value, first = 0, second = 0, third = 0) =>
            run('tmm_jet_from_derivatives', [value, first, second, third].flatMap(asComplex)),
        jetDerivatives: jet => run('tmm_jet_derivatives', jetValues(jet)),
        jetAdd: (left, right) => run('tmm_jet_add', jetValues(left), jetValues(right)),
        jetSubtract: (left, right) => run('tmm_jet_subtract', jetValues(left), jetValues(right)),
        jetScale: (jet, scalar) => run('tmm_jet_scale', jetValues(jet), [scalar]),
        jetMultiply: (left, right) => run('tmm_jet_multiply', jetValues(left), jetValues(right)),
        jetReciprocal: jet => run('tmm_jet_reciprocal', jetValues(jet)),
        jetDivide: (numerator, denominator) =>
            run('tmm_jet_divide', jetValues(numerator), jetValues(denominator)),
        jetSqrt: jet => run('tmm_jet_sqrt', jetValues(jet)),
        jetExp: jet => run('tmm_jet_exp', jetValues(jet)),
        jetLog: jet => run('tmm_jet_log', jetValues(jet)),
        jetPower: (jet, exponent) => run('tmm_jet_power', jetValues(jet), [exponent]),
        jetSinCos: jet => {
            const sine = runOnOne(wasm, 'tmm_jet_sin_cos', [jetValues(jet)], [8, 8]);
            const both = wasm._view(sine, 16);
            return { sine: readJet(both, 0), cosine: readJet(both, 8) };
        },
        jetSin: jet => run('tmm_jet_sin', jetValues(jet)),
        jetCos: jet => run('tmm_jet_cos', jetValues(jet)),
        jetCompose: (value, derivatives, xJet) => run('tmm_jet_compose', [value],
            [0, 1, 2].map(order => derivatives[order] ?? 0), jetValues(xJet)),
        wavelengthOmegaJet: (wavelength, omega) => run('tmm_wavelength_omega_jet', [wavelength], [omega]),
        jetWithImaginaryPart: (realJet, imaginaryJet) =>
            run('tmm_jet_with_imaginary_part', jetValues(realJet), jetValues(imaginaryJet)),
        jetRealPart: jet => run('tmm_jet_real_part', jetValues(jet)),
        jetClampRealMinimum: (jet, minimum) => run('tmm_jet_clamp_real_minimum', jetValues(jet), [minimum]),
        jetClampImaginary: (jet, limit) => run('tmm_jet_clamp_imaginary', jetValues(jet), [limit]),
        jetIsFinite: jet => {
            const flag = runOnOne(wasm, 'tmm_jet_is_finite', [jetValues(jet)], [1]);
            return new Int32Array(wasm.memory.buffer, flag, 1)[0] !== 0;
        },
        omegaFromLambdaNm: lambda_nm =>
            wasm._view(runOnOne(wasm, 'tmm_omega_from_lambda_nm', [[lambda_nm]], [1]), 1)[0],
    };
}
