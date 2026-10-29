/**
 * The kernel's exports as a C or Python caller binds them.
 *
 *   • the grid calls of the Jacobian, Hessian, needle scan and phase give, at
 *     every wavelength, exactly what the one-wavelength calls give there
 *   • every grid call, given one angle for the whole grid or one per
 *     wavelength, gives exactly what one call per wavelength at that angle
 *     gives
 *   • every jet export, applied to a grid of jets, agrees with src/js/jet.js
 *     jet by jet: bit for bit where the arithmetic is plain, to round-off
 *     where it goes through exp, log, sin, cos or a square root, whose last
 *     digits differ between the C library and JavaScript's Math
 *   • every jet call of an instance takes the arguments its namesake in
 *     src/js/jet.js takes, defaults and real numbers for complex values
 *     included, and agrees with it in the same way
 *   • a jet whose terms are typed arrays gives what the same jet as plain
 *     lists gives, and anything but 4 terms of two in place of a jet throws
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/kernel_exports.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { omegaFromLambdaNm } from '../src/js/phase.js';
import * as J from '../src/js/jet.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');
if (!existsSync(WASM)) {
    console.log('SKIP : src/wasm/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}
const kernel = await instantiateTmmWasm(readFileSync(WASM));

// Round-off allowed for the jet exports that call into the math library,
// relative to the largest coefficient of the jet.
const LIBM_TOL = 1e-13;

// ── Kernel memory ────────────────────────────────────────────────────────────
// Each check places its inputs, calls, reads, and frees everything it placed.
// A view is made for every read, since memory growth detaches older ones.
function withArena(body) {
    const blocks = [];
    const place = bytes => {
        const ptr = kernel.malloc(Math.max(8, bytes));
        assert.ok(ptr, 'kernel memory for the inputs');
        blocks.push(ptr);
        return ptr;
    };
    const put = values => {
        const ptr = place(8 * values.length);
        new Float64Array(kernel.memory.buffer, ptr, values.length).set(values);
        return ptr;
    };
    const arena = {
        put,
        // Outputs start as a sentinel no kernel writes, so an entry a kernel
        // leaves unwritten cannot pass for a computed zero.
        out: count => put(new Float64Array(count).fill(7)),
        ints: count => place(4 * count),
        read: (ptr, count) => Array.from(new Float64Array(kernel.memory.buffer, ptr, count)),
        readInts: (ptr, count) => Array.from(new Int32Array(kernel.memory.buffer, ptr, count)),
    };
    try { return body(arena); } finally { for (const ptr of blocks) kernel.free(ptr); }
}

let checks = 0;
function same(actual, expected, what) {
    assert.equal(actual.length, expected.length, `${what}: length`);
    for (let i = 0; i < actual.length; i++) {
        checks++;
        if (!Object.is(actual[i], expected[i])) assert.fail(`${what}[${i}]: ${actual[i]} against ${expected[i]}`);
    }
}

// ── A dispersive stack ───────────────────────────────────────────────────────
// Cauchy indices n = A + B/λ² (B in nm²) with a constant k: two dielectrics, a
// weakly absorbing layer and a thin metal, on a dispersive glass.
const LAMBDAS = [400, 470, 532, 600, 700, 850, 1064];
const THETAS = [0, 12, 25, 38, 52, 64, 78];
const nLam = LAMBDAS.length;
const MATERIALS = [
    { A: 2.2, B: 4.0e4, k: 0 },
    { A: 1.44, B: 4.0e3, k: 0 },
    { A: 2.0, B: 1.0e4, k: 0.05 },
    { A: 0.2, B: 0, k: 3.0 },
];
const THICK = [62, 105, 41, 8];
const N = THICK.length;
const INCIDENT = { A: 1, B: 0, k: 0 };
const SUBSTRATE = { A: 1.52, B: 4.0e3, k: 0 };
const CANDIDATES = [{ A: 1.38, B: 0, k: 0 }, { A: 2.1, B: 2.0e4, k: 0.001 }];
const GROWING = { A: 1.9, B: 1.5e4, k: 0.002 };

const index = (m, lambda) => [m.A + m.B / (lambda * lambda), m.k];
const indexJet = (m, lambda) => {
    const lambdaJet = J.wavelengthOmegaJet(lambda, omegaFromLambdaNm(lambda));
    const inverseSquare = J.jetReciprocal(J.jetMultiply(lambdaJet, lambdaJet));
    return J.jetAdd(J.jetConstant(m.A, m.k), J.jetScale(inverseSquare, m.B));
};
const sinJet = theta => J.jetFromDerivatives(Math.sin(theta * Math.PI / 180), 1e-4, -2e-6, 3e-8);

const perLambda = (m, f) => LAMBDAS.flatMap(lambda => f(m, lambda).flat());
const N0 = perLambda(INCIDENT, index), NS = perLambda(SUBSTRATE, index);
const MAT = MATERIALS.flatMap(m => perLambda(m, index));
const CAND = CANDIDATES.flatMap(m => perLambda(m, index));
const N0_JETS = perLambda(INCIDENT, indexJet), NS_JETS = perLambda(SUBSTRATE, indexJet);
const MAT_JETS = MATERIALS.flatMap(m => perLambda(m, indexJet));
const OMEGAS = LAMBDAS.map(omegaFromLambdaNm);
const SIN_JETS = THETAS.flatMap(theta => sinJet(theta).flat());
const at = (flat, li, width) => flat.slice(li * width, li * width + width);

// Arguments in the kernel's order: the grid, then its media and layers; one
// wavelength's, then its media and layers.
const gridAt = (a, thetas) => [a.put(LAMBDAS), nLam, a.put(thetas), thetas.length];
const gridStack = a => [a.put(N0), a.put(NS), a.put(MAT), a.put(THICK), N];
const gridJetStack = a => [a.put(N0_JETS), a.put(NS_JETS), a.put(MAT_JETS), a.put(THICK), N];
const pointStack = (a, li) => [...at(N0, li, 2), ...at(NS, li, 2),
    a.put(MATERIALS.flatMap((m, k) => [...index(m, LAMBDAS[li]), THICK[k]])), N];
const pointJetStack = (a, li) => [a.put(at(N0_JETS, li, 8)), a.put(at(NS_JETS, li, 8)),
    a.put(MATERIALS.flatMap(m => indexJet(m, LAMBDAS[li]).flat())), a.put(THICK), N];

// Calls a grid export once and its one-wavelength export at every wavelength;
// each output of the grid, wavelength-major, must hold the point's output.
function gridMatchesPoints(what, widths, grid, point) {
    withArena(a => {
        const outs = widths.map(width => a.out(nLam * width));
        grid(a, outs);
        const whole = outs.map((ptr, i) => a.read(ptr, nLam * widths[i]));
        for (let li = 0; li < nLam; li++) {
            const single = widths.map(width => a.out(width));
            point(a, li, single);
            single.forEach((ptr, i) => same(a.read(ptr, widths[i]), at(whole[i], li, widths[i]),
                `${what} λ${li} output ${i}`));
        }
    });
}

// ── R, T, A over a grid, one angle and one per wavelength ───────────────────
for (const thetas of [THETAS, [30]]) {
    withArena(a => {
        const outs = Array.from({ length: 6 }, () => a.out(nLam));
        kernel._tmm_rta_spectrum(...gridAt(a, thetas), ...gridStack(a), ...outs);
        const grid = outs.map(ptr => a.read(ptr, nLam));
        for (let li = 0; li < nLam * 2; li++) {
            const [lambdaIndex, pol] = [li >> 1, li & 1];
            const out = a.out(3);
            kernel._tmm_rta(LAMBDAS[lambdaIndex], thetas[thetas.length > 1 ? lambdaIndex : 0], pol,
                ...pointStack(a, lambdaIndex), out);
            same(a.read(out, 3), grid.slice(3 * pol, 3 * pol + 3).map(q => q[lambdaIndex]),
                `tmm_rta_spectrum θ×${thetas.length} λ${lambdaIndex} pol ${pol}`);
        }
    });
}

// ── Derivative and phase grids against their one-wavelength calls ──────────
// Each over one angle per wavelength and over one angle for the whole grid.
const FRACS = [0, 0.37, 1];
function derivativeAndPhaseGrids(thetas, pol) {
    const thetaAt = li => thetas[thetas.length > 1 ? li : 0];
    const tag = `θ×${thetas.length} pol ${pol}`;
    gridMatchesPoints(`tmm_jacobian_spectrum ${tag}`, [N, N, N, 3],
        (a, outs) => kernel._tmm_jacobian_spectrum(...gridAt(a, thetas), pol, ...gridStack(a), ...outs),
        (a, li, outs) => kernel._tmm_jacobian(LAMBDAS[li], thetaAt(li), pol, ...pointStack(a, li), ...outs));

    gridMatchesPoints(`tmm_hessian_spectrum ${tag}`, [N, N, N, N * N, N * N, N * N, 3],
        (a, outs) => kernel._tmm_hessian_spectrum(...gridAt(a, thetas), pol, ...gridStack(a), ...outs),
        (a, li, outs) => kernel._tmm_hessian(LAMBDAS[li], thetaAt(li), pol, ...pointStack(a, li), ...outs));

    for (const fracs of [FRACS, []]) {
        const nCand = CANDIDATES.length, nFrac = fracs.length;
        gridMatchesPoints(`tmm_needle_scan_spectrum ${tag} nFrac ${nFrac}`,
            [3, (N + 1) * nCand * 3, N * nFrac * nCand * 3],
            (a, outs) => kernel._tmm_needle_scan_spectrum(...gridAt(a, thetas), pol, ...gridStack(a),
                a.put(CAND), nCand, a.put(fracs), nFrac, ...outs),
            (a, li, outs) => kernel._tmm_needle_scan(LAMBDAS[li], thetaAt(li), pol, ...pointStack(a, li),
                a.put(CANDIDATES.flatMap(m => index(m, LAMBDAS[li]))), nCand, a.put(fracs), nFrac, ...outs));
    }

    // The incident sine as a jet replaces the angle when it is given.
    for (const sines of [null, SIN_JETS]) {
        const sineAt = (a, li) => sines ? a.put(at(sines, li, 8)) : 0;
        const label = `${tag}${sines ? ' sine jet' : ''}`;
        gridMatchesPoints(`tmm_phase_spectrum ${label}`, [10],
            (a, outs) => kernel._tmm_phase_spectrum(...gridAt(a, thetas), pol, ...gridJetStack(a),
                a.put(OMEGAS), sines ? a.put(sines) : 0, ...outs),
            (a, li, outs) => kernel._tmm_phase(LAMBDAS[li], thetaAt(li), pol, ...pointJetStack(a, li),
                OMEGAS[li], sineAt(a, li), ...outs));

        gridMatchesPoints(`tmm_phase_jacobian_spectrum ${label}`, [10, 10 * N],
            (a, outs) => kernel._tmm_phase_jacobian_spectrum(...gridAt(a, thetas), pol, ...gridJetStack(a),
                a.put(OMEGAS), sines ? a.put(sines) : 0, ...outs),
            (a, li, outs) => kernel._tmm_phase_jacobian(LAMBDAS[li], thetaAt(li), pol, ...pointJetStack(a, li),
                OMEGAS[li], sineAt(a, li), ...outs));
    }
}
for (const thetas of [THETAS, [30]]) {
    for (const pol of [0, 1]) derivativeAndPhaseGrids(thetas, pol);
}

// ── Deposition spectra and the growing stack, one angle per wavelength ──────
// The grid against one-wavelength grids: outputs are [step][λ] for the
// deposition spectra and one value per wavelength for a growing-stack sample.
const pointGrid = (a, li) => [a.put([LAMBDAS[li]]), 1, a.put([THETAS[li]]), 1,
    a.put(at(N0, li, 2)), a.put(at(NS, li, 2)),
    a.put(MATERIALS.flatMap(m => index(m, LAMBDAS[li]))), a.put(THICK), N];

withArena(a => {
    const outs = Array.from({ length: 6 }, () => a.out(N * nLam));
    assert.equal(kernel._tmm_deposition_spectra(...gridAt(a, THETAS), ...gridStack(a), ...outs), 1);
    const grid = outs.map(ptr => a.read(ptr, N * nLam));
    for (let li = 0; li < nLam; li++) {
        const single = Array.from({ length: 6 }, () => a.out(N));
        assert.equal(kernel._tmm_deposition_spectra(...pointGrid(a, li), ...single), 1);
        single.forEach((ptr, q) => same(a.read(ptr, N), grid[q].filter((_, i) => i % nLam === li),
            `tmm_deposition_spectra λ${li} output ${q}`));
    }
});

withArena(a => {
    const sample = (handle, count, d) => {
        const outs = Array.from({ length: 6 }, () => a.out(count));
        assert.equal(kernel._tmm_growing_stack_sample(handle, d, ...outs), 1);
        return outs.map(ptr => a.read(ptr, count));
    };
    const grid = kernel._tmm_growing_stack_create(...gridAt(a, THETAS), ...gridStack(a));
    assert.ok(grid, 'growing stack over the grid');
    kernel._tmm_growing_stack_set_top(grid, a.put(perLambda(GROWING, index)));
    for (let li = 0; li < nLam; li++) {
        const single = kernel._tmm_growing_stack_create(...pointGrid(a, li));
        assert.ok(single, 'growing stack at one wavelength');
        kernel._tmm_growing_stack_set_top(single, a.put(index(GROWING, LAMBDAS[li])));
        for (const d of [0, 37.5]) {
            const whole = sample(grid, nLam, d);
            sample(single, 1, d).forEach((values, q) => same(values, [whole[q][li]],
                `tmm_growing_stack λ${li} d ${d} output ${q}`));
        }
        kernel._tmm_growing_stack_free(single);
    }
    kernel._tmm_growing_stack_free(grid);
});

// ── Jet arithmetic over a grid ───────────────────────────────────────────────
let seed = 20261029;
const random = () => {
    seed = (seed * 1664525 + 1013904223) >>> 0;
    return seed / 2 ** 32;
};
const between = (low, high) => low + (high - low) * random();
const COUNT = 48;
const randomJets = (reLow, reHigh) => Array.from({ length: COUNT }, () =>
    Array.from({ length: 4 }, (_, order) => order === 0
        ? [between(reLow, reHigh), between(-1, 1)]
        : [between(-1, 1), between(-1, 1)]));
const numbers = (low, high) => Array.from({ length: COUNT }, () => between(low, high));

const anyJets = randomJets(-2, 2);
const otherJets = randomJets(-2, 2);
const rightHalfJets = randomJets(0.2, 2);     // away from the branch cut of log and root
const awayFromZero = randomJets(0.5, 2);
// Every fifth jet NaN, infinite or exactly at the clamp minimum used below.
const ODD_VALUES = [null, [0, 0, NaN], [1, 1, Infinity], [0, 0, 0.25], null];
const unusualJets = anyJets.map((jet, i) => {
    const copy = jet.map(pair => [...pair]);
    const odd = ODD_VALUES[i % 5];
    if (odd) copy[odd[0]][odd[1]] = odd[2];
    return copy;
});

const jetsInput = jets => ({ flat: jets.flatMap(jet => jet.flat()), item: i => jets[i] });
const numbersInput = values => ({ flat: values, item: i => values[i] });
const triplesInput = triples => ({ flat: triples.flat(), item: i => triples[i] });

// Each case: the C export, its inputs, the JS reference for one item, and
// whether the two must agree bit for bit.
const JET_CASES = [
    ['tmm_jet_constant', [numbersInput(numbers(-3, 3)), numbersInput(numbers(-1, 1))], J.jetConstant, true],
    ['tmm_jet_from_derivatives', [jetsInput(anyJets)], d => J.jetFromDerivatives(...d), true],
    ['tmm_jet_derivatives', [jetsInput(anyJets)], J.jetDerivatives, true],
    ['tmm_jet_add', [jetsInput(anyJets), jetsInput(otherJets)], J.jetAdd, true],
    ['tmm_jet_subtract', [jetsInput(anyJets), jetsInput(otherJets)], J.jetSubtract, true],
    ['tmm_jet_scale', [jetsInput(anyJets), numbersInput(numbers(-4, 4))], J.jetScale, true],
    ['tmm_jet_multiply', [jetsInput(anyJets), jetsInput(otherJets)], J.jetMultiply, true],
    ['tmm_jet_reciprocal', [jetsInput(awayFromZero)], J.jetReciprocal, true],
    ['tmm_jet_divide', [jetsInput(anyJets), jetsInput(awayFromZero)], J.jetDivide, true],
    ['tmm_jet_sqrt', [jetsInput(rightHalfJets)], J.jetSqrt, false],
    ['tmm_jet_exp', [jetsInput(anyJets)], J.jetExp, false],
    ['tmm_jet_log', [jetsInput(rightHalfJets)], J.jetLog, false],
    ['tmm_jet_power', [jetsInput(rightHalfJets), numbersInput(numbers(-2.5, 2.5))], J.jetPower, false],
    ['tmm_jet_sin', [jetsInput(anyJets)], J.jetSin, false],
    ['tmm_jet_cos', [jetsInput(anyJets)], J.jetCos, false],
    ['tmm_jet_compose', [numbersInput(numbers(-2, 2)),
        triplesInput(Array.from({ length: COUNT }, () => [between(-3, 3), between(-3, 3), between(-3, 3)])),
        jetsInput(anyJets)],
        J.jetCompose, true],
    ['tmm_wavelength_omega_jet', [numbersInput(numbers(200, 2500)), numbersInput(numbers(0.5, 9))],
        J.wavelengthOmegaJet, true],
    ['tmm_jet_with_imaginary_part', [jetsInput(anyJets), jetsInput(otherJets)], J.jetWithImaginaryPart, true],
    ['tmm_jet_real_part', [jetsInput(anyJets)], J.jetRealPart, true],
    ['tmm_jet_clamp_real_minimum', [jetsInput(unusualJets), numbersInput(Array(COUNT).fill(0.25))],
        J.jetClampRealMinimum, true],
    ['tmm_jet_clamp_imaginary', [jetsInput(anyJets), numbersInput(numbers(0.1, 1.2))], J.jetClampImaginary, true],
];

let worstLibm = 0;
function compareJet(actual, expected, what, exact) {
    if (exact) return same(actual, expected, what);
    const scale = Math.max(...expected.map(Math.abs));
    actual.forEach((value, i) => {
        checks++;
        const relative = Math.abs(value - expected[i]) / scale;
        worstLibm = Math.max(worstLibm, relative);
        assert.ok(relative <= LIBM_TOL, `${what}[${i}]: ${value} against ${expected[i]}`);
    });
}

for (const [name, inputs, reference, exact] of JET_CASES) {
    withArena(a => {
        const out = a.out(8 * COUNT);
        kernel[`_${name}`](...inputs.map(input => a.put(input.flat)), COUNT, out);
        const result = a.read(out, 8 * COUNT);
        for (let i = 0; i < COUNT; i++) {
            compareJet(at(result, i, 8), reference(...inputs.map(input => input.item(i))).flat(),
                `${name} jet ${i}`, exact);
        }
    });
}

withArena(a => {
    const sine = a.out(8 * COUNT), cosine = a.out(8 * COUNT);
    kernel._tmm_jet_sin_cos(a.put(jetsInput(anyJets).flat), COUNT, sine, cosine);
    const [s, c] = [a.read(sine, 8 * COUNT), a.read(cosine, 8 * COUNT)];
    anyJets.forEach((jet, i) => {
        const both = J.jetSinCos(jet);
        compareJet(at(s, i, 8), both.sine.flat(), `tmm_jet_sin_cos sine ${i}`, false);
        compareJet(at(c, i, 8), both.cosine.flat(), `tmm_jet_sin_cos cosine ${i}`, false);
    });
});

withArena(a => {
    const out = a.ints(COUNT);
    kernel._tmm_jet_is_finite(a.put(jetsInput(unusualJets).flat), COUNT, out);
    same(a.readInts(out, COUNT), unusualJets.map(jet => Number(J.jetIsFinite(jet))), 'tmm_jet_is_finite');
});

withArena(a => {
    const out = a.out(nLam);
    kernel._tmm_omega_from_lambda_nm(a.put(LAMBDAS), nLam, out);
    same(a.read(out, nLam), OMEGAS, 'tmm_omega_from_lambda_nm');
});

// ── The jet calls of an instance ─────────────────────────────────────────────
// Each call with the arguments of its namesake in jet.js, short argument lists
// and real numbers in place of complex ones included.
const A = [[1.7, 0.02], [-0.03, 0.001], [0.004, -0.0002], [-0.0005, 0.00001]];
const B = [[2.3, -0.1], [0.01, 0.02], [0.003, 0.0001], [0.0002, -0.00003]];
const NAN_JET = [[NaN, 0], [0, 0], [0, 0], [0, 0]];
const INSTANCE_CASES = [
    ['jetConstant', [1.5, 0.2]], ['jetConstant', [1.5]],
    ['jetFromDerivatives', [[1.5, 0.1], 0.2, [0.3, -0.1], 4]], ['jetFromDerivatives', [1.5]],
    ['jetDerivatives', [A]], ['jetAdd', [A, B]], ['jetSubtract', [A, B]], ['jetScale', [A, 2.5]],
    ['jetMultiply', [A, B]], ['jetReciprocal', [A]], ['jetDivide', [A, B]],
    ['jetSqrt', [A]], ['jetExp', [A]], ['jetLog', [A]], ['jetPower', [A, -1.7]],
    ['jetSinCos', [A]], ['jetSin', [A]], ['jetCos', [A]],
    ['jetCompose', [2.1, [0.5, -0.2, 0.05], A]], ['jetCompose', [2.1, [0.5], A]],
    ['wavelengthOmegaJet', [550, 3.42]],
    ['jetWithImaginaryPart', [A, B]], ['jetRealPart', [A]],
    ['jetClampRealMinimum', [A, 2]], ['jetClampRealMinimum', [A, 1]], ['jetClampImaginary', [A, 0.01]],
    ['jetIsFinite', [A]], ['jetIsFinite', [NAN_JET]],
    ['omegaFromLambdaNm', [550]],
];
const THROUGH_LIBM = new Set(['jetSqrt', 'jetExp', 'jetLog', 'jetPower', 'jetSinCos', 'jetSin', 'jetCos']);
const REFERENCE = { ...J, omegaFromLambdaNm };

// A jet, a sine and cosine pair, a flag or a number, as a list of values.
const valuesOf = result => (typeof result === 'object' ? Object.values(result).flat(2) : [result]);
for (const [name, args] of INSTANCE_CASES) {
    compareJet(valuesOf(kernel[name](...args)), valuesOf(REFERENCE[name](...args)),
        `${name}(${JSON.stringify(args)}) on an instance`, !THROUGH_LIBM.has(name));
}
const checked = new Set(INSTANCE_CASES.map(([name]) => name));
for (const name of Object.keys(REFERENCE)) {
    if (typeof kernel[name] === 'function') assert.ok(checked.has(name), `${name} on an instance is checked`);
}

const JETS = [A, B, NAN_JET];
const typedTerms = arg => (JETS.includes(arg) ? arg.map(term => Float64Array.from(term)) : arg);
for (const [name, args] of INSTANCE_CASES) {
    if (!args.some(arg => JETS.includes(arg))) continue;
    same(valuesOf(kernel[name](...args.map(typedTerms))), valuesOf(kernel[name](...args)),
        `${name} with typed terms on an instance`);
}
for (const notAJet of [A.slice(0, 3), [1.7, 0.02], A.map(term => term[0]), null]) {
    assert.throws(() => kernel.jetSqrt(notAJet), TypeError, `jetSqrt(${JSON.stringify(notAJet)})`);
    assert.throws(() => kernel.jetAdd(A, notAJet), TypeError, `jetAdd(A, ${JSON.stringify(notAJet)})`);
}

console.log(`${checks} values checked; worst jet difference through the math library ${worstLibm.toExponential(2)}`);
console.log('PASS : kernel exports');
