/**
 * The jet arithmetic against closed-form derivatives.
 *
 * The phase calls take every GD, GDD and TOD from jets, so each operation of
 * src/js/jet.js is checked here on its own: fed the jet of the variable at a
 * point z0, an operation must return f(z0), f'(z0), f''(z0) and f'''(z0) of
 * the function it stands for, written out by hand below. The points are
 * complex, as refractive indices are, with one real point for the calls that
 * take real arguments. The WebAssembly jets are held to these by
 * tests/kernel_exports.mjs.
 *
 * Also the unit helpers: C_NM_PER_FS is the SI speed of light in nm/fs, and
 * omegaFromLambdaNm(λ) = 2πc/λ.
 *
 *   node tests/jets.mjs
 */
import assert from 'node:assert/strict';
import {
    JET_ORDER, jetConstant, jetFromDerivatives, jetDerivatives, jetAdd, jetSubtract, jetScale,
    jetMultiply, jetReciprocal, jetDivide, jetSqrt, jetExp, jetLog, jetPower, jetSinCos, jetSin, jetCos,
    jetCompose, wavelengthOmegaJet, jetWithImaginaryPart, jetRealPart, jetClampRealMinimum,
    jetClampImaginary, jetIsFinite, C_NM_PER_FS, omegaFromLambdaNm,
} from '../src/js/index.js';

// Complex arithmetic of this test's own, so the reference shares no code
// with the arithmetic under test.
const add = (a, b) => [a[0] + b[0], a[1] + b[1]];
const mul = (a, b) => [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]];
const scale = (a, s) => [a[0] * s, a[1] * s];
const div = (a, b) => { const q = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / q, (a[1] * b[0] - a[0] * b[1]) / q]; };
const exp = z => [Math.exp(z[0]) * Math.cos(z[1]), Math.exp(z[0]) * Math.sin(z[1])];
const log = z => [Math.log(Math.hypot(z[0], z[1])), Math.atan2(z[1], z[0])];
const pow = (z, p) => exp(scale(log(z), p));
const sin = z => [Math.sin(z[0]) * Math.cosh(z[1]), Math.cos(z[0]) * Math.sinh(z[1])];
const cos = z => [Math.cos(z[0]) * Math.cosh(z[1]), -Math.sin(z[0]) * Math.sinh(z[1])];
const neg = z => scale(z, -1);

const TOLERANCE = 1e-13;
let checks = 0;
function derivativesAre(jet, expected, label) {
    jetDerivatives(jet).forEach((got, order) => {
        const want = Array.isArray(expected[order]) ? expected[order] : [expected[order], 0];
        const size = Math.max(1, Math.hypot(want[0], want[1]));
        const off = Math.hypot(got[0] - want[0], got[1] - want[1]) / size;
        checks++;
        assert.ok(off <= TOLERANCE, `${label}, derivative ${order}: got ${got}, expected ${want}`);
    });
}

assert.equal(JET_ORDER, 3, 'jets carry three derivatives');

for (const z0 of [[1.3, 0.4], [2.1, -0.7]]) {
    const x = jetFromDerivatives(z0, 1);
    const at = `at ${z0[0]}${z0[1] < 0 ? '' : '+'}${z0[1]}i`;
    const z2 = mul(z0, z0), z3 = mul(z2, z0), z4 = mul(z3, z0);

    derivativesAre(jetSubtract(jetAdd(jetScale(x, 3), jetConstant(1)), jetConstant(0.5, 2)),
        [add(scale(z0, 3), [0.5, -2]), 3, 0, 0], `3z + 1 − (0.5 + 2i) ${at}`);
    derivativesAre(jetMultiply(x, jetMultiply(x, x)), [z3, scale(z2, 3), scale(z0, 6), 6], `z³ ${at}`);
    derivativesAre(jetReciprocal(x), [div([1, 0], z0), neg(div([1, 0], z2)), div([2, 0], z3), div([-6, 0], z4)], `1/z ${at}`);
    {
        // (z + 2)/(z − 0.5) = 1 + 2.5/(z − 0.5)
        const w = add(z0, [-0.5, 0]);
        const w2 = mul(w, w), w3 = mul(w2, w), w4 = mul(w3, w);
        derivativesAre(jetDivide(jetAdd(x, jetConstant(2)), jetSubtract(x, jetConstant(0.5))),
            [add([1, 0], div([2.5, 0], w)), div([-2.5, 0], w2), div([5, 0], w3), div([-15, 0], w4)],
            `(z + 2)/(z − 0.5) ${at}`);
    }
    derivativesAre(jetSqrt(x), [pow(z0, 0.5), scale(pow(z0, -0.5), 0.5), scale(pow(z0, -1.5), -0.25),
        scale(pow(z0, -2.5), 0.375)], `√z ${at}`);
    derivativesAre(jetExp(x), [exp(z0), exp(z0), exp(z0), exp(z0)], `exp z ${at}`);
    derivativesAre(jetLog(x), [log(z0), div([1, 0], z0), neg(div([1, 0], z2)), div([2, 0], z3)], `log z ${at}`);
    {
        const p = 2.5;
        derivativesAre(jetPower(x, p), [pow(z0, p), scale(pow(z0, p - 1), p), scale(pow(z0, p - 2), p * (p - 1)),
            scale(pow(z0, p - 3), p * (p - 1) * (p - 2))], `z^2.5 ${at}`);
    }
    const sinCos = jetSinCos(x);
    derivativesAre(sinCos.sine, [sin(z0), cos(z0), neg(sin(z0)), neg(cos(z0))], `sin z ${at}`);
    derivativesAre(sinCos.cosine, [cos(z0), neg(sin(z0)), neg(cos(z0)), sin(z0)], `cos z ${at}`);
    derivativesAre(jetSin(x), [sin(z0), cos(z0), neg(sin(z0)), neg(cos(z0))], `jetSin ${at}`);
    derivativesAre(jetCos(x), [cos(z0), neg(sin(z0)), neg(cos(z0)), sin(z0)], `jetCos ${at}`);
}

// The chain rule, real: f(y) = y³ given by its derivatives at y0 = e^x0, so
// f(e^x) = e^{3x}, whose n-th derivative is 3^n e^{3x}.
{
    const x0 = 0.7;
    const y0 = Math.exp(x0);
    const inner = jetExp(jetFromDerivatives(x0, 1));
    const e3 = Math.exp(3 * x0);
    derivativesAre(jetCompose(y0 ** 3, [3 * y0 ** 2, 6 * y0, 6], inner), [e3, 3 * e3, 9 * e3, 27 * e3], 'y³ of e^x');
}

// λ(ω) = 2πc/ω: λ' = −λ/ω, λ'' = 2λ/ω², λ''' = −6λ/ω³.
{
    const lambda = 632.8;
    const omega = omegaFromLambdaNm(lambda);
    derivativesAre(wavelengthOmegaJet(lambda, omega),
        [lambda, -lambda / omega, 2 * lambda / omega ** 2, -6 * lambda / omega ** 3], 'λ(ω)');
    assert.equal(C_NM_PER_FS, 299.792458, 'c in nm/fs, exact in SI');
    assert.ok(Math.abs(omega - 2 * Math.PI * 299.792458 / lambda) <= 1e-15 * omega, 'ω = 2πc/λ');
}

// The helpers that move parts between jets, and the guards.
{
    const real = jetFromDerivatives(1.5, 2, 3, 4);
    const imaginary = jetFromDerivatives(0.25, -1, 0.5, 6);
    const both = jetWithImaginaryPart(real, imaginary);
    derivativesAre(both, [[1.5, 0.25], [2, -1], [3, 0.5], [4, 6]], 'real and imaginary jets combined');
    derivativesAre(jetRealPart(both), [1.5, 2, 3, 4], 'the real part of that');
    assert.deepEqual(jetClampRealMinimum(real, 2), jetConstant(2), 'a floor that bites flattens the jet');
    assert.equal(jetClampRealMinimum(real, 1), real, 'one that does not leaves it');
    const held = jetClampImaginary(jetFromDerivatives([1, 80], [2, 3]), 50);
    assert.deepEqual(held, [[1, 50], [2, 0], [0, 0], [0, 0]], 'an imaginary part past the limit is held there');
    assert.ok(jetIsFinite(both), 'a finite jet is finite');
    assert.ok(!jetIsFinite(jetFromDerivatives(1, Number.NaN)), 'a NaN anywhere is not');
}

console.log(`${checks} derivatives against closed forms`);
console.log('PASS : jet arithmetic');
