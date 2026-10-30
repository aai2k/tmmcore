/**
 * The calls of the package under one set of names. Each takes a stack, a
 * wavelength in nm or a list of them, and options { angle, pol }:
 *
 *   rta, rtaSpectrum                       R, T and A
 *   jacobian, jacobianSpectrum             with their thickness derivatives
 *   hessian, hessianSpectrum               and the second derivatives
 *   phase, phaseSpectrum                   phase, GD, GDD and TOD of r and t
 *   phaseJacobian, phaseJacobianSpectrum   with their thickness derivatives
 *
 * The stack is { incident, substrate, layers: [{ n, k, d }] } (stack.js).
 * `angle` is the angle of incidence in degrees, 0 when left out; a Spectrum
 * call also takes one angle per wavelength. `pol` is 's', 'p' or 'avg' in the
 * R/T and derivative calls, 'avg' when left out, which is the mean of the s and
 * p results; 's' or 'p' in the phase calls, where it must be given. An option
 * a call does not take throws.
 *
 * stackCalls() builds the calls on a set of kernels: the package index on the
 * JavaScript ones, each WebAssembly instance on its own. The two take the same
 * arguments and return the same shapes, every list of numbers a Float64Array.
 */

import {
    readWavelength, readWavelengths, readAngle, readRtaPol, readPhasePol, pointJetStack, gridJetStack,
} from './stack.js';

const OPTIONS = ['angle', 'pol'];
const FIRST = ['dRdd', 'dTdd', 'dAdd'];
const SECOND = ['d2Rdd', 'd2Tdd', 'd2Add'];
const PHASE_DERIVATIVES = ['dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared'];

// The angle and polarization of a call; nLam is the grid size of a Spectrum
// call, null at one wavelength.
function readOptions(options, nLam, readPol) {
    if (options === null || typeof options !== 'object') throw new TypeError('options must be { angle, pol }');
    const unknown = Object.keys(options).find(key => !OPTIONS.includes(key));
    if (unknown !== undefined) throw new TypeError(`unknown option '${unknown}'; the options are angle and pol`);
    return { angle: readAngle(options.angle, nLam), pol: readPol(options.pol) };
}

// The mean of the s and p results of one call, entry by entry, in the shape
// of either.
function mean(s, p) {
    if (typeof s === 'number') return (s + p) / 2;
    if (ArrayBuffer.isView(s)) {
        const out = new Float64Array(s.length);
        for (let i = 0; i < s.length; i++) out[i] = (s[i] + p[i]) / 2;
        return out;
    }
    if (Array.isArray(s)) return s.map((value, i) => mean(value, p[i]));
    return Object.fromEntries(Object.keys(s).map(key => [key, mean(s[key], p[key])]));
}

const polarized = (pol, run) => (pol === 'avg' ? mean(run('s'), run('p')) : run(pol));

// Derivatives at one wavelength as Float64Arrays: the JavaScript kernels build
// plain arrays; the WebAssembly ones return copies out of kernel memory, kept
// as they are.
const asFloat64 = values => (values instanceof Float64Array ? values : Float64Array.from(values));

// A derivative result in one key order on module and instance: R, T, A, the
// first derivatives, the second when there are any, then N, with the
// derivatives at one wavelength as Float64Arrays.
function derivativeResult(result, atOneWavelength) {
    const typed = values => (atOneWavelength ? asFloat64(values) : values);
    const out = { R: result.R, T: result.T, A: result.A };
    for (const key of FIRST) out[key] = typed(result[key]);
    if (result.d2Rdd) for (const key of SECOND) out[key] = atOneWavelength ? result[key].map(typed) : result[key];
    out.N = result.N;
    return out;
}

function typedPhaseSide(side) {
    if (!side) return side;
    const out = { ...side };
    for (const key of PHASE_DERIVATIVES) out[key] = asFloat64(side[key]);
    return out;
}

/**
 * The calls on `kernels`, with pol 's' or 'p':
 *   rta, jacobian, hessian     (λ, θ, pol, kernels.pointStack(stack))
 *   phase, phaseJacobian       (λ, θ, pol, { n0Jet, nsJet, layers: [{ nJet, d }] })
 *   rtaSpectrum                (grid), s and p together
 *   jacobianSpectrum, hessianSpectrum, phaseSpectrum, phaseJacobianSpectrum   (grid, pol)
 * kernels.pointStack(stack) and kernels.gridStack(stack, nLam) read the stack
 * of an R/T call the way those kernels take it: stack.js's pointStack and
 * gridStack for the JavaScript ones. The phase calls read theirs with
 * pointJetStack and gridJetStack. The grid of a Spectrum call is
 * { lambdas, theta_deg } together with what kernels.gridStack, or
 * gridJetStack for the phase calls, returns.
 */
export function stackCalls(kernels) {
    const atOne = (readPol, readMedia, run) => (stack, wavelength, options = {}) => {
        const lambda = readWavelength(wavelength);
        const { angle, pol } = readOptions(options, null, readPol);
        const media = readMedia(stack);
        return polarized(pol, side => run(lambda, angle, side, media));
    };
    const gridOf = (stack, wavelengths, options, readPol, readMedia) => {
        const lambdas = readWavelengths(wavelengths);
        const { angle, pol } = readOptions(options, lambdas.length, readPol);
        return { pol, grid: { lambdas, theta_deg: angle, ...readMedia(stack, lambdas.length) } };
    };
    const overGrid = (readPol, readMedia, run) => (stack, wavelengths, options = {}) => {
        const { pol, grid } = gridOf(stack, wavelengths, options, readPol, readMedia);
        return polarized(pol, side => run(grid, side));
    };

    // R, T and A of the polarization asked, and of s and p.
    const rtaSpectrum = (stack, wavelengths, options = {}) => {
        const { pol, grid } = gridOf(stack, wavelengths, options, readRtaPol, kernels.gridStack);
        const both = kernels.rtaSpectrum(grid);
        const pick = (s, p) => (pol === 'avg' ? mean(s, p) : (pol === 's' ? s : p).slice());
        return { R: pick(both.Rs, both.Rp), T: pick(both.Ts, both.Tp), A: pick(both.As, both.Ap), ...both };
    };

    return {
        rta: atOne(readRtaPol, kernels.pointStack, kernels.rta),
        rtaSpectrum,
        jacobian: atOne(readRtaPol, kernels.pointStack,
            (...args) => derivativeResult(kernels.jacobian(...args), true)),
        jacobianSpectrum: overGrid(readRtaPol, kernels.gridStack,
            (...args) => derivativeResult(kernels.jacobianSpectrum(...args), false)),
        hessian: atOne(readRtaPol, kernels.pointStack,
            (...args) => derivativeResult(kernels.hessian(...args), true)),
        hessianSpectrum: overGrid(readRtaPol, kernels.gridStack,
            (...args) => derivativeResult(kernels.hessianSpectrum(...args), false)),
        phase: atOne(readPhasePol, pointJetStack, kernels.phase),
        phaseSpectrum: overGrid(readPhasePol, gridJetStack, kernels.phaseSpectrum),
        phaseJacobian: atOne(readPhasePol, pointJetStack, (...args) => {
            const { r, t } = kernels.phaseJacobian(...args);
            return { r: typedPhaseSide(r), t: typedPhaseSide(t) };
        }),
        phaseJacobianSpectrum: overGrid(readPhasePol, gridJetStack, kernels.phaseJacobianSpectrum),
    };
}
