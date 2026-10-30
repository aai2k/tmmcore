/**
 * The calls of the package under one set of names. Each takes a stack, a
 * wavelength in nm or a list of them, and options:
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
 * p results; 's' or 'p' in the phase calls, which may leave it out at normal
 * incidence, where the two are the same. In place of `angle` a phase call
 * takes `sine`, the sine of the angle of incidence as a jet in ω (jet.js), one
 * jet per wavelength in a Spectrum call: for an angle held fixed outside a
 * dispersive incident medium, so that the angle inside it disperses. An option
 * a call does not take throws.
 *
 * stackCalls() builds the calls on a set of kernels: the package index on the
 * JavaScript ones, each WebAssembly instance on its own. The two take the same
 * arguments and return the same shapes, every list of numbers a Float64Array.
 */

import {
    readWavelength, readWavelengths, readAngle, readRtaPol, readPhasePol, readSine, atNormalIncidence,
    pointJetStack, gridJetStack,
} from './stack.js';

const RTA_OPTIONS = ['angle', 'pol'];
const PHASE_OPTIONS = ['angle', 'sine', 'pol'];
const FIRST = ['dRdd', 'dTdd', 'dAdd'];
const SECOND = ['d2Rdd', 'd2Tdd', 'd2Add'];
const PHASE_DERIVATIVES = ['dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared'];

// 'a and b', 'a, b and c'.
const listed = keys => `${keys.slice(0, -1).join(', ')} and ${keys[keys.length - 1]}`;

function checkOptions(options, keys) {
    if (options === null || typeof options !== 'object') throw new TypeError(`options must be { ${keys.join(', ')} }`);
    const unknown = Object.keys(options).find(key => !keys.includes(key));
    if (unknown !== undefined) throw new TypeError(`unknown option '${unknown}'; the options are ${listed(keys)}`);
}

// The angle and polarization of an R/T call; nLam is the grid size of a
// Spectrum call, null at one wavelength.
function readRtaOptions(options, nLam) {
    checkOptions(options, RTA_OPTIONS);
    return { angle: readAngle(options.angle, nLam), pol: readRtaPol(options.pol) };
}

// The angle, or the incident sine as jets, and the polarization of a phase
// call. With the sine the kernels take no angle, and are given 0.
function readPhaseOptions(options, nLam) {
    checkOptions(options, PHASE_OPTIONS);
    if (options.sine === undefined) {
        const angle = readAngle(options.angle, nLam);
        return { angle, sine: null, pol: readPhasePol(options.pol, atNormalIncidence(angle, null, nLam)) };
    }
    if (options.angle !== undefined) throw new TypeError('a phase call takes angle or sine, not both');
    const sine = readSine(options.sine, nLam);
    return { angle: 0, sine, pol: readPhasePol(options.pol, atNormalIncidence(0, sine, nLam)) };
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
 *   phase, phaseJacobian       (λ, θ, pol, { n0Jet, nsJet, layers: [{ nJet, d }] }, sine)
 *   rtaSpectrum                (grid), s and p together
 *   jacobianSpectrum, hessianSpectrum, phaseSpectrum, phaseJacobianSpectrum   (grid, pol)
 * kernels.pointStack(stack) and kernels.gridStack(stack, nLam) read the stack
 * of an R/T call the way those kernels take it: stack.js's pointStack and
 * gridStack for the JavaScript ones. The phase calls read theirs with
 * pointJetStack and gridJetStack, and give the incident sine as a jet, or
 * null for the angle's own. The grid of a Spectrum call is
 * { lambdas, theta_deg } together with what kernels.gridStack, or
 * gridJetStack for the phase calls, returns, and sinJets, the sine per
 * wavelength, when a phase call is given one.
 */
export function stackCalls(kernels) {
    const atOne = (read, readMedia, run) => (stack, wavelength, options = {}) => {
        const lambda = readWavelength(wavelength);
        const { angle, sine, pol } = read(options, null);
        const media = readMedia(stack);
        return polarized(pol, side => run(lambda, angle, side, media, sine));
    };
    const gridOf = (stack, wavelengths, options, read, readMedia) => {
        const lambdas = readWavelengths(wavelengths);
        const { angle, sine, pol } = read(options, lambdas.length);
        const grid = { lambdas, theta_deg: angle, ...readMedia(stack, lambdas.length) };
        if (sine) grid.sinJets = sine;
        return { pol, grid };
    };
    const overGrid = (read, readMedia, run) => (stack, wavelengths, options = {}) => {
        const { pol, grid } = gridOf(stack, wavelengths, options, read, readMedia);
        return polarized(pol, side => run(grid, side));
    };

    // R, T and A of the polarization asked, and of s and p.
    const rtaSpectrum = (stack, wavelengths, options = {}) => {
        const { pol, grid } = gridOf(stack, wavelengths, options, readRtaOptions, kernels.gridStack);
        const both = kernels.rtaSpectrum(grid);
        const pick = (s, p) => (pol === 'avg' ? mean(s, p) : (pol === 's' ? s : p).slice());
        return { R: pick(both.Rs, both.Rp), T: pick(both.Ts, both.Tp), A: pick(both.As, both.Ap), ...both };
    };

    return {
        rta: atOne(readRtaOptions, kernels.pointStack, kernels.rta),
        rtaSpectrum,
        jacobian: atOne(readRtaOptions, kernels.pointStack,
            (...args) => derivativeResult(kernels.jacobian(...args), true)),
        jacobianSpectrum: overGrid(readRtaOptions, kernels.gridStack,
            (...args) => derivativeResult(kernels.jacobianSpectrum(...args), false)),
        hessian: atOne(readRtaOptions, kernels.pointStack,
            (...args) => derivativeResult(kernels.hessian(...args), true)),
        hessianSpectrum: overGrid(readRtaOptions, kernels.gridStack,
            (...args) => derivativeResult(kernels.hessianSpectrum(...args), false)),
        phase: atOne(readPhaseOptions, pointJetStack, kernels.phase),
        phaseSpectrum: overGrid(readPhaseOptions, gridJetStack, kernels.phaseSpectrum),
        phaseJacobian: atOne(readPhaseOptions, pointJetStack, (...args) => {
            const { r, t } = kernels.phaseJacobian(...args);
            return { r: typedPhaseSide(r), t: typedPhaseSide(t) };
        }),
        phaseJacobianSpectrum: overGrid(readPhaseOptions, gridJetStack, kernels.phaseJacobianSpectrum),
    };
}
