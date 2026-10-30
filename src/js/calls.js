/**
 * The calls of the package under one set of names. Each takes a stack, a
 * wavelength in nm or a list of them, and options:
 *
 *   rta, rtaSpectrum                       R, T and A
 *   jacobian, jacobianSpectrum             with their thickness derivatives
 *   hessian, hessianSpectrum               and the second derivatives
 *   needleScan, needleScanSpectrum         their derivatives for a thin layer inserted
 *   phase, phaseSpectrum                   phase, GD, GDD and TOD of r and t
 *   phaseJacobian, phaseJacobianSpectrum   with their thickness derivatives
 *   monitorCurve                           a layer growing on the stack, at one wavelength
 *   depositionSpectra                      the spectra after each layer is deposited
 *   growingStack                           a layer growing on the stack, over a grid
 *
 * The stack is { incident, substrate, layers: [{ n, k, d }] } (stack.js).
 * `angle` is the angle of incidence in degrees, 0 when left out; a Spectrum
 * call also takes one angle per wavelength. `pol` is 's', 'p' or 'avg' in the
 * R/T, derivative and needle calls, 'avg' when left out, which is the mean of
 * the s and p results; 's' or 'p' in the phase calls, which may leave it out at
 * normal incidence, where the two are the same. In place of `angle` a phase
 * call takes `sine`, the sine of the angle of incidence as a jet in ω (jet.js),
 * one jet per wavelength in a Spectrum call: for an angle held fixed outside a
 * dispersive incident medium, so that the angle inside it disperses. A needle
 * scan takes `candidates`, the materials it tries inserting, each an index
 * given as a medium is, a number or { n, k }, and `fractions`, where inside each
 * layer it also tries them. An option a call does not take throws.
 *
 * The monitoring calls take `angle` only, and return s and p together, with
 * the reflectance of the coated surface seen from the substrate side:
 * { Rs, Ts, Rp, Tp, Rrs, Rrp }. monitorCurve also takes `top`, the growing
 * layer's index, given as a medium is, and `thicknesses`, the thicknesses of
 * it to sample. depositionSpectra deposits the layers from the substrate
 * outward, the last one listed first.
 *
 * stackCalls() builds the calls on a set of kernels: the package index on the
 * JavaScript ones, each WebAssembly instance on its own. The two take the same
 * arguments and return the same shapes, every list of numbers a Float64Array.
 */

import {
    readWavelength, readWavelengths, readAngle, readRtaPol, readPhasePol, readSine, atNormalIncidence, readFractions,
    readThicknesses, readThickness,
} from './options.js';
import { readIndex, readCandidates, pairOf, pairsOver, pointJetStack, gridJetStack } from './stack.js';

const RTA_OPTIONS = ['angle', 'pol'];
const PHASE_OPTIONS = ['angle', 'sine', 'pol'];
const NEEDLE_OPTIONS = ['angle', 'pol', 'candidates', 'fractions'];
const MONITOR_OPTIONS = ['angle', 'top', 'thicknesses'];
const GROWING_OPTIONS = ['angle'];
const FIRST = ['dRdd', 'dTdd', 'dAdd'];
const SECOND = ['d2Rdd', 'd2Tdd', 'd2Add'];
const PHASE_DERIVATIVES = ['dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared'];
const GROWING = ['Rs', 'Ts', 'Rp', 'Tp', 'Rrs', 'Rrp'];

// 'a and b', 'a, b and c'.
const listed = keys => `${keys.slice(0, -1).join(', ')} and ${keys[keys.length - 1]}`;

function checkOptions(options, keys) {
    if (options === null || typeof options !== 'object') throw new TypeError(`options must be { ${keys.join(', ')} }`);
    const unknown = Object.keys(options).find(key => !keys.includes(key));
    if (unknown === undefined) return;
    const known = keys.length === 1 ? `the only option is ${keys[0]}` : `the options are ${listed(keys)}`;
    throw new TypeError(`unknown option '${unknown}'; ${known}`);
}

// The options of a call as its kernels take them: the angle, the polarization,
// and in `extra` what else the kernel takes after the stack. nLam is the grid
// size of a Spectrum call, null at one wavelength.

// The angle and polarization of an R/T call.
function readRtaOptions(options, nLam) {
    checkOptions(options, RTA_OPTIONS);
    return { angle: readAngle(options.angle, nLam), pol: readRtaPol(options.pol), extra: null };
}

// The angle, or the incident sine as jets, and the polarization of a phase
// call, with the sine, or null, in `extra`. With the sine the kernels take no
// angle, and are given 0.
function readPhaseOptions(options, nLam) {
    checkOptions(options, PHASE_OPTIONS);
    if (options.sine === undefined) {
        const angle = readAngle(options.angle, nLam);
        return { angle, pol: readPhasePol(options.pol, atNormalIncidence(angle, null, nLam)), extra: null };
    }
    if (options.angle !== undefined) throw new TypeError('a phase call takes angle or sine, not both');
    const sine = readSine(options.sine, nLam);
    return { angle: 0, pol: readPhasePol(options.pol, atNormalIncidence(0, sine, nLam)), extra: sine };
}

// The angle and polarization of a needle scan, with its candidates and split
// fractions in `extra`.
function readNeedleOptions(options, nLam) {
    checkOptions(options, NEEDLE_OPTIONS);
    return {
        angle: readAngle(options.angle, nLam),
        pol: readRtaPol(options.pol),
        extra: { candidates: readCandidates(options.candidates, nLam), fractions: readFractions(options.fractions) },
    };
}

// The angle of a monitoring curve, with the growing layer's ñ as [re, im] and
// the thicknesses to sample in `extra`.
function readMonitorOptions(options) {
    checkOptions(options, MONITOR_OPTIONS);
    return {
        angle: readAngle(options.angle),
        extra: { top: pairOf(readIndex(options.top, null, 'top')), thicknesses: readThicknesses(options.thicknesses) },
    };
}

// The angle of the spectra after each deposited layer, or of a growing stack.
function readGrowingOptions(options, nLam) {
    checkOptions(options, GROWING_OPTIONS);
    return { angle: readAngle(options.angle, nLam) };
}

// What a growing stack's sample() writes into when given `out`: the result of
// an earlier sample of the same stack.
function checkOut(out, nLam) {
    const fits = key => out[key] instanceof Float64Array && out[key].length === nLam;
    if (out !== null && (typeof out !== 'object' || !GROWING.every(fits))) {
        throw new TypeError(`out must be { ${GROWING.join(', ')} }, each a Float64Array of ${nLam}`);
    }
}

// A growing stack's handle, the growing layer given as a medium is and the
// thickness checked before they reach the kernel's.
function growingHandle(handle, nLam) {
    return {
        setTop(top) {
            handle.setTop(pairsOver(readIndex(top, nLam, 'top'), nLam));
        },
        sample(d, out = null) {
            const thickness = readThickness(d);
            checkOut(out, nLam);
            return handle.sample(thickness, out);
        },
        free() {
            handle.free();
        },
    };
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

// [dR, dT, dA] triples, in the kernel's order, as three arrays.
function splitTriples(values) {
    const n = values.length / 3;
    const out = { dR: new Float64Array(n), dT: new Float64Array(n), dA: new Float64Array(n) };
    for (let i = 0; i < n; i++) {
        out.dR[i] = values[3 * i];
        out.dT[i] = values[3 * i + 1];
        out.dA[i] = values[3 * i + 2];
    }
    return out;
}

// A needle scan over a grid, from the kernel's layout (needleScanSpectrumKernel
// in derivatives.js): R, T and A per wavelength, and gaps and intra each as
// { dR, dT, dA }, flat in the kernel's order.
function needleGrid({ base, gaps, intra }, N) {
    const { dR: R, dT: T, dA: A } = splitTriples(base);
    return { R, T, A, gaps: splitTriples(gaps), intra: splitTriples(intra), N };
}

/**
 * The calls on `kernels`, with pol 's' or 'p':
 *   rta, jacobian, hessian     (λ, θ, pol, kernels.pointStack(stack))
 *   needleScan                 (λ, θ, pol, kernels.pointStack(stack), { candidates, fractions })
 *   phase, phaseJacobian       (λ, θ, pol, { n0Jet, nsJet, layers: [{ nJet, d }] }, sine)
 *   rtaSpectrum                (grid), s and p together
 *   jacobianSpectrum, hessianSpectrum, phaseSpectrum, phaseJacobianSpectrum   (grid, pol)
 *   needleScanSpectrum         (grid, pol, { candidates, fractions })
 *   monitorCurve               (λ, θ, kernels.pointStack(stack), { top, thicknesses }), top as [re, im]
 *   depositionSpectra          (grid), the layers in the stack's order, for the kernel to reverse
 *   growingStack               (grid), returning { setTop, sample, free }, setTop given ñ per
 *                              wavelength as [re, im]
 * kernels.pointStack(stack) and kernels.gridStack(stack, nLam) read the stack
 * of an R/T call the way those kernels take it: stack.js's pointStack and
 * gridStack for the JavaScript ones. The needle scans read theirs the same
 * way, and take the candidates as stack.js's readCandidates gives them. The
 * phase calls read theirs with pointJetStack and gridJetStack, and give the
 * incident sine as a jet, or null for the angle's own. The grid of a Spectrum
 * call is { lambdas, theta_deg } together with what kernels.gridStack, or
 * gridJetStack for the phase calls, returns, and sinJets, the sine per
 * wavelength, when a phase call is given one.
 */
export function stackCalls(kernels) {
    const atOne = (read, readMedia, run) => (stack, wavelength, options = {}) => {
        const lambda = readWavelength(wavelength);
        const { angle, pol, extra } = read(options, null);
        const media = readMedia(stack);
        return polarized(pol, side => run(lambda, angle, side, media, extra));
    };
    const gridOf = (stack, wavelengths, options, read, readMedia) => {
        const lambdas = readWavelengths(wavelengths);
        const { angle, pol, extra } = read(options, lambdas.length);
        return { pol, extra, grid: { lambdas, theta_deg: angle, ...readMedia(stack, lambdas.length) } };
    };
    const overGrid = (read, readMedia, run) => (stack, wavelengths, options = {}) => {
        const { pol, extra, grid } = gridOf(stack, wavelengths, options, read, readMedia);
        return polarized(pol, side => run(grid, side, extra));
    };
    // A phase grid with the incident sine, when the call is given one.
    const phaseGrid = kernel => (grid, side, sine) => kernel(sine ? { ...grid, sinJets: sine } : grid, side);

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
        needleScan: atOne(readNeedleOptions, kernels.pointStack, kernels.needleScan),
        needleScanSpectrum: overGrid(readNeedleOptions, kernels.gridStack,
            (grid, side, needles) => needleGrid(kernels.needleScanSpectrum(grid, side, needles), grid.thick.length)),
        phase: atOne(readPhaseOptions, pointJetStack, kernels.phase),
        phaseSpectrum: overGrid(readPhaseOptions, gridJetStack, phaseGrid(kernels.phaseSpectrum)),
        phaseJacobian: atOne(readPhaseOptions, pointJetStack, (...args) => {
            const { r, t } = kernels.phaseJacobian(...args);
            return { r: typedPhaseSide(r), t: typedPhaseSide(t) };
        }),
        phaseJacobianSpectrum: overGrid(readPhaseOptions, gridJetStack, phaseGrid(kernels.phaseJacobianSpectrum)),
        monitorCurve: (stack, wavelength, options = {}) => {
            const lambda = readWavelength(wavelength);
            const { angle, extra } = readMonitorOptions(options);
            return kernels.monitorCurve(lambda, angle, kernels.pointStack(stack), extra);
        },
        depositionSpectra: (stack, wavelengths, options = {}) =>
            kernels.depositionSpectra(gridOf(stack, wavelengths, options, readGrowingOptions, kernels.gridStack).grid),
        growingStack: (stack, wavelengths, options = {}) => {
            const { grid } = gridOf(stack, wavelengths, options, readGrowingOptions, kernels.gridStack);
            return growingHandle(kernels.growingStack(grid), grid.lambdas.length);
        },
    };
}
