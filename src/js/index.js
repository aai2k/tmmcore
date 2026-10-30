/**
 * tmmcore : transfer-matrix method for multilayer thin films.
 *
 * Two interchangeable implementations of the same physics:
 *
 *   • a JavaScript reference implementation, always available, no dependencies
 *   • a C kernel compiled to WebAssembly, roughly an order of magnitude faster
 *
 * The JavaScript path works immediately on import. WebAssembly is opt-in: load
 * the `.wasm` bytes, instantiate, and call the same names on the instance.
 * Results agree to float64 round-off either way : see `tests/equivalence.mjs`.
 *
 * Every call takes a stack, a wavelength in nm or a list of them, and options
 * { angle, pol } (calls.js):
 *
 *   const stack = { incident: 1, substrate: 1.52,
 *                   layers: [{ n: 2.35, d: 58.5 }, { n: 1.46, d: 94.2 }] };
 *   rta(stack, 550).R                          R for natural light at 550 nm
 *   rtaSpectrum(stack, wavelengths, { angle: 45, pol: 's' }).R
 *
 * Conventions, which matter more than anything else here:
 *   ñ = n + ik              k ≥ 0 for absorbing media
 *   time factor exp(−iωt)
 *   wavelengths and thicknesses in nm, angles in degrees from normal
 *   complex numbers in the primitives and jets are [re, im] pairs
 *   layers are ordered from the incident medium toward the substrate
 *
 * The names of 0.4 keep their arguments and results through 0.5.x, print one
 * warning on first use, and are removed in 0.6.0.
 *
 * MIT licensed.
 */

import { tmm as tmmKernel, rtaSpectrumKernel } from './rta.js';
import {
    tmmThicknessJacobian as jacobianKernel,
    tmmThicknessHessian as hessianKernel,
    jacobianSpectrumKernel, hessianSpectrumKernel,
} from './derivatives.js';
import {
    tmmPhaseDispersion as phaseKernel,
    tmmPhaseThicknessJacobian as phaseJacobianKernel,
    tmmCoefficientJets as coefficientJetsKernel,
    tmmCoefficientThicknessJets as coefficientThicknessJetsKernel,
    coefficientPhaseDispersion as coefficientPhaseKernel,
    coefficientPhaseThicknessDerivatives as coefficientPhaseDerivativesKernel,
    phaseSpectrumKernel, phaseJacobianSpectrumKernel,
} from './phase.js';
import { stackCalls } from './calls.js';
import { pointStack, gridStack } from './stack.js';
import { deprecated } from './deprecated.js';

export const {
    rta, rtaSpectrum,
    jacobian, jacobianSpectrum,
    hessian, hessianSpectrum,
    phase, phaseSpectrum,
    phaseJacobian, phaseJacobianSpectrum,
} = stackCalls({
    pointStack,
    gridStack,
    rta: (lambda, theta, pol, { n0, ns, layers }) => tmmKernel(lambda, theta, pol, n0, ns, layers),
    jacobian: (lambda, theta, pol, { n0, ns, layers }) => jacobianKernel(lambda, theta, pol, n0, ns, layers),
    hessian: (lambda, theta, pol, { n0, ns, layers }) => hessianKernel(lambda, theta, pol, n0, ns, layers),
    phase: (lambda, theta, pol, { n0Jet, nsJet, layers }) =>
        phaseKernel(lambda, theta, pol, n0Jet, nsJet, layers),
    phaseJacobian: (lambda, theta, pol, { n0Jet, nsJet, layers }) =>
        phaseJacobianKernel(lambda, theta, pol, n0Jet, nsJet, layers),
    rtaSpectrum: rtaSpectrumKernel,
    jacobianSpectrum: jacobianSpectrumKernel,
    hessianSpectrum: hessianSpectrumKernel,
    phaseSpectrum: phaseSpectrumKernel,
    phaseJacobianSpectrum: phaseJacobianSpectrumKernel,
});

// The names of 0.4, with their old arguments: (λ, θ, pol, n0, ns, layers),
// pol 's' or 'p', each ñ as [re, im].
export const tmm = deprecated('tmm', 'rta', tmmKernel);
export const tmmThicknessJacobian = deprecated('tmmThicknessJacobian', 'jacobian', jacobianKernel);
export const tmmThicknessHessian = deprecated('tmmThicknessHessian', 'hessian', hessianKernel);
export const tmmPhaseDispersion = deprecated('tmmPhaseDispersion', 'phase', phaseKernel);
export const tmmPhaseThicknessJacobian =
    deprecated('tmmPhaseThicknessJacobian', 'phaseJacobian', phaseJacobianKernel);

// The complex r and t as jets, used inside the phase calls; no public name
// replaces them.
export const tmmCoefficientJets = deprecated('tmmCoefficientJets', null, coefficientJetsKernel);
export const tmmCoefficientThicknessJets =
    deprecated('tmmCoefficientThicknessJets', null, coefficientThicknessJetsKernel);
export const coefficientPhaseDispersion =
    deprecated('coefficientPhaseDispersion', null, coefficientPhaseKernel);
export const coefficientPhaseThicknessDerivatives =
    deprecated('coefficientPhaseThicknessDerivatives', null, coefficientPhaseDerivativesKernel);

// Exact analytic needle P-function, one wavelength
export { tmmNeedleScan } from './derivatives.js';

// Low-level primitives, for building variants on the same conventions
export {
    cadd, csub, cmul, cdiv, cabs2, cconj, csqrt, ccos, csin, creal, cimag,
    matmul, rescaleMatrix, snellCosTheta, incidentCosTheta, layerMatrix, layerLogScale, cmatvec,
} from './primitives.js';

// Unit helpers
export { C_NM_PER_FS, omegaFromLambdaNm } from './phase.js';

// Taylor-jet arithmetic. The phase calls take a dispersive index as a jet;
// these compose one from any dispersion formula you can write.
export * from './jet.js';

export * from '../wasm/index.js';
