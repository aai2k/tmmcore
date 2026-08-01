/**
 * tmmcore — transfer-matrix method for multilayer thin films.
 *
 * Two interchangeable implementations of the same physics:
 *
 *   • a JavaScript reference implementation, always available, no dependencies
 *   • a C kernel compiled to WebAssembly, roughly an order of magnitude faster
 *
 * The JavaScript path works immediately on import. WebAssembly is opt-in: load
 * the `.wasm` bytes, instantiate, and call the instance methods. Results agree
 * to float64 round-off either way — see `tests/equivalence.mjs`.
 *
 * Conventions, which matter more than anything else here:
 *   ñ = n + ik              k ≥ 0 for absorbing media
 *   time factor exp(−iωt)
 *   wavelengths and thicknesses in nm, angles in degrees from normal
 *   complex numbers are [re, im] pairs
 *   layers are ordered from the incident medium toward the substrate
 *
 * MIT licensed.
 */

export {
    // Spectral quantities
    tmm,
    // Exact analytic derivatives
    tmmThicknessJacobian,
    tmmThicknessHessian,
    tmmNeedleScan,
    // Low-level primitives, for building variants on the same conventions
    cadd, csub, cmul, cdiv, cabs2, cconj, csqrt, ccos, csin, creal, cimag,
    matmul, rescaleMatrix, snellCosTheta, layerMatrix, cmatvec,
} from './tmm.js';

export * from './tmmWasm.js';
