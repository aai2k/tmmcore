/**
 * The cases the equivalence suite runs, shared with the Python package's
 * reference data (python/tests/reference.mjs): absorbing, dispersive and
 * oblique-incidence stacks, the needle candidates, and the dispersive
 * materials of the phase checks.
 */

import {
    jetAdd, jetConstant, jetMultiply, jetReciprocal, jetScale, wavelengthOmegaJet,
} from '../src/js/jet.js';

// Indices are ñ = n + ik with k ≥ 0. Ag-like values give a genuinely absorbing,
// numerically demanding stack; metals are where implementations diverge.

export const n0 = [1, 0], ns = [1.52, 0];
export const H = [2.35, 0.0005], L = [1.46, 0], AG = [0.15, 3.2];

export const STACKS = {
    'AR 4-layer': [[H, 110], [L, 95], [H, 60], [L, 130]],
    'Ag 7-layer': [[H, 88], [AG, 22], [L, 140], [H, 70], [AG, 15], [L, 210], [H, 55]],
    'QW 21-layer': Array.from({ length: 21 }, (_, i) => (i % 2 ? [L, 94.2] : [H, 58.5])),
    'opaque Ag': [[L, 120], [AG, 400], [H, 80]],
};
export const LAMBDAS = [400, 550, 632.8, 1064];
export const ANGLES = [0, 15, 45, 60];
export const POLS = [['s', 0], ['p', 1]];
export const CANDIDATES = [H, L, AG];
// Needle positions inside each layer, as fractions of its thickness from the
// incident side; 0 and 1 coincide with the gaps either side of it.
export const INTRA_FRACS = [0, 0.37, 1];

// ── Phase dispersion ─────────────────────────────────────────────────────────
// The phase kernel takes refractive indices as Taylor jets in angular frequency.
// tmmcore owns no material models, so the tests build their own: a Cauchy
// dispersion composed through the jet arithmetic, which is exactly how a caller
// is expected to feed any formula they can write.

export function cauchyJet(lambda, omega, [A, B, C, k]) {
    const lambdaJet = wavelengthOmegaJet(lambda, omega);
    const inverseSquare = jetReciprocal(jetMultiply(lambdaJet, lambdaJet));
    return jetAdd(
        jetConstant(A, k),
        jetAdd(jetScale(inverseSquare, B),
               jetScale(jetMultiply(inverseSquare, inverseSquare), C)),
    );
}

// A, B (nm²), C (nm⁴), k : high index, low index, and a lossy layer.
export const DISPERSIVE = {
    H: [2.25, 2.4e4, 1.1e9, 0.0005],
    L: [1.44, 3.6e3, 2.0e8, 0],
    M: [1.85, 1.2e4, 0, 0.02],
    // Opaque: 20 µm of it saturates the imaginary-phase clamp, so each layer
    // multiplies the matrix magnitude by about cosh(50) and a stack of six
    // pushes the product past the rescale threshold.
    K: [2.0, 0, 0, 0.6],
};
export const PHASE_STACKS = {
    'chirped 12-layer': Array.from({ length: 12 }, (_, i) =>
        [i % 2 ? 'L' : 'H', 70 + 7.5 * i]),
    'lossy 5-layer': [['H', 120], ['M', 45], ['L', 180], ['M', 30], ['H', 95]],
    'QW 21-layer': Array.from({ length: 21 }, (_, i) => (i % 2 ? ['L', 94.2] : ['H', 58.5])),
    'with a zero layer': [['H', 110], ['L', 0], ['H', 88]],
    'with a negative layer': [['H', 110], ['L', -25], ['H', 88]],
};

// Kept out of the ordinary phase checks: this one exists to take the matrix
// product past the rescale threshold, and its saturated phases make TOD a
// difference of nearly equal terms far noisier than the tolerances tuned for
// physically ordinary stacks.
export const OVERFLOW_STACK = [['H', 110], ...Array.from({ length: 6 }, () => ['K', 20000]), ['L', 90]];

// The substrate of the phase checks, a Cauchy glass.
export const PHASE_SUBSTRATE = [1.5, 4.2e3, 0, 0];

// The wavelength grid of the batched phase checks.
export const PHASE_GRID = Array.from({ length: 40 }, (_, i) => 500 + i * 7.5);
