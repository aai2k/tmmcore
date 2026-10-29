/**
 * rta.js : reflectance, transmittance and absorptance of a stack, at one
 * wavelength and over a wavelength grid. Conventions are in primitives.js.
 */

import {
    cadd, csub, cmul, cdiv, cabs2, creal, matmul, cmatvec, rescaleMatrix, incidence, snellCosTheta,
    incidentCosTheta, layerLogScale, layerAdmittance, layerMatrix, substrateVector, transmittedFlux,
    layerAbsorptance,
} from './primitives.js';
import { checkLayerKeys, pointArgs } from './stack.js';

// ── Core TMM for one wavelength ───────────────────────────────────────────────

/**
 * @param {number}   lambda_nm  wavelength in nm
 * @param {number}   theta_deg  angle of incidence in degrees
 * @param {string}   pol        's' or 'p'
 * @param {[number,number]} n0  complex n of incident medium
 * @param {[number,number]} ns  complex n of substrate (exit medium)
 * @param {{ n:[number,number], d:number }[]} layers  each layer { n: [re,im], d: thickness_nm }
 * @returns {{ R:number, T:number, A:number }}
 */
export function tmm(lambda_nm, theta_deg, pol, n0, ns, layers) {
    checkLayerKeys(layers, 'n', 'd');
    const { sinTheta0, cosTheta0 } = incidence(theta_deg);

    // Admittance of incident medium
    const eta0 = layerAdmittance(n0, incidentCosTheta(n0, sinTheta0, cosTheta0), pol);

    // [1, ηs], or [0, 1] at the substrate's critical angle in p
    const substrate = substrateVector(ns, snellCosTheta(n0, sinTheta0, ns, cosTheta0), pol);

    // Build total transfer matrix M = M1 × M2 × ... × MN
    let M = [[  [1, 0], [0, 0]  ], [  [0, 0], [1, 0]  ]]; // identity
    let logScale = 0;

    for (const { n, d } of layers) {
        // Not positive, NaN included: the layer is absent.
        if (!(d > 0)) continue;
        const cosThetaJ = snellCosTheta(n0, sinTheta0, n, cosTheta0);
        const Mj = layerMatrix(n, d, lambda_nm, cosThetaJ, pol);
        M = matmul(M, Mj);
        logScale += layerLogScale(n, d, lambda_nm, cosThetaJ) + rescaleMatrix(M);
    }

    // [B, C]^T = M × substrate
    const [B, C] = cmatvec(M, substrate);

    // r = (η0 B - C) / (η0 B + C)
    const eta0B = cmul(eta0, B);
    const r = cdiv(csub(eta0B, C), cadd(eta0B, C));

    // t = 2 η0 bs / (η0 B + C)
    const t = cmul(cdiv(cmul([2, 0], eta0), cadd(eta0B, C)), substrate[0]);

    const R = cabs2(r);

    // T = Re(etaS) / Re(eta0) * |t|²
    const T = Math.max(0, transmittedFlux(substrate) / creal(eta0) * cabs2(t) * Math.exp(-2 * logScale));

    const A = Math.max(0, layerAbsorptance(R, T, eta0, r));

    return { R, T, A };
}

// ── Over a wavelength grid ───────────────────────────────────────────────────

/** R, T and A in s and p at every wavelength: { Rs, Ts, As, Rp, Tp, Ap }. */
export function rtaSpectrumKernel(grid) {
    const nLam = grid.lambdas.length;
    const out = {};
    for (const key of ['Rs', 'Ts', 'As', 'Rp', 'Tp', 'Ap']) out[key] = new Float64Array(nLam);
    for (let i = 0; i < nLam; i++) {
        const s = tmm(...pointArgs(grid, 's', i));
        const p = tmm(...pointArgs(grid, 'p', i));
        out.Rs[i] = s.R; out.Ts[i] = s.T; out.As[i] = s.A;
        out.Rp[i] = p.R; out.Tp[i] = p.T; out.Ap[i] = p.A;
    }
    return out;
}
