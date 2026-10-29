/**
 * phase.js : phase, GD, GDD and TOD, and their thickness derivatives, on a
 * WebAssembly instance, at one wavelength and over a wavelength grid, with
 * the methods of 0.4 that run them.
 */

import { checkLayerKeys } from '../js/stack.js';
import { omegaFromLambdaNm } from '../js/phase.js';
import { warnOnce } from '../js/deprecated.js';
import { onGrid, named, JETS } from './memory.js';
import { writeJet } from './jet.js';

const PHASE_KEYS = ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared'];
const PHASE_DERIVATIVES = ['dPhaseDeg', 'dGd', 'dGdd', 'dTod', 'dLogMagnitudeSquared'];

// The kernel writes NaN where the coefficient is exactly zero and the phase is
// undefined; the JS reference returns null there, so agree with it.
function readPhase(buf, offset) {
    const magnitudeSquared = buf[offset + 4];
    if (Number.isNaN(magnitudeSquared)) return null;
    const phaseRad = buf[offset];
    return {
        phaseRad,
        phaseDeg: phaseRad * 180 / Math.PI,
        gd: buf[offset + 1],
        gdd: buf[offset + 2],
        tod: buf[offset + 3],
        magnitudeSquared,
    };
}

export const phaseMethods = {
    /** True: every build since 0.5.0 carries the phase-dispersion kernel. */
    hasPhase() { return true; },

    /**
     * Phase, group delay, GDD and TOD at one wavelength : mirrors
     * tmmPhaseDispersion() in src/js/phase.js.
     *
     * @param {number[][]} n0Jet  incident-medium index jet, 4 × [re, im]
     * @param {number[][]} nsJet  substrate index jet
     * @param {{nJet:number[][], d:number}[]} layers
     * @param {{omega?:number, sinTheta0Jet?:number[][]}} [options]
     * @returns {{r: object|null, t: object|null}}
     */
    _phase(lambda_nm, theta_deg, polCode, { n0Jet, nsJet, layers }, options = {}) {
        checkLayerKeys(layers, 'nJet', 'd');
        const N = layers.length;
        const omega = options.omega ?? omegaFromLambdaNm(lambda_nm);
        const sinJet = options.sinTheta0Jet || null;
        // arena: layerJets[8N] | thick[N] | n0[8] | ns[8] | sin[8] | out[10]
        const oLay = 0, oThick = 8 * N, oN0 = oThick + N, oNs = oN0 + 8,
              oSin = oNs + 8, oOut = oSin + 8;
        const need = oOut + 10;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            writeJet(buf, oLay + 8 * i, layers[i].nJet);
            buf[oThick + i] = layers[i].d;
        }
        writeJet(buf, oN0, n0Jet);
        writeJet(buf, oNs, nsJet);
        if (sinJet) writeJet(buf, oSin, sinJet);
        const P = (off) => ptr + off * 8;
        this._tmm_phase(lambda_nm, theta_deg, polCode | 0,
            P(oN0), P(oNs), P(oLay), P(oThick), N, omega, sinJet ? P(oSin) : 0, P(oOut));
        const out = this._view(ptr, need);
        return { r: readPhase(out, oOut), t: readPhase(out, oOut + 5) };
    },

    /**
     * Deprecated, removed in 0.6.0: use phase(stack, wavelength, { angle, pol }).
     * Takes the jets as three arguments, n0Jet, nsJet and layers [{ nJet, d }],
     * then options { omega, sinTheta0Jet }.
     * @returns {{r: object|null, t: object|null}}
     */
    tmmPhaseOne(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options = {}) {
        warnOnce('tmmPhaseOne', 'phase');
        return this._phase(lambda_nm, theta_deg, polCode, { n0Jet, nsJet, layers }, options);
    },

    /**
     * Deprecated, removed in 0.6.0: use phaseSpectrum(stack, wavelengths,
     * { angle, pol }). Phase dispersion over a λ grid for one polarization.
     *
     * @param {number[]} lambdas
     * @param {number[][][]} n0Jets  index jet per λ
     * @param {number[][][]} nsJets  index jet per λ
     * @param {number[][][][]} layerJets  [layer][λ] = index jet
     * @param {number[]} thick  layer thicknesses (nm), length N
     * @param {{omegas?:number[], sinJets?:number[][][]}} [options]
     * @returns {{r: object, t: object}} each `{phaseRad, gd, gdd, tod,
     *   magnitudeSquared}` of Float64Array(nLam). Failed samples hold NaN.
     */
    tmmPhaseSpectrum(lambdas, n0Jets, nsJets, layerJets, thick, theta_deg, polCode, options = {}) {
        warnOnce('tmmPhaseSpectrum', 'phaseSpectrum');
        return phaseSpectrumWasm(this, {
            lambdas, n0Jets, nsJets, layerJets, thick, theta_deg: [Number(theta_deg)],
            omegas: options.omegas, sinJets: options.sinJets,
        }, polCode | 0);
    },

    /**
     * Phase dispersion plus exact thickness derivatives : mirrors
     * tmmPhaseThicknessJacobian(). Layers used AS-IS (index parity).
     *
     * @returns {{r, t}} each the phase quantities plus `dPhaseDeg`, `dGd`,
     *   `dGdd`, `dTod` and `dLogMagnitudeSquared` as Float64Array(N).
     */
    _phaseJacobian(lambda_nm, theta_deg, polCode, { n0Jet, nsJet, layers }, options = {}) {
        checkLayerKeys(layers, 'nJet', 'd');
        const N = layers.length;
        const M = Math.max(1, N);
        const omega = options.omega ?? omegaFromLambdaNm(lambda_nm);
        const sinJet = options.sinTheta0Jet || null;
        // arena: layerJets[8N] | thick[N] | n0[8] | ns[8] | sin[8] | out[10] | deriv[10M]
        const oLay = 0, oThick = 8 * N, oN0 = oThick + N, oNs = oN0 + 8,
              oSin = oNs + 8, oOut = oSin + 8, oDeriv = oOut + 10;
        const need = oDeriv + 10 * M;
        const ptr = this._scratch(need);
        const buf = this._view(ptr, need);
        for (let i = 0; i < N; i++) {
            writeJet(buf, oLay + 8 * i, layers[i].nJet);
            buf[oThick + i] = layers[i].d;
        }
        writeJet(buf, oN0, n0Jet);
        writeJet(buf, oNs, nsJet);
        if (sinJet) writeJet(buf, oSin, sinJet);
        const P = (off) => ptr + off * 8;
        this._tmm_phase_jacobian(lambda_nm, theta_deg, polCode | 0,
            P(oN0), P(oNs), P(oLay), P(oThick), N, omega, sinJet ? P(oSin) : 0,
            P(oOut), P(oDeriv));
        const out = this._view(ptr, need);
        const side = (phaseBase, derivBase) => {
            const base = readPhase(out, phaseBase);
            if (!base) return null;
            const take = (q) => out.slice(oDeriv + (derivBase + q) * N, oDeriv + (derivBase + q) * N + N);
            return {
                ...base,
                dPhaseDeg: take(0), dGd: take(1), dGdd: take(2), dTod: take(3),
                dLogMagnitudeSquared: take(4),
            };
        };
        return { r: side(oOut, 0), t: side(oOut + 5, 5) };
    },

    /**
     * Deprecated, removed in 0.6.0: use phaseJacobian(stack, wavelength,
     * { angle, pol }). Arguments as tmmPhaseOne.
     * @returns {{r, t}}
     */
    tmmPhaseJacobian(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options = {}) {
        warnOnce('tmmPhaseJacobian', 'phaseJacobian');
        return this._phaseJacobian(lambda_nm, theta_deg, polCode, { n0Jet, nsJet, layers }, options);
    },

    /** True: every build since 0.5.0 carries the batched phase Jacobian. */
    hasPhaseJacobianSpectrum() { return true; },

    /**
     * Deprecated, removed in 0.6.0: use phaseJacobianSpectrum(stack,
     * wavelengths, { angle, pol }). Arguments as tmmPhaseSpectrum.
     *
     * @returns {{r: object, t: object}} each the five Float64Array(nLam) of
     *   tmmPhaseSpectrum plus `dPhaseDeg`, `dGd`, `dGdd`, `dTod` and
     *   `dLogMagnitudeSquared` as Float64Array(nLam × N), the derivative for
     *   wavelength i and layer k at `[i * N + k]`.
     */
    tmmPhaseJacobianSpectrum(lambdas, n0Jets, nsJets, layerJets, thick, theta_deg, polCode, options = {}) {
        warnOnce('tmmPhaseJacobianSpectrum', 'phaseJacobianSpectrum');
        return phaseJacobianSpectrumWasm(this, {
            lambdas, n0Jets, nsJets, layerJets, thick, theta_deg: [Number(theta_deg)],
            omegas: options.omegas, sinJets: options.sinJets,
        }, polCode | 0);
    },
};

// The angular frequencies and incident sines that follow the stack in the
// phase exports, read by index; no sines means the angle's own.
function phaseTail(grid, put, putTable) {
    const nLam = grid.lambdas.length;
    return [
        put(Float64Array.from({ length: nLam },
            (_, i) => (grid.omegas ? grid.omegas[i] : omegaFromLambdaNm(grid.lambdas[i])))),
        grid.sinJets ? putTable([grid.sinJets], 1, JETS) : 0,
    ];
}

// One side of a phase grid. Per wavelength the kernel writes the five
// quantities of r, then those of t, and with the derivatives a
// [side][quantity][layer] block of 10 × N.
function phaseSide(values, derivatives, side, N) {
    const nLam = values.length / 10;
    const quantities = named(PHASE_KEYS, PHASE_KEYS.map((_, j) =>
        Float64Array.from({ length: nLam }, (__, i) => values[10 * i + 5 * side + j])));
    if (!derivatives) return quantities;
    PHASE_DERIVATIVES.forEach((key, j) => {
        quantities[key] = new Float64Array(nLam * N);
        for (let i = 0; i < nLam; i++) {
            const from = 10 * N * i + (5 * side + j) * N;
            quantities[key].set(derivatives.subarray(from, from + N), i * N);
        }
    });
    return quantities;
}

/** Phase, GD, GDD and TOD over the grid for one polarization, as phaseSpectrumKernel. */
export function phaseSpectrumWasm(wasm, grid, polCode) {
    const nLam = grid.lambdas.length;
    return onGrid(wasm, grid, JETS, ({ lead, put, putTable, out, read }) => {
        const values = out(10 * nLam);
        wasm._tmm_phase_spectrum(...lead(polCode), ...phaseTail(grid, put, putTable), values);
        const all = read(values, 10 * nLam);
        return { r: phaseSide(all, null, 0, 0), t: phaseSide(all, null, 1, 0) };
    });
}

/** The phase grid with its thickness derivatives, as phaseJacobianSpectrumKernel. */
export function phaseJacobianSpectrumWasm(wasm, grid, polCode) {
    const nLam = grid.lambdas.length;
    const N = grid.thick.length;
    return onGrid(wasm, grid, JETS, ({ lead, put, putTable, out, read }) => {
        const values = out(10 * nLam);
        const derivatives = out(10 * N * nLam);
        wasm._tmm_phase_jacobian_spectrum(...lead(polCode), ...phaseTail(grid, put, putTable), values, derivatives);
        const all = read(values, 10 * nLam);
        const allDerivatives = read(derivatives, 10 * N * nLam);
        return { r: phaseSide(all, allDerivatives, 0, N), t: phaseSide(all, allDerivatives, 1, N) };
    });
}
