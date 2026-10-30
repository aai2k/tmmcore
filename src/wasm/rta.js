/**
 * rta.js : R, T and A on a WebAssembly instance, at one wavelength and over a
 * wavelength grid, with the methods of 0.4 that run them.
 */

import { checkLayerKeys } from '../js/stack.js';
import { warnOnce } from '../js/deprecated.js';
import { onGrid, named, PAIRS, writePairLayers } from './memory.js';

const RTA_KEYS = ['Rs', 'Ts', 'As', 'Rp', 'Tp', 'Ap'];

export const rtaMethods = {
    /**
     * R, T and A at one (λ, θ, pol), pol 0 for s and 1 for p : mirrors tmm()
     * in src/js/rta.js. The layers are written by `writeLayers` (memory.js)
     * and their keys are checked by the caller.
     * @returns {{R:number,T:number,A:number}}
     */
    _rta(lambda_nm, theta_deg, polCode, { n0, ns, layers }, writeLayers) {
        const N = layers.length;
        const need = 3 * N + 3;                  // layers [0..3N) + out [3N..3N+3)
        const ptr = this._scratch(need);         // may grow→detach; view created AFTER
        const buf = this._view(ptr, need);
        writeLayers(buf, layers);
        const outPtr = ptr + 3 * N * 8;
        this._tmm_rta(lambda_nm, theta_deg, polCode | 0,
            n0[0], n0[1], ns[0], ns[1], ptr, N, outPtr);
        return { R: buf[3 * N], T: buf[3 * N + 1], A: buf[3 * N + 2] };
    },

    /**
     * Deprecated, removed in 0.6.0: use rta(stack, wavelength, { angle, pol }).
     * R, T and A at one (λ, θ, pol), pol 0 for s and 1 for p.
     * @returns {{R:number,T:number,A:number}}
     */
    tmmOne(lambda_nm, theta_deg, polCode, n0, ns, layers) {
        warnOnce('tmmOne', 'rta');
        checkLayerKeys(layers, 'n', 'd');
        return this._rta(lambda_nm, theta_deg, polCode, { n0, ns, layers }, writePairLayers);
    },

    /**
     * Deprecated, removed in 0.6.0: use rtaSpectrum(stack, wavelengths,
     * { angle, pol }). R, T and A over a λ grid for both polarizations.
     * @param {number[]} lambdas
     * @param {[number,number][]} n0List  incident ñ per λ
     * @param {[number,number][]} nsList  substrate ñ per λ
     * @param {[number,number][][]} layerNK  [layer][λ] = ñ
     * @param {number[]} thick   layer thicknesses (nm), length N
     * @param {number}   theta_deg
     * @returns {{Rs,Ts,As,Rp,Tp,Ap}} each a Float64Array(nLam)
     */
    tmmSpectrum(lambdas, n0List, nsList, layerNK, thick, theta_deg) {
        warnOnce('tmmSpectrum', 'rtaSpectrum');
        return rtaSpectrumWasm(this, { lambdas, n0List, nsList, layerNK, thick, theta_deg: [Number(theta_deg)] }, PAIRS);
    },
};

/**
 * R, T and A in s and p over the grid: { Rs, Ts, As, Rp, Tp, Ap }. `indices`
 * says how the grid holds its indices (memory.js): PAIRS or STACK.
 */
export function rtaSpectrumWasm(wasm, grid, indices) {
    const nLam = grid.lambdas.length;
    return onGrid(wasm, grid, indices, ({ lead, out, read }) => {
        const outs = RTA_KEYS.map(() => out(nLam));
        wasm._tmm_rta_spectrum(...lead(), ...outs);
        return named(RTA_KEYS, outs.map(ptr => read(ptr, nLam)));
    });
}
