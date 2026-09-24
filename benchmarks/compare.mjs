/**
 * Accuracy against an independent implementation. Needs nothing but Node.
 *
 * Feeds tmmcore the same inputs that were fed to Steven Byrnes' `tmm` package
 * (a different author, a different language) and reports the disagreement.
 *
 * The inputs live in cases.json: wavelength grid, per-layer complex indices and
 * thicknesses, all precomputed. No dispersion evaluation, no material lookup,
 * no unit conversion happens on either side, so the only thing that differs
 * between the two implementations is the transfer-matrix mathematics itself.
 * The cases run at normal incidence and at 30° and 60°, where s and p part.
 *
 * R and T are compared from tmm(), and the complex reflection coefficients from
 * the phase kernel, which reports r as |r|² and its phase. r_s agrees as it
 * stands. r_p is compared with its sign reversed: Byrnes follows the Fresnel
 * convention, r_p = −r_s at normal incidence, and tmmcore the tilted
 * admittances n/cosθ, r_p = r_s there (Macleod 5th ed., §9.7.1, Eq. 9.36).
 *
 * Byrnes' outputs are committed in reference_byrnes.json so this runs with no
 * Python installed. If you would rather not trust a committed file, which is
 * a reasonable instinct, regenerate it; see benchmarks/README.md.
 *
 *   node benchmarks/compare.mjs
 */

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { tmmPhaseDispersion } from '../src/phase.js';
import { jetConstant } from '../src/taylorJet.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const read = f => JSON.parse(readFileSync(join(HERE, f), 'utf8'));

const cases = read('cases.json');
const ref = read('reference_byrnes.json');

const WASM_PATH = join(HERE, '..', 'src', 'tmm_kernel.wasm');
const kernel = existsSync(WASM_PATH)
    ? await instantiateTmmWasm(readFileSync(WASM_PATH))
    : null;

console.log(`Reference: ${ref.source.package} ${ref.source.version}, ${ref.source.author}`);
console.log(`           ${ref.source.url}`);
console.log(`Convention: ${ref.source.convention}\n`);

// r from the phase kernel's |r|² and phase, where phaseRad = −arg r.
const coefficient = side => {
    const magnitude = Math.sqrt(side.magnitudeSquared);
    return [magnitude * Math.cos(side.phaseRad), -magnitude * Math.sin(side.phaseRad)];
};
const distance = ([a, b], [c, d]) => Math.hypot(a - c, b - d);

const rows = [];
let worstAll = 0;
let worstR = 0;

for (const c of cases.cases) {
    const r = ref.cases[c.name];
    if (!r) continue;

    const N = c.nLayers;
    let worstJS = 0, worstWA = 0, worstCoefficient = 0;

    for (let i = 0; i < c.lambdas.length; i++) {
        const lam = c.lambdas[i];
        const n0 = c.n0[i], ns = c.ns[i];
        // layerNK is [layer][lambda] = [re, im]
        const layers = Array.from({ length: N }, (_, k) => ({
            n: c.layerNK[k][i], d: c.thick[k],
        }));
        const jetLayers = layers.map(({ n, d }) => ({ nJet: jetConstant(n[0], n[1]), d }));

        for (const [pol, code, Rref, Tref, rref, sign] of [
            ['s', 0, r.Rs[i], r.Ts[i], r.rs[i], 1],
            ['p', 1, r.Rp[i], r.Tp[i], r.rp[i], -1],
        ]) {
            const js = tmm(lam, c.theta_deg, pol, n0, ns, layers);
            worstJS = Math.max(worstJS, Math.abs(js.R - Rref), Math.abs(js.T - Tref));
            const expected = [sign * rref[0], sign * rref[1]];
            const phase = tmmPhaseDispersion(lam, c.theta_deg, pol,
                jetConstant(n0[0], n0[1]), jetConstant(ns[0], ns[1]), jetLayers);
            worstCoefficient = Math.max(worstCoefficient, distance(coefficient(phase.r), expected));

            if (kernel) {
                const wa = kernel.tmmOne(lam, c.theta_deg, code, n0, ns, layers);
                worstWA = Math.max(worstWA, Math.abs(wa.R - Rref), Math.abs(wa.T - Tref));
                const waPhase = kernel.tmmPhaseOne(lam, c.theta_deg, code,
                    jetConstant(n0[0], n0[1]), jetConstant(ns[0], ns[1]), jetLayers);
                worstCoefficient = Math.max(worstCoefficient,
                    distance(coefficient(waPhase.r), expected));
            }
        }
    }

    worstAll = Math.max(worstAll, worstJS, worstWA);
    worstR = Math.max(worstR, worstCoefficient);
    rows.push({
        name: c.name, N, nLam: c.lambdas.length, theta: c.theta_deg,
        js: worstJS, wasm: worstWA, r: worstCoefficient,
    });
}

const pad = (s, n) => String(s).padEnd(n);
const lpad = (s, n) => String(s).padStart(n);
const e = v => v.toExponential(1);

console.log(`${pad('case', 18)}${lpad('angle', 6)}${lpad('layers', 7)}${lpad('points', 8)}`
    + `${lpad('max |Δ| JS', 13)}${lpad('max |Δ| WASM', 14)}${lpad('max |Δr|', 10)}`);
console.log('-'.repeat(76));
for (const r of rows) {
    console.log(pad(r.name, 18) + lpad(`${r.theta}°`, 6) + lpad(r.N, 7) + lpad(r.nLam, 8) +
                lpad(e(r.js), 13) + lpad(kernel ? e(r.wasm) : 'not built', 14) + lpad(e(r.r), 10));
}

const total = rows.reduce((a, r) => a + r.nLam * 4, 0);
console.log(`\n${total} values of R and T compared across ${rows.length} cases, both polarizations,`);
console.log(`and the complex r at each of their ${total / 2} pairs of wavelength and polarization.`);
console.log(`Worst disagreement with an independently written implementation: ${e(worstAll)} in R and T, `
    + `${e(worstR)} in r`);

// Double precision carries ~16 significant digits; accumulated round-off over a
// 40-layer matrix product lands near 1e-14. Anything above this is a real
// difference in the mathematics, not arithmetic noise.
const LIMIT = 1e-12;
if (!kernel) console.log('\nNote: tmm_kernel.wasm not built; JavaScript path only.');
const pass = worstAll <= LIMIT && worstR <= LIMIT;
console.log(pass
    ? '\nPASS : agreement is at the level of float64 accumulation noise.'
    : `\nFAIL : disagreement exceeds ${LIMIT}.`);
process.exit(pass ? 0 : 1);
