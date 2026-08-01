/**
 * Accuracy against an independent implementation — needs nothing but Node.
 *
 * Feeds tmmcore the same inputs that were fed to Steven Byrnes' `tmm` package
 * (a different author, a different language) and reports the disagreement.
 *
 * The inputs live in cases.json: wavelength grid, per-layer complex indices and
 * thicknesses, all precomputed. No dispersion evaluation, no material lookup,
 * no unit conversion happens on either side, so the only thing that differs
 * between the two implementations is the transfer-matrix mathematics itself.
 *
 * Byrnes' outputs are committed in reference_byrnes.json so this runs with no
 * Python installed. If you would rather not trust a committed file — a
 * reasonable instinct — regenerate it; see benchmarks/README.md.
 *
 *   node benchmarks/compare.mjs
 */

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const read = f => JSON.parse(readFileSync(join(HERE, f), 'utf8'));

const cases = read('cases.json');
const ref = read('reference_byrnes.json');

const WASM_PATH = join(HERE, '..', 'src', 'tmm_kernel.wasm');
const kernel = existsSync(WASM_PATH)
    ? await instantiateTmmWasm(readFileSync(WASM_PATH))
    : null;

console.log(`Reference: ${ref.source.package} ${ref.source.version} — ${ref.source.author}`);
console.log(`           ${ref.source.url}`);
console.log(`Convention: ${ref.source.convention}\n`);

const rows = [];
let worstAll = 0;

for (const c of cases.cases) {
    const r = ref.cases[c.name];
    if (!r) continue;

    const N = c.nLayers;
    let worstJS = 0, worstWA = 0;

    for (let i = 0; i < c.lambdas.length; i++) {
        const lam = c.lambdas[i];
        const n0 = c.n0[i], ns = c.ns[i];
        // layerNK is [layer][lambda] = [re, im]
        const layers = Array.from({ length: N }, (_, k) => ({
            n: c.layerNK[k][i], d: c.thick[k],
        }));

        for (const [pol, code, Rref, Tref] of [
            ['s', 0, r.Rs[i], r.Ts[i]],
            ['p', 1, r.Rp[i], r.Tp[i]],
        ]) {
            const js = tmm(lam, c.theta_deg, pol, n0, ns, layers);
            worstJS = Math.max(worstJS, Math.abs(js.R - Rref), Math.abs(js.T - Tref));

            if (kernel) {
                const wa = kernel.tmmOne(lam, c.theta_deg, code, n0, ns, layers);
                worstWA = Math.max(worstWA, Math.abs(wa.R - Rref), Math.abs(wa.T - Tref));
            }
        }
    }

    worstAll = Math.max(worstAll, worstJS, worstWA);
    rows.push({ name: c.name, N, nLam: c.lambdas.length, js: worstJS, wasm: worstWA });
}

const pad = (s, n) => String(s).padEnd(n);
const lpad = (s, n) => String(s).padStart(n);
const e = v => v.toExponential(1);

console.log(`${pad('case', 14)}${lpad('layers', 7)}${lpad('points', 8)}${lpad('max |Δ| JS', 13)}${lpad('max |Δ| WASM', 14)}`);
console.log('-'.repeat(56));
for (const r of rows) {
    console.log(pad(r.name, 14) + lpad(r.N, 7) + lpad(r.nLam, 8) +
                lpad(e(r.js), 13) + lpad(kernel ? e(r.wasm) : 'not built', 14));
}

const total = rows.reduce((a, r) => a + r.nLam * 4, 0);
console.log(`\n${total} values compared across ${rows.length} cases, both polarizations.`);
console.log(`Worst disagreement with an independently written implementation: ${e(worstAll)}`);

// Double precision carries ~16 significant digits; accumulated round-off over a
// 40-layer matrix product lands near 1e-14. Anything above this is a real
// difference in the mathematics, not arithmetic noise.
const LIMIT = 1e-12;
if (!kernel) console.log('\nNote: tmm_kernel.wasm not built — JavaScript path only.');
console.log(worstAll <= LIMIT
    ? '\nPASS — agreement is at the level of float64 accumulation noise.'
    : `\nFAIL — disagreement exceeds ${LIMIT}.`);
process.exit(worstAll <= LIMIT ? 0 : 1);
