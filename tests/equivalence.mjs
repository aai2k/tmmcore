/**
 * JavaScript ⇆ WebAssembly equivalence.
 *
 * The C kernel is a line-by-line port of the JavaScript implementation. This
 * test drives both with identical inputs across absorbing, dispersive and
 * oblique-incidence cases in s and p polarization, and checks every returned
 * quantity — R, T, A, the thickness Jacobian, the thickness Hessian, and the
 * needle P-function.
 *
 * Agreement is not bit-exact by design. The only divergence is libm: the WASM
 * build uses musl's sin/cos/exp/atan2, the JavaScript engine uses its own, and
 * they differ at roughly 1 ULP. The tolerances below sit far above that noise
 * and far below any tolerance that matters physically.
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/equivalence.mjs
 */

import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import {
    tmm, tmmThicknessJacobian, tmmThicknessHessian, tmmNeedleScan,
} from '../src/tmm.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');

if (!existsSync(WASM)) {
    console.log('SKIP — src/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}

const ABS_RTA = 1e-9;      // R, T, A live in [0,1]
const REL_DERIV = 1e-7;    // analytic derivatives span many orders of magnitude
const ABS_DERIV = 1e-12;

// ── Cases ────────────────────────────────────────────────────────────────────
// Indices are ñ = n + ik with k ≥ 0. Ag-like values give a genuinely absorbing,
// numerically demanding stack; metals are where implementations diverge.

const n0 = [1, 0], ns = [1.52, 0];
const H = [2.35, 0.0005], L = [1.46, 0], AG = [0.15, 3.2];

const STACKS = {
    'AR 4-layer': [[H, 110], [L, 95], [H, 60], [L, 130]],
    'Ag 7-layer': [[H, 88], [AG, 22], [L, 140], [H, 70], [AG, 15], [L, 210], [H, 55]],
    'QW 21-layer': Array.from({ length: 21 }, (_, i) => (i % 2 ? [L, 94.2] : [H, 58.5])),
    'opaque Ag': [[L, 120], [AG, 400], [H, 80]],
};
const LAMBDAS = [400, 550, 632.8, 1064];
const ANGLES = [0, 15, 45, 60];
const POLS = [['s', 0], ['p', 1]];
const CANDIDATES = [H, L, AG];

const mk = rows => rows.map(([n, d]) => ({ n, d }));

// ── Comparison ───────────────────────────────────────────────────────────────

let checks = 0, failures = 0;
const worst = { rta: 0, derivAbs: 0, derivRel: 0 };

function near(js, wa, kind, label) {
    checks++;
    const diff = Math.abs(js - wa);
    if (kind === 'rta') {
        if (diff > worst.rta) worst.rta = diff;
        if (diff <= ABS_RTA) return;
    } else {
        const rel = diff / Math.max(Math.abs(js), Math.abs(wa), Number.MIN_VALUE);
        if (diff > worst.derivAbs) worst.derivAbs = diff;
        // A relative figure is only meaningful once the absolute difference is
        // above the noise floor; below it, rel divides one rounding error by
        // another and reports nonsense.
        if (diff > ABS_DERIV && rel > worst.derivRel) worst.derivRel = rel;
        if (diff <= ABS_DERIV || rel <= REL_DERIV) return;
    }
    failures++;
    if (failures <= 10) console.log(`  FAIL ${label}: js=${js} wasm=${wa} Δ=${diff}`);
}

const nearAll = (js, wa, kind, label) => {
    for (let i = 0; i < js.length; i++) near(js[i], wa[i], kind, `${label}[${i}]`);
};

// ── Run ──────────────────────────────────────────────────────────────────────

const k = await instantiateTmmWasm(readFileSync(WASM));
if (!k) { console.log('FAIL — instantiation returned null'); process.exit(1); }

for (const [stackName, rows] of Object.entries(STACKS)) {
    const layers = mk(rows);
    for (const lam of LAMBDAS) {
        for (const th of ANGLES) {
            for (const [pol, code] of POLS) {
                const at = `${stackName} λ=${lam} θ=${th} ${pol}`;

                const a = tmm(lam, th, pol, n0, ns, layers);
                const b = k.tmmOne(lam, th, code, n0, ns, layers);
                near(a.R, b.R, 'rta', `${at} R`);
                near(a.T, b.T, 'rta', `${at} T`);
                near(a.A, b.A, 'rta', `${at} A`);

                const ja = tmmThicknessJacobian(lam, th, pol, n0, ns, layers);
                const jb = k.tmmJacobian(lam, th, code, n0, ns, layers);
                nearAll(ja.dRdd, jb.dRdd, 'deriv', `${at} dR/dd`);
                nearAll(ja.dTdd, jb.dTdd, 'deriv', `${at} dT/dd`);
                nearAll(ja.dAdd, jb.dAdd, 'deriv', `${at} dA/dd`);

                if (k.hasHessian()) {
                    const ha = tmmThicknessHessian(lam, th, pol, n0, ns, layers);
                    const hb = k.tmmHessian(lam, th, code, n0, ns, layers);
                    for (let i = 0; i < ha.N; i++) {
                        nearAll(ha.d2Rdd[i], hb.d2Rdd[i], 'deriv', `${at} d²R[${i}]`);
                        nearAll(ha.d2Tdd[i], hb.d2Tdd[i], 'deriv', `${at} d²T[${i}]`);
                        nearAll(ha.d2Add[i], hb.d2Add[i], 'deriv', `${at} d²A[${i}]`);
                    }
                }

                const na = tmmNeedleScan(lam, th, pol, n0, ns, layers, CANDIDATES);
                const nb = k.tmmNeedleScan(lam, th, code, n0, ns, layers, CANDIDATES);
                for (let pos = 0; pos < na.gaps.length; pos++) {
                    for (let c = 0; c < CANDIDATES.length; c++) {
                        const g = na.gaps[pos][c], h = nb.gaps[pos][c];
                        near(g.dR, h.dR, 'deriv', `${at} needle[${pos}][${c}] dR`);
                        near(g.dT, h.dT, 'deriv', `${at} needle[${pos}][${c}] dT`);
                        near(g.dA, h.dA, 'deriv', `${at} needle[${pos}][${c}] dA`);
                    }
                }
            }
        }
    }
}

console.log(`\n${checks} comparisons across ${Object.keys(STACKS).length} stacks, ` +
            `${LAMBDAS.length} wavelengths, ${ANGLES.length} angles, s and p.`);
console.log(`worst |Δ| on R/T/A     : ${worst.rta.toExponential(2)}  (tolerance ${ABS_RTA})`);
console.log(`worst |Δ| on derivatives: ${worst.derivAbs.toExponential(2)}  (tolerance ${ABS_DERIV} abs / ${REL_DERIV} rel)`);
console.log(`worst relative, above the ${ABS_DERIV} floor: ` +
            (worst.derivRel === 0 ? 'none exceeded the floor' : worst.derivRel.toExponential(2)));
console.log(failures === 0 ? '\nPASS — JavaScript and WebAssembly agree.'
                           : `\nFAIL — ${failures} comparisons outside tolerance.`);
process.exit(failures === 0 ? 0 : 1);
