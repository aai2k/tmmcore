/**
 * The JavaScript ports of the kernel's grid calls against the kernel itself:
 * the needle scan over a grid, the monitoring curve, the spectra after each
 * deposited layer, the growing top layer, and the phase grids given their own
 * angular frequencies and incident sines. The ports run the C's arithmetic,
 * so the two agree to the round-off of the math library, as in
 * equivalence.mjs. A grid of typed arrays gives the ports, bit for bit, what
 * the same grid of plain lists gives.
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/js_grid_kernels.mjs
 */

import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import * as J from '../src/js/jet.js';
import { needleScanSpectrumKernel } from '../src/js/derivatives.js';
import { omegaFromLambdaNm, phaseSpectrumKernel, phaseJacobianSpectrumKernel } from '../src/js/phase.js';
import { monitorCurveKernel, depositionSpectraKernel, growingStackKernel } from '../src/js/growing.js';
import { needleScanSpectrumWasm } from '../src/wasm/derivatives.js';
import { PAIRS, writePairLayers } from '../src/wasm/memory.js';
import { phaseSpectrumWasm, phaseJacobianSpectrumWasm } from '../src/wasm/phase.js';
import { monitorCurveWasm, depositionSpectraWasm, growingStackWasm } from '../src/wasm/growing.js';
import { instantiateTmmWasm } from '../src/wasm/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');
if (!existsSync(WASM)) {
    console.log('SKIP : src/wasm/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}
const wasm = await instantiateTmmWasm(readFileSync(WASM));

const ABS_RT = 1e-9;       // R and T live in [0,1]
const REL_DERIV = 1e-7;    // derivatives and phase quantities span many orders of magnitude
const ABS_DERIV = 1e-12;

let checks = 0;
const worst = { rt: 0, rel: 0 };

function near(js, wa, absolute, what) {
    checks++;
    if (Number.isNaN(js) && Number.isNaN(wa)) return;
    const diff = Math.abs(js - wa);
    if (absolute) {
        worst.rt = Math.max(worst.rt, diff);
        assert.ok(diff <= ABS_RT, `${what}: ${js} against ${wa}`);
        return;
    }
    const rel = diff / Math.max(Math.abs(js), Math.abs(wa), Number.MIN_VALUE);
    if (diff > ABS_DERIV) worst.rel = Math.max(worst.rel, rel);
    assert.ok(diff <= ABS_DERIV || rel <= REL_DERIV, `${what}: ${js} against ${wa}`);
}

// Every array of `js` against the same array of `wa`.
function nearAll(js, wa, absolute, what) {
    assert.deepEqual(Object.keys(wa).sort(), Object.keys(js).sort(), `${what}: keys`);
    for (const [key, values] of Object.entries(js)) {
        if (!ArrayBuffer.isView(values)) {
            nearAll(values, wa[key], absolute, `${what}.${key}`);
            continue;
        }
        assert.equal(wa[key].length, values.length, `${what}.${key}: length`);
        values.forEach((value, i) => near(value, wa[key][i], absolute, `${what}.${key}[${i}]`));
    }
}

// ── A dispersive, absorbing design over a grid ───────────────────────────────
// Cauchy n = A + B/λ² (B in nm²) with a constant k.
const LAMBDAS = [400, 470, 550, 700, 1064];
const ANGLES = [0, 18, 37, 58, 76];
const nLam = LAMBDAS.length;
const cauchy = (A, B, k) => ({ A, B, k });
const nkAt = (m, lambda) => [m.A + m.B / (lambda * lambda), m.k];
const jetAt = (m, lambda) => {
    const lambdaJet = J.wavelengthOmegaJet(lambda, omegaFromLambdaNm(lambda));
    return J.jetAdd(J.jetConstant(m.A, m.k), J.jetScale(J.jetReciprocal(J.jetMultiply(lambdaJet, lambdaJet)), m.B));
};

const AIR = cauchy(1, 0, 0);
const GLASS = cauchy(1.5, 4200, 0);
const HIGH = cauchy(2.2, 4.0e4, 2e-4);
const LOW = cauchy(1.44, 4.0e3, 0);
const METAL = cauchy(0.15, 0, 3.2);
const LAYERS = [[HIGH, 62], [LOW, 105], [METAL, 22], [HIGH, 0], [LOW, 140], [HIGH, 41]];
const N = LAYERS.length;

const gridOf = theta_deg => ({
    lambdas: LAMBDAS,
    theta_deg,
    n0List: LAMBDAS.map(lambda => nkAt(AIR, lambda)),
    nsList: LAMBDAS.map(lambda => nkAt(GLASS, lambda)),
    layerNK: LAYERS.map(([m]) => LAMBDAS.map(lambda => nkAt(m, lambda))),
    thick: LAYERS.map(([, d]) => d),
});
const jetGridOf = theta_deg => ({
    lambdas: LAMBDAS,
    theta_deg,
    n0Jets: LAMBDAS.map(lambda => jetAt(AIR, lambda)),
    nsJets: LAMBDAS.map(lambda => jetAt(GLASS, lambda)),
    layerJets: LAYERS.map(([m]) => LAMBDAS.map(lambda => jetAt(m, lambda))),
    thick: LAYERS.map(([, d]) => d),
});
const THETAS = [0, 45, ANGLES];
const thetaTag = theta => (typeof theta === 'number' ? `θ ${theta}` : 'θ per λ');

// ── Needle scan over a grid ──────────────────────────────────────────────────
const CANDIDATES = [cauchy(1.38, 0, 0), cauchy(2.1, 2.0e4, 0.001), METAL];
const candNK = CANDIDATES.map(m => LAMBDAS.map(lambda => nkAt(m, lambda)));
for (const theta of THETAS) {
    for (const fracs of [[0, 0.37, 1], []]) {
        for (const [pol, code] of [['s', 0], ['p', 1]]) {
            const tag = `needle scan ${thetaTag(theta)} pol ${pol} fracs ${fracs.length}`;
            const js = needleScanSpectrumKernel(gridOf(theta), pol, { candNK, fracs });
            const wa = needleScanSpectrumWasm(wasm, gridOf(theta), code, { candNK, fracs }, PAIRS);
            assert.equal(js.base.length, 3 * nLam, `${tag} base length`);
            nearAll({ base: js.base }, { base: wa.base }, true, tag);
            nearAll({ gaps: js.gaps, intra: js.intra }, { gaps: wa.gaps, intra: wa.intra }, false, tag);
        }
    }
}

// ── Phase grids with ω and the incident sine given ───────────────────────────
// The sine as a jet stands for an angle held fixed outside a dispersive
// incident medium: sin θ_outside / n_incident(ω).
const OMEGAS = LAMBDAS.map(lambda => 1.0005 * omegaFromLambdaNm(lambda));
const SIN_JETS = LAMBDAS.map((lambda, i) =>
    J.jetDivide(J.jetConstant(Math.sin(ANGLES[i] * Math.PI / 180)), jetAt(cauchy(1.02, 900, 0), lambda)));
for (const [pol, code] of [['s', 0], ['p', 1]]) {
    for (const extra of [{ omegas: OMEGAS }, { sinJets: SIN_JETS }, { omegas: OMEGAS, sinJets: SIN_JETS }]) {
        const grid = { ...jetGridOf(ANGLES), ...extra };
        const tag = `pol ${pol} ${Object.keys(extra).join(' and ')}`;
        nearAll(phaseSpectrumKernel(grid, pol), phaseSpectrumWasm(wasm, grid, code), false, `phase grid ${tag}`);
        nearAll(phaseJacobianSpectrumKernel(grid, pol), phaseJacobianSpectrumWasm(wasm, grid, code), false,
            `phase Jacobian grid ${tag}`);
    }
}

// ── Monitoring curve ─────────────────────────────────────────────────────────
// The finished stack outermost first; zero thicknesses are skipped and
// d ≤ 0 samples the finished stack alone.
const SAMPLES = [-5, 0, 0.5, 12, 58.5, 130, 275];
for (const lambda of [400, 550, 1064]) {
    for (const theta of [0, 37, 76]) {
        for (const top of [HIGH, METAL]) {
            const args = {
                lambda_nm: lambda, theta_deg: theta, n0: nkAt(AIR, lambda), ns: nkAt(GLASS, lambda),
                base: LAYERS.map(([m, d]) => ({ n: nkAt(m, lambda), d })), top: nkAt(top, lambda),
                thicknesses: SAMPLES,
            };
            nearAll(monitorCurveKernel(args), monitorCurveWasm(wasm, args, writePairLayers), true,
                `monitor curve λ ${lambda} θ ${theta} top ${top.A}`);
        }
    }
}
const bareGlass = {
    lambda_nm: 550, theta_deg: 30, n0: [1, 0], ns: [1.52, 0], base: [], top: [2.35, 0], thicknesses: SAMPLES,
};
nearAll(monitorCurveKernel(bareGlass), monitorCurveWasm(wasm, bareGlass, writePairLayers), true,
    'monitor curve on bare glass');

// ── Spectra after each deposited layer ───────────────────────────────────────
for (const theta of THETAS) {
    nearAll(depositionSpectraKernel(gridOf(theta)), depositionSpectraWasm(wasm, gridOf(theta), PAIRS), true,
        `deposition spectra ${thetaTag(theta)}`);
}
const empty = { ...gridOf(30), layerNK: [], thick: [] };
nearAll(depositionSpectraKernel(empty), depositionSpectraWasm(wasm, empty, PAIRS), true,
    'deposition spectra, no layers');

// ── Growing top layer ────────────────────────────────────────────────────────
for (const theta of THETAS) {
    const tag = `growing stack ${thetaTag(theta)}`;
    const js = growingStackKernel(gridOf(theta));
    const wa = growingStackWasm(wasm, gridOf(theta), PAIRS);
    nearAll(js.sample(0), wa.sample(0), true, `${tag} d 0 before setTop`);
    assert.throws(() => js.sample(10), /sampled before setTop/, `${tag} JS sample before setTop`);
    assert.throws(() => wa.sample(10), /sampled before setTop/, `${tag} WASM sample before setTop`);
    for (const top of [HIGH, METAL]) {
        const ngList = LAMBDAS.map(lambda => nkAt(top, lambda));
        js.setTop(ngList);
        wa.setTop(ngList);
        for (const d of [-3, 0, 0.5, 37.5, 140]) nearAll(js.sample(d), wa.sample(d), true, `${tag} top ${top.A} d ${d}`);
    }
    const jsOut = js.sample(0);
    const waOut = wa.sample(0);
    assert.equal(js.sample(25, jsOut), jsOut, `${tag} JS writes into out`);
    assert.equal(wa.sample(25, waOut), waOut, `${tag} WASM writes into out`);
    nearAll(jsOut, waOut, true, `${tag} written into out`);
    js.free();
    wa.free();
    for (const [name, handle] of [['JS', js], ['WASM', wa]]) {
        assert.throws(() => handle.sample(0), /already freed/, `${tag} ${name} sample after free`);
        assert.throws(() => handle.setTop([]), /already freed/, `${tag} ${name} setTop after free`);
        handle.free();
    }
}

// ── A grid of typed arrays ───────────────────────────────────────────────────
// Wavelengths, angles, thicknesses, ñ pairs and angular frequencies.
const typed = values => Float64Array.from(values);
const typedGridOf = theta => {
    const grid = gridOf(typeof theta === 'number' ? theta : typed(theta));
    return {
        ...grid, lambdas: typed(grid.lambdas), thick: typed(grid.thick),
        n0List: grid.n0List.map(typed), nsList: grid.nsList.map(typed),
        layerNK: grid.layerNK.map(row => row.map(typed)),
    };
};
const FRACS = [0, 0.37, 1];
const HIGH_NK = LAMBDAS.map(lambda => nkAt(HIGH, lambda));
for (const theta of THETAS) {
    const tag = `typed grid ${thetaTag(theta)}`;
    assert.deepStrictEqual(
        needleScanSpectrumKernel(typedGridOf(theta), 'p', { candNK: candNK.map(row => row.map(typed)), fracs: typed(FRACS) }),
        needleScanSpectrumKernel(gridOf(theta), 'p', { candNK, fracs: FRACS }), `needle scan ${tag}`);
    assert.deepStrictEqual(depositionSpectraKernel(typedGridOf(theta)), depositionSpectraKernel(gridOf(theta)),
        `deposition spectra ${tag}`);
    const plain = growingStackKernel(gridOf(theta));
    const typedStack = growingStackKernel(typedGridOf(theta));
    plain.setTop(HIGH_NK);
    typedStack.setTop(HIGH_NK.map(typed));
    assert.deepStrictEqual(typedStack.sample(37.5), plain.sample(37.5), `growing stack ${tag}`);
}
const plainJetGrid = { ...jetGridOf(ANGLES), omegas: OMEGAS };
const typedJetGrid = {
    ...plainJetGrid, lambdas: typed(LAMBDAS), theta_deg: typed(ANGLES), thick: typed(plainJetGrid.thick),
    omegas: typed(OMEGAS),
};
assert.deepStrictEqual(phaseSpectrumKernel(typedJetGrid, 's'), phaseSpectrumKernel(plainJetGrid, 's'),
    'phase grid, typed grid');
assert.deepStrictEqual(phaseJacobianSpectrumKernel(typedJetGrid, 'p'), phaseJacobianSpectrumKernel(plainJetGrid, 'p'),
    'phase Jacobian grid, typed grid');

console.log(`${checks} values checked; worst R/T difference ${worst.rt.toExponential(2)}, `
    + `worst relative ${worst.rel.toExponential(2)}`);
console.log('PASS : JavaScript grid kernels');
