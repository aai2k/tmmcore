/**
 * JavaScript ⇆ WebAssembly equivalence.
 *
 * The C kernel is a line-by-line port of the JavaScript implementation. This
 * test drives both with identical inputs across absorbing, dispersive and
 * oblique-incidence cases in s and p polarization, and checks every returned
 * quantity : R, T, A, the thickness Jacobian, the thickness Hessian, the
 * needle P-function, and the phase quantities with their thickness derivatives.
 *
 * Agreement is not bit-exact by design. The only divergence is libm: the WASM
 * build uses musl's sin/cos/exp/atan2, the JavaScript engine uses its own, and
 * they differ at roughly 1 ULP. The tolerances below sit far above that noise
 * and far below any tolerance that matters physically.
 *
 * The phase quantities amplify that noise, because each derivative order is a
 * difference of nearly equal terms: the coefficient agrees to a few ULP, GD to
 * tens, GDD and TOD to a few hundred. Still 1e-11 relative in the worst case,
 * which is orders below the precision of any measured n and k.
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
import {
    omegaFromLambdaNm, tmmPhaseDispersion, tmmPhaseThicknessJacobian,
} from '../src/phase.js';
import {
    jetAdd, jetConstant, jetDivide, jetMultiply, jetReciprocal, jetScale,
    wavelengthOmegaJet,
} from '../src/taylorJet.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');

if (!existsSync(WASM)) {
    console.log('SKIP : src/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
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
// Needle positions inside each layer, as fractions of its thickness from the
// incident side; 0 and 1 coincide with the gaps either side of it.
const INTRA_FRACS = [0, 0.37, 1];

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

// The needle output inside every layer, JavaScript against WebAssembly.
function nearIntra(na, nb, at) {
    const entries = na.intra.flatMap((row, layer) => row.flatMap(({ frac, perCand }, f) =>
        perCand.map((g, c) => [`${at} intra[${layer}] at ${frac} [${c}]`, g, nb.intra[layer][f].perCand[c]])));
    for (const [where, g, h] of entries) {
        for (const q of ['dR', 'dT', 'dA']) near(g[q], h[q], 'deriv', `${where} ${q}`);
    }
}

// ── Run ──────────────────────────────────────────────────────────────────────

const k = await instantiateTmmWasm(readFileSync(WASM));
if (!k) { console.log('FAIL : instantiation returned null'); process.exit(1); }

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

                const na = tmmNeedleScan(lam, th, pol, n0, ns, layers, CANDIDATES, INTRA_FRACS);
                const nb = k.tmmNeedleScan(lam, th, code, n0, ns, layers, CANDIDATES, INTRA_FRACS);
                for (let pos = 0; pos < na.gaps.length; pos++) {
                    for (let c = 0; c < CANDIDATES.length; c++) {
                        const g = na.gaps[pos][c], h = nb.gaps[pos][c];
                        near(g.dR, h.dR, 'deriv', `${at} needle[${pos}][${c}] dR`);
                        near(g.dT, h.dT, 'deriv', `${at} needle[${pos}][${c}] dT`);
                        near(g.dA, h.dA, 'deriv', `${at} needle[${pos}][${c}] dA`);
                    }
                }
                nearIntra(na, nb, at);
            }
        }
    }
}

// ── Phase dispersion ─────────────────────────────────────────────────────────
// The phase kernel takes refractive indices as Taylor jets in angular frequency.
// tmmcore owns no material models, so the test builds its own: a Cauchy
// dispersion composed through the jet arithmetic, which is exactly how a caller
// is expected to feed any formula they can write.

function cauchyJet(lambda, omega, A, B, C, k) {
    const lambdaJet = wavelengthOmegaJet(lambda, omega);
    const inverseSquare = jetReciprocal(jetMultiply(lambdaJet, lambdaJet));
    return jetAdd(
        jetConstant(A, k),
        jetAdd(jetScale(inverseSquare, B),
               jetScale(jetMultiply(inverseSquare, inverseSquare), C)),
    );
}

// A, B (nm²), C (nm⁴), k : high index, low index, and a lossy layer.
const DISPERSIVE = {
    H: [2.25, 2.4e4, 1.1e9, 0.0005],
    L: [1.44, 3.6e3, 2.0e8, 0],
    M: [1.85, 1.2e4, 0, 0.02],
    // Opaque: 20 µm of it saturates the imaginary-phase clamp, so each layer
    // multiplies the matrix magnitude by about cosh(50) and a stack of six
    // pushes the product past the rescale threshold.
    K: [2.0, 0, 0, 0.6],
};
const PHASE_STACKS = {
    'chirped 12-layer': Array.from({ length: 12 }, (_, i) =>
        [i % 2 ? 'L' : 'H', 70 + 7.5 * i]),
    'lossy 5-layer': [['H', 120], ['M', 45], ['L', 180], ['M', 30], ['H', 95]],
    'QW 21-layer': Array.from({ length: 21 }, (_, i) => (i % 2 ? ['L', 94.2] : ['H', 58.5])),
    'with a zero layer': [['H', 110], ['L', 0], ['H', 88]],
    'with a negative layer': [['H', 110], ['L', -25], ['H', 88]],
};

// Kept out of the loop above: this one exists to take the matrix product past
// the rescale threshold, and its saturated phases make TOD a difference of
// nearly equal terms far noisier than the tolerances tuned for physically
// ordinary stacks.
const OVERFLOW_STACK = [['H', 110], ...Array.from({ length: 6 }, () => ['K', 20000]), ['L', 90]];

let phaseChecks = 0;
const worstPhase = { abs: 0, rel: 0 };

// `floor` replaces ABS_DERIV where a quantity's own cancellation sets a higher
// noise level than round-off in the result.
function nearPhase(js, wa, label, floor = ABS_DERIV) {
    phaseChecks++;
    if (js === null && wa === null) return;
    if (js === null || wa === null) {
        failures++;
        if (failures <= 10) console.log(`  FAIL ${label}: js=${js} wasm=${wa} (validity differs)`);
        return;
    }
    const diff = Math.abs(js - wa);
    const rel = diff / Math.max(Math.abs(js), Math.abs(wa), Number.MIN_VALUE);
    if (diff > worstPhase.abs) worstPhase.abs = diff;
    if (diff > floor && rel > worstPhase.rel) worstPhase.rel = rel;
    if (diff <= floor || rel <= REL_DERIV) return;
    failures++;
    if (failures <= 10) console.log(`  FAIL ${label}: js=${js} wasm=${wa} Δ=${diff}`);
}

const PHASE_KEYS = ['phaseRad', 'phaseDeg', 'gd', 'gdd', 'tod', 'magnitudeSquared'];

function comparePhaseSide(a, b, label) {
    if (a === null || b === null) { nearPhase(a, b, label); return; }
    for (const key of PHASE_KEYS) nearPhase(a[key], b[key], `${label}.${key}`);
}

if (!k.hasPhase()) {
    console.log('\nNOTE : this .wasm predates the phase kernel; phase checks skipped.');
} else {
    for (const [stackName, rows] of Object.entries(PHASE_STACKS)) {
        for (const lam of LAMBDAS) {
            const omega = omegaFromLambdaNm(lam);
            const jetFor = name => cauchyJet(lam, omega, ...DISPERSIVE[name]);
            const layers = rows.map(([name, d]) => ({ nJet: jetFor(name), d }));
            const n0Jet = jetConstant(1);
            const nsJet = cauchyJet(lam, omega, 1.5, 4.2e3, 0, 0);
            for (const th of ANGLES) {
                for (const [pol, code] of POLS) {
                    const at = `${stackName} λ=${lam} θ=${th} ${pol}`;

                    const a = tmmPhaseDispersion(lam, th, pol, n0Jet, nsJet, layers);
                    const b = k.tmmPhaseOne(lam, th, code, n0Jet, nsJet, layers);
                    comparePhaseSide(a.r, b.r, `${at} r`);
                    comparePhaseSide(a.t, b.t, `${at} t`);

                    const ja = tmmPhaseThicknessJacobian(lam, th, pol, n0Jet, nsJet, layers);
                    const jb = k.tmmPhaseJacobian(lam, th, code, n0Jet, nsJet, layers);
                    for (const s of ['r', 't']) {
                        comparePhaseSide(ja[s], jb[s], `${at} ${s} (jac)`);
                        if (!ja[s] || !jb[s]) continue;
                        for (const q of ['dPhaseDeg', 'dGd', 'dGdd', 'dTod',
                            'dLogMagnitudeSquared']) {
                            for (let i = 0; i < layers.length; i++) {
                                nearPhase(ja[s][q][i], jb[s][q][i], `${at} ${s}.${q}[${i}]`);
                            }
                        }
                    }

                    // An embedded stack at a fixed external angle: the internal
                    // sine disperses, so it arrives as a jet of its own.
                    const externalSine = Math.sin(th * Math.PI / 180);
                    const sinTheta0Jet = jetDivide(
                        jetScale(jetConstant(1), externalSine), nsJet);
                    const options = { sinTheta0Jet };
                    const ea = tmmPhaseDispersion(lam, th, pol, nsJet, n0Jet, layers, options);
                    const eb = k.tmmPhaseOne(lam, th, code, nsJet, n0Jet, layers, options);
                    comparePhaseSide(ea.r, eb.r, `${at} embedded r`);
                    comparePhaseSide(ea.t, eb.t, `${at} embedded t`);
                }
            }
        }
    }

    // Batched grid against the point evaluator, which is where a layout mistake
    // in the [layer][λ] packing would show up.
    const gridLambdas = Array.from({ length: 40 }, (_, i) => 500 + i * 7.5);
    const gridRows = PHASE_STACKS['chirped 12-layer'];
    const thick = gridRows.map(([, d]) => d);
    const omegas = gridLambdas.map(omegaFromLambdaNm);
    const jetAt = (name, i) => cauchyJet(gridLambdas[i], omegas[i], ...DISPERSIVE[name]);
    const n0Jets = gridLambdas.map(() => jetConstant(1));
    const nsJets = gridLambdas.map((lam, i) => cauchyJet(lam, omegas[i], 1.5, 4.2e3, 0, 0));
    const layerJets = gridRows.map(([name]) => gridLambdas.map((_, i) => jetAt(name, i)));
    for (const [pol, code] of POLS) {
        const batch = k.tmmPhaseSpectrum(gridLambdas, n0Jets, nsJets, layerJets, thick, 22.5, code);
        for (let i = 0; i < gridLambdas.length; i++) {
            const layers = gridRows.map(([name, d], j) => ({ nJet: layerJets[j][i], d }));
            const point = tmmPhaseDispersion(gridLambdas[i], 22.5, pol,
                n0Jets[i], nsJets[i], layers);
            for (const s of ['r', 't']) {
                for (const key of ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared']) {
                    nearPhase(point[s][key], batch[s][key][i],
                        `batch ${pol} λ=${gridLambdas[i]} ${s}.${key}`);
                }
            }
        }
    }

    // The batched Jacobian against the point Jacobian, on a grid and on a stack
    // holding a zero-thickness layer whose derivative slot must survive.
    if (!k.hasPhaseJacobianSpectrum()) {
        console.log('\nNOTE : this .wasm predates the batched phase Jacobian; its checks skipped.');
    } else {
        for (const stackName of ['chirped 12-layer', 'with a zero layer']) {
            const rows = PHASE_STACKS[stackName];
            const stackThick = rows.map(([, d]) => d);
            const stackJets = rows.map(([name]) => gridLambdas.map((_, i) => jetAt(name, i)));
            for (const [pol, code] of POLS) {
                const batch = k.tmmPhaseJacobianSpectrum(
                    gridLambdas, n0Jets, nsJets, stackJets, stackThick, 22.5, code);
                for (let i = 0; i < gridLambdas.length; i++) {
                    const layers = rows.map(([, d], j) => ({ nJet: stackJets[j][i], d }));
                    const point = tmmPhaseThicknessJacobian(gridLambdas[i], 22.5, pol,
                        n0Jets[i], nsJets[i], layers);
                    for (const s of ['r', 't']) {
                        const at = `batched jacobian ${stackName} ${pol} λ=${gridLambdas[i]} ${s}`;
                        for (const key of ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared']) {
                            nearPhase(point[s][key], batch[s][key][i], `${at}.${key}`);
                        }
                        for (const q of ['dPhaseDeg', 'dGd', 'dGdd', 'dTod',
                            'dLogMagnitudeSquared']) {
                            for (let j = 0; j < layers.length; j++) {
                                nearPhase(point[s][q][j], batch[s][q][i * layers.length + j],
                                    `${at}.${q}[${j}]`);
                            }
                        }
                    }
                }
            }
        }
    }

    // The rescaled path. The matrix product of this stack runs past the rescale
    // threshold. The Jacobian carries binary exponents through its prefix and
    // suffix products, the point evaluator rescales its one product, and the
    // two reach the coefficient by different roundings. The derivatives must
    // come back, and the point call, the kernel and the batched kernel must
    // agree on them.
    //
    // The q-th frequency order of each coefficient is a difference of terms of
    // order GD^q, GD the group delay through the whole stack, about 800 fs
    // here; its thickness derivative, of terms of order |Q|·GD^q with Q = dδ/dd
    // of the opaque layers. That cancellation, not round-off in the result,
    // sets how closely two roundings of the same quantity can agree, and the
    // derivatives of the opaque layers are zero up to exactly that noise.
    {
        const lam = 550;
        const omega = omegaFromLambdaNm(lam);
        const n0Jet = jetConstant(1);
        const nsJet = cauchyJet(lam, omega, 1.5, 4.2e3, 0, 0);
        const rows = OVERFLOW_STACK;
        const stackThick = rows.map(([, d]) => d);
        const layers = rows.map(([name, d]) => ({
            nJet: cauchyJet(lam, omega, ...DISPERSIVE[name]), d,
        }));
        const point = tmmPhaseDispersion(lam, 0, 's', n0Jet, nsJet, layers);
        const jacobian = tmmPhaseThicknessJacobian(lam, 0, 's', n0Jet, nsJet, layers);
        const kernel = k.tmmPhaseJacobian(lam, 0, 0, n0Jet, nsJet, layers);
        const groupDelay = Math.abs(point.t.gd);
        const opaque = DISPERSIVE.K;
        const phasePerNm = 2 * Math.PI * Math.hypot(opaque[0], opaque[3]) / lam;
        const cancellation = 64 * Number.EPSILON;
        const BASE = { phaseRad: 0, phaseDeg: 0, gd: 1, gdd: 2, tod: 3, magnitudeSquared: 0 };
        const DERIVATIVES = { dPhaseDeg: 0, dGd: 1, dGdd: 2, dTod: 3, dLogMagnitudeSquared: 0 };
        const baseFloor = order => Math.max(ABS_DERIV, cancellation * groupDelay ** order);
        const derivativeFloor = order =>
            Math.max(ABS_DERIV, cancellation * phasePerNm * groupDelay ** order);
        // [side, quantity, frequency order, layer] for every derivative entry.
        const entries = ['r', 't'].flatMap(side => Object.entries(DERIVATIVES).flatMap(
            ([quantity, order]) => layers.map((_, j) => [side, quantity, order, j])));
        for (const side of ['r', 't']) {
            for (const [key, order] of Object.entries(BASE)) {
                nearPhase(point[side][key], jacobian[side][key], `rescaled ${side} base.${key}`,
                    baseFloor(order));
                nearPhase(jacobian[side][key], kernel[side][key], `rescaled ${side} kernel.${key}`,
                    baseFloor(order));
            }
        }
        for (const [side, quantity, order, j] of entries) {
            nearPhase(jacobian[side][quantity][j], kernel[side][quantity][j],
                `rescaled ${side}.${quantity}[${j}]`, derivativeFloor(order));
        }

        if (k.hasPhaseJacobianSpectrum()) {
            const batch = k.tmmPhaseJacobianSpectrum(
                [lam], [n0Jet], [nsJet],
                rows.map(([name]) => [cauchyJet(lam, omega, ...DISPERSIVE[name])]),
                stackThick, 0, 0);
            for (const side of ['r', 't']) {
                nearPhase(kernel[side].magnitudeSquared, batch[side].magnitudeSquared[0],
                    `rescaled batched ${side}.magnitudeSquared`);
            }
            for (const [side, quantity, order, j] of entries) {
                nearPhase(kernel[side][quantity][j], batch[side][quantity][j],
                    `rescaled batched ${side}.${quantity}[${j}]`, derivativeFloor(order));
            }
        }
    }

    // The point path skips negative thickness. The Jacobian must report the
    // same base coefficient and reserve a zero derivative at that layer index.
    {
        const lam = 550;
        const omega = omegaFromLambdaNm(lam);
        const n0Jet = jetConstant(1);
        const nsJet = cauchyJet(lam, omega, 1.5, 4.2e3, 0, 0);
        const layers = PHASE_STACKS['with a negative layer']
            .map(([name, d]) => ({ nJet: cauchyJet(lam, omega, ...DISPERSIVE[name]), d }));
        const point = tmmPhaseDispersion(lam, 0, 's', n0Jet, nsJet, layers);
        const jacobian = tmmPhaseThicknessJacobian(lam, 0, 's', n0Jet, nsJet, layers);
        for (const side of ['r', 't']) {
            comparePhaseSide(point[side], jacobian[side], `negative thickness ${side} base`);
            for (const quantity of ['dPhaseDeg', 'dGd', 'dGdd', 'dTod',
                'dLogMagnitudeSquared']) {
                nearPhase(jacobian[side][quantity][1], 0,
                    `negative thickness ${side}.${quantity}[1]`);
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
if (phaseChecks) {
    console.log(`\n${phaseChecks} phase comparisons (phase, GD, GDD, TOD and their ` +
                `thickness derivatives, point and batched).`);
    console.log(`worst |Δ| on phase quantities: ${worstPhase.abs.toExponential(2)}`);
    console.log(`worst relative, above the ${ABS_DERIV} floor: ` +
                (worstPhase.rel === 0 ? 'none exceeded the floor' : worstPhase.rel.toExponential(2)));
}
console.log(failures === 0 ? '\nPASS : JavaScript and WebAssembly agree.'
                           : `\nFAIL : ${failures} comparisons outside tolerance.`);
process.exit(failures === 0 ? 0 : 1);
