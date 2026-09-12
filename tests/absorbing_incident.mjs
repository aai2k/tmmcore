/**
 * An absorbing incident medium at oblique incidence.
 *
 * Every wave in the stack shares the real transverse invariant n0 sinθ0, with
 * n0 the real part of the incident index. With that, the only way R + T
 * departs from 1 over lossless layers is the interference of the incident and
 * reflected waves inside the absorbing medium, which the transmittance
 * definition of Macleod's Eq. 2.83 carries:
 *
 *     1 − R − T = −2 (Im η0 / Re η0) Im(r)      (in this module's sign convention)
 *
 * exact at every angle and polarization, and second order in k0 for a bare
 * interface, where Im(r) is itself of order k0. Carrying the complex index
 * into the invariant instead makes the incident wave's amplitude vary along
 * the interface, energy flows sideways inside lossless layers, and R + T
 * exceeds 1 by a further term linear in k0 that grows with the angle and with
 * the stack's resonance.
 *
 *  1. The identity holds to round-off for lossless stacks at 0, 30 and 60
 *     degrees, both polarizations, incident k from 1e-6 to 1e-2, and the
 *     complex invariant would have broken it.
 *  2. A bare interface between an absorbing medium and its conjugate at
 *     normal incidence gives T = 1 + k0²/n0² and R = k0²/n0² (Eq. 2.84).
 *  3. A transparent incident medium, and normal incidence, are unchanged bit
 *     for bit.
 *  4. The WebAssembly kernel agrees with the JavaScript for an absorbing
 *     incident medium: R/T/A, the thickness Jacobian, the needle scan, the
 *     phase quantities and their thickness Jacobian, with and without a
 *     dispersive incident sine. Skipped when the kernel has not been built.
 *
 *   node tests/absorbing_incident.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import {
    tmm, tmmThicknessJacobian, tmmNeedleScan,
    cadd, csub, cmul, cdiv, csqrt, cabs2, matmul, layerMatrix, rescaleMatrix,
    snellCosTheta, incidentCosTheta,
} from '../src/tmm.js';
import { tmmPhaseDispersion, tmmPhaseThicknessJacobian } from '../src/phase.js';
import { jetConstant, jetDivide, jetScale } from '../src/taylorJet.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));

const near = (actual, expected, tolerance, message) => assert.ok(
    Math.abs(actual - expected) <= tolerance,
    `${message}: got ${actual}, expected ${expected}, tolerance ${tolerance}`);

// The kernel's own coefficients, from its exported primitives, so r and η0
// are available; `cosineOf` chooses the incident-side geometry.
function coefficients(lam, thetaDeg, pol, n0, ns, layers, cosineOf) {
    const sin0 = [Math.sin(thetaDeg * Math.PI / 180), 0];
    const eta = (n, cos) => (pol === 's' ? cmul(n, cos) : cdiv(n, cos));
    const eta0 = eta(n0, cosineOf.incident(n0, sin0));
    const etaS = eta(ns, cosineOf.medium(n0, sin0, ns));
    let M = [[[1, 0], [0, 0]], [[0, 0], [1, 0]]];
    let logScale = 0;
    for (const { n, d } of layers) {
        M = matmul(M, layerMatrix(n, d, lam, cosineOf.medium(n0, sin0, n), pol));
        logScale += rescaleMatrix(M);
    }
    const B = cadd(M[0][0], cmul(M[0][1], etaS));
    const C = cadd(M[1][0], cmul(M[1][1], etaS));
    const den = cadd(cmul(eta0, B), C);
    const r = cdiv(csub(cmul(eta0, B), C), den);
    const t = cdiv(cmul([2, 0], eta0), den);
    return { r, eta0, R: cabs2(r), T: etaS[0] / eta0[0] * cabs2(t) * Math.exp(-2 * logScale) };
}
const plainCosine = sin0 => csqrt(csub([1, 0], cmul(sin0, sin0)));
// The kernel's geometry: real invariant everywhere.
const realInvariant = { incident: incidentCosTheta, medium: snellCosTheta };
// The former geometry: the complex index carried into the invariant.
const complexInvariant = {
    incident: (n0, sin0) => plainCosine(sin0),
    medium: (n0, sin0, nj) => {
        const s = cdiv(cmul(n0, sin0), nj);
        return csqrt(csub([1, 0], cmul(s, s)));
    },
};

const H = [2.3, 0], L = [1.46, 0], GLASS = [1.52, 0];
const STACKS = {
    'AR 4-layer': [[H, 110], [L, 95], [H, 60], [L, 130]],
    'QW 21-layer': Array.from({ length: 21 }, (_, i) => (i % 2 ? [L, 94.2] : [H, 58.5])),
};
const mk = rows => rows.map(([n, d]) => ({ n, d }));

// ── 1. The energy identity ───────────────────────────────────────────────────

let worstGap = 0;
let formerWorst = 0;
for (const [name, rows] of Object.entries(STACKS)) {
    const layers = mk(rows);
    for (const k0 of [1e-6, 1e-4, 1e-2]) {
        const n0 = [1.5, k0];
        for (const theta of [0, 30, 60]) {
            for (const pol of ['s', 'p']) {
                const at = `${name} k0=${k0} ${theta}° ${pol}`;
                const { r, eta0, R, T } = coefficients(550, theta, pol, n0, GLASS, layers, realInvariant);
                const kernel = tmm(550, theta, pol, n0, GLASS, layers);
                near(kernel.R, R, 1e-15, `${at}: kernel R against the primitives`);
                near(kernel.T, T, 1e-15, `${at}: kernel T against the primitives`);
                const interference = -2 * (eta0[1] / eta0[0]) * r[1];
                const gap = Math.abs(1 - R - T - interference);
                worstGap = Math.max(worstGap, gap);
                assert.ok(gap <= 1e-14,
                    `${at}: 1 − R − T = ${1 - R - T}, interference term ${interference}, gap ${gap}`);
                if (theta > 0) {
                    const former = coefficients(550, theta, pol, n0, GLASS, layers, complexInvariant);
                    const formerTerm = -2 * (former.eta0[1] / former.eta0[0]) * former.r[1];
                    formerWorst = Math.max(formerWorst, Math.abs(1 - former.R - former.T - formerTerm) / k0);
                }
            }
        }
    }
}
assert.ok(formerWorst > 0.1, `the complex invariant breaks the identity by a term linear in k0 (${formerWorst} per unit k0)`);

// ── 2. The bare interface of Eq. 2.82 ───────────────────────────────────────

for (const [n, k] of [[1.5, 0.01], [1.5, 0.1], [2, 0.3]]) {
    const { R, T } = tmm(550, 0, 's', [n, k], [n, -k], []);
    near(T, 1 + (k / n) ** 2, 1e-15, `T into the conjugate admittance, n0 = ${n} + ${k}i`);
    near(R, (k / n) ** 2, 1e-15, `R from the conjugate admittance, n0 = ${n} + ${k}i`);
}

// ── 3. Nothing moves for a transparent incident medium or at normal incidence

{
    const sin0 = [Math.sin(Math.PI / 4), 0];
    assert.deepEqual(incidentCosTheta([1.5, 0], sin0), plainCosine(sin0),
        'a transparent incident medium keeps the plain cosine');
    assert.deepEqual(snellCosTheta([1.5, 0], sin0, [2.3, 0.001]),
        complexInvariant.medium([1.5, 0], sin0, [2.3, 0.001]),
        'and the layers keep the former arithmetic exactly');
    assert.deepEqual(incidentCosTheta([1.5, 0.02], [0, 0]), [1, 0],
        'an absorbing incident medium at normal incidence has cosθ0 = 1 exactly');
    const layers = mk(STACKS['AR 4-layer']);
    const a = tmm(550, 0, 'p', [1.5, 0.02], GLASS, layers);
    const b = coefficients(550, 0, 'p', [1.5, 0.02], GLASS, layers, complexInvariant);
    assert.equal(a.R, b.R, 'normal incidence into an absorbing medium is unchanged');
}

// ── 4. The WebAssembly kernel ───────────────────────────────────────────────

const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');
let wasmChecks = 0;
if (existsSync(WASM)) {
    const k = await instantiateTmmWasm(readFileSync(WASM));
    assert.ok(k, 'the kernel instantiates');
    const close = (js, wa, label) => {
        wasmChecks++;
        const diff = Math.abs(js - wa);
        const rel = diff / Math.max(Math.abs(js), Math.abs(wa), Number.MIN_VALUE);
        assert.ok(diff <= 1e-12 || rel <= 1e-7, `${label}: js=${js} wasm=${wa}`);
    };
    const n0 = [1.5, 0.02];
    const layers = mk([[[2.3, 0.0005], 110], [L, 95], [[0.15, 3.2], 18], [L, 130], [H, 60]]);
    for (const theta of [0, 30, 60]) {
        for (const [pol, code] of [['s', 0], ['p', 1]]) {
            const at = `θ=${theta} ${pol}`;
            const a = tmm(550, theta, pol, n0, GLASS, layers);
            const b = k.tmmOne(550, theta, code, n0, GLASS, layers);
            for (const q of ['R', 'T', 'A']) close(a[q], b[q], `${at} ${q}`);
            const ja = tmmThicknessJacobian(550, theta, pol, n0, GLASS, layers);
            const jb = k.tmmJacobian(550, theta, code, n0, GLASS, layers);
            for (const q of ['dRdd', 'dTdd', 'dAdd']) {
                for (let i = 0; i < layers.length; i++) close(ja[q][i], jb[q][i], `${at} ${q}[${i}]`);
            }
            const na = tmmNeedleScan(550, theta, pol, n0, GLASS, layers, [H, L]);
            const nb = k.tmmNeedleScan(550, theta, code, n0, GLASS, layers, [H, L]);
            for (let pos = 0; pos < na.gaps.length; pos++) {
                for (let c = 0; c < 2; c++) {
                    for (const q of ['dR', 'dT', 'dA']) close(na.gaps[pos][c][q], nb.gaps[pos][c][q], `${at} needle[${pos}][${c}] ${q}`);
                }
            }
            if (k.hasPhase()) {
                const n0Jet = jetConstant(1.5, 0.02);
                const nsJet = jetConstant(1.52);
                const jetLayers = layers.map(({ n, d }) => ({ nJet: jetConstant(n[0], n[1]), d }));
                // Once at the external angle, once embedded behind a dispersive
                // face so the incident sine arrives as a jet.
                const sinTheta0Jet = jetDivide(jetScale(jetConstant(1), Math.sin(theta * Math.PI / 180)), n0Jet);
                for (const [label, options] of [['', {}], [' embedded', { sinTheta0Jet }]]) {
                    const pa = tmmPhaseDispersion(550, theta, pol, n0Jet, nsJet, jetLayers, options);
                    const pb = k.tmmPhaseOne(550, theta, code, n0Jet, nsJet, jetLayers, options);
                    const ja = tmmPhaseThicknessJacobian(550, theta, pol, n0Jet, nsJet, jetLayers, options);
                    const jb = k.tmmPhaseJacobian(550, theta, code, n0Jet, nsJet, jetLayers, options);
                    for (const side of ['r', 't']) {
                        for (const q of ['phaseRad', 'gd', 'gdd', 'tod', 'magnitudeSquared']) {
                            close(pa[side][q], pb[side][q], `${at}${label} phase ${side}.${q}`);
                        }
                        for (const q of ['dPhaseDeg', 'dGd', 'dGdd', 'dTod']) {
                            for (let i = 0; i < layers.length; i++) {
                                close(ja[side][q][i], jb[side][q][i], `${at}${label} phase jacobian ${side}.${q}[${i}]`);
                            }
                        }
                    }
                }
            }
        }
    }
} else {
    console.log('NOTE : src/tmm_kernel.wasm not built; the WebAssembly checks were skipped.');
}

console.log(`worst gap in the energy identity: ${worstGap.toExponential(2)}`);
console.log(`the former complex invariant broke it by ${formerWorst.toExponential(2)} per unit k0`);
console.log(`WebAssembly comparisons: ${wasmChecks}`);
console.log('PASS : absorbing incident medium');
