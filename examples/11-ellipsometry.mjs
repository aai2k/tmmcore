/**
 * Ellipsometric Ψ and Δ of an oxide on silicon, checked against closed forms.
 *
 * An ellipsometer measures the ratio of the p and s reflection coefficients,
 * ρ_p/ρ_s = tan Ψ·e^{iΔ}. As a transparent film grows on a known substrate the
 * pair (Ψ, Δ) runs round a closed loop and returns to the bare-substrate point
 * after one period of thickness. Here the film is SiO2 on silicon, at 632.8 nm
 * and 70° in air, from zero thickness through one period.
 *
 * Δ follows Macleod, Thin-Film Optical Filters 5th ed., §9.7.1, Eqs. 9.35 and 9.36:
 *
 *     tan Ψ = |ρ_p/ρ_s|,    Δ = φ_p − φ_s ± π  in reflection,
 *
 * with φ the phase of ρ in his convention, ñ = n − ik and exp(+iωt), and ρ_p from
 * the tilted admittance n/cos θ, as tmmcore's is. The ±π stands in for the flip
 * of the reflected p direction that ellipsometry makes and thin-film work does
 * not, so that Δ is 180° at normal incidence. Δ is printed in 0-360°.
 *
 * tmmcore's r is the complex conjugate of Macleod's ρ, and phase() reports
 * phaseRad as minus the argument of r, which is φ with Macleod's sign. The
 * conjugation is therefore already made, and Δ = phaseRad_p − phaseRad_s + π.
 *
 * Indices at 632.8 nm:
 *   silicon   Aspnes and Studna, Phys. Rev. B 27, 985 (1983), interpolated
 *             linearly in wavelength between their points at 2.0 and 1.9 eV;
 *   oxide     fused silica, the Sellmeier formula of Malitson,
 *             J. Opt. Soc. Am. 55, 1205 (1965).
 *
 * Two checks against closed forms:
 *
 *   1. At zero thickness the stack is the bare substrate, with single-interface
 *      coefficients ρ = (η0 − η)/(η0 + η) (§2.6, Eq. 2.112). For the substrate
 *      η_s = √(n² − k² − sin²θ0 − 2ink) in the fourth quadrant and
 *      η_p = (n − ik)²/η_s (Eqs. 2.116 and 2.117); for air η_s = cos θ0 and
 *      η_p = 1/cos θ0 (Eq. 2.113).
 *
 *   2. The phase thickness of a film is δ = 2πnd·cos θ/λ (§9.2, Eq. 9.2). Adding
 *      π to δ changes only the sign of the film's characteristic matrix
 *      (Eq. 2.113), which leaves ρ = (η0B − C)/(η0B + C) as it was; at δ = π the
 *      matrix is −I and the film is absent (§3.1). So for a film without loss,
 *      Ψ and Δ repeat with period d = λ/(2√(n² − sin²θ0)), from
 *      n cos θ = √(n² − sin²θ0) by Snell's law in air (Eq. 9.5).
 *
 * Writes the plot to docs/img/ellipsometry.svg.
 *
 *   node examples/11-ellipsometry.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { phase } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const LAMBDA = 632.8;          // nm
const ANGLE = 70;              // degrees, in air

// Silicon: the two tabulated points either side of 632.8 nm, [λ in nm, n, k].
const SI_BELOW = [619.9, 3.906, 0.022];
const SI_ABOVE = [652.5, 3.847, 0.016];
const fraction = (LAMBDA - SI_BELOW[0]) / (SI_ABOVE[0] - SI_BELOW[0]);
const nSi = SI_BELOW[1] + fraction * (SI_ABOVE[1] - SI_BELOW[1]);
const kSi = SI_BELOW[2] + fraction * (SI_ABOVE[2] - SI_BELOW[2]);

// Fused silica: n² − 1 = Σ Bᵢλ²/(λ² − Cᵢ²), λ in µm.
const MALITSON = [[0.6961663, 0.0684043], [0.4079426, 0.1162414], [0.8974794, 9.896161]];
const um = LAMBDA / 1000;
const nOx = Math.sqrt(1 + MALITSON.reduce((sum, [B, C]) => sum + B * um * um / (um * um - C * C), 0));

const sin0 = Math.sin(ANGLE * Math.PI / 180);
const cos0 = Math.cos(ANGLE * Math.PI / 180);
const PERIOD = LAMBDA / (2 * Math.sqrt(nOx * nOx - sin0 * sin0));

const toDegrees = rad => rad * 180 / Math.PI;
const wrap360 = deg => ((deg % 360) + 360) % 360;

// Ψ and Δ in degrees for an oxide of thickness d (nm).
function ellipsometry(d) {
    const stack = { incident: 1, substrate: { n: nSi, k: kSi }, layers: [{ n: nOx, d }] };
    const s = phase(stack, LAMBDA, { pol: 's', angle: ANGLE }).r;
    const p = phase(stack, LAMBDA, { pol: 'p', angle: ANGLE }).r;
    return {
        psi: toDegrees(Math.atan2(Math.sqrt(p.magnitudeSquared), Math.sqrt(s.magnitudeSquared))),
        delta: wrap360(toDegrees(p.phaseRad - s.phaseRad + Math.PI)),
    };
}

// ── Check 1: the bare substrate, in Macleod's convention ─────────────────────
// Complex numbers as [re, im]. sqrt takes the principal root, which puts the
// root of n² − k² − sin²θ0 − 2ink, whose imaginary part is negative, in the
// fourth quadrant as Eq. 2.116 asks.
const mul = ([a, b], [c, d]) => [a * c - b * d, a * d + b * c];
const div = ([a, b], [c, d]) => { const m = c * c + d * d; return [(a * c + b * d) / m, (b * c - a * d) / m]; };
const sqrt = ([a, b]) => {
    const r = Math.hypot(a, b);
    const re = Math.sqrt((r + a) / 2);
    return [re, b < 0 ? -Math.sqrt((r - a) / 2) : Math.sqrt((r - a) / 2)];
};
const fresnel = (eta0, eta) => div([eta0[0] - eta[0], eta0[1] - eta[1]], [eta0[0] + eta[0], eta0[1] + eta[1]]);

const N = [nSi, -kSi];                                   // n − ik
const etaS = sqrt([nSi * nSi - kSi * kSi - sin0 * sin0, -2 * nSi * kSi]);
const etaP = div(mul(N, N), etaS);
const rhoS = fresnel([cos0, 0], etaS);
const rhoP = fresnel([1 / cos0, 0], etaP);
const bare = {
    psi: toDegrees(Math.atan2(Math.hypot(...rhoP), Math.hypot(...rhoS))),
    delta: wrap360(toDegrees(Math.atan2(rhoP[1], rhoP[0]) - Math.atan2(rhoS[1], rhoS[0]) + Math.PI)),
};

// ── Check 2: one period on, at zero thickness and at a thickness inside the loop
const atZero = ellipsometry(0);
const atPeriod = ellipsometry(PERIOD);
const PROBE = 100;             // nm
const atProbe = ellipsometry(PROBE);
const afterPeriod = ellipsometry(PROBE + PERIOD);

console.log(`SiO2 on Si at ${LAMBDA} nm, ${ANGLE}° in air`);
console.log(`  Si    ñ = ${nSi.toFixed(4)} + ${kSi.toFixed(4)}i   SiO2  n = ${nOx.toFixed(5)}`);
console.log(`  one period of oxide, λ/(2√(n² − sin²θ)) = ${PERIOD.toFixed(2)} nm\n`);

console.log('   d (nm)      Ψ (°)      Δ (°)');
for (let d = 0; d <= PERIOD; d += 20) {
    const { psi, delta } = ellipsometry(d);
    console.log(`  ${String(d).padStart(7)}  ${psi.toFixed(3).padStart(9)}  ${delta.toFixed(3).padStart(9)}`);
}

const row = (label, a, b, digits) =>
    console.log(`  ${label.padEnd(26)}${a.toFixed(digits).padStart(16)}${b.toFixed(digits).padStart(18)}`);
const gap = (label, a, b) => console.log(`  ${label.padEnd(26)}` +
    `${Math.abs(a.psi - b.psi).toExponential(1).padStart(16)}${Math.abs(a.delta - b.delta).toExponential(1).padStart(18)}`);
console.log(`\n${'checks'.padEnd(28)}${'Ψ (°)'.padStart(16)}${'Δ (°)'.padStart(18)}`);
row('bare Si, tmmcore', atZero.psi, atZero.delta, 12);
row('bare Si, closed form', bare.psi, bare.delta, 12);
gap('difference', atZero, bare);
row('one period of oxide', atPeriod.psi, atPeriod.delta, 12);
gap('difference from bare Si', atPeriod, atZero);
row(`${PROBE} nm of oxide`, atProbe.psi, atProbe.delta, 12);
row(`${PROBE} nm plus one period`, afterPeriod.psi, afterPeriod.delta, 12);
gap('difference', afterPeriod, atProbe);

// The trajectory for the figure. Near 142 nm Ψ comes close to 90°, |r_s| is
// small, and Δ moves by more than 90° within 3 nm, so the steps are fine.
const psis = [], deltas = [];
const STEPS = 1500;
for (let i = 0; i <= STEPS; i++) {
    const { psi, delta } = ellipsometry((i / STEPS) * PERIOD);
    psis.push(psi); deltas.push(delta);
}

// The Ψ axis spans 0-90°, the whole range Ψ can take.
mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'ellipsometry.svg'), linePlot({
    title: `SiO2 on Si through one period of oxide, ${LAMBDA} nm at ${ANGLE}°`,
    xLabel: 'Ψ (degrees)',
    yLabel: 'Δ (degrees)',
    xMin: 0, xMax: 90, yMin: 50, yMax: 300,
    series: [{ x: psis, y: deltas }],
}));
console.log('\nwrote docs/img/ellipsometry.svg');
