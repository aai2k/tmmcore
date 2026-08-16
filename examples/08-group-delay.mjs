/**
 * Group delay, GDD and TOD, checked against a closed form.
 *
 * A slab whose surrounding media share its own index has no interfaces at all,
 * so nothing reflects and the transmission coefficient is a pure phase,
 * t = exp(iδ) with δ = n(ω)·ω·d/c. Its derivatives are then exactly the bulk
 * propagation result:
 *
 *     GD  = δ'   = (d/c)(n + ω n')
 *     GDD = δ''  = (d/c)(2n' + ω n'')
 *     TOD = δ''' = (d/c)(3n'' + ω n''')
 *
 * That gives a check on the phase kernel that does not go through another
 * transfer-matrix calculation. A Cauchy index makes the closed form tidy: with
 * n(λ) = A + B/λ² and λ = 2πc/ω, the index is exactly quadratic in ω,
 * n(ω) = A + βω² with β = B/(2πc)², so n' = 2βω, n'' = 2β and n''' = 0.
 *
 * The second half shows what the kernel is actually for: reading GDD off a
 * chirped mirror, and getting its thickness gradient in the same call.
 *
 *   node examples/08-group-delay.mjs
 */

import {
    C_NM_PER_FS, omegaFromLambdaNm, tmmPhaseDispersion, tmmPhaseThicknessJacobian,
    jetAdd, jetConstant, jetMultiply, jetReciprocal, jetScale, wavelengthOmegaJet,
} from '../src/index.js';

// ── Building an index jet ────────────────────────────────────────────────────
// n(λ) = A + B/λ², composed through the jet arithmetic. Nothing here is
// differentiated by hand: wavelengthOmegaJet carries the λ(ω) derivatives and
// every operation after it propagates them.

function cauchyJet(lambda_nm, omega, A, B, k = 0) {
    const lambda = wavelengthOmegaJet(lambda_nm, omega);
    const inverseSquare = jetReciprocal(jetMultiply(lambda, lambda));
    return jetAdd(jetConstant(A, k), jetScale(inverseSquare, B));
}

// ── Part 1: the matched slab ─────────────────────────────────────────────────

const A = 1.45, B = 3.6e3;          // fused-silica-ish Cauchy, B in nm²
const thickness = 1000;             // 1 µm, in nm
const lambda = 800;
const omega = omegaFromLambdaNm(lambda);

const nJet = cauchyJet(lambda, omega, A, B);
const { t } = tmmPhaseDispersion(lambda, 0, 's', nJet, nJet, [{ nJet, d: thickness }]);

const beta = B / (2 * Math.PI * C_NM_PER_FS) ** 2;
const overC = thickness / C_NM_PER_FS;
const exact = {
    gd: overC * (A + 3 * beta * omega * omega),
    gdd: overC * (6 * beta * omega),
    tod: overC * (6 * beta),
};

console.log(`matched slab, ${thickness / 1000} µm of n(λ) = ${A} + ${B}/λ², at ${lambda} nm\n`);
console.log('            tmmcore            closed form         relative');
for (const [key, unit] of [['gd', 'fs'], ['gdd', 'fs²'], ['tod', 'fs³']]) {
    const relative = Math.abs(t[key] - exact[key]) / Math.abs(exact[key]);
    console.log(`  ${key.toUpperCase().padEnd(4)} ${t[key].toPrecision(15).padStart(20)}` +
        ` ${exact[key].toPrecision(15).padStart(20)}   ${relative.toExponential(1)}  ${unit}`);
}

// |t|² is 1 because there is nothing to reflect from.
console.log(`\n  |t|² = ${t.magnitudeSquared.toPrecision(17)}  (no interfaces, so nothing is lost)`);

// How far that holds. TOD is a difference of terms of order (group delay × ω)³,
// so the cancellation in it worsens as the square of the thickness. Coatings sit
// at the top of this table; the bottom is millimetres of glass, where a
// dedicated bulk-propagation formula is the better tool.
console.log('\nhow the closed form and the kernel part company as the slab thickens:\n');
console.log('      d        GD (fs)    relative error on TOD');
for (const micron of [1, 10, 100, 1000]) {
    const d = micron * 1000;
    const one = tmmPhaseDispersion(lambda, 0, 's', nJet, nJet, [{ nJet, d }]).t;
    const reference = (d / C_NM_PER_FS) * 6 * beta;
    const label = micron < 1000 ? `${micron} µm` : '1 mm';
    console.log(`  ${label.padStart(7)}  ${((d / C_NM_PER_FS) * (A + 3 * beta * omega * omega))
        .toFixed(1).padStart(10)}    ${(Math.abs(one.tod - reference) / reference).toExponential(1)}`);
}

// ── Part 2: a chirped mirror ─────────────────────────────────────────────────
// Layer thicknesses ramp along the stack, so different wavelengths turn around
// at different depths and the reflected red is delayed relative to the blue.
// That is negative GDD, which is what compresses a pulse.

const H = l => cauchyJet(l, omegaFromLambdaNm(l), 2.21, 2.4e4);
const L = l => cauchyJet(l, omegaFromLambdaNm(l), 1.44, 3.6e3);
const air = jetConstant(1);
const substrate = l => cauchyJet(l, omegaFromLambdaNm(l), 1.5, 4.2e3);

const stack = Array.from({ length: 24 }, (_, i) => ({
    material: i % 2 ? L : H,
    d: (i % 2 ? 132 : 85) * (1 + 0.011 * i),
}));
const layersAt = l => stack.map(layer => ({ nJet: layer.material(l), d: layer.d }));

console.log('\n24-layer chirped mirror, normal incidence\n');
console.log('   λ (nm)      R        GD (fs)   GDD (fs²)   TOD (fs³)');
for (const l of [700, 750, 800, 850, 900]) {
    const { r } = tmmPhaseDispersion(l, 0, 's', air, substrate(l), layersAt(l));
    console.log(`  ${String(l).padStart(6)}  ${r.magnitudeSquared.toFixed(4)}` +
        `  ${r.gd.toFixed(2).padStart(9)}  ${r.gdd.toFixed(1).padStart(9)}` +
        `  ${r.tod.toFixed(0).padStart(10)}`);
}

// The thickness gradient of GDD comes out of the same matrix product, so a
// design step costs one evaluation rather than one per layer.
const { r } = tmmPhaseThicknessJacobian(800, 0, 's', air, substrate(800), layersAt(800));
const ranked = [...r.dGdd.keys()]
    .sort((a, b) => Math.abs(r.dGdd[b]) - Math.abs(r.dGdd[a]))
    .slice(0, 5);

console.log('\nlayers that move GDD at 800 nm the most, from ∂GDD/∂d:\n');
console.log('   layer   d (nm)    ∂GDD/∂d (fs²/nm)');
for (const i of ranked) {
    console.log(`  ${String(i).padStart(6)}  ${stack[i].d.toFixed(1).padStart(7)}` +
        `  ${r.dGdd[i].toFixed(3).padStart(17)}`);
}

// A finite difference on one layer, as a check that the gradient means what it
// says. Small enough to be linear, large enough to stay off the rounding floor.
const probe = ranked[0];
const step = 0.01;
const shifted = layersAt(800);
shifted[probe] = { ...shifted[probe], d: shifted[probe].d + step };
const moved = tmmPhaseDispersion(800, 0, 's', air, substrate(800), shifted).r;
const difference = (moved.gdd - r.gdd) / step;

console.log(`\n  analytic  ∂GDD/∂d[${probe}] = ${r.dGdd[probe].toPrecision(10)}`);
console.log(`  difference quotient      = ${difference.toPrecision(10)}   (h = ${step} nm)`);
