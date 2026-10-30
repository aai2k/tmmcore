/**
 * tmmcore against Essential Macleod.
 *
 * tests/macleod/*.json hold Essential Macleod's results for ten designs built
 * from constant-index materials of this test's own, written through its COM
 * server by tests/macleod/export.ps1. Each file names the program version,
 * the media and the layers, and for every angle and polarization a wavelength
 * grid of the quantities below. tmmcore computes the same designs and every
 * number is compared:
 *
 *  1. R, T and A in s and p (rtaSpectrum), and for natural light (pol 'avg').
 *  2. The phase of r and t, and the GD, GDD and TOD of each (phaseSpectrum).
 *  3. Psi and Delta of r and of t.
 *  4. The reflectance from the substrate side, which Essential Macleod
 *     defines at the angle θ0 taken in the substrate: rtaSpectrum of the
 *     reversed stack, and at normal incidence monitorCurve's Rrs and Rrp.
 *  5. The admittance loci, the reflectance with the coating cut at a depth:
 *     monitorCurve, growingStack and depositionSpectra of the cut coating.
 *
 * Every item runs on the JavaScript and on the WebAssembly kernel; the
 * WebAssembly checks are skipped when the kernel has not been built.
 *
 * The exports write indices as n − ik with k ≥ 0, the numbers tmmcore takes
 * as n + ik, and phases in degrees with Macleod's sign, which tmmcore's
 * phaseRad carries. Psi is atan sqrt(Rp/Rs) for r and atan sqrt(Tp/Ts) for t.
 * Delta is in the Azzam-Bashara sign: phase(r_p) − phase(r_s) + 180° for r
 * and phase(t_p) − phase(t_s) for t. Essential Macleod writes the phase of t
 * as 0 where T < 1e-28, so the phase of t is compared only above that.
 *
 * Tolerances. R, T, A and the locus reflectance agree to TOL.fraction and
 * phases, Psi and Delta to TOL.degrees. GD, GDD and TOD agree to TOL.delay
 * times the largest value of the scan, or of 1 fs, 1 fs² or 1 fs³ where that
 * is larger: a delay that is zero in exact arithmetic, such as the reflection
 * GDD of an opaque metal of constant index, comes out as round-off of about
 * 1e-11 in both programs. Each tolerance is the worst difference over the ten
 * designs, rounded up to a power of ten, times ten.
 *
 * Two places take another reference:
 *
 *  - Essential Macleod's reflection GDD and TOD lose accuracy where r nears
 *    zero. The air gap below its critical angle (frustrated-tir at 30°) passes
 *    through r = 0 at 529 nm: at 530 nm, |r|² = 4.7e-7 in p, its TOD is off
 *    the closed form below by 1.9e-6 fs³, where tmmcore's is off by 3.6e-11.
 *    That scan's reflection delays are checked against the closed form.
 *  - An absorbing incident medium. Essential Macleod takes the angles from
 *    the real part of the incident index, as tmmcore does, but it also takes
 *    the real part a of the incident tilted admittance η0, which leaves the
 *    incident medium's absorption out of its R and T. tmmcore keeps the
 *    complex η0 of Macleod's Eq. 2.83, so the two differ by a term linear in
 *    k0. The coating is the same in both: its surface admittance
 *    Y = η0(1 − r)/(1 + r), from tmmcore's r, gives Essential Macleod's r as
 *    (a − Y)/(a + Y), and its t in the same ratio, and those are compared.
 *    This covers the two absorbing-incident designs and the substrate-side
 *    reflectance of the absorbing substrate. For the GD, GDD and TOD of the
 *    two designs the conversion is carried through the frequency jets of r
 *    and t, on the JavaScript kernel alone, since the WebAssembly kernel
 *    returns no jets.
 *
 *   node tests/macleod.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import * as tmmcore from '../src/js/index.js';
import {
    C_NM_PER_FS, cadd, csub, cmul, cdiv, csqrt, cabs2, omegaFromLambdaNm, wavelengthOmegaJet,
    jetConstant, jetAdd, jetSubtract, jetMultiply, jetDivide, jetScale, jetRealPart,
} from '../src/js/index.js';
import { tmmCoefficientJets, coefficientPhaseDispersion } from '../src/js/phase.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const DIR = join(HERE, 'macleod');
const WASM = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');

const TOL = { fraction: 1e-12, degrees: 1e-10, delay: 1e-9 };
const T_PHASE_FLOOR = 1e-28;
const DEG = 180 / Math.PI;
const POLS = ['s', 'p'];
const FORWARD = { s: 'Rs', p: 'Rp' };
const REVERSE = { s: 'Rrs', p: 'Rrp' };
const DELAYS = [['gd', 'Gd', 1], ['gdd', 'Gdd', 2], ['tod', 'Tod', 3]];

const CASES = readdirSync(DIR).filter(f => f.endsWith('.json')).sort()
    .map(f => JSON.parse(readFileSync(join(DIR, f), 'utf8')));

const worst = {};
let comparisons = 0;
function check(group, diff, tolerance, at) {
    comparisons++;
    const size = Math.abs(diff);
    assert.ok(size <= tolerance, `${at}: ${group} off by ${diff}, tolerance ${tolerance}`);
    if (!(size <= (worst[group]?.size ?? -1))) worst[group] = { size, at };
}
const wrapDegrees = d => ((d % 360) + 540) % 360 - 180;

const medium = m => ({ n: m.n[0], k: m.n[1] });
const stackOf = c => ({
    incident: medium(c.incident),
    substrate: medium(c.substrate),
    layers: c.layers.map(l => ({ ...medium(l), d: l.d })),
});
const reversed = stack => ({ incident: stack.substrate, substrate: stack.incident, layers: [...stack.layers].reverse() });

// ── Essential Macleod's view of an absorbing incident medium ─────────────────
// In Macleod's convention, N0 = n0 − ik0, η0 = sqrt(N0² − n0² sin²θ0) for s
// and N0²/that for p. With a = Re η0 and Y = η0(1 − r)/(1 + r):
//   r' = (a − Y)/(a + Y) = (a(1 + r) − η0(1 − r)) / (a(1 + r) + η0(1 − r))
//   t'/t = 2a / (a(1 + r) + η0(1 − r)),   T' = T |t'/t|²
function incidentAdmittance({ n, k }, angle, pol) {
    const N0 = [n, -k];
    const sine = n * Math.sin(angle / DEG);
    const N0sq = cmul(N0, N0);
    const cosine = csqrt(csub(N0sq, [sine * sine, 0]));
    return pol === 's' ? cosine : cdiv(N0sq, cosine);
}

function macleodView(medium0, angle, pol, rt, ph) {
    if (medium0.k === 0) {
        return { R: rt.R, T: rt.T, A: rt.A, rPhase: ph.r.phaseRad, tPhase: ph.t.phaseRad, delays: ph };
    }
    const eta0 = incidentAdmittance(medium0, angle, pol);
    const a = [eta0[0], 0];
    const view = { R: [], T: [], A: rt.A, rPhase: [], tPhase: [], delays: null };
    rt.R.forEach((_, i) => {
        const size = Math.sqrt(ph.r.magnitudeSquared[i]);
        const r = [size * Math.cos(ph.r.phaseRad[i]), size * Math.sin(ph.r.phaseRad[i])];
        const onePlus = cadd([1, 0], r);
        const oneMinus = csub([1, 0], r);
        const denominator = cadd(cmul(a, onePlus), cmul(eta0, oneMinus));
        const rView = cdiv(csub(cmul(a, onePlus), cmul(eta0, oneMinus)), denominator);
        const tRatio = cdiv([2 * a[0], 0], denominator);
        view.R.push(cabs2(rView));
        view.T.push(rt.T[i] * cabs2(tRatio));
        view.rPhase.push(Math.atan2(rView[1], rView[0]));
        view.tPhase.push(ph.t.phaseRad[i] + Math.atan2(tRatio[1], tRatio[0]));
    });
    return view;
}

// The same conversion carried through the frequency jets of r and t, which
// gives GD, GDD and TOD of Essential Macleod's r and t. The jets are the
// JavaScript kernel's; tests/absorbing_incident.mjs holds the WebAssembly
// phase calls to it under an absorbing incident medium.
function macleodDelaysAt(stack, lambda, angle, pol) {
    const index = m => jetConstant(m.n, m.k);
    const jets = tmmCoefficientJets({
        wavelengthJet: wavelengthOmegaJet(lambda, omegaFromLambdaNm(lambda)),
        thetaDeg: angle,
        polarization: pol,
        incidentIndexJet: index(stack.incident),
        substrateIndexJet: index(stack.substrate),
        layers: stack.layers.map(l => ({ indexJet: index(l), thicknessNm: l.d })),
    });
    const r = jets.reflection, eta0 = jets.incidentEta;
    const a = jetRealPart(eta0);
    const one = jetConstant(1);
    const aOnePlus = jetMultiply(a, jetAdd(one, r));
    const etaOneMinus = jetMultiply(eta0, jetSubtract(one, r));
    const denominator = jetAdd(aOnePlus, etaOneMinus);
    return {
        r: coefficientPhaseDispersion(jetDivide(jetSubtract(aOnePlus, etaOneMinus), denominator)),
        t: coefficientPhaseDispersion(jetDivide(jetScale(jetMultiply(a, jets.transmission), 2), denominator),
            jets.transmissionLogScale),
    };
}

// ── Reflection delays of one film between identical media ────────────────────
// One layer of admittance η1 between media of admittance η: the matrix of
// Macleod's Eq. 2.113 gives
//   r = i(η² − η1²) sin δ / (2ηη1 cos δ + i(η² + η1²) sin δ),
// which is r01 (1 − e^{−2iδ}) / (1 − r01² e^{−2iδ}), r01 = (η − η1)/(η + η1).
// With x = 2δ = aω, a = 2 n1 d cosθ1 / c, the phase of 1 − e^{−ix} is
// −x/2 + const, so φ = const − x/2 − arg(1 − u), u = r01² e^{−ix}. With hn the
// n-th x-derivative of log(1 − u),
//   h1 = i u/(1 − u),  h2 = u/(1 − u)²,  h3 = −i u (1 + u)/(1 − u)³,
// GD = −dφ/dω = a/2 + a Im h1, GDD = a² Im h2, TOD = a³ Im h3.
function filmDelays(c, scan, pol) {
    const nOuter = c.incident.n[0], nFilm = c.layers[0].n[0];
    const cosOuter = Math.cos(scan.angle / DEG);
    const sineFilm = nOuter * Math.sin(scan.angle / DEG) / nFilm;
    const cosFilm = Math.sqrt(1 - sineFilm * sineFilm);
    const tilted = pol === 's' ? (n, cos) => n * cos : (n, cos) => n / cos;
    const etaOuter = tilted(nOuter, cosOuter), etaFilm = tilted(nFilm, cosFilm);
    const r01 = (etaOuter - etaFilm) / (etaOuter + etaFilm);
    const a = 2 * nFilm * c.layers[0].d * cosFilm / C_NM_PER_FS;
    return scan.wavelengths.map(lambda => {
        const x = a * 2 * Math.PI * C_NM_PER_FS / lambda;
        const u = [r01 * r01 * Math.cos(x), -r01 * r01 * Math.sin(x)];
        const w = csub([1, 0], u);
        const h1 = cmul([0, 1], cdiv(u, w));
        const h2 = cdiv(u, cmul(w, w));
        const h3 = cmul([0, -1], cdiv(cmul(u, cadd([1, 0], u)), cmul(w, cmul(w, w))));
        return { gd: a / 2 + a * h1[1], gdd: a * a * h2[1], tod: a * a * a * h3[1] };
    });
}
const exactReflectionDelays = (c, scan) => c.case === 'frustrated-tir' && scan.angle === 30;

// ── 1 to 4: spectra ──────────────────────────────────────────────────────────

function compareScan(api, backend, c, scan) {
    const stack = stackOf(c);
    const L = scan.wavelengths;
    const views = {};
    for (const pol of POLS) {
        const options = { angle: scan.angle, pol };
        const rt = api.rtaSpectrum(stack, L, options);
        const ph = api.phaseSpectrum(stack, L, options);
        const view = macleodView(stack.incident, scan.angle, pol, rt, ph);
        views[pol] = view;
        const em = scan[pol];
        const back = reversed(stack);
        const backView = macleodView(back.incident, scan.angle, pol,
            api.rtaSpectrum(back, L, options), api.phaseSpectrum(back, L, options));
        L.forEach((lambda, i) => {
            const at = `${backend} ${c.case} ${scan.angle}° ${pol} ${lambda} nm`;
            for (const q of ['R', 'T', 'A']) check('R, T, A', view[q][i] - em[q][i], TOL.fraction, `${at} ${q}`);
            check('substrate-side R', backView.R[i] - em.Rback[i], TOL.fraction, `${at} Rback`);
            check('phase', wrapDegrees(view.rPhase[i] * DEG - em.rPhaseDeg[i]), TOL.degrees, `${at} r phase`);
            if (em.T[i] > T_PHASE_FLOOR) {
                check('phase', wrapDegrees(view.tPhase[i] * DEG - em.tPhaseDeg[i]), TOL.degrees, `${at} t phase`);
            }
        });
        if (view.delays) compareDelays(view.delays, c, scan, pol, backend);
    }
    L.forEach((lambda, i) => {
        const at = `${backend} ${c.case} ${scan.angle}° ${lambda} nm`;
        const { s, p } = views;
        const psiR = Math.atan(Math.sqrt(p.R[i] / s.R[i])) * DEG;
        const psiT = Math.atan(Math.sqrt(p.T[i] / s.T[i])) * DEG;
        check('Psi, Delta', psiR - scan.s.rPsiDeg[i], TOL.degrees, `${at} Psi of r`);
        check('Psi, Delta', wrapDegrees((p.rPhase[i] - s.rPhase[i]) * DEG + 180 - scan.s.rDeltaDeg[i]), TOL.degrees, `${at} Delta of r`);
        if (scan.s.T[i] > T_PHASE_FLOOR) {
            check('Psi, Delta', psiT - scan.s.tPsiDeg[i], TOL.degrees, `${at} Psi of t`);
            check('Psi, Delta', wrapDegrees((p.tPhase[i] - s.tPhase[i]) * DEG - scan.s.tDeltaDeg[i]), TOL.degrees, `${at} Delta of t`);
        }
    });
    if (stack.incident.k === 0) {
        const natural = api.rtaSpectrum(stack, L, { angle: scan.angle });
        L.forEach((lambda, i) => {
            for (const q of ['R', 'T', 'A']) {
                check('R, T, A', natural[q][i] - (scan.s[q][i] + scan.p[q][i]) / 2, TOL.fraction,
                    `${backend} ${c.case} ${scan.angle}° ${lambda} nm ${q}, pol 'avg'`);
            }
        });
    }
}

function compareDelays(ph, c, scan, pol, backend) {
    const em = scan[pol];
    const L = scan.wavelengths;
    const film = exactReflectionDelays(c, scan) ? filmDelays(c, scan, pol) : null;
    for (const side of ['r', 't']) {
        for (const [q, key, power] of DELAYS) {
            const exact = side === 'r' && film !== null;
            const reference = exact ? film.map(delays => delays[q]) : em[`${side}${key}`];
            const scale = Math.max(1, ...reference.map(Math.abs));
            L.forEach((lambda, i) => {
                const at = `${backend} ${c.case} ${scan.angle}° ${pol} ${lambda} nm ${side} ${q}${exact ? ' (closed form)' : ''}`;
                check(`GD, GDD, TOD`, (ph[side][q][i] - reference[i]) / scale, TOL.delay, `${at}, fs^${power} over ${scale}`);
            });
        }
    }
}

function compareAbsorbingIncidentDelays(c) {
    const stack = stackOf(c);
    for (const [scan, pol] of c.scans.flatMap(scan => POLS.map(pol => [scan, pol]))) {
        const points = scan.wavelengths.map(lambda => macleodDelaysAt(stack, lambda, scan.angle, pol));
        const series = side => Object.fromEntries(DELAYS.map(([q]) => [q, points.map(p => p[side][q])]));
        compareDelays({ r: series('r'), t: series('t') }, c, scan, pol, 'js jets');
    }
}

// ── 5: the admittance loci ───────────────────────────────────────────────────
// Depth z is measured from the outer surface. With layer j spanning
// [top_j, top_j + d_j], the coating cut at z is layers j+1.. finished and
// layer j grown to top_j + d_j − z.

function cutAt(layers, z) {
    let top = 0;
    for (let j = 0; j < layers.length; j++) {
        if (z < top + layers[j].d) return { j, grown: top + layers[j].d - z };
        top += layers[j].d;
    }
    return { j: layers.length - 1, grown: 0 };
}

const indexOf = ({ n, k }) => ({ n, k });

function compareLoci(api, backend, c) {
    const stack = stackOf(c);
    for (const locus of c.loci) {
        const { wavelength: lambda, angle, pol } = locus;
        const key = FORWARD[pol];
        locus.depth.forEach((z, k) => {
            const at = `${backend} ${c.case} locus ${lambda} nm ${angle}° ${pol} at depth ${z} nm`;
            const expected = locus.Rpct[k] / 100;
            const { j, grown } = cutAt(stack.layers, z);
            const finished = { ...stack, layers: stack.layers.slice(j + 1) };
            const top = indexOf(stack.layers[j]);
            const curve = api.monitorCurve(finished, lambda, { angle, top, thicknesses: [grown] });
            check('locus R', curve[key][0] - expected, TOL.fraction, `${at}, monitorCurve`);
            const growing = api.growingStack(finished, [lambda], { angle });
            growing.setTop(top);
            check('locus R', growing.sample(grown)[key][0] - expected, TOL.fraction, `${at}, growingStack`);
            growing.free();
            const cut = { ...stack, layers: [{ ...top, d: grown }, ...finished.layers] };
            const steps = api.depositionSpectra(cut, [lambda], { angle });
            check('locus R', steps[key][cut.layers.length - 1] - expected, TOL.fraction, `${at}, depositionSpectra`);
        });
    }
}

// Normal incidence: the substrate-side reflectance at the same invariant,
// which is the one Essential Macleod reports when θ0 = 0. An absorbing
// substrate is the incident medium of that pass, where the two programs take
// different admittances (above), so it is left to the reversed stack.
function compareReverseAtNormalIncidence(api, backend, c) {
    if (c.substrate.n[1] !== 0) return;
    const stack = stackOf(c);
    const [outer, ...below] = stack.layers;
    const options = { top: indexOf(outer), thicknesses: [outer.d] };
    for (const scan of c.scans.filter(s => s.angle === 0)) {
        scan.wavelengths.forEach((lambda, i) => {
            const curve = api.monitorCurve({ ...stack, layers: below }, lambda, options);
            for (const pol of POLS) {
                check('substrate-side R', curve[REVERSE[pol]][0] - scan[pol].Rback[i], TOL.fraction,
                    `${backend} ${c.case} 0° ${pol} ${lambda} nm, monitorCurve Rr`);
            }
        });
    }
}

// ── Run ──────────────────────────────────────────────────────────────────────

const backends = [['js', tmmcore]];
if (existsSync(WASM)) {
    backends.push(['wasm', await tmmcore.instantiateTmmWasm(readFileSync(WASM))]);
} else {
    console.log('NOTE : src/wasm/tmm_kernel.wasm not built; the WebAssembly checks were skipped.');
}

for (const c of CASES) {
    for (const [backend, api] of backends) {
        for (const scan of c.scans) compareScan(api, backend, c, scan);
        compareLoci(api, backend, c);
        compareReverseAtNormalIncidence(api, backend, c);
    }
}
CASES.filter(c => c.incident.n[1] !== 0).forEach(compareAbsorbingIncidentDelays);

console.log(`${[...new Set(CASES.map(c => c.program))].join(', ')}: ${CASES.length} designs, ${comparisons} comparisons`);
for (const [group, { size, at }] of Object.entries(worst)) {
    console.log(`  worst ${group.padEnd(18)} ${size.toExponential(1)}  (${at})`);
}
console.log('PASS : Essential Macleod');
