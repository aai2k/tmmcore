/**
 * A narrowband filter tilted from 0° to 40°.
 *
 * A single-cavity (Fabry-Perot) filter of the high-index-cavity type in
 * Macleod's Figure 8.8 (Thin-Film Optical Filters 5th ed., §8.2.2):
 * quarter-wave reflectors (HL)⁴ and (LH)⁴ either side of a half-wave HH spacer,
 * first order. Tilting it moves the passband to shorter wavelengths, because the
 * phase thickness 2*pi*n*d*cos(theta)/lambda of every layer falls as the angle
 * inside the layer grows.
 *
 * For a cavity of index n* between reflectors whose phase does not change, the
 * peak moves to (§8.2.5.1 "Simple Tilts in Collimated Light")
 *
 *     lambda(theta0) = lambda0 cos(theta),  theta = arcsin(sin(theta0) / n*)
 *                    = lambda0 sqrt(1 - sin²(theta0) / n*²)
 *
 * In a real filter the reflector phase changes with angle too. Pidgeon and
 * Smith (J. Opt. Soc. Am. 54, 1459, 1964) showed the peak still moves by that
 * form, with an effective index n* between nL and nH. Macleod gives n* for a
 * high-index cavity of order m as Eq. 8.35 (§8.2.5.2), and for the
 * first order as Eq. 8.36:
 *
 *     n* = sqrt(nH nL)
 *
 * The derivation is a small-angle expansion about normal incidence that takes
 * both admittances as n, so it gives one n* for s and p. The table
 * prints how far each polarization strays from it.
 *
 * The exact check is at 0°. Every layer is then a quarter or half wave at
 * lambda0, and a half wave is an absentee layer (§3.1): taking out HH
 * leaves LL, taking that out leaves HH, and so on until no layer is left. At
 * lambda0 the filter transmits exactly what the bare substrate does,
 * T = 4 n0 ns / (n0 + ns)² (Eq. 2.53).
 *
 * Indices are non-dispersive, so the tilt is the only effect in view.
 *
 * Writes the plot to docs/img/tilted-filter.svg.
 *
 *   node examples/06-tilted-filter.mjs
 */

import { writeFileSync, mkdirSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { rta, rtaSpectrum } from '../src/js/index.js';
import { linePlot } from './_plot.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const OUT = join(HERE, '..', 'docs', 'img');

const nH = 2.35;               // high index, TiO2-like
const nL = 1.46;               // low index, SiO2-like
const n0 = 1;                  // air
const ns = 1.52;               // glass

const LAMBDA0 = 550;           // nm, reference wavelength of the quarter waves
const PAIRS = 4;               // (HL)^4 HH (LH)^4  ->  17 layers

const H = { n: nH, d: LAMBDA0 / (4 * nH) };
const L = { n: nL, d: LAMBDA0 / (4 * nL) };
const spacer = { n: nH, d: LAMBDA0 / (2 * nH) };

const layers = [];
for (let i = 0; i < PAIRS; i++) layers.push(H, L);
layers.push(spacer);
for (let i = 0; i < PAIRS; i++) layers.push(L, H);
const filter = { incident: n0, substrate: ns, layers };

const nStar = Math.sqrt(nH * nL);   // Macleod Eq. 8.36
const shifted = theta0 => {
    const s = Math.sin((theta0 * Math.PI) / 180) / nStar;
    return LAMBDA0 * Math.sqrt(1 - s * s);
};

// A scan wide enough to hold the 0° to 40° peaks, fine against the passband.
const STEP = 0.05;             // nm
const SCAN = Array.from({ length: 1201 }, (_, i) => 500 + i * STEP);

// The peak: the best scan sample, then a golden-section search on rta between
// its two neighbours, which holds only that one maximum.
function peak(angle, pol) {
    const { T } = rtaSpectrum(filter, SCAN, { angle, pol });
    let best = 0;
    for (let i = 1; i < T.length; i++) if (T[i] > T[best]) best = i;
    const at = lam => rta(filter, lam, { angle, pol }).T;
    const golden = (Math.sqrt(5) - 1) / 2;
    let a = SCAN[best] - STEP, b = SCAN[best] + STEP;
    let c = b - golden * (b - a), d = a + golden * (b - a);
    let fc = at(c), fd = at(d);
    while (b - a > 1e-6) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - golden * (b - a); fc = at(c); }
        else { a = c; c = d; fc = fd; d = a + golden * (b - a); fd = at(d); }
    }
    return (a + b) / 2;
}

// Where T falls to half its peak value, by bisection between the peak and a
// point well into the reflector's stop band.
function halfPoint(top, far) {
    const half = rta(filter, top).T / 2;
    let a = top, b = far;
    while (Math.abs(b - a) > 1e-6) {
        const m = (a + b) / 2;
        if (rta(filter, m).T > half) a = m; else b = m;
    }
    return (a + b) / 2;
}

const top = peak(0, 's');
const width = halfPoint(top, top + 20) - halfPoint(top, top - 20);

const exact = (4 * n0 * ns) / (n0 + ns) ** 2;
const atLambda0 = rta(filter, LAMBDA0).T;

console.log(`(HL)^${PAIRS} HH (LH)^${PAIRS} on glass, nH = ${nH}, nL = ${nL}, lambda0 = ${LAMBDA0} nm, ` +
    `halfwidth ${width.toFixed(2)} nm at 0°`);
console.log(`effective index n* = sqrt(nH nL) = ${nStar.toFixed(4)}\n`);

console.log('at 0° the layers cancel in half-wave pairs, and the filter transmits as bare glass');
console.log(`  T at ${LAMBDA0} nm, tmmcore     ${atLambda0.toPrecision(17)}`);
console.log(`  4 n0 ns / (n0 + ns)²     ${exact.toPrecision(17)}`);
console.log(`  difference               ${Math.abs(atLambda0 - exact).toExponential(1)}\n`);

// Adding 0 turns -0 into 0, so a difference that rounds to zero prints unsigned.
const nm = (v, width) => String((Math.round(v * 1000) / 1000 + 0).toFixed(3)).padStart(width);
console.log('peak wavelength, nm');
console.log('  AOI' + 's'.padStart(10) + 'p'.padStart(10) + 'formula'.padStart(10) +
    's - formula'.padStart(14) + 'p - formula'.padStart(14));
for (const theta of [0, 10, 20, 30, 40]) {
    const s = peak(theta, 's'), p = peak(theta, 'p'), f = shifted(theta);
    console.log(`  ${String(theta).padStart(2)}°` + nm(s, 10) + nm(p, 10) + nm(f, 10) + nm(s - f, 14) + nm(p - f, 14));
}

// The window runs past 560 nm so the legend sits over the stop band, not a peak.
const plotGrid = Array.from({ length: 751 }, (_, i) => 500 + i * 0.1);
const curve = angle => Array.from(rtaSpectrum(filter, plotGrid, { angle }).T, T => T * 100);

mkdirSync(OUT, { recursive: true });
writeFileSync(join(OUT, 'tilted-filter.svg'), linePlot({
    title: `Single-cavity filter (HL)⁴ HH (LH)⁴ tilted, unpolarized, λ₀ = ${LAMBDA0} nm`,
    xLabel: 'wavelength (nm)',
    yLabel: 'transmittance (%)',
    yMin: 0,
    yMax: 100,
    series: [
        { label: '0°', x: plotGrid, y: curve(0) },
        { label: '20°', x: plotGrid, y: curve(20) },
        { label: '40°', x: plotGrid, y: curve(40) },
    ],
}));
console.log('\nwrote docs/img/tilted-filter.svg');
