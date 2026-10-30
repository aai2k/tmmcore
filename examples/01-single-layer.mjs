/**
 * A single quarter-wave layer, checked against the closed-form result.
 *
 * For one layer of index n1 at quarter-wave optical thickness on a substrate ns,
 * at normal incidence, reflectance has an exact solution:
 *
 *     R = ((n0·ns − n1²) / (n0·ns + n1²))²
 *
 * This is the single-layer reflectance of Macleod, Thin-Film Optical Filters
 * 5th ed., Eq. 4.5, at a quarter wave.
 *
 *   node examples/01-single-layer.mjs
 */

import assert from 'node:assert/strict';
import { rta } from '../src/js/index.js';

const n0 = 1.0;      // air
const ns = 1.52;     // glass
const n1 = 1.38;     // MgF2
const lambda0 = 550; // nm

const d = lambda0 / (4 * n1);   // quarter wave at lambda0

const { R, T, A } = rta({ incident: n0, substrate: ns, layers: [{ n: n1, d }] }, lambda0);

const num = n0 * ns - n1 * n1;
const den = n0 * ns + n1 * n1;
const exact = (num / den) ** 2;

console.log(`quarter-wave thickness   ${d.toFixed(4)} nm`);
console.log(`R from tmmcore           ${R.toPrecision(17)}`);
console.log(`R closed form            ${exact.toPrecision(17)}`);
console.log(`difference               ${Math.abs(R - exact).toExponential(2)}`);
// The bound is the difference above, rounded up to a power of ten, times ten.
assert.ok(Math.abs(R - exact) <= 1e-15, 'R matches the closed form');
console.log(`R + T + A                ${(R + T + A).toPrecision(17)}`);

// Without the coating, the bare substrate reflects the Fresnel amount.
const bare = rta({ incident: n0, substrate: ns, layers: [] }, lambda0);
console.log(`\nbare glass R             ${(bare.R * 100).toFixed(2)} %`);
console.log(`coated R                 ${(R * 100).toFixed(2)} %`);
