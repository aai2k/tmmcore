/**
 * A single quarter-wave layer, checked against the closed-form result.
 *
 * For one layer of index n1 at quarter-wave optical thickness on a substrate ns,
 * at normal incidence, reflectance has an exact solution:
 *
 *     R = ((n0·ns − n1²) / (n0·ns + n1²))²
 *
 * This is the standard single-layer antireflection result — Macleod, Thin-Film
 * Optical Filters 5th ed., §3.2. It is the simplest case where a TMM
 * implementation can be checked against something other than another TMM.
 *
 *   node examples/01-single-layer.mjs
 */

import { tmm } from '../src/index.js';

const n0 = 1.0;      // air
const ns = 1.52;     // glass
const n1 = 1.38;     // MgF2
const lambda0 = 550; // nm

const d = lambda0 / (4 * n1);   // quarter wave at lambda0

const { R, T, A } = tmm(lambda0, 0, 's', [n0, 0], [ns, 0], [{ n: [n1, 0], d }]);

const num = n0 * ns - n1 * n1;
const den = n0 * ns + n1 * n1;
const exact = (num / den) ** 2;

console.log(`quarter-wave thickness   ${d.toFixed(4)} nm`);
console.log(`R from tmmcore           ${R.toPrecision(17)}`);
console.log(`R closed form            ${exact.toPrecision(17)}`);
console.log(`difference               ${Math.abs(R - exact).toExponential(2)}`);
console.log(`R + T + A                ${(R + T + A).toPrecision(17)}`);

// Without the coating, the bare substrate reflects the Fresnel amount.
const bare = tmm(lambda0, 0, 's', [n0, 0], [ns, 0], []);
console.log(`\nbare glass R             ${(bare.R * 100).toFixed(2)} %`);
console.log(`coated R                 ${(R * 100).toFixed(2)} %`);
