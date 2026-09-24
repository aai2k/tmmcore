/**
 * The WebAssembly kernel with its heap exhausted.
 *
 * Every exported kernel that allocates working state must notice when it gets
 * none. In wasm32 a null pointer is address 0 of linear memory, so a kernel
 * that works through one does not trap: it writes over whatever lives there
 * and hands back numbers. Here a fresh instance has its heap taken up, block by
 * halving block, until not one more byte can be had; every kernel must then
 * fill its outputs with NaN, and once the heap is returned the same instance
 * must compute correctly again.
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/wasm_allocation.mjs
 */
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { tmm } from '../src/tmm.js';
import { omegaFromLambdaNm } from '../src/phase.js';
import { instantiateTmmWasm } from '../src/tmmWasm.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'tmm_kernel.wasm');
if (!existsSync(WASM)) {
    console.log('SKIP : src/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}

const kernel = await instantiateTmmWasm(readFileSync(WASM));
const N = 3, nLam = 2, nCand = 2, nFrac = 1;
const INDICES = [2.35, 1.46, 2.35], THICK = [110, 95, 60];

// Inputs and outputs are placed before the heap is taken up.
function place(values) {
    const ptr = kernel.malloc(8 * values.length);
    assert.ok(ptr, 'room for the inputs');
    new Float64Array(kernel.memory.buffer, ptr, values.length).set(values);
    return ptr;
}
const jet = re => [re, 0, 0, 0, 0, 0, 0, 0];
const layers = place(INDICES.flatMap((n, k) => [n, 0, THICK[k]]));
const thick = place(THICK);
const lambdas = place([500, 600]);
const omegas = place([omegaFromLambdaNm(500), omegaFromLambdaNm(600)]);
const n0Pairs = place([1, 0, 1, 0]);
const nsPairs = place([1.52, 0, 1.52, 0]);
const matNK = place(INDICES.flatMap(n => [n, 0, n, 0]));
const candidates = place([2.35, 0, 1.46, 0]);
const fracs = place([0.5]);
const n0Jet = place(jet(1));
const nsJet = place(jet(1.52));
const layerJets = place(INDICES.flatMap(jet));
const n0Jets = place([...jet(1), ...jet(1)]);
const nsJets = place([...jet(1.52), ...jet(1.52)]);
const matJets = place(INDICES.flatMap(n => [...jet(n), ...jet(n)]));
const OMEGA = omegaFromLambdaNm(550);

// How many entries each kernel writes into each of its outputs, and the call.
const KERNELS = {
    tmm_spectrum: {
        sizes: Array(6).fill(nLam),
        call: out => kernel._tmm_spectrum(lambdas, nLam, n0Pairs, nsPairs, matNK, thick, N, 30, ...out),
    },
    tmm_jacobian: {
        sizes: [N, N, N, 3],
        call: out => kernel._tmm_jacobian(550, 30, 0, 1, 0, 1.52, 0, layers, N, ...out),
    },
    tmm_needle_scan: {
        sizes: [3, (N + 1) * nCand * 3, N * nFrac * nCand * 3],
        call: out => kernel._tmm_needle_scan(550, 30, 0, 1, 0, 1.52, 0, layers, N,
            candidates, nCand, fracs, nFrac, ...out),
    },
    tmm_hessian: {
        sizes: [N, N, N, N * N, N * N, N * N, 3],
        call: out => kernel._tmm_hessian(550, 30, 0, 1, 0, 1.52, 0, layers, N, ...out),
    },
    tmm_phase_one: {
        sizes: [10],
        call: out => kernel._tmm_phase_one(550, OMEGA, 30, 0, n0Jet, nsJet, layerJets, thick, N, 0, ...out),
    },
    tmm_phase_spectrum: {
        sizes: [10 * nLam],
        call: out => kernel._tmm_phase_spectrum(lambdas, omegas, nLam, n0Jets, nsJets, matJets, thick, N,
            30, 0, 0, ...out),
    },
    tmm_phase_jacobian: {
        sizes: [10, 10 * N],
        call: out => kernel._tmm_phase_jacobian(550, OMEGA, 30, 0, n0Jet, nsJet, layerJets, thick, N, 0,
            ...out),
    },
    tmm_phase_jacobian_spectrum: {
        sizes: [10 * nLam, 10 * N * nLam],
        call: out => kernel._tmm_phase_jacobian_spectrum(lambdas, omegas, nLam, n0Jets, nsJets, matJets,
            thick, N, 30, 0, 0, ...out),
    },
};
// Outputs start as a sentinel that is neither NaN nor a plausible result.
for (const entry of Object.values(KERNELS)) entry.out = entry.sizes.map(size => place(Array(size).fill(7)));
const read = (ptr, size) => Array.from(new Float64Array(kernel.memory.buffer, ptr, size));
const outputs = entry => entry.out.map((ptr, i) => read(ptr, entry.sizes[i]));

const hogs = [];
for (let size = 2 ** 30; size >= 1;) {
    const ptr = kernel.malloc(size);
    if (ptr) hogs.push(ptr);
    else size = Math.floor(size / 2);
}

let checked = 0;
for (const [name, entry] of Object.entries(KERNELS)) {
    entry.call(entry.out);
    const values = outputs(entry).flat();
    checked += values.length;
    assert.ok(values.every(Number.isNaN), `${name} without memory wrote ${values.filter(v => !Number.isNaN(v))}`);
}
for (const ptr of hogs) kernel.free(ptr);

// With the heap back, the same instance computes as before.
for (const [name, entry] of Object.entries(KERNELS)) {
    entry.call(entry.out);
    assert.ok(outputs(entry).flat().every(Number.isFinite), `${name} with memory back`);
}
const [R] = read(KERNELS.tmm_jacobian.out[3], 3);
const reference = tmm(550, 30, 's', [1, 0], [1.52, 0], INDICES.map((n, k) => ({ n: [n, 0], d: THICK[k] })));
assert.ok(Math.abs(R - reference.R) <= 1e-12, `R after the heap is returned: ${R} against ${reference.R}`);

console.log(`heap exhausted at ${kernel.memory.buffer.byteLength / 2 ** 20} MiB of linear memory`);
console.log(`${checked} outputs of ${Object.keys(KERNELS).length} kernels checked NaN without memory`);
console.log('PASS : allocation failure');
