/**
 * One set of calls in C, in the JavaScript module and on a WebAssembly
 * instance. Every kernel export of the C kernel (src/c) names a call: the
 * export without its tmm_ prefix, in camelCase (tmm_rta_spectrum is
 * rtaSpectrum).
 *
 *   • each export's call is on the module and on an instance, with the same
 *     parameter names, except the exports in NOT_CALLS; calls.js builds the
 *     stack calls of both, so the names compared are those of the jet calls
 *   • each call an instance holds as its own has its export; the methods it
 *     inherits, the names of 0.4 among them, are in old_names.mjs
 *   • an export in NOT_CALLS that the module starts to offer must leave the
 *     list
 *
 * Skips cleanly if tmm_kernel.wasm has not been built.
 *
 *   node tests/parity.mjs
 */

import assert from 'node:assert/strict';
import { readFileSync, readdirSync, existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import * as tmmcore from '../src/js/index.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const WASM = join(HERE, '..', 'src', 'wasm', 'tmm_kernel.wasm');
if (!existsSync(WASM)) {
    console.log('SKIP : src/wasm/tmm_kernel.wasm not built. Build it with: npm run build:wasm');
    process.exit(0);
}
const wasm = await tmmcore.instantiateTmmWasm(readFileSync(WASM));

// Exports the module does not offer under their own names: the growing-layer
// kernels, reached through an instance's monitorCurve, depositionSpectra and
// growingEval.
const NOT_CALLS = [
    'monitorCurve', 'depositionSpectra',
    'growingStackCreate', 'growingStackSetTop', 'growingStackSample', 'growingStackFree',
];

const C = join(HERE, '..', 'src', 'c');
const source = readdirSync(C).filter(file => file.endsWith('.c'))
    .map(file => readFileSync(join(C, file), 'utf8')).join('\n');
const exports = [...source.matchAll(/TMM_EXPORT\("(tmm_[a-z0-9_]+)"\)/g)].map(match => match[1]);
assert.ok(exports.length > 0, 'kernel exports found');
const camel = name => name.replace(/^tmm_/, '').replace(/_([a-z0-9])/g, (_, c) => c.toUpperCase());

// The parameter names of a function, read from its source; the package's
// default values hold no parentheses or commas.
function parameters(fn) {
    const text = fn.toString();
    const list = text.match(/^([A-Za-z_$][\w$]*)\s*=>/)?.[1] ?? text.match(/\(([^)]*)\)/)[1];
    return list.split(',').map(part => part.split('=')[0].trim()).filter(Boolean);
}

let compared = 0;
for (const name of exports.map(camel)) {
    if (NOT_CALLS.includes(name)) {
        assert.equal(typeof tmmcore[name], 'undefined', `${name} is on the module now: take it off NOT_CALLS`);
        continue;
    }
    assert.equal(typeof tmmcore[name], 'function', `${name} on the module`);
    assert.equal(typeof wasm[name], 'function', `${name} on an instance`);
    assert.deepEqual(parameters(wasm[name]), parameters(tmmcore[name]), `${name}: parameters`);
    compared++;
}

const calls = new Set(exports.map(camel));
const ownCalls = Object.keys(wasm).filter(key =>
    typeof wasm[key] === 'function' && !key.startsWith('_') && key !== 'malloc' && key !== 'free');
for (const name of ownCalls) assert.ok(calls.has(name), `${name} on an instance has no kernel export`);

console.log(`${exports.length} kernel exports; ${compared} calls the same on module and instance, `
    + `${NOT_CALLS.length} not offered under their own names`);
console.log('PASS : parity');
