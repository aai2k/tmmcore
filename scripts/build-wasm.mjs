/**
 * Build src/wasm/tmm_kernel.wasm with Emscripten, then check that the module
 * uses no WebAssembly feature Chromium 108 lacks.
 *
 * Chromium 108 is the engine of Electron 22, the last Electron that runs on
 * Windows 7. A newer Emscripten can start emitting newer instructions without
 * any change to the C, so the build validates the output with wasm-opt under
 * the features listed below, all of which Chromium 108 shipped, and fails if it
 * does not validate. The module is built under build/ and replaces
 * src/wasm/tmm_kernel.wasm only once it passes.
 *
 * The exported functions are the `export_name` attributes in the sources under
 * src/c; only the allocator is listed here. The sources are compiled together
 * with link-time optimization, so a helper inlines across files.
 *
 *   npm run build:wasm
 *
 * Needs emcc on PATH (activate emsdk first, or use src/build.ps1 on Windows,
 * which finds and activates it).
 */
import { spawnSync } from 'node:child_process';
import { existsSync, mkdirSync, renameSync, rmSync, statSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const ROOT = join(dirname(fileURLToPath(import.meta.url)), '..');
const BUILT = 'build/tmm_kernel.wasm';
const OUTPUT = 'src/wasm/tmm_kernel.wasm';
const WINDOWS = process.platform === 'win32';

const SOURCES = ['rta.c', 'growing.c', 'derivatives.c', 'phase.c', 'jet.c'].map(file => `src/c/${file}`);

const EMCC_ARGS = [
    ...SOURCES, '-O3', '-flto', '--no-entry',
    '-sSTANDALONE_WASM=1', '-sALLOW_MEMORY_GROWTH=1',
    '-sMIN_CHROME_VERSION=108',
    '-sEXPORTED_FUNCTIONS=_malloc,_free',
    '-o', BUILT,
];

// WebAssembly features past the first release that Chromium 108 ships: sign
// extension and mutable globals (Chrome 74), non-trapping float-to-int and
// bulk memory (75), multi-value (85), SIMD (91), reference types (96). Threads
// (74) and exception handling (95) are left out, since this C needs neither:
// the check is stricter than Chromium 108, never looser. Exception handling
// stays out on purpose, as wasm-opt's flag for it also admits the newer
// exnref form, which Chromium 108 cannot load.
const CHROMIUM_108_FEATURES = [
    '--mvp-features',
    '--enable-sign-ext', '--enable-mutable-globals',
    '--enable-nontrapping-float-to-int', '--enable-bulk-memory',
    '--enable-multivalue', '--enable-simd', '--enable-reference-types',
];

// emcc and em-config are batch files on Windows, which only a shell can run.
// The arguments hold no spaces or shell metacharacters, so one command line
// passes them through unchanged.
function runTool(name, args, options = {}) {
    return spawnSync([name, ...args].join(' '), {
        cwd: ROOT, shell: true, encoding: 'utf8', ...options,
    });
}

// A module that failed is removed, so it can neither be loaded nor committed.
function fail(message) {
    rmSync(join(ROOT, BUILT), { force: true });
    console.error(`build-wasm: ${message}`);
    process.exit(1);
}

mkdirSync(join(ROOT, dirname(BUILT)), { recursive: true });
const built = runTool('emcc', EMCC_ARGS, { stdio: 'inherit' });
if (built.error || built.status !== 0) {
    fail(`emcc failed (exit ${built.status}). Is emsdk activated? A prebuilt ${OUTPUT} is ` +
        'committed, so emcc is needed only after a change to the C under src/c.');
}

const config = runTool('em-config', ['BINARYEN_ROOT']);
const binaryenRoot = (config.stdout || '').trim();
const wasmOpt = join(binaryenRoot, 'bin', WINDOWS ? 'wasm-opt.exe' : 'wasm-opt');
if (config.status !== 0 || !binaryenRoot || !existsSync(wasmOpt)) {
    fail(`wasm-opt not found under BINARYEN_ROOT "${binaryenRoot}"; cannot check the module, ` +
        `so ${OUTPUT} is left as it was.`);
}

// No passes and no output file: wasm-opt only parses and validates.
const check = spawnSync(wasmOpt, [...CHROMIUM_108_FEATURES, join(ROOT, BUILT)], { encoding: 'utf8' });
if (check.error || check.status !== 0) {
    console.error(check.stderr || check.error?.message);
    fail(`the module uses a WebAssembly feature Chromium 108 lacks, which the Windows 7 build ` +
        `could not load; ${OUTPUT} is left as it was.`);
}

renameSync(join(ROOT, BUILT), join(ROOT, OUTPUT));
console.log(`Built ${OUTPUT} (${statSync(join(ROOT, OUTPUT)).size} bytes), valid for Chromium 108.`);
