/**
 * The C kernel built to WebAssembly. Load it with instantiateTmmWasm() or
 * initTmmWasmFromUrl(), then call the package's calls on the instance under
 * the same names, with the same arguments and results: rta, rtaSpectrum,
 * jacobian, jacobianSpectrum, hessian, hessianSpectrum, needleScan,
 * needleScanSpectrum, phase, phaseSpectrum, phaseJacobian,
 * phaseJacobianSpectrum (src/js/calls.js), and the jet arithmetic of
 * src/js/jet.js with omegaFromLambdaNm. The methods of 0.4 keep their
 * arguments through 0.5.x.
 *
 * The bindings are split by area as the C sources are: rta.js, growing.js,
 * derivatives.js, phase.js and jet.js, on the kernel memory of memory.js.
 * loader.js loads the module and wraps an instance.
 */

export {
    TmmWasmInstance, instantiateTmmWasm, initTmmWasmFromUrl, setTmmWasmEnabled,
    __setTmmWasmInstanceForTest, initTmmWasmMainThread, getTmmWasmBytesForWorker,
    noteTmmWasmBytes, awaitTmmWasmReady, isTmmWasmEnabled, isTmmWasmReady,
    tmmWasmActive, getTmmWasm,
} from './loader.js';
