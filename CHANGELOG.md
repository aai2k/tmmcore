# Changelog

## [0.5.0] - unreleased

### Added
- A Python package, `pip install tmmcore`, on the same C kernel: the R/T, derivative, needle, phase and jet calls under snake_case names, taking numpy arrays.
- One call shape for the R/T, derivative, needle and phase calls: `name(stack, wavelength or wavelengths, { angle, pol })`, with the stack given once as `{ incident, substrate, layers: [{ n, k, d }] }`. The needle scan takes its candidate materials and split fractions as the options `candidates` and `fractions`.
- Grid calls in JavaScript for R/T, phase and the phase Jacobian, which only the WebAssembly instance had.
- Grid versions of the Jacobian and the Hessian, in JavaScript, on the WebAssembly instance and in C.
- `needleScanSpectrum`, the needle scan over a grid, in JavaScript and on the WebAssembly instance.
- `pol: 'avg'`, the mean of s and p, the default in the R/T, derivative and needle calls.
- Grid calls take one angle or one per wavelength.
- The phase calls take a plain number as a constant index; a jet is needed only for dispersion.
- The phase calls may leave out `pol` at normal incidence, where s and p are the same.
- The jet arithmetic on the WebAssembly instance and in C, each C call over a whole grid of jets.
- C exports `tmm_needle_scan_spectrum` and `tmm_omega_from_lambda_nm`.
- `tmmcore/wasm`, the WebAssembly loader as its own entry point.

### Changed
- The calls have one set of names, the same in JavaScript and on the WebAssembly instance. The 0.4 names keep their arguments and results through 0.5.x, warn once on first use, and are removed in 0.6.0. [The names of 0.4](docs/api.md#the-names-of-04) maps each to its 0.5 name.
- `tmmCoefficientJets`, `tmmCoefficientThicknessJets`, `coefficientPhaseDispersion` and `coefficientPhaseThicknessDerivatives` warn once and are removed in 0.6.0, with no replacement.
- The 0.5 phase calls take the incident sine as `sine`, which the 0.4 phase names take as `sinTheta0Jet` and `sinJets`. They take no `omega`; the 0.4 phase names still do.
- A layer missing its index or its thickness throws a `TypeError` naming the layer. 0.4 skipped a layer without a thickness as if it were absent, and failed with an unrelated error on one without an index.
- C: `tmm_one`, `tmm_spectrum`, `tmm_phase_one` and `tmm_growing_eval_*` are renamed `tmm_rta`, `tmm_rta_spectrum`, `tmm_phase` and `tmm_growing_stack_*`, and the 0.4 names are gone.
- C: every call takes the angle after the wavelength, and the grid calls take it as an array of one or one per wavelength; the phase calls take ω after the stack.
- The sources are in one folder per language, `src/c`, `src/js`, `src/wasm` and `python`, each split into `rta`, `growing`, `derivatives`, `jet` and `phase`; the 0.4 import paths `tmmcore/tmm.js`, `tmmcore/phase.js`, `tmmcore/taylorJet.js`, `tmmcore/tmmWasm.js`, `tmmcore/tmm_kernel.wasm` and `tmmcore/tmm_kernel.c` still resolve.

### Fixed
- `tmmNeedleScan` in JavaScript returned NaN in place of each split result when the split fractions were a `Float64Array`.
- `depositionSpectra` on the WebAssembly instance returned one stray value per output for an empty stack.

## [0.4.2] - 2026-09-24

### Added
- `layerLogScale`, the log of the factor `layerMatrix` leaves out of an opaque layer, for callers who build their own matrix product and need T from it.
- `transmissionLogScale` in the results of `tmmCoefficientJets` and `tmmCoefficientThicknessJets`, and an optional `logScale` argument to `coefficientPhaseDispersion`, so the true |t| can be read from the coefficient jets.

### Fixed
- T through an opaque layer, such as a thick metal or a thick layer beyond the critical angle, stopped falling as the layer got thicker and came out far too large; every call that reports T or |t|² now gives the true value.
- The thickness derivatives of T from `tmmThicknessJacobian` and `tmmThicknessHessian` (`tmmJacobian` and `tmmHessian` in WebAssembly), and `t.dLogMagnitudeSquared` from `tmmPhaseThicknessJacobian`, `tmmPhaseJacobian` and `tmmPhaseJacobianSpectrum`, ignored an opaque layer's absorption; they now give the rate at which T really falls as that layer gets thicker.
- `tmmNeedleScan` gave a wrong `dT` for a needle inside an opaque layer; it is now right.

## [0.4.1] - 2026-09-24

### Added
- `snellCosTheta` and `incidentCosTheta` take an optional `cosTheta0`, which keeps the cosine in a medium close to the incident index accurate at grazing incidence.

### Changed
- `tmmThicknessJacobian` and `tmmThicknessHessian` give zero derivatives for a layer of negative thickness, as for an absent layer, where they used to treat it as a layer of zero thickness.
- `tmmPhaseThicknessJacobian`, `tmmPhaseJacobian`, `tmmPhaseJacobianSpectrum` and `tmmCoefficientThicknessJets` return thickness derivatives for opaque stacks, where they used to return `null` (NaN in the batched call).
- `tmmCoefficientThicknessJets` no longer reads a `rescaleThreshold` option.

### Fixed
- With an absorbing incident medium, A was 1 - R - T and came out nonzero even over lossless layers; `tmm`, `tmmSpectrum` and the derivative functions now return what the layers absorb, so R + T + A no longer sums to one in that case.
- Near 90° incidence results lost accuracy, and at exactly 90° `tmm` returned NaN; every call now holds up to and including 90°, where R is one.
- At the exact critical angle of a layer, or of the substrate in p, `tmm`, the derivative functions and the monitoring calls returned NaN; they now return the finite limit, and a substrate at its critical angle reflects everything.
- `tmmThicknessJacobian` and `tmmThicknessHessian` (`tmmJacobian` and `tmmHessian` in WebAssembly) returned wrong derivatives, sometimes NaN, for the thickness of an opaque layer such as a thick metal; they are now the derivatives of what `tmm` returns.
- `tmmThicknessJacobian`, `tmmThicknessHessian` and `tmmNeedleScan` returned all-zero derivatives at wavelengths where the stack as a whole was strongly opaque; they now return the true values.
- `tmmNeedleScan` gave wrong `intra` gradients for needles inside an opaque layer, too large, or nonzero where the true value is zero; they are now right.
- A layer of NaN thickness made `tmm` and the derivative functions return NaN; it is now left out, like a layer of negative thickness.
- When the WebAssembly heap ran out, some instance methods returned numbers computed from memory they never got; they now return NaN for that call, and the phase methods return `null`.

## [0.4.0] - 2026-09-13

### Added
- `dLogMagnitudeSquared` in the results of `tmmPhaseThicknessJacobian` and `tmmPhaseJacobian`: the derivative of ln|r|² and ln|t|² with respect to each layer thickness. `coefficientPhaseThicknessDerivatives` returns it as `logMagnitudeSquared`.
- `tmmPhaseJacobianSpectrum` on the WebAssembly instance (`tmm_phase_jacobian_spectrum` in C): the phase Jacobian at every wavelength of a grid in one call, with `hasPhaseJacobianSpectrum()` to check the loaded build has it.

### Changed
- In C, `tmm_phase_jacobian` writes 10 × N derivatives instead of 8 × N, five per coefficient, so callers must allocate the larger buffer, and the transmission block now starts at 5 × N.

## [0.3.1] - 2026-09-12

### Added
- `incidentCosTheta`, the cosine of the angle in the incident medium itself.
- `jetRealPart`, the real part of a complex jet.

### Changed
- `snellCosTheta` takes the real part of the incident index, so its result differs from 0.3.0 when that index is complex.

### Fixed
- With an absorbing incident medium at oblique incidence, R, T, A, their derivatives, the phase functions and the monitoring calls were wrong; they are now right.

## [0.3.0] - 2026-09-01

### Added
- `monitorCurve` on the WebAssembly instance (`tmm_monitor_curve` in C): R and T in s and p of one layer growing on a finished stack, at one wavelength, for a list of its thicknesses, with the reflectance of the coated surface seen from the substrate side.
- `depositionSpectra` (`tmm_deposition_spectra` in C): the spectrum after every layer of a deposition run in one call, layers in deposition order, with the same outputs as `monitorCurve`.
- `growingEval` (`tmm_growing_eval_create`, `tmm_growing_eval_set_top`, `tmm_growing_eval_sample` and `tmm_growing_eval_free` in C): a handle that holds a finished stack over a wavelength grid and, through `setTop`, `sample` and `free`, returns the spectrum for any thickness of the layer growing on it.
- `hasGrowingKernels()` and `hasGrowingEval()`, to check the loaded build has these.

## [0.2.0] - 2026-08-16

### Added
- `tmmPhaseDispersion`: phase, group delay, GDD and TOD of r and t at one wavelength, with refractive indices given as Taylor jets in angular frequency; options `omega` to set the time unit and `sinTheta0Jet` for a stack in a dispersive incident medium.
- `tmmPhaseThicknessJacobian`: the same quantities with their exact derivatives with respect to every layer thickness.
- `tmmCoefficientJets`, `tmmCoefficientThicknessJets`, `coefficientPhaseDispersion` and `coefficientPhaseThicknessDerivatives`, for callers who want the r and t jets themselves.
- Taylor-jet arithmetic for building index jets from any dispersion formula: `jetConstant`, `jetFromDerivatives`, `jetDerivatives`, `jetAdd`, `jetSubtract`, `jetMultiply`, `jetDivide`, `jetReciprocal`, `jetScale`, `jetSqrt`, `jetExp`, `jetLog`, `jetPower`, `jetSin`, `jetCos`, `jetSinCos`, `jetCompose`, `wavelengthOmegaJet`, `jetClampImaginary`, `jetClampRealMinimum`, `jetWithImaginaryPart`, `jetIsFinite` and `JET_ORDER`.
- `C_NM_PER_FS` and `omegaFromLambdaNm`, for results in fs, fs² and fs³.
- `tmmPhaseOne`, `tmmPhaseJacobian` and `tmmPhaseSpectrum` on the WebAssembly instance (`tmm_phase_one`, `tmm_phase_jacobian` and `tmm_phase_spectrum` in C), with `hasPhase()` to check the loaded build has them.

## [0.1.0] - 2026-08-01

### Added
- `tmm`: reflectance, transmittance and absorptance of a multilayer at one wavelength, angle and polarization, with complex indices for the layers, the incident medium and the substrate.
- `tmmThicknessJacobian`: exact first derivatives of R, T and A with respect to every layer thickness, returned with R, T and A.
- `tmmThicknessHessian`: exact second derivatives of R, T and A with respect to every pair of layer thicknesses.
- `tmmNeedleScan`: the needle P-function, the gradient of R, T and A for inserting a thin layer of each candidate material at every interface and, if asked, at chosen fractions inside each layer.
- Complex-number and 2×2 matrix primitives for building variants on the same conventions: `cadd`, `csub`, `cmul`, `cdiv`, `cabs2`, `cconj`, `csqrt`, `ccos`, `csin`, `creal`, `cimag`, `matmul`, `cmatvec`, `snellCosTheta`, `layerMatrix` and `rescaleMatrix`.
- The same calculations as C99 source needing only libm: `tmm_one`, `tmm_spectrum`, `tmm_jacobian`, `tmm_hessian` and `tmm_needle_scan`.
- A WebAssembly build of the C, loaded with `instantiateTmmWasm(bytes)` or `initTmmWasmFromUrl(url)`, whose instance has `tmmOne`, `tmmJacobian`, `tmmHessian`, `tmmNeedleScan` and `hasHessian()`.
- `tmmSpectrum` on the WebAssembly instance: a whole wavelength grid in both polarizations in one call, with an index per layer per wavelength.
- Loader state and worker helpers: `getTmmWasm`, `setTmmWasmEnabled`, `isTmmWasmEnabled`, `isTmmWasmReady`, `tmmWasmActive`, `initTmmWasmMainThread`, `getTmmWasmBytesForWorker`, `noteTmmWasmBytes` and `awaitTmmWasmReady`.
