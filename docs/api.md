# API reference

All functions take the same leading arguments:

| Argument | Type | Meaning |
|---|---|---|
| `lambda_nm` | number | Wavelength in nanometres |
| `theta_deg` | number | Angle of incidence, degrees from normal |
| `pol` | `'s'` \| `'p'` | Polarization |
| `n0` | `[re, im]` | Incident medium index, $\tilde n = n + ik$ |
| `ns` | `[re, im]` | Substrate index |
| `layers` | `{ n: [re, im], d }[]` | Stack, incident medium first |

Thicknesses are in nanometres. The sign and phase conventions are in [Getting started](getting-started.md#conventions).

---

## Spectral quantities

### `tmm(lambda_nm, theta_deg, pol, n0, ns, layers)`

Reflectance, transmittance and absorptance.

```js
const { R, T, A } = tmm(550, 0, 's', [1, 0], [1.52, 0], layers);
```

**Returns** `{ R, T, A }`, each a number in $[0, 1]$.

Absorptance is computed as $A = 1 - R - T$, so the three always sum to one by construction. That makes `R + T + A` useless as a self-check; use agreement with a closed form or another implementation instead.

Layers of zero or negative thickness are skipped.

!!! note "Very thick absorbing layers"

    The phase thickness of a strongly absorbing layer grows without bound, and $\cosh$ of it overflows to infinity, which would poison the whole matrix product with `NaN`. The imaginary phase is clamped at a magnitude where the layer is already optically opaque, with single-pass transmittance below $10^{-43}$, so the result is exact to machine precision and the matrix stays finite. For non-absorbing or thin layers the clamp never engages.

---

## Analytic derivatives

Each derivative comes from the same characteristic-matrix product that produced the spectrum, computed exactly rather than by finite differences or automatic differentiation.

### `tmmThicknessJacobian(lambda_nm, theta_deg, pol, n0, ns, layers)`

First derivatives of R, T and A with respect to every layer thickness.

```js
const { R, dRdd } = tmmThicknessJacobian(550, 0, 's', [1, 0], [1.52, 0], layers);
// dRdd[j] = ∂R/∂d_j, in units of 1/nm
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, N }` where the three derivative arrays have length `N`, the layer count.

The spectrum comes back with the derivatives, so a Gauss–Newton or Levenberg–Marquardt step costs one evaluation rather than `N + 1`. See the [refinement example](examples.md#refinement-with-the-analytic-jacobian).

### `tmmThicknessHessian(lambda_nm, theta_deg, pol, n0, ns, layers)`

Second derivatives, for Newton-type methods and curvature analysis.

```js
const { d2Rdd } = tmmThicknessHessian(550, 0, 's', [1, 0], [1.52, 0], layers);
// d2Rdd[i][j] = ∂²R/∂d_i∂d_j
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N }`. The second derivatives are `N × N` symmetric arrays; first derivatives and the spectrum come along at no extra cost.

Cost is $O(N^2)$ small-matrix operations per evaluation.

### `tmmNeedleScan(lambda_nm, theta_deg, pol, n0, ns, layers, candidateNs, intraFracs?)`

The needle P-function: the gradient of each spectral quantity with respect to inserting an infinitesimally thin layer of a candidate material, evaluated at every position in the stack at once.

```js
const scan = tmmNeedleScan(550, 0, 's', [1, 0], [1.52, 0], layers,
                           [[2.35, 0], [1.46, 0]]);
// scan.gaps[position][candidate] = { dR, dT, dA }
```

**Arguments**

| | |
|---|---|
| `candidateNs` | `[re, im][]`, the materials that could be inserted |
| `intraFracs` | optional `number[]`, fractional positions for splitting a host layer |

**Returns** `{ R, T, A, gaps, intra, N }`.

`gaps` has `N + 1` entries, one per interface: index 0 is before the first layer, index `N` is against the substrate. Each holds one `{ dR, dT, dA }` per candidate. Insert where the merit-function gradient is negative.

This is the $d \to 0$ limit of Sullivan's numerical pre/post method, the analytic P-function of Tikhonravov et al. Unlike the numerical form it needs no trial thickness and no second spectrum evaluation.

!!! tip "A free consistency check"

    A needle of the host layer's own index, inserted at that layer's boundary,
    is just a thicker layer. So `gaps[0][c].dR` for a candidate matching layer
    0's index equals `dRdd[0]` from the Jacobian. The two code paths are derived
    independently and agree to round-off.

Unlike the other functions, this one does **not** skip zero-thickness layers, so gap indices line up with your design array exactly.

---

## Phase dispersion

Phase, group delay, group delay dispersion and third-order dispersion, for chirped mirrors and any other design where the shape of a pulse matters.

These take refractive indices as **Taylor jets** rather than plain complex numbers, because the answer depends on how the index varies with frequency and not only on its value. See [Refractive index as a jet](#refractive-index-as-a-jet) below.

### `tmmPhaseDispersion(lambda_nm, theta_deg, pol, n0Jet, nsJet, layers, options?)`

```js
const { r, t } = tmmPhaseDispersion(800, 0, 's', n0Jet, nsJet, layers);
// r.gdd = group delay dispersion on reflection, fs²
```

`layers` is `{ nJet, d }[]` here, `d` in nanometres.

**Returns** `{ r, t }`, one per coefficient, each `{ phaseRad, phaseDeg, gd, gdd, tod, magnitudeSquared }`, or `null` where that coefficient is exactly zero and the phase is undefined.

| `options` | |
|---|---|
| `omega` | Angular frequency. Sets the time unit of the result; defaults to $2\pi c/\lambda$ in rad/fs. |
| `sinTheta0Jet` | Sine of the incident angle as a jet, for a stack embedded in a dispersive medium at a fixed external angle. Defaults to the constant $\sin\theta_0$. |

!!! note "Sampling does not change these numbers"

    GD, GDD and TOD here are pointwise logarithmic derivatives of the coefficient, computed from analytic frequency derivatives of the characteristic matrix. No phase unwrapping and no wavelength stencil takes part, so a value at one wavelength does not depend on any neighbouring sample. Plot on whatever grid suits you; only the drawn polyline changes.

    That is the practical difference from differencing an unwrapped phase curve, where the answer moves with the step, and a step that is too coarse can lock onto a wrong value that looks stable as you refine it.

### `tmmPhaseThicknessJacobian(lambda_nm, theta_deg, pol, n0Jet, nsJet, layers, options?)`

The same quantities plus their exact derivatives with respect to every layer thickness, for gradient-based dispersion design.

```js
const { r } = tmmPhaseThicknessJacobian(800, 0, 's', n0Jet, nsJet, layers);
// r.dGdd[j] = ∂GDD/∂d_j, in fs²/nm
```

**Returns** `{ r, t }`, each the phase quantities plus `dPhaseDeg`, `dGd`, `dGdd` and `dTod`, arrays of length `layers.length`.

Frequency remains the Taylor variable throughout, so every thickness derivative is itself a third-order frequency jet and all four quantities come out of one matrix product. Fitting a GDD target costs one evaluation per step rather than `N + 1`, the same bargain `tmmThicknessJacobian` offers for reflectance.

Zero-thickness layers are retained so their derivatives and indices line up with your design array. Negative or non-finite thicknesses are skipped, matching the point evaluator, and their derivative entries are zero.

The derivative arrays are `null` if the matrix product overflowed, since the prefix/suffix decomposition cannot carry a rescaling; the phase quantities themselves are still returned.

### Refractive index as a jet

A jet is a truncated Taylor series, `[f, f', f''/2!, f'''/3!]`, each entry a complex `[re, im]` pair. The phase functions differentiate with respect to angular frequency, so an index jet means $[\tilde n,\ d\tilde n/d\omega,\ (d^2\tilde n/d\omega^2)/2,\ (d^3\tilde n/d\omega^3)/6]$.

tmmcore ships no material models on purpose: your dispersion data is yours. What it exports instead is the [jet arithmetic](#taylor-jet-arithmetic), so you compose the jet from whatever formula you already use. A Cauchy index, written straight through:

```js
import { wavelengthOmegaJet, jetReciprocal, jetMultiply,
         jetScale, jetAdd, jetConstant, omegaFromLambdaNm } from 'tmmcore';

// n(λ) = A + B/λ² + C/λ⁴, with a constant k
function cauchyJet(lambda_nm, A, B, C, k) {
    const omega = omegaFromLambdaNm(lambda_nm);
    const lambda = wavelengthOmegaJet(lambda_nm, omega);      // λ(ω)
    const inverseSquare = jetReciprocal(jetMultiply(lambda, lambda));
    return jetAdd(jetConstant(A, k),
        jetAdd(jetScale(inverseSquare, B),
               jetScale(jetMultiply(inverseSquare, inverseSquare), C)));
}
```

Nothing there differentiates anything by hand. `wavelengthOmegaJet` carries the $\lambda(\omega)$ derivatives and every operation after it propagates them. Sellmeier, Lorentz oscillators, and any composition of them work the same way.

For a non-dispersive medium, `jetConstant(n, k)` is the whole story.

If your data is a table rather than a formula, interpolate it with something at least three times differentiable and use `jetCompose`, which takes a value and its first three derivatives and chains them onto the $\omega$ jet. An interpolant that is only $C^1$, such as PCHIP, gives a GDD and TOD that jump at every knot. That is a real property of the interpolant rather than a bug, and it is worth showing your users rather than smoothing over.

### Units

The differentiation variable is angular frequency, and you choose its unit by choosing `omega`. GD comes back in the reciprocal of that unit, GDD in its square, TOD in its cube. The default uses `C_NM_PER_FS` with wavelengths in nm, giving **fs, fs² and fs³**.

| | |
|---|---|
| `C_NM_PER_FS` | Speed of light in vacuum, 299.792458 nm/fs |
| `omegaFromLambdaNm(lambda_nm)` | Angular frequency in rad/fs |

### Sign convention

tmmcore uses $\tilde n = n + ik$ with $e^{-i\omega t}$, the complex conjugate of Macleod's convention. Macleod defines each reported order as minus the derivative of physical phase, so all three quantities equal derivatives of the raw transfer-matrix phase computed here, and the numbers agree with his.

Positive GDD means the red end of the pulse arrives first.

---

## Taylor-jet arithmetic

Exported so you can build index jets from your own dispersion model, as above. Same status as the complex primitives: lower level, and correspondingly less stable across versions.

Jets are arrays of four complex `[re, im]` pairs.

| | |
|---|---|
| `jetConstant(re, im?)` | A constant: value with zero derivatives |
| `jetFromDerivatives(f, f', f'', f''')` | Build from plain derivatives, applying the factorials |
| `jetDerivatives(jet)` | Read back as `[f, f', f'', f''']` |
| `jetAdd`, `jetSubtract`, `jetMultiply`, `jetDivide`, `jetReciprocal` | Arithmetic |
| `jetScale(jet, s)` | Multiply every order by a real scalar |
| `jetSqrt`, `jetExp`, `jetLog`, `jetPower`, `jetSin`, `jetCos`, `jetSinCos` | Elementary functions |
| `jetCompose(value, [f', f'', f'''], xJet)` | Chain rule for a function known only through its derivatives |
| `wavelengthOmegaJet(lambda, omega)` | $\lambda(\omega)$, the usual starting point |
| `jetClampImaginary`, `jetClampRealMinimum`, `jetWithImaginaryPart`, `jetIsFinite` | Guards and helpers |
| `JET_ORDER` | 3 |

`wavelengthOmegaJet` needs no value for $c$: given both $\lambda$ and $\omega$, every derivative follows from $\lambda' = -\lambda/\omega$. Its output is in whatever unit $\lambda$ was given in.

---

## WebAssembly

Load the bytes, instantiate, then call methods on the returned instance. See [Getting started](getting-started.md#turning-on-webassembly) for setup.

Polarization is an integer here: `0` for s, `1` for p.

| Method | Mirrors |
|---|---|
| `tmmOne(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmm` |
| `tmmJacobian(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmmThicknessJacobian` |
| `tmmHessian(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmmThicknessHessian` |
| `tmmNeedleScan(lambda_nm, theta_deg, polCode, n0, ns, layers, candidateNs, intraFracs?)` | `tmmNeedleScan` |
| `tmmPhaseOne(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options?)` | `tmmPhaseDispersion` |
| `tmmPhaseJacobian(lambda_nm, theta_deg, polCode, n0Jet, nsJet, layers, options?)` | `tmmPhaseThicknessJacobian` |
| `hasHessian()`, `hasPhase()` | Whether this build carries those kernels |

### `tmmSpectrum(lambdas, n0List, nsList, layerNK, thick, theta_deg)`

Evaluates an entire wavelength grid for **both** polarizations in a single call. That is where the speed advantage comes from: it amortizes the JavaScript-to-WebAssembly boundary crossing, which otherwise dominates.

| Argument | Shape |
|---|---|
| `lambdas` | `number[]` |
| `n0List`, `nsList` | `[re, im][]`, one per wavelength |
| `layerNK` | `[re, im][][]`, indexed `[layer][wavelength]` |
| `thick` | `number[]`, length `N` |

**Returns** `{ Rs, Ts, As, Rp, Tp, Ap }`, each a `Float64Array` of grid length.

### `tmmPhaseSpectrum(lambdas, n0Jets, nsJets, layerJets, thick, theta_deg, polCode, options?)`

The same amortization for the phase kernel, and the call worth reaching for first: this kernel is roughly an order of magnitude dearer per sample than the plain spectrum, so a per-wavelength boundary crossing hurts correspondingly more.

| Argument | Shape |
|---|---|
| `lambdas` | `number[]` |
| `n0Jets`, `nsJets` | jet per wavelength |
| `layerJets` | indexed `[layer][wavelength]` |
| `thick` | `number[]`, length `N` |

| `options` | |
|---|---|
| `omegas` | `number[]`, one per wavelength. Defaults to $2\pi c/\lambda$ in rad/fs. |
| `sinJets` | jet per wavelength, as `sinTheta0Jet` above |

**Returns** `{ r, t }`, each `{ phaseRad, gd, gdd, tod, magnitudeSquared }` of `Float64Array` at grid length. Samples where the coefficient vanishes hold `NaN`.

Unlike `tmmSpectrum` this takes one polarization rather than doing both, because the cost per sample is high enough that a caller at normal incidence, where s and p are identical, should not pay twice for the same numbers. Call it twice and average if you want the unpolarized answer.

### Loader functions

| | |
|---|---|
| `instantiateTmmWasm(bytes)` | Instantiate from an `ArrayBuffer`/`Uint8Array`. Async, returns the instance. |
| `initTmmWasmFromUrl(url)` | Fetch and instantiate. Resolves `true`/`false`; **does not throw** on failure. |
| `getTmmWasm()` | The current instance, or `null`. |
| `setTmmWasmEnabled(on)` | Feature flag, default off. |
| `isTmmWasmReady()`, `isTmmWasmEnabled()`, `tmmWasmActive()` | State queries. |

Instances are not shared between threads, since there is no shared memory. Each worker instantiates its own from the same bytes.

---

## Low-level primitives

Also exported, for building transfer-matrix variants on the same conventions. These are lower level than the functions above and correspondingly less stable across versions.

Complex numbers are `[re, im]` pairs throughout.

| | |
|---|---|
| `cadd`, `csub`, `cmul`, `cdiv` | Complex arithmetic |
| `cabs2`, `cconj`, `csqrt`, `ccos`, `csin` | Modulus squared, conjugate, roots, trigonometry |
| `creal`, `cimag` | Parts |
| `matmul(A, B)` | 2×2 complex matrix product |
| `cmatvec(M, v)` | 2×2 matrix times 2-vector |
| `snellCosTheta(n0, sinTheta0, nj)` | Complex $\cos\theta$ in a medium |
| `layerMatrix(nj, dj_nm, lambda_nm, cosTheta_j, pol)` | Characteristic matrix of one layer |
| `rescaleMatrix(M)` | Rescales in place past an overflow threshold; returns the accumulated log scale |

`rescaleMatrix` is what keeps opaque stacks finite. A common real factor cancels from reflectance but not from transmittance, so callers must carry the returned log scale and apply $e^{-2\,\text{logScale}}$ to T.
