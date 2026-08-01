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

Thicknesses are in nanometres. The sign and phase conventions are in
[Getting started](getting-started.md#conventions).

---

## Spectral quantities

### `tmm(lambda_nm, theta_deg, pol, n0, ns, layers)`

Reflectance, transmittance and absorptance.

```js
const { R, T, A } = tmm(550, 0, 's', [1, 0], [1.52, 0], layers);
```

**Returns** `{ R, T, A }`, each a number in $[0, 1]$.

Absorptance is computed as $A = 1 - R - T$, so the three always sum to one by
construction. That makes `R + T + A` useless as a self-check; use agreement with
a closed form or another implementation instead.

Layers of zero or negative thickness are skipped.

!!! note "Very thick absorbing layers"

    The phase thickness of a strongly absorbing layer grows without bound, and
    $\cosh$ of it overflows to infinity, which would poison the whole matrix
    product with `NaN`. The imaginary phase is clamped at a magnitude where the
    layer is already optically opaque, with single-pass transmittance below
    $10^{-43}$, so the result is exact to machine precision and the matrix stays
    finite. For non-absorbing or thin layers the clamp never engages.

---

## Analytic derivatives

Each derivative comes from the same characteristic-matrix product that produced
the spectrum, computed exactly rather than by finite differences or automatic
differentiation.

### `tmmThicknessJacobian(lambda_nm, theta_deg, pol, n0, ns, layers)`

First derivatives of R, T and A with respect to every layer thickness.

```js
const { R, dRdd } = tmmThicknessJacobian(550, 0, 's', [1, 0], [1.52, 0], layers);
// dRdd[j] = ∂R/∂d_j, in units of 1/nm
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, N }` where the three derivative arrays
have length `N`, the layer count.

The spectrum comes back with the derivatives, so a Gauss–Newton or
Levenberg–Marquardt step costs one evaluation rather than `N + 1`. See the
[refinement example](examples.md#refinement-with-the-analytic-jacobian).

### `tmmThicknessHessian(lambda_nm, theta_deg, pol, n0, ns, layers)`

Second derivatives, for Newton-type methods and curvature analysis.

```js
const { d2Rdd } = tmmThicknessHessian(550, 0, 's', [1, 0], [1.52, 0], layers);
// d2Rdd[i][j] = ∂²R/∂d_i∂d_j
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N }`. The second
derivatives are `N × N` symmetric arrays; first derivatives and the spectrum
come along at no extra cost.

Cost is $O(N^2)$ small-matrix operations per evaluation.

### `tmmNeedleScan(lambda_nm, theta_deg, pol, n0, ns, layers, candidateNs, intraFracs?)`

The needle P-function: the gradient of each spectral quantity with respect to
inserting an infinitesimally thin layer of a candidate material, evaluated at
every position in the stack at once.

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

`gaps` has `N + 1` entries, one per interface: index 0 is before the first
layer, index `N` is against the substrate. Each holds one `{ dR, dT, dA }` per
candidate. Insert where the merit-function gradient is negative.

This is the $d \to 0$ limit of Sullivan's numerical pre/post method, the
analytic P-function of Tikhonravov et al. Unlike the numerical form it needs no
trial thickness and no second spectrum evaluation.

!!! tip "A free consistency check"

    A needle of the host layer's own index, inserted at that layer's boundary,
    is just a thicker layer. So `gaps[0][c].dR` for a candidate matching layer
    0's index equals `dRdd[0]` from the Jacobian. The two code paths are derived
    independently and agree to round-off.

Unlike the other functions, this one does **not** skip zero-thickness layers, so
gap indices line up with your design array exactly.

---

## WebAssembly

Load the bytes, instantiate, then call methods on the returned instance. See
[Getting started](getting-started.md#turning-on-webassembly) for setup.

Polarization is an integer here: `0` for s, `1` for p.

| Method | Mirrors |
|---|---|
| `tmmOne(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmm` |
| `tmmJacobian(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmmThicknessJacobian` |
| `tmmHessian(lambda_nm, theta_deg, polCode, n0, ns, layers)` | `tmmThicknessHessian` |
| `tmmNeedleScan(lambda_nm, theta_deg, polCode, n0, ns, layers, candidateNs, intraFracs?)` | `tmmNeedleScan` |
| `hasHessian()` | Returns whether this build carries the Hessian kernel |

### `tmmSpectrum(lambdas, n0List, nsList, layerNK, thick, theta_deg)`

Evaluates an entire wavelength grid for **both** polarizations in a single call.
That is where the speed advantage comes from: it amortizes the
JavaScript-to-WebAssembly boundary crossing, which otherwise dominates.

| Argument | Shape |
|---|---|
| `lambdas` | `number[]` |
| `n0List`, `nsList` | `[re, im][]`, one per wavelength |
| `layerNK` | `[re, im][][]`, indexed `[layer][wavelength]` |
| `thick` | `number[]`, length `N` |

**Returns** `{ Rs, Ts, As, Rp, Tp, Ap }`, each a `Float64Array` of grid length.

### Loader functions

| | |
|---|---|
| `instantiateTmmWasm(bytes)` | Instantiate from an `ArrayBuffer`/`Uint8Array`. Async, returns the instance. |
| `initTmmWasmFromUrl(url)` | Fetch and instantiate. Resolves `true`/`false`; **does not throw** on failure. |
| `getTmmWasm()` | The current instance, or `null`. |
| `setTmmWasmEnabled(on)` | Feature flag, default off. |
| `isTmmWasmReady()`, `isTmmWasmEnabled()`, `tmmWasmActive()` | State queries. |

Instances are not shared between threads, since there is no shared memory. Each
worker instantiates its own from the same bytes.

---

## Low-level primitives

Also exported, for building transfer-matrix variants on the same conventions.
These are lower level than the functions above and correspondingly less stable
across versions.

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

`rescaleMatrix` is what keeps opaque stacks finite. A common real factor cancels
from reflectance but not from transmittance, so callers must carry the returned
log scale and apply $e^{-2\,\text{logScale}}$ to T.
