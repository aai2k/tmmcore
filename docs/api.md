# API reference

The R/T, derivative, needle, phase and monitoring calls take a stack, a wavelength in nm or, in a `Spectrum` call, `depositionSpectra` and `growingStack`, a list of them, and an options object:

```js
import { rta, rtaSpectrum } from 'tmmcore';

const stack = {
    incident: 1,
    substrate: 1.52,
    layers: [{ n: 2.35, d: 58.5 }, { n: 1.46, d: 94.2 }],   // d in nm
};
const wavelengths = [450, 500, 550, 600, 650];

rta(stack, 550).R;                                            // natural light, normal incidence
rtaSpectrum(stack, wavelengths, { angle: 45, pol: 's' }).R;   // one value per wavelength
```

The [WebAssembly instance](#webassembly) has the same calls, and the [Python package](python.md) has them under snake_case names. The sign and phase conventions are in [Getting started](getting-started.md#conventions).

| Call | Returns |
|---|---|
| [`rta`, `rtaSpectrum`](#r-t-and-a) | R, T and A |
| [`jacobian`, `jacobianSpectrum`](#thickness-derivatives) | R, T and A, and their derivatives in each layer's thickness |
| [`hessian`, `hessianSpectrum`](#thickness-derivatives) | the same, and the second derivatives |
| [`needleScan`, `needleScanSpectrum`](#needle-scan) | the derivatives of R, T and A for inserting a thin layer anywhere in the stack |
| [`phase`, `phaseSpectrum`](#phase-gd-gdd-and-tod) | phase, GD, GDD and TOD of r and t |
| [`phaseJacobian`, `phaseJacobianSpectrum`](#phase-gd-gdd-and-tod) | the same, and their thickness derivatives |
| [`jetConstant` … `jetIsFinite`](#jet-arithmetic) | Taylor-jet arithmetic, which builds the dispersive indices the phase calls take |
| [`monitorCurve`, `depositionSpectra`, `growingStack`](#monitoring-and-deposition) | R and T in s and p while the layers are deposited |

## Arguments

### The stack

`{ incident, substrate, layers }`, the layers ordered from the incident medium toward the substrate.

| | |
|---|---|
| `layers` | `{ n, k, d }[]`, with `d` the thickness in nm |
| `incident`, `substrate` | a number `n`, or `{ n, k }` |
| `n`, `k` | a number, or in a `Spectrum` call a list with one value per wavelength |

$\tilde n = n + ik$, with $k \ge 0$ for absorbing media. `k` may be left out and is then 0. A medium given as a bare list is refused, since an `[re, im]` pair would read as two wavelengths. A layer of zero, negative or NaN thickness adds nothing to the spectrum.

In the phase calls an index is a number, which is a constant index with no dispersion, or a [jet](#indices-as-jets) with `k` inside it: one jet at one wavelength, a list with one jet per wavelength in a `Spectrum` call. A list of plain numbers per wavelength is refused there, since the group delay computed from it would leave out the material's dispersion.

### Options

| | |
|---|---|
| `angle` | Angle of incidence in degrees from the normal, 0 when left out. A `Spectrum` call also takes a list with one angle per wavelength. |
| `sine` | Phase calls only, in place of `angle`: the sine of the angle of incidence as a [jet](#indices-as-jets), a list with one jet per wavelength in a `Spectrum` call. For [an angle fixed outside the incident medium](#an-angle-fixed-outside-the-incident-medium). |
| `pol` | In the R/T, derivative and needle calls `'s'`, `'p'` or `'avg'`, the mean of the two, which is the default. In the phase calls `'s'` or `'p'`; it may be left out at normal incidence, where the two are the same. The monitoring calls take none and return both. |
| `candidates` | Needle calls only, and required there: the materials that could be inserted, a list of indices, each a number or `{ n, k }` as a medium is. |
| `fractions` | Needle calls only: positions inside each layer, as fractions of its thickness from the incident side, where the needle is also tried. None when left out. |
| `top` | `monitorCurve` only, and required there: the index of the growing layer, a number or `{ n, k }` as a medium is. |
| `thicknesses` | `monitorCurve` only, and required there: the thicknesses of the growing layer to sample, in nm. |

### Results over a grid

Each value is a `Float64Array` with one entry per wavelength. `N`, where a result has it, is the number of layers. A derivative over a grid is flat and wavelength-major: wavelength `i`, layer `k` at `[i * N + k]`.

### Errors

A wrong argument throws a `TypeError` that names it, such as `layer 0 needs {n, d}`, `pol must be 's', 'p' or 'avg'` or `unknown option 'polarization'; the options are angle and pol`.

## R, T and A

### `rta(stack, wavelength, options?)`

```js
const { R, T, A } = rta(stack, 550, { pol: 's' });
```

**Returns** `{ R, T, A }`: reflectance, transmittance into the substrate, and the absorptance of the layers, each a fraction of the incident irradiance.

$A$ is the net irradiance entering the front of the stack less the irradiance leaving it into the substrate (Macleod, Eq. 2.120), so it is what the layers absorb. With a transparent incident medium $R + T + A = 1$ by construction, so the sum is no check on the calculation.

With an absorbing incident medium they do not sum to one. The incident and reflected waves interfere in the irradiance of the incident medium, and

$$R + T + A = 1 + 2\,\frac{\mathrm{Im}\,\eta_0}{\mathrm{Re}\,\eta_0}\,\mathrm{Im}(r)$$

in tmmcore's sign convention. $T$ can then exceed one: Macleod finds $T = 1 + k_0^2/n_0^2$ for a surface whose admittance is the complex conjugate of the incident one (Eq. 2.84). See [Validation](validation.md#an-absorbing-incident-medium).

### `rtaSpectrum(stack, wavelengths, options?)`

R, T and A over a wavelength grid, s and p in one pass.

```js
const { R, Rs, Rp } = rtaSpectrum(stack, wavelengths, { angle: 45 });
```

**Returns** `{ R, T, A, Rs, Ts, As, Rp, Tp, Ap }`: R, T and A for the `pol` asked, and the s and p values they come from.

## Thickness derivatives

The derivatives come from the same characteristic-matrix product as the spectrum, exactly, not by finite differences. They are the derivatives of what `rta` computes, [opaque layers](#opaque-layers) included, and [Validation](validation.md#derivatives-against-finite-differences) checks them against finite differences.

They line up with `stack.layers` as given. A layer of zero thickness keeps its derivative, the rate at which the spectrum moves as that layer starts to grow. A layer of negative or NaN thickness is absent, and its derivatives are zero. With `pol: 'avg'` the derivatives are the mean of the s and p ones.

### `jacobian(stack, wavelength, options?)`

```js
import { jacobian } from 'tmmcore';

const { R, dRdd } = jacobian(stack, 550, { pol: 's' });   // dRdd[j] = ∂R/∂d_j, in 1/nm
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, N }`: R, T and A, and their derivatives as `Float64Array`s of length `N`.

The spectrum comes with its derivatives, so a Gauss-Newton or Levenberg-Marquardt step costs one evaluation rather than $N + 1$. See the [refinement example](examples.md#refinement-with-the-analytic-jacobian).

### `jacobianSpectrum(stack, wavelengths, options?)`

The Jacobian at every wavelength of a grid. **Returns** R, T and A over the grid, `dRdd`, `dTdd` and `dAdd` with nLam × N values each, and `N`.

### `hessian(stack, wavelength, options?)`

```js
import { hessian } from 'tmmcore';

const { d2Rdd } = hessian(stack, 550, { pol: 's' });   // d2Rdd[i][j] = ∂²R/∂d_i∂d_j
```

**Returns** `{ R, T, A, dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, N }`: the fields of `jacobian`, and the second derivatives as `N` rows of `N`, symmetric. The cost is $O(N^2)$ small-matrix operations. See the [tolerance example](examples.md#what-thickness-errors-cost-from-the-exact-hessian).

### `hessianSpectrum(stack, wavelengths, options?)`

The Hessian at every wavelength of a grid. **Returns** the fields of `jacobianSpectrum`, and `d2Rdd`, `d2Tdd` and `d2Add` with nLam × N × N values each, wavelength `i` and layers `j`, `k` at `[(i * N + j) * N + k]`.

## Needle scan

The derivatives of R, T and A for inserting a layer of zero thickness of a candidate material, at every interface of the stack at once and, if asked, inside its layers. This is the analytic P-function of Tikhonravov et al. (*Appl. Opt.* **35**, 5493, 1996), the $d \to 0$ limit of the numerical method of Sullivan and Dobrowolski (*Appl. Opt.* **35**, 5484, 1996), so it needs no trial thickness and no second evaluation.

### `needleScan(stack, wavelength, options)`

```js
import { needleScan } from 'tmmcore';

const scan = needleScan(stack, 550, { pol: 's', candidates: [2.35, 1.46] });
// scan.gaps[position][candidate] = { dR, dT, dA }
```

**Returns** `{ R, T, A, gaps, intra, N }`.

- `gaps` has $N + 1$ entries, one per interface: 0 is in front of the first layer and `N` against the substrate. Each holds one `{ dR, dT, dA }` per candidate, per nm of needle. Insert where the derivative of your merit function is negative.
- `intra[k][i]` is `{ frac, perCand }` for a needle inside layer `k` at fraction `fractions[i]`, with `perCand` shaped like an entry of `gaps`. Fractions 0 and 1 give the interfaces either side of the layer.

The positions count the layers as given, zero-thickness ones included. A needle of a layer's own material at its boundary is that layer made thicker, so `gaps[0][c].dR` for a candidate with layer 0's index equals `dRdd[0]` from `jacobian`; the two are computed independently. With `pol: 'avg'` the derivatives are the mean of the s and p ones. See the [needle example](examples.md#the-needle-p-function-where-the-next-layer-goes).

### `needleScanSpectrum(stack, wavelengths, options)`

The needle scan at every wavelength of a grid, where a candidate's `n` and `k` may also be one value per wavelength. **Returns** R, T and A over the grid, `gaps` and `intra` each as `{ dR, dT, dA }`, and `N`. Each list is flat and wavelength-major: wavelength `i`, interface `p` and candidate `c` at `[(i * (N + 1) + p) * nCand + c]` in `gaps`, and wavelength `i`, layer `k`, fraction `f` and candidate `c` at `[((i * N + k) * nFrac + f) * nCand + c]` in `intra`.

## Phase, GD, GDD and TOD

Phase, group delay, group delay dispersion and third-order dispersion of r and t. They depend on how each index changes with frequency, so these calls take an index as a number, which has no dispersion, or as a [jet](#indices-as-jets) that carries it.

GD, GDD and TOD are exact frequency derivatives at the wavelength asked, from the frequency derivatives of the characteristic matrix. No phase unwrapping and no wavelength stencil is involved, so a value does not depend on the grid around it.

### `phase(stack, wavelength, options?)`

```js
import { phase } from 'tmmcore';

const { r, t } = phase(stack, 800, { pol: 's' });   // r.gdd: GDD on reflection, in fs²
```

**Returns** `{ r, t }`, each `{ phaseRad, phaseDeg, gd, gdd, tod, magnitudeSquared }`, or `null` where that coefficient is exactly zero and has no phase. Units are in [Units and signs](#units-and-signs).

Both are also `null` at the exact critical angle of a layer, and at the substrate's critical angle in p. At the substrate's critical angle in s the phase and `magnitudeSquared` come back, and `gd`, `gdd` and `tod` are `NaN`.

### `phaseSpectrum(stack, wavelengths, options?)`

`phase` at every wavelength of a grid. **Returns** `{ r, t }`, each `{ phaseRad, gd, gdd, tod, magnitudeSquared }` with one value per wavelength, `NaN` where the coefficient vanishes.

### `phaseJacobian(stack, wavelength, options?)`

The phase quantities and their exact derivatives in every layer's thickness.

```js
import { phaseJacobian } from 'tmmcore';

const { r } = phaseJacobian(stack, 800, { pol: 's' });   // r.dGdd[j] = ∂GDD/∂d_j, in fs²/nm
```

**Returns** `{ r, t }`, each the fields of `phase` and `dPhaseDeg`, `dGd`, `dGdd`, `dTod` and `dLogMagnitudeSquared`, `Float64Array`s of length `N`. `dLogMagnitudeSquared` is the derivative of ln |r|² or ln |t|², the relative change in intensity. The derivatives line up with `stack.layers` as in the [thickness derivatives](#thickness-derivatives), and a GDD fit costs one evaluation per step rather than $N + 1$.

### `phaseJacobianSpectrum(stack, wavelengths, options?)`

`phaseJacobian` at every wavelength of a grid, for fitting measured phase quantities such as group delay or ellipsometric Ψ and Δ. **Returns** `{ r, t }`, each the five arrays of `phaseSpectrum` and the five derivatives with nLam × N values each.

### Indices as jets

A jet is a truncated Taylor series in angular frequency, `[f, f', f''/2!, f'''/3!]`, each entry a complex `[re, im]` pair. An index jet is $[\tilde n,\ d\tilde n/d\omega,\ (d^2\tilde n/d\omega^2)/2,\ (d^3\tilde n/d\omega^3)/6]$.

tmmcore has no material models. You build the jet from your own dispersion formula with the [jet arithmetic](#jet-arithmetic). A Cauchy index:

```js
import { phaseSpectrum, wavelengthOmegaJet, jetReciprocal, jetMultiply,
         jetScale, jetAdd, jetConstant, omegaFromLambdaNm } from 'tmmcore';

// n(λ) = A + B/λ² + C/λ⁴ with a constant k, λ in nm
function cauchyJet(lambda_nm, A, B, C, k) {
    const lambda = wavelengthOmegaJet(lambda_nm, omegaFromLambdaNm(lambda_nm));   // λ(ω)
    const inverseSquare = jetReciprocal(jetMultiply(lambda, lambda));
    return jetAdd(jetConstant(A, k),
        jetAdd(jetScale(inverseSquare, B),
               jetScale(jetMultiply(inverseSquare, inverseSquare), C)));
}

const grid = [700, 750, 800, 850, 900];
const glass = { n: grid.map(l => cauchyJet(l, 1.5, 4200, 0, 0)) };
const layers = [{ n: grid.map(l => cauchyJet(l, 2.21, 2.4e4, 0, 0)), d: 85 }];
const { r } = phaseSpectrum({ incident: 1, substrate: glass, layers }, grid, { pol: 's' });
```

Nothing there is differentiated by hand: `wavelengthOmegaJet` carries the derivatives of λ(ω), and every operation after it carries them on. Sellmeier and Lorentz formulas work the same way.

For tabulated data, interpolate with something at least three times differentiable and pass the value and its first three derivatives to `jetCompose`. A PCHIP interpolant is only $C^1$, and its GDD and TOD jump at every knot.

### An angle fixed outside the incident medium

Light that reaches a coating through its own substrate, at an angle fixed in air, meets the coating at an angle that changes with frequency: by Snell's law $\sin\theta_0 = \sin\theta_{\text{air}}/n_0(\omega)$ in the glass, $n_0$ the real part of its index. An `angle` set to the inside angle at each wavelength gets the phase right, but GD, GDD and TOD then leave out how that angle moves. Give the sine as a jet instead, with `cauchyJet` from above:

```js
import { jetDivide } from 'tmmcore';

// The two-layer coating seen from inside the glass, air beyond it, 45° in air.
const lams = [700, 750, 800, 850, 900];
const n0 = lams.map(l => cauchyJet(l, 1.5, 4000, 0, 0));
const sine = n0.map(n => jetDivide(jetConstant(Math.sin(Math.PI / 4)), n));
const back = { incident: { n: n0 }, substrate: 1, layers: [{ n: 1.46, d: 94.2 }, { n: 2.35, d: 58.5 }] };
const reflected = phaseSpectrum(back, lams, { pol: 'p', sine }).r;
```

At 800 nm this gives a GD of 1.3265 fs and a GDD of −0.132 fs² on reflection. The inside angle held at its 800 nm value gives 1.3247 fs and −0.136 fs².

### Units and signs

The phase calls differentiate with respect to the angular frequency $\omega = 2\pi c/\lambda$ in rad/fs, with wavelengths in nm, so GD is in fs, GDD in fs² and TOD in fs³.

| | |
|---|---|
| `C_NM_PER_FS` | the speed of light in vacuum, 299.792458 nm/fs |
| `omegaFromLambdaNm(lambda_nm)` | angular frequency in rad/fs |

tmmcore uses $\tilde n = n + ik$ with $e^{-i\omega t}$, the complex conjugate of Macleod's convention. `phaseRad` and `phaseDeg` are the phase $\phi$ in Macleod's convention, the negative of the argument of $r$ or $t$ as tmmcore computes them, and GD, GDD and TOD are $-d\phi/d\omega$, $-d^2\phi/d\omega^2$ and $-d^3\phi/d\omega^3$, as he defines them. So all of them equal his. Positive GDD means the red end of a pulse arrives first.

$r_p$ has the sign of thin-film practice, $r_p = r_s$ at normal incidence. [Conventions](getting-started.md#conventions) explains the 180° this puts into an ellipsometric $\Delta$, and the [ellipsometry example](examples.md#ellipsometry-of-an-oxide-on-silicon) computes one.

## Jet arithmetic

Jets are arrays of four complex `[re, im]` pairs. These calls are lower level than the ones above and may change between versions.

| | |
|---|---|
| `jetConstant(re, im?)` | a constant: the value, with zero derivatives |
| `jetFromDerivatives(f, f', f'', f''')` | a jet from plain derivatives, applying the factorials |
| `jetDerivatives(jet)` | the plain derivatives `[f, f', f'', f''']` back |
| `jetAdd`, `jetSubtract`, `jetMultiply`, `jetDivide`, `jetReciprocal` | arithmetic |
| `jetScale(jet, s)` | every order times a real number |
| `jetSqrt`, `jetExp`, `jetLog`, `jetPower`, `jetSin`, `jetCos`, `jetSinCos` | elementary functions |
| `jetCompose(value, [f', f'', f'''], xJet)` | the chain rule, for a function known through its derivatives |
| `wavelengthOmegaJet(lambda, omega)` | λ as a jet in ω, the usual starting point |
| `jetRealPart(jet)` | the real part, order by order |
| `jetWithImaginaryPart(realJet, imaginaryJet)` | a real jet and an imaginary jet combined into one complex jet |
| `jetClampImaginary`, `jetClampRealMinimum`, `jetIsFinite` | guards |
| `JET_ORDER` | 3 |

`wavelengthOmegaJet` needs no value of $c$: given λ and ω, every derivative follows from $\lambda' = -\lambda/\omega$, in the unit λ was given in.

## Monitoring and deposition

What an optical monitor reads while a coating is deposited. These calls take `angle` and no `pol`, and return s and p together: `{ Rs, Ts, Rp, Tp, Rrs, Rrp }`, the reflectance and transmittance of the coated surface, and its reflectance seen from the substrate side. A caller that treats the substrate as a thick slab combines all three itself; for a semi-infinite substrate, $A = 1 - R - T$.

### `monitorCurve(stack, wavelength, options)`

The signal of one layer growing on a finished stack, at one wavelength, for each thickness of it in `thicknesses`. `stack` is what is already deposited, and `top` the index of the growing layer.

```js
import { monitorCurve } from 'tmmcore';

const finished = { incident: 1, substrate: 1.52, layers: [{ n: 1.46, d: 94.2 }] };
const { Rs } = monitorCurve(finished, 550, { top: 2.35, thicknesses: [0, 20, 40, 58.5] });
```

The finished stack's matrix is built once, and each thickness then costs one 2×2 product, so a whole monitoring curve is one call: the incremental control algorithm of Tikhonravov and Trubetskov, *Appl. Opt.* **44**, 6877 (2005).

**Returns** `{ Rs, Ts, Rp, Tp, Rrs, Rrp }`, each a `Float64Array` with one value per thickness.

### `depositionSpectra(stack, wavelengths, options?)`

The spectrum after each layer of a deposition run, in one call. The layers are deposited from the substrate outward, the last one listed first. Each is added to the running product and the grid evaluated after it, so all N spectra together cost about what the last one costs alone. `n`, `k` and `angle` may be one per wavelength, as in a `Spectrum` call.

**Returns** `{ Rs, Ts, Rp, Tp, Rrs, Rrp }`, each a `Float64Array` of N × nLam values: after `k + 1` layers, the last `k + 1` of `stack.layers`, wavelength `i` at `[k * nLam + i]`. A layer of zero thickness repeats the spectrum before it.

### `growingStack(stack, wavelengths, options?)`

A layer growing on a finished stack over a whole wavelength grid, one thickness per call: what a broadband monitor reads at each scan. The finished stack is multiplied out once for every wavelength and polarization and kept, so each sample costs one layer matrix and one 2×2 product per wavelength and polarization. `n`, `k` and `angle` may be one per wavelength, as in a `Spectrum` call. **Returns a handle**:

| | |
|---|---|
| `setTop(top)` | sets the growing layer's index, given as a medium is, `n` and `k` each a number or one per wavelength. Call it before the first `sample` with `d > 0`, and again to change the material |
| `sample(d, out?)` | `{ Rs, Ts, Rp, Tp, Rrs, Rrp }` over the grid for the growing layer at thickness `d` nm; `d <= 0` gives the finished stack. Pass the previous result as `out` to reuse its arrays. Throws if called with `d > 0` before `setTop` |
| `free()` | releases what the handle holds. On the WebAssembly instance that is kernel memory: a dropped handle is reclaimed by a finalizer in time, but calling `free()` keeps a long run's memory flat. `setTop` and `sample` throw after it |

## WebAssembly

The C kernel compiled to WebAssembly, with the same calls. [Getting started](getting-started.md#turning-on-webassembly) shows how to load it in Node and in a browser.

```js
import { instantiateTmmWasm } from 'tmmcore';

const wasm = await instantiateTmmWasm(bytes);   // bytes: the contents of tmm_kernel.wasm
const { R } = wasm.rtaSpectrum(stack, wavelengths);
```

The grid calls are where the speed is: they cross from JavaScript into WebAssembly per grid, not per wavelength. The loader is also its own entry point, `tmmcore/wasm`, for code that needs only the WebAssembly side.

An instance is not shared between threads, since there is no shared memory. Each worker instantiates its own from the same bytes.

### The instance

| Method | |
|---|---|
| `rta`, `rtaSpectrum`, `jacobian`, `jacobianSpectrum`, `hessian`, `hessianSpectrum`, `needleScan`, `needleScanSpectrum`, `phase`, `phaseSpectrum`, `phaseJacobian`, `phaseJacobianSpectrum`, `monitorCurve`, `depositionSpectra`, `growingStack` | the calls above, with the same arguments and results |
| `jetConstant` … `jetIsFinite`, `wavelengthOmegaJet`, `omegaFromLambdaNm` | the jet arithmetic |
| `hasHessian()`, `hasPhase()`, `hasPhaseJacobianSpectrum()`, `hasGrowingKernels()`, `hasGrowingEval()` | always `true`: a module missing a kernel is refused when it loads |

If the WebAssembly heap is exhausted, a call that cannot allocate memory for its arguments and results throws `Error('tmmWasm: malloc failed')`; the one-wavelength calls reuse theirs from call to call and need more only for a larger stack. A kernel that cannot allocate its own working state writes `NaN` to every output instead: numbers and arrays come back `NaN`, and `phase` and `phaseJacobian` return `null` for both coefficients. `depositionSpectra` and `growingStack` throw.

### Loading

| | |
|---|---|
| `instantiateTmmWasm(bytes)` | instantiates from an `ArrayBuffer` or `Uint8Array`, and resolves to the instance |
| `initTmmWasmFromUrl(url)` | fetches and instantiates. Resolves `true` or `false` and does not throw; a second call returns the first call's promise |
| `getTmmWasm()` | the last instance loaded, or `null`. A failed `initTmmWasmFromUrl` also leaves it `null` |
| `isTmmWasmReady()` | whether an instance is loaded |
| `setTmmWasmEnabled(on)`, `isTmmWasmEnabled()` | a flag for whether to use the instance, off at first |
| `tmmWasmActive()` | whether the flag is on and an instance is loaded |
| `TmmWasmInstance` | the class of the instance |

For workers that cannot fetch the file themselves, the main thread hands them its bytes:

| | |
|---|---|
| `initTmmWasmMainThread(bytes, enabled)` | main thread: keeps the bytes for workers, then instantiates and turns the flag on if `enabled`, or turns it off if not. Resolves `true` when the flag is on |
| `getTmmWasmBytesForWorker()` | main thread: the bytes to send a worker, or `null` while the flag is off |
| `noteTmmWasmBytes(bytes)` | worker: starts instantiating from the bytes it was sent and turns its flag on; does nothing once an instance is loaded |
| `awaitTmmWasmReady()` | worker: resolves when that instantiation ends, `true` if an instance is loaded |

## Edge cases

### Grazing incidence and critical angles

The angle of incidence runs up to and including 90°. $\cos\theta_0$ is taken directly rather than as $\sqrt{1 - \sin^2\theta_0}$, which keeps nothing of it near 90°, and a layer of the incident index keeps that cosine too. At 90° exactly $R$ is one and $T$ zero to rounding.

At a medium's critical angle to the last bit, $\cos\theta$ in it is exactly zero and the layer matrix as written is $0/0$. Its limit is finite and tmmcore uses it: the layer's matrix becomes linear in its thickness, and a substrate at its critical angle reflects everything. The derivative calls take the same limit. The phase calls do not; see [`phase`](#phase-gd-gdd-and-tod).

### Opaque layers

The imaginary part of a layer's phase thickness grows with its thickness, and $\cosh$ of it overflows, which would turn the whole matrix product into `NaN`. So the imaginary phase is held at 50 when the matrix is built. Past that point the held matrix is the true one divided by a real factor, $e^{|\mathrm{Im}\,\delta| - 50}$, to within $e^{-100}$ relative. The factor cancels from $r$ and from every phase, and tmmcore carries it as a log scale for $t$, so $T$ is right however opaque the layer, until it falls below about $10^{-308}$, the range of a double, and reads 0.

The hold depends on the phase, not on $k$. A thick absorbing layer reaches it, and so does a lossless one beyond the critical angle, where the wave in it is evanescent: an air gap of 10 µm under a glass prism at 70° is well past it.

The derivatives include the factor, so $T$ falls with an opaque layer's thickness at the true rate, and $R$ does not move. In a deep stack the partial products outgrow the range of a double; each carries a power-of-two scale that cancels exactly from every derivative, so no wavelength drops to zero. A needle inside an opaque layer sees the field that is there: near the surface its derivatives are those of an ordinary layer, and deep inside $dR$ is zero while $dT$ keeps its value relative to $T$.

## The names of 0.4

0.5.x keeps the names of 0.4 with their arguments and results. Each prints one warning on its first call, and 0.6.0 removes them.

| 0.4, module | 0.4, WebAssembly instance | 0.5 |
|---|---|---|
| `tmm` | `tmmOne` | `rta` |
| | `tmmSpectrum` | `rtaSpectrum` |
| `tmmThicknessJacobian` | `tmmJacobian` | `jacobian` |
| `tmmThicknessHessian` | `tmmHessian` | `hessian` |
| `tmmNeedleScan` | `tmmNeedleScan` | `needleScan` |
| `tmmPhaseDispersion` | `tmmPhaseOne` | `phase` |
| | `tmmPhaseSpectrum` | `phaseSpectrum` |
| `tmmPhaseThicknessJacobian` | `tmmPhaseJacobian` | `phaseJacobian` |
| | `tmmPhaseJacobianSpectrum` | `phaseJacobianSpectrum` |
| | `monitorCurve`, `depositionSpectra`, with the arguments of 0.4 | `monitorCurve`, `depositionSpectra` |
| | `growingEval` | `growingStack` |
| `tmmCoefficientJets`, `tmmCoefficientThicknessJets`, `coefficientPhaseDispersion`, `coefficientPhaseThicknessDerivatives` | | none |

The one-wavelength names take `(lambda_nm, theta_deg, pol, n0, ns, layers)`, each index as `[re, im]` and the layers as `{ n, d }[]`, or in the phase names `{ nJet, d }[]` with the media as jets; the instance takes `pol` as `0` or `1`. `tmmNeedleScan` takes after them the candidates, each as `[re, im]`, and the fractions. The [0.4.2 API reference](https://github.com/aai2k/tmmcore/blob/0a432b7/docs/api.md) has them in full.

The instance's `monitorCurve` and `depositionSpectra` tell the two forms apart by the first argument: a stack in 0.5, a wavelength or a list of them in 0.4. `depositionSpectra` and `growingEval` in 0.4 take the wavelengths, the media and the layers as lists of `[re, im]`; `depositionSpectra` takes the layers in the order they are deposited, the first deposited first, where the 0.5 call takes the stack's order.

The 0.4 phase names take the incident sine as the option `sinTheta0Jet`, and over a grid `sinJets`, where the 0.5 calls take `sine`. The 0.5 calls take no `omega`, which sets the time unit; the 0.4 phase names still do.

## Low-level primitives

Exported for building transfer-matrix variants on the same conventions. These are lower level than the calls above and may change between versions. Complex numbers are `[re, im]` pairs.

| | |
|---|---|
| `cadd`, `csub`, `cmul`, `cdiv` | complex arithmetic |
| `cabs2`, `cconj`, `csqrt`, `ccos`, `csin` | modulus squared, conjugate, square root, cosine and sine |
| `creal`, `cimag` | the parts |
| `matmul(A, B)` | 2×2 complex matrix product |
| `cmatvec(M, v)` | 2×2 matrix times a 2-vector |
| `snellCosTheta(n0, sinTheta0, nj, cosTheta0?)` | complex $\cos\theta$ in a medium, from the real invariant $\mathrm{Re}(n_0)\sin\theta_0$. Pass `cosTheta0` as `[cos θ0, 0]`, as `rta` does, and a medium with an index close to the incident one keeps its cosine at grazing incidence |
| `incidentCosTheta(n0, sinTheta0, cosTheta0?)` | $\cos\theta_0$ of the incident medium itself: its own cosine when it is transparent, from the same invariant when it absorbs |
| `layerMatrix(nj, dj_nm, lambda_nm, cosTheta_j, pol)` | the characteristic matrix of one layer; its limit where `cosTheta_j` is exactly zero, and past the [opaque-layer hold](#opaque-layers) the held matrix |
| `layerLogScale(nj, dj_nm, lambda_nm, cosTheta_j)` | the log of the factor `layerMatrix` leaves out past the hold, 0 short of it |
| `rescaleMatrix(M)` | rescales in place past an overflow threshold and returns the log scale it took out |

`layerLogScale` and `rescaleMatrix` keep opaque stacks finite and their transmittance right. Both are real factors that cancel from reflectance but not from transmittance, so a caller building its own product adds both to one log scale and multiplies T by $e^{-2\,\text{logScale}}$, as `rta` does. A field computed from a partial product carries its own log scale the same way.
