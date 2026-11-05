# Validation

Every check on this page runs from a clone with Node, except the native C build, which needs a C compiler, and the Python package's tests, which need Python. The outputs shown are from tmmcore 0.5.0 on 5 November 2026, and yours should match them to the last digit or two. `npm test` runs every test file below, `npm run compare` the comparison with Byrnes' `tmm`, and `npm run examples` the examples; CI runs all three on Linux, macOS and Windows.

Each check answers a different question, and none of them covers the others. Two implementations can agree perfectly and both be wrong; matching a closed form in one special case says nothing about the general one.

| Check | Against | Command | Result |
|---|---|---|---|
| [Essential Macleod](#against-essential-macleod) | a commercial design program, on ten designs | `node tests/macleod.mjs` | worst $1.4\times10^{-14}$ in $R$, $T$, $A$ |
| [Byrnes' `tmm`](#against-byrnes-tmm) | an open-source implementation by another author | `npm run compare` | worst $8.6\times10^{-14}$ in $R$ and $T$ |
| [A quarter-wave layer](#closed-forms) | its closed form | `node examples/01-single-layer.mjs` | $1.4\times10^{-17}$ in $R$ |
| [GD, GDD and TOD](#group-delay-gdd-and-tod) | the closed form for a slab with no interfaces | `node examples/12-group-delay.mjs` | worst $2.0\times10^{-13}$ relative, in TOD |
| [Grazing incidence, critical angle](#grazing-incidence-and-the-critical-angle) | Fresnel's equations, and the limit of the layer matrix | `node tests/edge_cases.mjs` | passes |
| [Jet arithmetic](#the-jet-arithmetic) | derivatives written out by hand | `node tests/jets.mjs` | 112 derivatives pass |
| [Absorbing incident medium](#an-absorbing-incident-medium) | the energy identity | `node tests/absorbing_incident.mjs` | worst gap $1.6\times10^{-15}$ |
| [Opaque layers](#transmittance-through-opaque-layers) | the Airy recursion | `node tests/opaque_transmittance.mjs` | 750 checks pass |
| [Derivatives](#derivatives-against-finite-differences) | finite differences of `rta()` and `phase()` | `node tests/derivatives_fd.mjs` | worst at 0.24 of its tolerance |
| [JavaScript and WebAssembly](#javascript-against-webassembly) | each other | `node tests/equivalence.mjs` | worst $4.4\times10^{-16}$ in $R$, $T$, $A$ |
| [Heap exhausted](#with-the-heap-exhausted) | the rule for a failed allocation: `NaN` outputs, or a status of 0 | `node tests/wasm_allocation.mjs` | 486 outputs of 13 kernels pass |
| [The C, built natively](#the-c-built-natively) | compiler warnings, as errors | CI, with GCC, Clang and MSVC | no warnings |

## Against Essential Macleod

Essential Macleod, from Thin Film Center, computed ten designs made for this check. Every material is a constant $n$ and $k$ of the test's own, so no material data from either program enters. Its results are committed, so the check needs no copy of the program; `tests/macleod/export.ps1` writes them again through its COM automation server on a machine that has it.

```bash
node tests/macleod.mjs
```

```
Essential Macleod 11.09.0605: 10 designs, 71142 comparisons
  worst R, T, A            1.4e-14  (js grazing 89° p 690 nm R)
  worst substrate-side R   2.1e-14  (js metal-mirror 70° p 510 nm Rback)
  worst phase              3.4e-12  (js frustrated-tir 30° p 530 nm r phase)
  worst Psi, Delta         2.8e-12  (js frustrated-tir 30° 530 nm Delta of r)
  worst GD, GDD, TOD       3.6e-11  (js frustrated-tir 30° p 530 nm r tod (closed form), fs^3 over 1)
  worst locus R            1.0e-15  (js metal-mirror locus 560 nm 45° s at depth 42.5 nm, monitorCurve)
PASS : Essential Macleod
```

Each design is scanned across the visible at two or three angles, in s and p:

| Design | What it reaches |
|---|---|
| nine quarter-wave layers at 560 nm | a stop band |
| a low-index layer over a metal | a mirror out to 70° |
| four layers, two of them absorbing | absorption inside the stack |
| two layers on an absorbing substrate | transmittance into an absorbing medium, and reflectance seen from inside it |
| three layers under an incident medium with $k = 0.001$, and again with $k = 0.016$ | the absorbing incident medium |
| a 400 nm air gap between glasses | frustrated total reflection at 45° and 60°, and a reflection zero at 30° |
| one layer at 80°, 89° and 89.9° | grazing incidence |
| 1 µm of metal | the opaque-layer clamp |
| a 20 µm layer, 500 to 600 nm in 2.5 nm steps | a phase that turns fast with wavelength, and long group delay |

In each scan the test compares $R$, $T$ and $A$, and for natural light too; the phase of $r$ and of $t$, with the GD, GDD and TOD of each; Psi and Delta of $r$ and of $t$; and the reflectance from the substrate side. For five of the designs Essential Macleod also gives monitoring curves, the reflectance at one wavelength with the coating cut at each depth, which `monitorCurve`, `growingStack` and `depositionSpectra` must reproduce. Everything runs on the JavaScript and the WebAssembly kernel, except the delays under an absorbing incident medium, which only the JavaScript kernel gives.

| Quantity | As Essential Macleod reports it | Compared with |
|---|---|---|
| index | $n - ik$ | $n + ik$ |
| phase of $r$ and $t$ | degrees, with Macleod's sign | `phaseRad`, which carries the same sign |
| Psi of $r$ and of $t$ | $\arctan\sqrt{R_p/R_s}$ and $\arctan\sqrt{T_p/T_s}$ | `rtaSpectrum` in s and p |
| Delta of $r$ | $\varphi_p - \varphi_s + 180^\circ$ | `phaseSpectrum` in s and p |
| Delta of $t$ | $\varphi_p - \varphi_s$ | `phaseSpectrum` in s and p |
| phase of $t$ where $T < 10^{-28}$ | written as 0 | not compared |
| reflectance from the substrate side | at the angle of incidence taken in the substrate | `rtaSpectrum` of the reversed stack |

$R$, $T$ and $A$ are held to $10^{-12}$, phases, Psi and Delta to $10^{-10}$ degrees, and GD, GDD and TOD to $10^{-9}$ of the largest value in the scan, or of 1 fs, fs² or fs³ where that is larger. Each bound is the worst difference over the ten designs, rounded up to a power of ten, times ten.

Two places take another reference. Essential Macleod's reflection GDD and TOD lose digits where $r$ nears zero: the air gap at 30° passes through $r = 0$ near 529 nm, and at 530 nm its TOD is off by $2\times10^{-6}$ fs³. The reflection delays of that scan are held instead to the closed form for one film between identical media, from the characteristic matrix of Macleod's Eq. 2.113. Under an absorbing incident medium the two programs define $R$ and $T$ differently, as described [below](#an-absorbing-incident-medium): tmmcore's $r$ and $t$ are carried to Essential Macleod's through the surface admittance, and their GD, GDD and TOD through the frequency jets of $r$ and $t$.

The Python package runs the same comparison in `python/tests/test_macleod.py`, except the delays under an absorbing incident medium, which need calls the package does not expose.

## Against Byrnes' `tmm`

The same inputs given to [Steven Byrnes' `tmm`](https://github.com/sbyrnes321/tmm)[^byrnes]: MIT, pure Python, sharing no code and no author with tmmcore.

It uses the **same** $\tilde n = n + ik$ convention and time factor as tmmcore, so values transfer verbatim with no conjugation. There is no material-data confound: the wavelength grid, the complex indices and the thicknesses are precomputed into a shared file that both read. Only the mathematics differs.

The cases run at normal incidence, and at 30° and 60°, where s and p part and the p admittance $n/\cos\theta$ is under test. R and T come from `rta()`, and the complex reflection coefficient from `phase()`. $r_s$ is compared as it stands and $r_p$ with its sign reversed, since Byrnes takes the Fresnel sign for p (see [Conventions](getting-started.md#conventions)), so the comparison also checks that statement about the sign.

Byrnes' outputs are committed, so this needs no Python:

```bash
npm run compare
```

```
case               angle layers  points   max |Δ| JS  max |Δ| WASM  max |Δr|
----------------------------------------------------------------------------
AR4/g71               0°      4      71      2.1e-15       2.1e-15   4.6e-16
HR21/g71              0°     21      71      1.2e-14       1.2e-14   6.0e-15
AG7/g71               0°      7      71      1.4e-15       1.6e-15   1.4e-15
BIG40/g71             0°     40      71      1.7e-14       1.7e-14   2.3e-14
AR4/g701              0°      4     701      2.6e-15       2.6e-15   7.5e-16
HR21/g701             0°     21     701      2.4e-14       2.4e-14   1.4e-14
AG7/g701              0°      7     701      3.1e-15       3.1e-15   1.8e-15
BIG40/g701            0°     40     701      8.6e-14       8.6e-14   5.6e-14
AR4/g71/30deg        30°      4      71      2.3e-15       2.1e-15   1.0e-15
HR21/g71/30deg       30°     21      71      1.7e-14       1.7e-14   1.5e-14
AG7/g71/30deg        30°      7      71      2.0e-15       1.9e-15   2.6e-15
BIG40/g71/30deg      30°     40      71      5.4e-14       5.4e-14   3.9e-14
AR4/g71/60deg        60°      4      71      2.4e-15       2.4e-15   1.0e-15
HR21/g71/60deg       60°     21      71      1.5e-14       1.5e-14   2.8e-14
AG7/g71/60deg        60°      7      71      2.6e-15       2.6e-15   2.0e-15
BIG40/g71/60deg      60°     40      71      8.1e-14       8.1e-14   5.0e-14

14624 values of R and T compared across 16 cases, both polarizations,
and the complex r at each of their 7312 pairs of wavelength and polarization.
Worst disagreement with an independently written implementation: 8.6e-14 in R and T, 5.6e-14 in r
```

The worst case is the forty-layer stack, where round-off accumulates through the longest matrix product. The oblique cases are held to the same $10^{-12}$ as the rest. Regenerating the reference file rather than trusting the committed one takes one `pip install`, of `tmm` and `numpy`; see [`benchmarks/README.md`](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md).

[Comparison with other packages](comparison.md) times tmmcore against Byrnes' `tmm` and three further Python packages, and checks $R_s$ of one of them, `tmm_faster`, against Byrnes'.

[^byrnes]: S. J. Byrnes, *Multilayer optical calculations*, arXiv:1603.02720.

## Closed forms

For a single quarter-wave layer at normal incidence, reflectance has an exact solution, Macleod's Eq. 4.5 at a quarter wave:

$$R = \left(\frac{n_0 n_s - n_1^2}{n_0 n_s + n_1^2}\right)^2$$

```bash
node examples/01-single-layer.mjs
```

```
R from tmmcore           0.012600790214630274
R closed form            0.012600790214630288
difference               1.39e-17
```

Below double-precision epsilon, $2.2\times10^{-16}$. The example fails past $10^{-15}$, and CI runs every example.

### Group delay, GDD and TOD

A slab whose surrounding media share its own index has no interfaces, so nothing reflects and the transmission coefficient is a pure phase, $t = e^{i\delta}$ with $\delta = n(\omega)\,\omega d/c$. The three phase quantities then reduce exactly to bulk propagation:

$$\mathrm{GD} = \frac{d}{c}(n + \omega n'), \qquad
  \mathrm{GDD} = \frac{d}{c}(2n' + \omega n''), \qquad
  \mathrm{TOD} = \frac{d}{c}(3n'' + \omega n''')$$

A Cauchy index is exactly quadratic in $\omega$, which makes the right-hand sides elementary.

```bash
node examples/12-group-delay.mjs
```

```
matched slab, 1 µm of n(λ) = 1.45 + 3600/λ², at 800 nm

            tmmcore            closed form         relative
  GD       4.89296832143789     4.89296832143789   1.8e-16  fs
  GDD    0.0478126142151489   0.0478126142151501   2.6e-14  fs²
  TOD    0.0203063517881787   0.0203063517881746   2.0e-13  fs³
```

The example fails if any of the three drifts past its agreement here, rounded up to a power of ten, times ten.

!!! warning "GDD and TOD lose precision on thick elements"

    GDD is a difference of terms of order $\mathrm{GD}^2$, and TOD of terms of order $\mathrm{GD}^3$. On a thick element those terms are far larger than GDD and TOD themselves, so both lose digits as the element thickens, TOD the faster. The same command goes on to print, for thicker slabs of the same material:

    | Thickness | Group delay | Relative error on GDD | Relative error on TOD |
    |---|---|---|---|
    | 1 µm | 4.9 fs | $3\times10^{-14}$ | $2\times10^{-13}$ |
    | 10 µm | 48.9 fs | $3\times10^{-13}$ | $4\times10^{-10}$ |
    | 100 µm | 489 fs | $8\times10^{-12}$ | $4\times10^{-8}$ |
    | 1 mm | 4893 fs | $1\times10^{-11}$ | $7\times10^{-7}$ |

    Coatings live at the top of that table and are unaffected. For the dispersion of millimetres of glass, differentiate the propagation phase directly instead of asking a transfer matrix for it: the closed forms above are the whole calculation and they do not cancel.

    GD is not such a difference, and does not lose digits this way.

### Grazing incidence and the critical angle

Two angles where a formula that is right everywhere else loses everything. Near 90°, $\sqrt{1 - \sin^2\theta_0}$ keeps none of $\cos\theta_0$: computed that way, $1 - R$ of a bare interface is wrong by $7\times10^{-4}$ relative at 89.99999°, and $R$ is `NaN` in p at 90°. At a medium's critical angle to the last bit, $\cos\theta$ in it is exactly zero and the layer matrix as written is $0/0$.

```bash
node tests/edge_cases.mjs
```

```
exact critical angle found at θ0 = 30.000000000000004°
WebAssembly comparisons: 772
PASS : edge cases
```

At grazing incidence the check is Fresnel's equations for a bare interface, $T = 4\eta_0\eta_s/(\eta_0 + \eta_s)^2$, with the exact cosine; the test holds $T$ to $10^{-12}$ relative from 89.9° to 90°, and an air layer under air changes nothing. The critical angle is searched for among the doubles near 30°, from glass of index 2 into air, and $R$ and $T$ of an air gap there match the limit of the layer matrix in closed form, and the doubles either side of it, which the matrix reaches smoothly. A substrate at its critical angle reflects everything. The same file checks that a negative or NaN thickness leaves a layer out of every evaluator, with zero derivative in both Jacobians, and that $r_p = r_s$ at normal incidence.

### The jet arithmetic

Every GD, GDD and TOD is read off Taylor jets in frequency, carried through the transfer matrix by the operations of `src/js/jet.js`. Each operation is checked on its own: fed the jet of the variable at a complex point, it must return the value and the first three derivatives of the function it stands for, written out by hand, to $10^{-13}$ of the value or of 1, whichever is larger.

```bash
node tests/jets.mjs
```

```
112 derivatives against closed forms
PASS : jet arithmetic
```

`tests/kernel_exports.mjs` holds the WebAssembly jets to these.

## An absorbing incident medium

Light arriving from inside an absorbing medium, a cemented cube or an immersion liquid, is the one case where $R + T$ over lossless layers is not 1, and the size of the departure is fixed by the definitions. Every wave in the stack shares the real invariant $n_0 \sin\theta_0$ of Macleod's Eq. 2.85, with $n_0$ the real part of the incident index as Essential Macleod takes it, and the transmittance carries the incident admittance the way Macleod's Eq. 2.83 does. What is left is the interference of the incident and reflected waves inside the absorbing medium:

$$1 - R - T = -2\,\frac{\mathrm{Im}\,\eta_0}{\mathrm{Re}\,\eta_0}\,\mathrm{Im}(r)$$

in tmmcore's sign convention, exact at every angle and polarization. On a bare interface $\mathrm{Im}(r)$ is itself of order $k_0$. Into the complex conjugate of the incident admittance, Macleod's matched surface (Eq. 2.84), $T = 1 + k_0^2/n_0^2$ as he finds, and $R = k_0^2/n_0^2$.

```bash
node tests/absorbing_incident.mjs
```

```
worst gap in the energy identity: 1.56e-15
the former complex invariant broke it by 5.10e+1 per unit k0
WebAssembly comparisons: 924
PASS : absorbing incident medium
```

The second line shows what the test guards against: carrying the complex index into Snell's invariant makes the incident wave's amplitude vary along the interface, energy flows sideways inside lossless layers, and $R + T$ exceeds 1 by a further term linear in $k_0$ that grows with the angle and with the stack's resonance.

The same test checks $A$ against what the layers absorb: the net irradiance $\tfrac12\mathrm{Re}(BC^*)$ entering the front of the stack less $\tfrac12\mathrm{Re}(\eta_s)$ leaving it (Macleod Eq. 2.120), per unit incident irradiance, computed from $[B, C]$ rather than from $r$. It holds $A$ to within $10^{-14}$ of that on a stack with a silver layer, and below $10^{-14}$ on the lossless 21-layer stack.

Essential Macleod reports $R$ and $T$ differently here. It takes the angles from the same real invariant, but works with the real part $a$ of the incident tilted admittance $\eta_0$, so over lossless layers its $R + T$ is exactly 1. The coating is the same in both: from the surface admittance $Y = \eta_0(1 - r)/(1 + r)$ of tmmcore's $r$, Essential Macleod's reflection coefficient is $(a - Y)/(a + Y)$, and the [comparison with it](#against-essential-macleod) matches its values to round-off that way.

## Transmittance through opaque layers

Past the imaginary-phase clamp a layer's matrix is held short by a real factor, which tmmcore carries as a log scale for $t$ (see [Opaque layers](api.md#opaque-layers)). The check is the Airy recursion, which builds $r$ and $t$ outward from the substrate from interface coefficients and propagation factors $e^{i\delta}$, never forms $\cosh$, and shares no code with the characteristic matrix. Every evaluator that reports $T$ or $|t|$ is held to it: `rta()`, $T$ from the Jacobian, Hessian and needle scan, $|t|^2$ and the phase of $t$ from the phase kernel, and in WebAssembly the same plus the batched spectrum and the three growing-stack kernels.

```bash
node tests/opaque_transmittance.mjs
```

```
750 checks against the Airy recursion, WebAssembly included.
without the held factor, T would be off by up to 2.7e+123
PASS : transmittance past the opaque-layer bound.
```

The stacks are 1 µm and 2.5 µm of aluminium at 550 nm, a 10 µm air gap in glass at 45° and 80°, and a dielectric stack as a control. The test holds $T$ to $10^{-12}$ relative everywhere, down to $T$ below $10^{-167}$.

Without the held factor, $T$ stops near $10^{-43}$ for each clamped layer whatever its thickness. $R$ is the same either way, since the factor cancels from it.

## Derivatives against finite differences

The JavaScript and WebAssembly derivative kernels share their formulas, so agreeing with each other proves nothing about the formulas. This check holds every derivative to finite differences of `rta()` and `phase()`, which share nothing with the kernels but the layer matrix: the thickness Jacobian, every entry of the Hessian, the needle P-function at every gap and inside every layer, and the thickness derivatives of phase, GD, GDD, TOD and $\ln|c|^2$ for both $r$ and $t$. A needle of the material of the layer beside it or around it is that layer made thicker, so there the needle scan must also equal the Jacobian. Through a clamped layer $T$ is below $10^{-43}$, where an absolute tolerance sees nothing, so every derivative of $T$ is also held relative to $T$.

```bash
node tests/derivatives_fd.mjs
```

```
23564 comparisons against finite differences and the port, 9342 of them WebAssembly.
worst difference as a fraction of its tolerance: 2.42e-1, six clamped layers θ=30 λ=550 p r dGdd[3]
PASS : analytic derivatives match finite differences.
```

The stacks are the ones where derivative kernels go wrong while R and T stay right:

| Stack | What it reaches |
|---|---|
| 700 nm and 1500 nm of Al between silica, 550 nm | the imaginary-phase clamp, just past it and far past it |
| 3 µm and 12 µm of a layer with $k = 5$, 532 nm | the clamp, where the unclamped phase would overflow $\cosh$ |
| 10 µm air gap under a glass prism at 70° | the clamp with $k = 0$, on an evanescent wave |
| 16 pairs of ZnS and Ge at 579 nm, where Ge absorbs | partial products past $10^{77}$, beyond which a plain complex division by $(\eta_0 B + C)^2$ overflows |
| an air gap between glasses of index 2, at its exact critical angle | a layer matrix that is $0/0$ as written, and a needle of air there |
| three layers with Al, under an incident medium with $k = 0.01$ | the absorptance's own term in $\mathrm{Im}\,r$, which $dR$ and $dT$ do not carry |
| four dielectric layers | nothing, as a control |
| one and six 20 µm opaque layers, phase kernel | the clamp, and a matrix product past the rescale threshold |
| nine dielectric layers, phase kernel, from air and from an incident medium with $k = 0.016$ | TOD, which passes $10^5$ fs³ at a resonance |

Differences use five-point central stencils in the layer thickness, with a step of 0.01 nm, or 0.25 nm for the phase quantities through the opaque layers, where GDD carries enough noise to need the wider step. Their TOD carries more noise than a difference at any step resolves, so its thickness derivatives there are left to the JavaScript ⇆ WebAssembly check below. A needle cannot have negative thickness, so its differences are one-sided, third order. The worst case above is that GDD noise.

## JavaScript against WebAssembly

The C kernel is a line-by-line port of the JavaScript. They are driven with identical inputs across absorbing, dispersive and oblique-incidence cases, in both polarizations, and every returned quantity is compared: R, T, A, the thickness Jacobian, the thickness Hessian, the needle P-function at gaps and inside layers, and the phase quantities with their thickness derivatives, point by point and over grids.

```bash
node tests/equivalence.mjs
```

```
94656 comparisons across 4 stacks, 4 wavelengths, 4 angles, s and p.
worst |Δ| on R/T/A     : 4.44e-16  (tolerance 1e-9)
worst |Δ| on derivatives: 2.22e-16  (tolerance 1e-12 abs / 1e-7 rel)
worst relative, above the 1e-12 floor: none exceeded the floor

34448 phase comparisons (phase, GD, GDD, TOD and their thickness derivatives, point and batched).
worst |Δ| on phase quantities: 5.52e-6
worst relative, above the 1e-12 floor: 1.05e-11

PASS : JavaScript and WebAssembly agree.
```

The two do not agree bit for bit. The only divergence is libm: the WebAssembly build uses musl's `sin`/`cos`/`exp`/`atan2`, the JavaScript engine uses its own, and they differ at roughly one unit in the last place.

The phase quantities amplify that noise, because every derivative order is a difference of nearly equal terms, which is why their worst relative difference above is larger than that of R, T and A.

One stack is held to a looser bound: six 20 µm opaque layers, whose product passes the rescale threshold. There the $q$-th frequency order is a difference of terms of order $\mathrm{GD}^q$, with GD the group delay through the whole stack, about 800 fs, so two roundings of it agree to about $64\,\varepsilon\,\mathrm{GD}^q$ and no better. The derivatives of $r$ with respect to the opaque layers are zero to exactly that noise; those of $t$ are not, since $t$ falls and turns with each opaque layer's thickness.

The test skips cleanly if `tmm_kernel.wasm` has not been built.

The Python package runs the same C, compiled natively, and its tests hold every call to the WebAssembly build on these cases; see [Python](python.md#agreement-with-the-npm-package).

### With the heap exhausted

In wasm32 a null pointer is address 0 of linear memory, so a kernel that works through a failed allocation does not trap: it writes over whatever lives there and returns numbers. This test takes up the whole heap of a fresh instance, block by halving block. Every kernel that allocates working state must then fill its outputs with `NaN`, or, if it returns a status, return 0 and leave its outputs as they were. Once the heap is returned, the same instance must compute correctly again.

```bash
node tests/wasm_allocation.mjs
```

```
heap exhausted at 2048 MiB of linear memory
486 outputs of 13 kernels checked without memory
PASS : allocation failure
```

## The C, built natively

The kernel is C99 with no dependencies beyond libm, so it can be checked outside WebAssembly:

```bash
cc -std=c99 -pedantic -Wall -Wextra -O2 -c src/c/tmm_kernel.c
```

It compiles without warnings under GCC 16.1.0. GitHub Actions compiles the same sources on Linux with GCC and on macOS with Clang, with `-std=c99 -pedantic -Wall -Wextra -Werror`, and on Windows with MSVC, with `/W3 /WX`.

The numbers of the native build are checked through the Python package, which is this C compiled for each platform: its tests compare every call with the WebAssembly build, as [above](#javascript-against-webassembly), and run the comparison with [Essential Macleod](#against-essential-macleod).

## What is not tested

- **Coated stacks at oblique incidence against a closed-form value of $R$ or $T$.** The oblique closed forms here are a bare interface near grazing, an air gap at its critical angle, and the reflection delays of an air gap below it. Stacks at oblique incidence are checked against Essential Macleod, Byrnes' implementation and the energy identity above.
- **tmmcore's own $R$ and $T$ under an absorbing incident medium against another implementation.** Essential Macleod defines them differently there, so the comparison with it checks $r$ and $t$ through the surface admittance, not that definition, and Byrnes' `tmm` asserts that $n_0 \sin\theta_0$ is real, so it refuses a complex $n_0$ at a real angle of incidence. The definition rests on the energy identity, on the absorptance from the net irradiance, and on the JavaScript ⇆ WebAssembly agreement.
- **The thickness derivatives of TOD through an opaque layer against finite differences.** There TOD carries more round-off than a difference resolves, so those derivatives rest on the JavaScript ⇆ WebAssembly agreement.
- **The phase kernel at an exact critical angle.** It returns `null` there, or `NaN` group delay and its derivatives at a substrate in s, rather than a limit; see the [API reference](api.md#phase-gd-gdd-and-tod).
- **Extreme parameter ranges.** Very large layer counts and indices far outside the optical range are exercised by neither the equivalence suite nor the cross-library comparison.

If tmmcore disagrees with something you trust, please open an issue with the inputs.
