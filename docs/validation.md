# Validation

Everything below runs from a clone. The numbers quoted are the output on a developer machine, and yours should match to the last digit or two.

There are four independent levels of check, each answering a different question.

| Level | Question it answers |
|---|---|
| Closed-form oracles | Does it match the *equations*? |
| Finite differences | Are the derivatives derivatives of *what `tmm()` computes*? |
| JavaScript ⇆ WebAssembly | Do the two implementations match *each other*? |
| Cross-library | Does it match an *independent author's* implementation? |

None of these subsumes the others. Two implementations can agree perfectly and both be wrong; matching a closed form in one special case says nothing about the general one.

---

## Closed-form oracles

For a single quarter-wave layer at normal incidence, reflectance has an exact solution (Macleod §3.2):

$$R = \left(\frac{n_0 n_s - n_1^2}{n_0 n_s + n_1^2}\right)^2$$

```bash
node examples/01-single-layer.mjs
```

```
R from tmmcore           0.012600790214630274
R closed form            0.012600790214630288
difference               1.39e-17
```

Below double-precision epsilon ($2.2\times10^{-16}$).

### Group delay, GDD and TOD

A slab whose surrounding media share its own index has no interfaces, so nothing reflects and the transmission coefficient is a pure phase, $t = e^{i\delta}$ with $\delta = n(\omega)\,\omega d/c$. The three phase quantities then reduce exactly to bulk propagation:

$$\mathrm{GD} = \frac{d}{c}(n + \omega n'), \qquad
  \mathrm{GDD} = \frac{d}{c}(2n' + \omega n''), \qquad
  \mathrm{TOD} = \frac{d}{c}(3n'' + \omega n''')$$

A Cauchy index is exactly quadratic in $\omega$, which makes the right-hand sides elementary. This checks the phase kernel against something that is not another transfer-matrix calculation.

```bash
node examples/08-group-delay.mjs
```

```
matched slab, 1 µm of n(λ) = 1.45 + 3600/λ², at 800 nm

            tmmcore            closed form         relative
  GD       4.89296832143789     4.89296832143789   1.8e-16  fs
  GDD    0.0478126142151489   0.0478126142151501   2.6e-14  fs²
  TOD    0.0203063517881787   0.0203063517881746   2.0e-13  fs³
```

!!! warning "TOD loses precision on thick elements"

    TOD is a difference of terms of order $(\mathrm{GD}\cdot\omega)^3$, so the
    cancellation in it grows as the square of the thickness. Measured on the same
    slab:

    | Thickness | Group delay | Relative error on TOD |
    |---|---|---|
    | 1 µm | 4.9 fs | $2\times10^{-13}$ |
    | 10 µm | 48.9 fs | $4\times10^{-10}$ |
    | 100 µm | 489 fs | $4\times10^{-8}$ |
    | 1 mm | 4893 fs | $7\times10^{-7}$ |

    Coatings live at the top of that table and are unaffected. If you want the
    dispersion of millimetres of glass, differentiate the propagation phase
    directly instead of asking a transfer matrix for it: the closed forms above
    are the whole calculation and they do not cancel.

    GD and GDD do not suffer this; GD is at machine precision throughout.

### An absorbing incident medium

Light arriving from inside an absorbing medium, a cemented cube or an immersion liquid, is the one case where $R + T$ over lossless layers is not 1, and the size of the departure is fixed by the definitions. Every wave in the stack shares the real invariant $n_0 \sin\theta_0$, with $n_0$ the real part of the incident index (Macleod §10.2), and the transmittance carries the incident admittance the way Macleod's Eq. 2.83 does. What is left is the interference of the incident and reflected waves inside the absorbing medium:

$$1 - R - T = -2\,\frac{\mathrm{Im}\,\eta_0}{\mathrm{Re}\,\eta_0}\,\mathrm{Im}(r)$$

in tmmcore's sign convention, exact at every angle and polarization. On a bare interface $\mathrm{Im}(r)$ is itself of order $k_0$, and the excess reduces to Macleod's Eq. 2.84: into the conjugate admittance, $T = 1 + k_0^2/n_0^2$ and $R = k_0^2/n_0^2$.

```bash
node tests/absorbing_incident.mjs
```

```
worst gap in the energy identity: 1.56e-15
the former complex invariant broke it by 5.10e+1 per unit k0
WebAssembly comparisons: 924
PASS : absorbing incident medium
```

The second line is what the test guards against. Carrying the complex index into Snell's invariant makes the incident wave's amplitude vary along the interface, energy flows sideways inside lossless layers, and $R + T$ exceeds 1 by a further term linear in $k_0$ that grows with the angle and with the stack's resonance. Versions before 0.3.1 did that.

The same test holds $A$ to what the layers absorb: the net irradiance $\tfrac12\mathrm{Re}(BC^*)$ entering the front of the stack less $\tfrac12\mathrm{Re}(\eta_s)$ leaving it (Macleod Eq. 2.120), per unit incident irradiance, computed from $[B, C]$ rather than from $r$. It agrees to $10^{-14}$ on a stack with a silver layer, and $A$ stays below $10^{-14}$ on the lossless 21-layer stack. Up to 0.4.0, $A$ was $1 - R - T$ clamped at zero, which carried the interference term: 0.017 on that lossless stack, and $R + T + A = 1.047$ where the clamp held it at zero.

### Grazing incidence and the critical angle

Two angles where a formula that is right everywhere else loses everything. Near 90°, $\sqrt{1 - \sin^2\theta_0}$ keeps none of $\cos\theta_0$: up to 0.4.0, $1 - R$ of a bare interface was wrong by $7\times10^{-4}$ relative at 89.99999°, and $R$ was `NaN` in p at 90°. At a medium's critical angle to the last bit, $\cos\theta$ in it is exactly zero and the layer matrix is $0/0$, so $R$ and $T$ were `NaN` in both polarizations.

```bash
node tests/edge_cases.mjs
```

```
exact critical angle found at θ0 = 30.000000000000004°
WebAssembly comparisons: 772
PASS : edge cases
```

At grazing incidence the check is Fresnel's equations for a bare interface, $T = 4\eta_0\eta_s/(\eta_0 + \eta_s)^2$, with the exact cosine; tmmcore matches them to $10^{-12}$ relative in $T$ from 89.9° to 90°, and an air layer under air changes nothing. The critical angle is searched for among the doubles near 30°, from glass of index 2 into air, and $R$ and $T$ of an air gap there match the limit of the layer matrix in closed form, and the doubles either side of it, which the matrix reaches smoothly. A substrate at its critical angle reflects everything. The same file checks that a negative or NaN thickness leaves a layer out of every evaluator, with zero derivative in both Jacobians, and that $r_p = r_s$ at normal incidence.

---

## Derivatives against finite differences

The JavaScript and WebAssembly derivative kernels share their formulas, so agreeing with each other proves nothing about the formulas. This check holds every derivative to finite differences of `tmm()` and `tmmPhaseDispersion()`, which share nothing with the kernels but the layer matrix: the thickness Jacobian, every entry of the Hessian, the needle P-function at every gap and inside every layer, and the thickness derivatives of phase, GD, GDD and $\ln|r|^2$.

```bash
node tests/derivatives_fd.mjs
```

```
15668 comparisons against finite differences and the port, 9342 of them WebAssembly.
worst difference as a fraction of its tolerance: 2.43e-1, six clamped layers λ=550 p dGdd[3]
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

Differences use five-point central stencils in the layer thickness, with a step of 0.01 nm, or 0.25 nm for the phase quantities, where GDD carries enough cancellation to need the wider step. A needle cannot have negative thickness, so its differences are one-sided, third order. The worst case above is that GDD noise, at a quarter of its tolerance.

Up to 0.4.0 these kernels differentiated a different matrix from the one in the product past the clamp, and returned exactly zero at any wavelength where the partial products passed $10^{77}$. On the stacks above that gave derivatives of $10^{38}$ or `NaN` where the true value is zero, zero where it is $10^{-2}$, and needle gradients inside a metal an order of magnitude too large, or large where the true value is zero. The phase Jacobian returned `null` past the rescale threshold.

---

## JavaScript against WebAssembly

The C kernel is a line-by-line port of the JavaScript. They are driven with identical inputs across absorbing, dispersive and oblique-incidence cases, in both polarizations, and every returned quantity is compared: R, T, A, the thickness Jacobian, the thickness Hessian, the needle P-function at gaps and inside layers, and the phase quantities with their thickness derivatives, point by point and batched.

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

`npm test` runs this together with the other test files on this page; `npm run compare` runs the cross-library check below.

Agreement is not bit-exact by design. The only divergence is libm: the WebAssembly build uses musl's `sin`/`cos`/`exp`/`atan2`, the JavaScript engine uses its own, and they differ at roughly one unit in the last place. The observed disagreement sits seven orders of magnitude inside the tolerance.

The phase quantities amplify that noise, because every derivative order is a difference of nearly equal terms. Measured per quantity on dispersive stacks:

| Quantity | Worst JS ⇆ WASM difference |
|---|---|
| $\lvert r\rvert^2$ | 4 ulp |
| phase | 2 ulp |
| GD | 36 ulp |
| GDD | 325 ulp |
| TOD | 293 ulp |

About a decade per order of differentiation, which is what cancellation costs. Even the worst case is $10^{-11}$ relative, orders below the precision of any measured $n$ and $k$.

One stack is held to a looser bound: six 20 µm opaque layers, whose product passes the rescale threshold. There the $q$-th frequency order is a difference of terms of order $\mathrm{GD}^q$, with GD the group delay through the whole stack, about 800 fs, so two roundings of it agree to about $64\,\varepsilon\,\mathrm{GD}^q$ and no better. The derivatives with respect to the opaque layers are zero to exactly that noise.

The test skips cleanly if `tmm_kernel.wasm` has not been built.

### With the heap exhausted

In wasm32 a null pointer is address 0 of linear memory, so a kernel that works through a failed allocation does not trap. It writes over whatever lives there and returns numbers. This test takes up the whole heap of a fresh instance, block by halving block, then requires every kernel that allocates working state to write `NaN` to all its outputs, and the same instance to compute correctly once the heap is returned.

```bash
node tests/wasm_allocation.mjs
```

```
heap exhausted at 2048 MiB of linear memory
258 outputs of 8 kernels checked NaN without memory
PASS : allocation failure
```

Up to 0.4.0 the derivative and phase kernels did not check their allocations, and returned numbers here.

---

## The C, built natively

The kernel is C99 with no dependencies beyond libm, so it can be checked outside WebAssembly entirely:

```bash
cc -std=c99 -pedantic -Wall -Wextra -O2 -c src/tmm_kernel.c
```

Verified warning-free under **GCC 16.1.0**, with the native build reproducing the closed-form quarter-wave result to $1.4\times10^{-17}$, and s and p identical at normal incidence.

!!! success "Portable across three native toolchains"

    GitHub Actions compiles the same source on Linux with GCC, on macOS with
    Clang, and on Windows with MSVC. GCC and Clang use `-std=c99 -pedantic
    -Wall -Wextra -Werror`; MSVC uses `/W3 /WX`. All three builds pass.

---

## Against an independent implementation

The same inputs fed to [Steven Byrnes' `tmm`](https://github.com/sbyrnes321/tmm): MIT, peer-reviewed[^byrnes], pure Python, sharing no code and no author with tmmcore.

It uses the **same** $\tilde n = n + ik$ convention and time factor as tmmcore, so values transfer verbatim with no conjugation. There is no material-data confound: the wavelength grid, the complex indices and the thicknesses are precomputed into a shared file that every implementation reads. Only the mathematics differs.

The cases run at normal incidence, and at 30° and 60°, where s and p part and the p admittance $n/\cos\theta$ is under test. R and T come from `tmm()`, and the complex reflection coefficient from the phase kernel. $r_s$ is compared as it stands and $r_p$ with its sign reversed, since Byrnes takes the Fresnel sign for p (see [Conventions](getting-started.md#conventions)); the comparison is what holds that statement about the sign.

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

The worst case is the forty-layer stack, where round-off accumulates through the longest matrix product. The oblique cases are held to the same $10^{-12}$ as the rest. Regenerating the reference file rather than trusting the committed one takes two `pip install`s; see [`benchmarks/README.md`](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md).

See [Comparison with other packages](comparison.md) for the full tables, including three further libraries.

[^byrnes]: S. J. Byrnes, *Multilayer optical calculations*, arXiv:1603.02720.

---

## What is not tested

- **Coated stacks at oblique incidence against a closed-form value.** The
  oblique closed forms here are a bare interface near grazing and an air gap at
  its critical angle. Stacks at oblique incidence are checked against Byrnes'
  implementation and against the energy identity above, not against an
  analytic $R$ or $T$.
- **An absorbing incident medium against another implementation.** Byrnes'
  `tmm` asserts that $n_0 \sin\theta_0$ is real, so it refuses a complex $n_0$
  at a real angle of incidence. That case rests on the identity, on the
  absorptance from the net irradiance, and on the JavaScript ⇆ WebAssembly
  agreement.
- **The phase kernel at an exact critical angle.** It returns `null` there, or
  `NaN` group delay and its derivatives at a substrate in s, rather than a
  limit; see the [API reference](api.md#phase-dispersion).
- **Extreme parameter ranges.** Very large layer counts and indices far outside
  the optical range are exercised by neither the equivalence suite nor the
  cross-library comparison.

If you hit a case where tmmcore disagrees with something you trust, that is a useful bug report. Please open an issue with the inputs.
