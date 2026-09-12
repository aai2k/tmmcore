# Validation

Everything below runs from a clone. The numbers quoted are the output on a developer machine, and yours should match to the last digit or two.

There are three independent levels of check, each answering a different question.

| Level | Question it answers |
|---|---|
| Closed-form oracles | Does it match the *equations*? |
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
worst gap in the energy identity: 1.13e-15
the former complex invariant broke it by 5.10e+1 per unit k0
WebAssembly comparisons: 924
PASS : absorbing incident medium
```

The second line is what the test guards against. Carrying the complex index into Snell's invariant makes the incident wave's amplitude vary along the interface, energy flows sideways inside lossless layers, and $R + T$ exceeds 1 by a further term linear in $k_0$ that grows with the angle and with the stack's resonance. Versions before 0.3.1 did that.

---

## JavaScript against WebAssembly

The C kernel is a line-by-line port of the JavaScript. They are driven with identical inputs across absorbing, dispersive and oblique-incidence cases, in both polarizations, and every returned quantity is compared: R, T, A, the thickness Jacobian, the thickness Hessian, the needle P-function, and the phase quantities with their thickness derivatives, point by point and batched.

```bash
npm test
```

```
64416 comparisons across 4 stacks, 4 wavelengths, 4 angles, s and p.
worst |Δ| on R/T/A     : 4.44e-16  (tolerance 1e-9)
worst |Δ| on derivatives: 5.55e-17  (tolerance 1e-12 abs / 1e-7 rel)

15904 phase comparisons (phase, GD, GDD, TOD and their thickness derivatives,
point and batched).
worst |Δ| on phase quantities: 5.52e-6
worst relative, above the 1e-12 floor: 1.05e-11

PASS — JavaScript and WebAssembly agree.
```

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

The test skips cleanly if `tmm_kernel.wasm` has not been built.

---

## The C, built natively

The kernel is C99 with no dependencies beyond libm, so it can be checked outside WebAssembly entirely:

```bash
cc -std=c99 -pedantic -Wall -Wextra -O2 -c src/tmm_kernel.c
```

Verified warning-free under **GCC 16.1.0**, with the native build reproducing the closed-form quarter-wave result to $1.4\times10^{-17}$, `R + T + A = 1` exactly, and s and p identical at normal incidence.

!!! success "Portable across three native toolchains"

    GitHub Actions compiles the same source on Linux with GCC, on macOS with
    Clang, and on Windows with MSVC. GCC and Clang use `-std=c99 -pedantic
    -Wall -Wextra -Werror`; MSVC uses `/W3 /WX`. All three builds pass.

---

## Against an independent implementation

The same inputs fed to [Steven Byrnes' `tmm`](https://github.com/sbyrnes321/tmm): MIT, peer-reviewed[^byrnes], pure Python, sharing no code and no author with tmmcore.

It uses the **same** $\tilde n = n + ik$ convention as tmmcore, so values transfer verbatim with no conjugation. There is no material-data confound: the wavelength grid, the complex indices and the thicknesses are precomputed into a shared file that every implementation reads. Only the mathematics differs.

Byrnes' outputs are committed, so this needs no Python:

```bash
npm run compare
```

```
case           layers  points   max |Δ| JS  max |Δ| WASM
--------------------------------------------------------
AR4/g71             4      71      2.1e-15       2.1e-15
HR21/g71           21      71      1.2e-14       1.2e-14
AG7/g71             7      71      1.4e-15       1.6e-15
BIG40/g71          40      71      1.7e-14       1.7e-14
AR4/g701            4     701      2.6e-15       2.6e-15
HR21/g701          21     701      2.4e-14       2.4e-14
AG7/g701            7     701      3.1e-15       3.1e-15
BIG40/g701         40     701      8.6e-14       8.6e-14

12352 values compared across 8 cases, both polarizations.
Worst disagreement with an independently written implementation: 8.6e-14
```

The worst case is the forty-layer stack, where round-off accumulates through the longest matrix product. Regenerating the reference file rather than trusting the committed one takes two `pip install`s; see [`benchmarks/README.md`](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md).

See [Comparison with other packages](comparison.md) for the full tables, including three further libraries.

[^byrnes]: S. J. Byrnes, *Multilayer optical calculations*, arXiv:1603.02720.

---

## What is not tested

- **Non-normal incidence against a closed-form value.** The oblique cases are
  checked against other implementations and against the energy identity above,
  not against an analytic $R$ or $T$.
- **An absorbing incident medium against another implementation.** Byrnes'
  `tmm` asserts that $n_0 \sin\theta_0$ is real, so it refuses a complex $n_0$
  at a real angle of incidence. That case rests on the identity and on the
  JavaScript ⇆ WebAssembly agreement alone.
- **Extreme parameter ranges.** Very large layer counts, indices far outside the
  optical range, and grazing incidence are exercised by neither the equivalence
  suite nor the cross-library comparison.

If you hit a case where tmmcore disagrees with something you trust, that is a useful bug report. Please open an issue with the inputs.
