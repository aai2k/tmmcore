# Validation

Everything below runs from a clone. The numbers quoted are the output on a
developer machine, and yours should match to the last digit or two.

There are three independent levels of check, each answering a different
question.

| Level | Question it answers |
|---|---|
| Closed-form oracles | Does it match the *equations*? |
| JavaScript ⇆ WebAssembly | Do the two implementations match *each other*? |
| Cross-library | Does it match an *independent author's* implementation? |

None of these subsumes the others. Two implementations can agree perfectly and
both be wrong; matching a closed form in one special case says nothing about the
general one.

---

## Closed-form oracles

For a single quarter-wave layer at normal incidence, reflectance has an exact
solution (Macleod §3.2):

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

---

## JavaScript against WebAssembly

The C kernel is a line-by-line port of the JavaScript. They are driven with
identical inputs across absorbing, dispersive and oblique-incidence cases, in
both polarizations, and every returned quantity is compared: R, T, A, the
thickness Jacobian, the thickness Hessian, and the needle P-function.

```bash
npm test
```

```
64416 comparisons across 4 stacks, 4 wavelengths, 4 angles, s and p.
worst |Δ| on R/T/A     : 4.44e-16  (tolerance 1e-9)
worst |Δ| on derivatives: 5.55e-17  (tolerance 1e-12 abs / 1e-7 rel)

PASS — JavaScript and WebAssembly agree.
```

Agreement is not bit-exact by design. The only divergence is libm: the
WebAssembly build uses musl's `sin`/`cos`/`exp`/`atan2`, the JavaScript engine
uses its own, and they differ at roughly one unit in the last place. The
observed disagreement sits seven orders of magnitude inside the tolerance.

The test skips cleanly if `tmm_kernel.wasm` has not been built.

---

## The C, built natively

The kernel is C99 with no dependencies beyond libm, so it can be checked outside
WebAssembly entirely:

```bash
cc -std=c99 -pedantic -Wall -Wextra -O2 -c src/tmm_kernel.c
```

Verified warning-free under **GCC 16.1.0**, with the native build reproducing
the closed-form quarter-wave result to $1.4\times10^{-17}$, `R + T + A = 1`
exactly, and s and p identical at normal incidence.

!!! success "Portable across three native toolchains"

    GitHub Actions compiles the same source on Linux with GCC, on macOS with
    Clang, and on Windows with MSVC. GCC and Clang use `-std=c99 -pedantic
    -Wall -Wextra -Werror`; MSVC uses `/W3 /WX`. All three builds pass.

---

## Against an independent implementation

The same inputs fed to [Steven Byrnes'
`tmm`](https://github.com/sbyrnes321/tmm): MIT, peer-reviewed[^byrnes], pure
Python, sharing no code and no author with tmmcore.

It uses the **same** $\tilde n = n + ik$ convention as tmmcore, so values
transfer verbatim with no conjugation. There is no material-data confound: the
wavelength grid, the complex indices and the thicknesses are precomputed into a
shared file that every implementation reads. Only the mathematics differs.

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

The worst case is the forty-layer stack, where round-off accumulates through the
longest matrix product. Regenerating the reference file rather than trusting the
committed one takes two `pip install`s; see
[`benchmarks/README.md`](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md).

See [Comparison with other packages](comparison.md) for the full tables,
including three further libraries.

[^byrnes]: S. J. Byrnes, *Multilayer optical calculations*, arXiv:1603.02720.

---

## What is not tested

- **Non-normal incidence against a closed form.** The oblique cases are checked
  against other implementations, not against an analytic result.
- **Extreme parameter ranges.** Very large layer counts, indices far outside the
  optical range, and grazing incidence are exercised by neither the equivalence
  suite nor the cross-library comparison.

If you hit a case where tmmcore disagrees with something you trust, that is a
useful bug report. Please open an issue with the inputs.
