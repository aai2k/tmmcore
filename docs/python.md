# Python

The Python package carries the same C kernel as the npm package, compiled for your platform and called through ctypes. The calls take the same stack and return the same results as those in the [API reference](api.md), under snake_case names, with numpy arrays where the JavaScript has lists. This page covers the Python names, types and shapes; the API reference describes each call.

## Install

```bash
pip install tmmcore
```

Wheels are built for Windows (x86-64), macOS (Intel and Apple silicon) and Linux (x86-64 and ARM64, glibc), for Python 3.9 and later. The only dependency is numpy.

## The calls

| JavaScript | Python |
|---|---|
| `rta`, `rtaSpectrum` | `rta`, `rta_spectrum` |
| `jacobian`, `jacobianSpectrum` | `jacobian`, `jacobian_spectrum` |
| `hessian`, `hessianSpectrum` | `hessian`, `hessian_spectrum` |
| `phase`, `phaseSpectrum` | `phase`, `phase_spectrum` |
| `phaseJacobian`, `phaseJacobianSpectrum` | `phase_jacobian`, `phase_jacobian_spectrum` |
| `jetConstant` … `jetIsFinite`, `wavelengthOmegaJet`, `omegaFromLambdaNm` | `jet_constant` … `jet_is_finite`, `wavelength_omega_jet`, `omega_from_lambda_nm` |

Each call takes the stack, then a wavelength in nm or, in a `_spectrum` call, an array of them, and `angle` and `pol` as keyword arguments:

```python
import numpy as np
import tmmcore

stack = {"incident": 1.0, "substrate": 1.52,
         "layers": [{"n": 2.35, "d": 58.5}, {"n": 1.46, "d": 94.2}]}

tmmcore.rta(stack, 550)["R"]                                              # natural light
tmmcore.rta_spectrum(stack, np.linspace(400, 800, 401), angle=45, pol="s")["R"]
```

The phase calls also take `sine` in place of `angle`: a jet, or in a `_spectrum` call an array of shape (nLam, 4), for [an angle fixed outside the incident medium](api.md#an-angle-fixed-outside-the-incident-medium).

Not in the Python package yet: the needle scan, the monitoring calls (monitoring curve, deposition spectra, growing stack), the low-level primitives, and the constants `C_NM_PER_FS` and `JET_ORDER`.

## The stack

A dict, the layers ordered from the incident medium toward the substrate:

```python
stack = {"incident": 1.0,
         "substrate": {"n": 1.52},
         "layers": [{"n": 2.35, "k": 0.0005, "d": 58.5},
                    {"n": 1.46, "d": 94.2}]}   # d in nm
```

- A medium is a number or a dict with `n` and `k`. Each `n` and `k` is a number, or in a `_spectrum` call an array with one value per wavelength. `k` may be left out and is then 0.
- `n` may also be complex, `n + ik`, with `k` left out: `{"n": 2.35 + 0.0005j, "d": 58.5}`.
- In the phase calls an index is a number, which is a constant index with no dispersion, or a jet (below). An array of plain indices per wavelength is refused there, since group delay computed from it would leave out the material's dispersion.

A wrong value raises `TypeError` with a message naming it, such as `layer 1 needs {n, d}`. An unknown keyword raises Python's own `TypeError`.

## Results

A dict. At one wavelength the values are numbers, each first derivative an array with one value per layer and each second derivative an (N, N) array; over a grid each gains a leading wavelength axis. `N`, where a result has it, is the number of layers. The R, T and A keys and the derivative keys are those of the JavaScript: `R`, `T`, `A`, `Rs` … `Ap`, `dRdd`, `d2Rdd` and the rest, and `N`. The phase keys are in snake_case:

| JavaScript | Python |
|---|---|
| `phaseRad`, `phaseDeg` | `phase_rad`, `phase_deg` |
| `gd`, `gdd`, `tod` | `gd`, `gdd`, `tod` |
| `magnitudeSquared` | `magnitude_squared` |
| `dPhaseDeg`, `dGd`, `dGdd`, `dTod`, `dLogMagnitudeSquared` | `d_phase_deg`, `d_gd`, `d_gdd`, `d_tod`, `d_log_magnitude_squared` |

Derivatives over a grid come as arrays indexed by wavelength and then layer, where the JavaScript returns them flat:

| | JavaScript | Python |
|---|---|---|
| `dRdd` of `jacobian_spectrum`, `d_gdd` of `phase_jacobian_spectrum` | nLam × N values, `[i * N + k]` | shape (nLam, N), `[i, k]` |
| `d2Rdd` of `hessian_spectrum` | nLam × N × N values, `[(i * N + j) * N + k]` | shape (nLam, N, N), `[i, j, k]` |

`phase` and `phase_jacobian` return `None` for a coefficient that is exactly zero, where the JavaScript returns `null`.

## Jets

A jet is a complex array of shape (4,), the Taylor terms $[f,\ f',\ f''/2!,\ f'''/3!]$ of a function of angular frequency, and jets over a grid are an array of shape (…, 4). Every jet call broadcasts over the leading axes as numpy does, so one call builds a material's jets for a whole grid. A Cauchy index:

```python
def cauchy(lam, a, b):
    """n(λ) = a + b/λ², λ in nm, as a jet in ω at each wavelength."""
    lam_jet = tmmcore.wavelength_omega_jet(lam, tmmcore.omega_from_lambda_nm(lam))
    inverse_square = tmmcore.jet_reciprocal(tmmcore.jet_multiply(lam_jet, lam_jet))
    return tmmcore.jet_add(tmmcore.jet_constant(a), tmmcore.jet_scale(inverse_square, b))
```

The 24-layer chirped mirror of the [group-delay example](examples.md#group-delay-against-a-closed-form), its layer thicknesses ramping along the stack:

```python
lam = np.array([700.0, 750.0, 800.0, 850.0, 900.0])
high, low = cauchy(lam, 2.21, 2.4e4), cauchy(lam, 1.44, 3.6e3)

layers = [{"n": low if i % 2 else high, "d": (132 if i % 2 else 85) * (1 + 0.011 * i)}
          for i in range(24)]
mirror = {"incident": 1.0, "substrate": {"n": cauchy(lam, 1.5, 4.2e3)}, "layers": layers}

r = tmmcore.phase_spectrum(mirror, lam, pol="s")["r"]
```

`r["magnitude_squared"]`, `r["gd"]`, `r["gdd"]` and `r["tod"]` are then R, GD, GDD and TOD at each wavelength, the table the [group-delay example](examples.md#group-delay-against-a-closed-form) prints.

## Agreement with the npm package

The package's tests run every call on the cases of the npm package's equivalence suite and compare each value with its WebAssembly build: R, T and A to within 1e-9, everything else to within 1e-7 of its size or 1e-12, whichever is larger. Both run the same C, so they differ only in the last digits, where the platform's math library rounds differently from the one compiled into the WebAssembly. These are the tolerances the npm package's tests hold its JavaScript kernels to against the same WebAssembly.

They also run the [comparison with Essential Macleod](validation.md#against-essential-macleod) on its committed results, to the same tolerances as the npm package.

## From source

The wheels build `src/c` with the platform's C compiler. To build and test a checkout, with Node for the reference results:

```bash
git clone https://github.com/aai2k/tmmcore
cd tmmcore
pip install ./python
node python/tests/reference.mjs
python -m unittest discover -s python/tests
```

`reference.mjs` runs the npm package on the equivalence cases and writes the results the tests compare with.
