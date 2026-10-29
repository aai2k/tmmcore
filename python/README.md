# tmmcore

**Transfer-matrix method for multilayer thin-film optics, with exact analytic derivatives.**

[Documentation](https://aai2k.github.io/tmmcore/) · [Python](https://aai2k.github.io/tmmcore/python/) · [Source](https://github.com/aai2k/tmmcore)

Takes a stack of layers and returns reflectance, transmittance and absorptance for absorbing and dispersive materials, at any angle of incidence, in s and p polarization. Alongside them it returns the exact thickness Jacobian and Hessian, and phase, group delay, GDD and third-order dispersion with their thickness derivatives, all computed analytically rather than by finite differences.

This is the Python package of tmmcore. The wheel carries the same C kernel as the npm package, compiled for your platform and called through ctypes. It needs numpy and nothing else.

## Install

```bash
pip install tmmcore
```

Wheels are built for Windows (x86-64), macOS (Intel and Apple silicon) and Linux (x86-64 and ARM64, glibc), for Python 3.9 and later.

## Use

```python
import tmmcore

# A quarter-wave MgF2 layer on glass, at 550 nm, normal incidence.
stack = {"incident": 1.0, "substrate": 1.52,
         "layers": [{"n": 1.38, "d": 550 / (4 * 1.38)}]}   # d in nm

print(tmmcore.rta(stack, 550)["R"])   # 0.012600790214630274
```

Layers run from the incident medium toward the substrate. In `rta_spectrum` any n or k can be one value per wavelength. A quarter-half-quarter antireflection coating on BK7-like glass (Macleod, *Thin-Film Optical Filters* 5th ed., ch. 4), with Cauchy indices:

```python
import numpy as np
import tmmcore

def cauchy(a, b):
    return lambda lam: a + b / lam**2

low, high, mid, glass = cauchy(1.377, 2000), cauchy(2.170, 8000), cauchy(1.680, 5000), cauchy(1.504, 4200)
lam = np.arange(400.0, 701.0, 2.0)   # nm

# Quarter, half and quarter wave at 510 nm, the layer facing air first.
layers = [{"n": m(lam), "d": waves * 510 / m(510)} for m, waves in ((low, 0.25), (high, 0.5), (mid, 0.25))]
coated = tmmcore.rta_spectrum({"incident": 1.0, "substrate": {"n": glass(lam)}, "layers": layers}, lam)

print(f"{coated['R'].mean():.2%}")   # 0.25%, against 4.25% for bare glass
```

## Calls

| Call | Returns |
|---|---|
| `rta`, `rta_spectrum` | R, T and A |
| `jacobian`, `jacobian_spectrum` | R, T and A, and their derivatives in each layer's thickness |
| `hessian`, `hessian_spectrum` | the same, and the second derivatives |
| `phase`, `phase_spectrum` | phase, GD, GDD and TOD of r and t |
| `phase_jacobian`, `phase_jacobian_spectrum` | the same, and their thickness derivatives |
| `jet_*`, `wavelength_omega_jet`, `omega_from_lambda_nm` | Taylor-jet arithmetic, which builds the dispersive indices the phase calls take |

Each call takes the stack, then one wavelength or, in a `_spectrum` call, an array of them. `angle` and `pol` are keyword arguments. `angle` is the angle of incidence, 0 when left out, and in a `_spectrum` call may be one angle per wavelength. `pol` is `"s"`, `"p"` or `"avg"` (the default, the mean of the two) in the R/T and derivative calls, and `"s"` or `"p"` in the phase calls, where it must be given. Results are dicts. At one wavelength the values are numbers, each first derivative an array with one value per layer and each second derivative an (N, N) array; over a grid each gains a leading wavelength axis. The phase calls return one dict for r and one for t, and at one wavelength None for a coefficient that is exactly zero.

## Conventions

| | |
|---|---|
| Refractive index | ñ = n + ik, with **k ≥ 0** for absorbing media; n may instead be complex, n + ik |
| Time factor | exp(−iωt), so a wave exp(i(kz − ωt)) decays for k > 0 |
| p-polarized amplitudes | from the tilted admittance n/cos θ, so r_p = r_s at normal incidence |
| Wavelength, thickness | nanometres |
| Angle of incidence | degrees from normal |
| Phase | radians, in Macleod's sign convention; GD, GDD and TOD in fs, fs² and fs³ |
| Layer order | incident medium first, substrate last |

## Agreement with the npm package

The package's tests run every call on the cases of the npm package's equivalence suite and compare each value with its WebAssembly build: R, T and A to within 1e-9, everything else to within 1e-7 of its size or 1e-12, whichever is larger. Both run the same C, so they differ only in the last digits, where the platform's math library rounds differently from the one compiled into the WebAssembly.

## Licence

[MIT](https://github.com/aai2k/tmmcore/blob/main/LICENSE) © Andrey Achapovsky
