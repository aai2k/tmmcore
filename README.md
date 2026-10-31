<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://raw.githubusercontent.com/aai2k/tmmcore/main/docs/img/banner-on-dark.png">
  <img alt="tmmcore"
       src="https://raw.githubusercontent.com/aai2k/tmmcore/main/docs/img/banner-on-light.png"
       width="320">
</picture>

**Transfer-matrix method for multilayer thin-film optics, with exact analytic derivatives.**

[![npm](https://img.shields.io/npm/v/tmmcore)](https://www.npmjs.com/package/tmmcore) [![PyPI](https://img.shields.io/pypi/v/tmmcore)](https://pypi.org/project/tmmcore/) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)

**[Documentation](https://aai2k.github.io/tmmcore/)** · **[Getting started](https://aai2k.github.io/tmmcore/getting-started/)** · **[API](https://aai2k.github.io/tmmcore/api/)** · **[Python](https://aai2k.github.io/tmmcore/python/)** · **[Validation](https://aai2k.github.io/tmmcore/validation/)**

</div>

Takes a stack of layers and returns reflectance, transmittance and absorptance for absorbing and dispersive materials, at any angle of incidence, in s and p polarization. Alongside the spectra it returns the exact thickness Jacobian, the exact thickness Hessian, and the needle-insertion P-function, computed analytically rather than by finite differences.

It also computes phase, group delay, GDD and third-order dispersion, by carrying the same matrix in third-order Taylor arithmetic. They are exact frequency derivatives, so they do not depend on the wavelength grid and there is no finite-difference step to tune, and their thickness derivatives come from the same call.

Ships as JavaScript, as C, as a WebAssembly build of the C, and as a Python package on the same C. The JavaScript has no dependencies and works on import. WebAssembly is opt-in and, over a wavelength grid, about an order of magnitude faster than the JavaScript.

## Install

```bash
npm install tmmcore
```

The WebAssembly binary ships prebuilt, so nothing is compiled at install time. For Python, `pip install tmmcore`; see [Python](https://aai2k.github.io/tmmcore/python/).

## Use

```js
import { rta } from 'tmmcore';

// A quarter-wave MgF2 layer on glass, at 550 nm, normal incidence.
const stack = {
    incident: 1.0,
    substrate: 1.52,
    layers: [{ n: 1.38, d: 550 / (4 * 1.38) }],   // thickness in nm
};

console.log(rta(stack, 550).R);   // 0.012600790214630274
```

Layers run from the incident medium toward the substrate. In `rtaSpectrum` any n or k can be one value per wavelength. A quarter-half-quarter antireflection coating on BK7-like glass (Macleod, *Thin-Film Optical Filters* 5th ed., ch. 4), with Cauchy indices:

```js
import { rtaSpectrum } from 'tmmcore';

const cauchy = (a, b) => lam => a + b / lam ** 2;
const low = cauchy(1.377, 2000), high = cauchy(2.170, 8000), mid = cauchy(1.680, 5000), glass = cauchy(1.504, 4200);
const lams = Array.from({ length: 151 }, (_, i) => 400 + 2 * i);   // nm

// Quarter, half and quarter wave at 510 nm, the layer facing air first.
const layers = [[low, 0.25], [high, 0.5], [mid, 0.25]]
    .map(([m, waves]) => ({ n: lams.map(m), d: waves * 510 / m(510) }));
const coated = rtaSpectrum({ incident: 1, substrate: { n: lams.map(glass) }, layers }, lams);

const mean = values => values.reduce((sum, v) => sum + v, 0) / values.length;
console.log(`${(100 * mean(coated.R)).toFixed(2)} %`);   // 0.25 %, against 4.25 % for bare glass
```

The R/T, derivative, needle and phase calls all take a stack, a wavelength or a list of them, and `{ angle, pol }`: the angle of incidence in degrees, 0 when left out, and `'s'`, `'p'` or `'avg'`, the mean of the two. `'avg'` is the default in the R/T, derivative and needle calls; the phase calls take `'s'` or `'p'`, and may leave it out at normal incidence. The needle calls also take `candidates`, the materials they try inserting. The [API reference](https://aai2k.github.io/tmmcore/api/) lists them.

## Conventions

Check these before comparing results with another code.

| | |
|---|---|
| Refractive index | ñ = n + i·k, with **k ≥ 0** for absorbing media |
| Time factor | exp(−iωt), so a wave exp(i(kz − ωt)) decays for k > 0 |
| p-polarized amplitudes | from the tilted admittance n/cos θ, so r_p = r_s at normal incidence, the opposite sign to Born and Wolf's Fresnel r_p |
| Wavelength, thickness | nanometres |
| Angle of incidence | degrees from normal |
| Complex numbers | `[re, im]` pairs, in jets and the primitives |
| Layer order | incident medium first, substrate last |

This is the complex conjugate of Macleod's convention. R, T and A are the same under conjugation; the phase calls report phase in Macleod's convention, as [Units and signs](https://aai2k.github.io/tmmcore/api/#units-and-signs) explains. [Getting started](https://aai2k.github.io/tmmcore/getting-started/#conventions) adds the conventions for an absorbing incident medium.

## Verify it yourself

Four commands, none needing anything but Node:

| Command | Checks | Worst difference |
|---|---|---|
| `node examples/01-single-layer.mjs` | R of a quarter-wave layer against its closed form (Macleod Eq. 4.5) | 1.4e-17 |
| `npm test` | the JavaScript against the C built to WebAssembly, value by value; every derivative against finite differences; grazing incidence and the exact critical angle against closed forms; and the Essential Macleod comparison below | 4.4e-16 in R, T, A and their derivatives; 1.1e-11 relative in phase, GD, GDD and TOD; JavaScript against WebAssembly |
| `npm run compare` | R, T and the complex r against [Steven Byrnes' `tmm`](https://github.com/sbyrnes321/tmm), written independently in Python, at 0°, 30° and 60° | 8.6e-14 |
| `node tests/macleod.mjs` | ten designs against Essential Macleod's results, committed here: R, T and A, the phase of r and t with their GD, GDD and TOD, Psi and Delta, and monitoring curves | 2.1e-14 in R, T and A, from either side; 3.4e-12° in phase |

The first checks the equations. The JavaScript ⇆ WebAssembly part of the second compares two implementations by the same author, so it catches porting bugs and nothing more. The third compares with another author's code under the same complex-index convention; its worst case is a forty-layer stack, where round-off builds up through the longest matrix product. The fourth compares with a commercial design program, on designs of constant-index materials made for the check, so no material data enters.

[Validation](https://aai2k.github.io/tmmcore/validation/) sets out what each check does and does not establish, and lists what is not tested. [Comparison](https://aai2k.github.io/tmmcore/comparison/) puts tmmcore beside four other TMM packages for accuracy and speed.

## Documentation

- [Getting started](https://aai2k.github.io/tmmcore/getting-started/)
- [API reference](https://aai2k.github.io/tmmcore/api/)
- [Python](https://aai2k.github.io/tmmcore/python/)
- [Examples](https://aai2k.github.io/tmmcore/examples/)
- [Validation](https://aai2k.github.io/tmmcore/validation/)
- [Comparison with other TMM packages](https://aai2k.github.io/tmmcore/comparison/)
- [Changelog](./CHANGELOG.md)

## Using the C directly

The C in `src/c` is C99 with no dependencies beyond libm. Copy the folder into a project, include `tmmcore.h`, which declares the calls, and compile `tmm_kernel.c`, which includes the other sources:

```bash
cc -std=c99 -O2 -c src/c/tmm_kernel.c
```

[C](https://aai2k.github.io/tmmcore/api/#c) in the API reference describes the arrays the calls take and write.

## Licence

[MIT](./LICENSE) © Andrey Achapovsky

Built for and used by [TFStudio](https://github.com/aai2k/TFStudio), an open-source optical coating design application.
