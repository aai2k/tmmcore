<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://raw.githubusercontent.com/aai2k/tmmcore/main/docs/img/banner-on-dark.png">
  <img alt="tmmcore"
       src="https://raw.githubusercontent.com/aai2k/tmmcore/main/docs/img/banner-on-light.png"
       width="640">
</picture>

**Transfer-matrix method for multilayer thin-film optics, with exact analytic derivatives.**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)

</div>

Takes a stack of layers and returns reflectance, transmittance and absorptance
for absorbing and dispersive materials, at any angle of incidence, in s and p
polarization. Alongside the spectra it returns the exact thickness Jacobian, the
exact thickness Hessian, and the needle-insertion P-function, computed
analytically rather than by finite differences.

Ships as JavaScript, as C, and as a WebAssembly build of the C. The JavaScript
has no dependencies and works on import. WebAssembly is opt-in and roughly an
order of magnitude faster.

**[Documentation](https://aai2k.github.io/tmmcore/)**

## Install

```bash
npm install tmmcore
```

The `.wasm` is prebuilt and included, so no Emscripten toolchain is required.

## Use

```js
import { tmm } from 'tmmcore';

// A quarter-wave MgF2 layer on glass, at 550 nm, normal incidence.
const { R, T, A } = tmm(
    550,            // wavelength, nm
    0,              // angle of incidence, degrees from normal
    's',            // polarization: 's' or 'p'
    [1.0, 0],       // incident medium, ñ = [n, k]
    [1.52, 0],      // substrate
    [{ n: [1.38, 0], d: 550 / (4 * 1.38) }]   // quarter wave, thickness in nm
);

console.log(R);   // 0.012600790214630274
```

Layers run from the incident medium toward the substrate.

## Conventions

Mismatched conventions are the most common cause of two TMM codes disagreeing,
so check these first.

| | |
|---|---|
| Refractive index | ñ = n + i·k, with **k ≥ 0** for absorbing media |
| Time factor | exp(−iωt), so a wave exp(i(kz − ωt)) decays for k > 0 |
| Wavelength, thickness | nanometres |
| Angle | degrees from normal |
| Complex numbers | `[re, im]` pairs |
| Layer order | incident medium → substrate |

This is the complex conjugate of Macleod's convention. R, T and A are identical
under conjugation; phase-sensitive quantities are not.

## Verify it yourself

Two commands, neither needing anything but Node:

```bash
npm test        # the JavaScript and the C agree, to float64 round-off
npm run compare # and both agree with an independent implementation
```

The first drives both implementations with identical inputs across absorbing,
dispersive and oblique-incidence cases and compares every returned quantity.
64,416 comparisons, worst disagreement 4.4e-16.

The second checks them against [Steven Byrnes'
`tmm`](https://github.com/sbyrnes321/tmm), written independently in Python under
the same complex-index convention, so only the mathematics is under test. 12,352
values, worst disagreement 8.6e-14, which is float64 accumulation noise over a
forty-layer matrix product.

Timing comparisons against four other packages are in the
[documentation](https://aai2k.github.io/tmmcore/comparison/), and need a Python
environment.

## Documentation

- [Getting started](https://aai2k.github.io/tmmcore/getting-started/)
- [API reference](https://aai2k.github.io/tmmcore/api/)
- [Examples](https://aai2k.github.io/tmmcore/examples/)
- [Validation](https://aai2k.github.io/tmmcore/validation/)
- [Comparison with other TMM packages](https://aai2k.github.io/tmmcore/comparison/)

## Using the C directly

`src/tmm_kernel.c` is C99 with no dependencies beyond libm. Drop it into a
project and compile:

```bash
cc -std=c99 -O2 -c src/tmm_kernel.c
```

## Licence

[MIT](./LICENSE) © Andrey Achapovsky

Built for and used by [TFStudio](https://github.com/aai2k/TFStudio), an
open-source optical coating design application.
