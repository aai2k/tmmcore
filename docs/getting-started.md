# Getting started

## Install

```bash
npm install tmmcore
```

The WebAssembly binary ships prebuilt, so nothing is compiled at install time. For Python, `pip install tmmcore`; see [Python](python.md).

## Conventions

Check these before comparing results with another code.

| | |
|---|---|
| Refractive index | $\tilde n = n + ik$, with $k \ge 0$ for absorbing media |
| Time factor | $\exp(-i\omega t)$, so a wave $\exp(i(kz - \omega t))$ decays for $k > 0$ |
| p-polarized amplitudes | from the tilted admittance $n/\cos\theta$, so $r_p = r_s$ at normal incidence |
| Wavelength, thickness | nanometres |
| Angle of incidence | degrees from normal |
| Complex numbers | `[re, im]` pairs, in jets and the primitives |
| Layer order | incident medium first, substrate last |
| Absorbing incident medium | every angle follows from the real invariant $n_0 \sin\theta_0$, with $n_0$ the real part of the index; the incident wave decays along the normal only, and $A$ is what the layers absorb, so $R + T + A$ is not one |

This is the complex conjugate of the convention in Macleod's *Thin-Film Optical Filters* ($\tilde n = n - ik$, $\exp(+i\omega t)$). Reflectance, transmittance and absorptance are the same under conjugation. Phase is not: the phase calls report it in Macleod's convention, as [Units and signs](api.md#units-and-signs) explains.

The p amplitudes follow thin-film practice, as Macleod's do, rather than the Fresnel coefficients of Born and Wolf. Byrnes' `tmm` takes the Fresnel sign, and the ellipsometric convention flips the reflected p direction to the same effect: $r_p = -r_s$ at normal incidence. Their $r_p$ is the negative of tmmcore's at every angle, and a $\Delta$ from $r_p/r_s$ differs by 180°. Macleod amends $\Delta$ by $\pm\pi$ in reflection for this reason (§9.7.1, Eq. 9.36). The [cross-library check](validation.md#against-an-independent-implementation) compares both coefficients with Byrnes' at 0°, 30° and 60°.

!!! warning "k is positive for loss"

    Some codes and datasets use $\tilde n = n - ik$. A negative $k$ given to tmmcore describes a **gain** medium, which gives transmittance above one rather than an error.

## Describing a stack

A stack is the incident medium, the substrate and the layers between them, ordered from the incident medium toward the substrate. A layer is `{ n, k, d }`, with `k` left out when it is 0.

```js
import { rta } from 'tmmcore';

const stack = {
    incident: 1,
    substrate: 1.52,
    layers: [
        { n: 1.38, d: 92.4 },    // outermost, facing the incident medium
        { n: 2.20, d: 115.9 },
        { n: 1.70, d: 75.0 },    // innermost, against the substrate
    ],
};

const { R, T, A } = rta(stack, 550);
```

The incident medium and the substrate are semi-infinite and have no thickness: each is a number, or `{ n, k }` when it absorbs.

tmmcore has no material data: give it the index at the wavelengths you evaluate. Over a grid, any `n` or `k` can be a list with one value per wavelength.

```js
import { rtaSpectrum } from 'tmmcore';

const wavelengths = [450, 500, 550, 600, 650];
const nSiO2 = lam => 1.46 + 3000 / (lam * lam);   // a Cauchy formula, lam in nm

const layers = [{ n: wavelengths.map(nSiO2), d: 94.2 }];
const { R } = rtaSpectrum({ incident: 1, substrate: 1.52, layers }, wavelengths);
```

## Polarization

`pol` is `'s'`, `'p'`, or `'avg'`, the mean of the two, which is natural light. The R/T and derivative calls take `'avg'` when it is left out:

```js
rta(stack, 550, { angle: 45 }).R;              // natural light
rta(stack, 550, { angle: 45, pol: 's' }).R;
```

`rtaSpectrum` returns `Rs` and `Rp` beside `R`, so one call gives all three.

The phase calls take `'s'` or `'p'`, and `pol` must be given.

## Turning on WebAssembly

The JavaScript works on import. WebAssembly is opt-in: load the bytes, instantiate, and call the same names on the instance.

=== "Node"

    ```js
    import { readFileSync } from 'node:fs';
    import { createRequire } from 'node:module';
    import { instantiateTmmWasm } from 'tmmcore';

    const require = createRequire(import.meta.url);
    const kernel = await instantiateTmmWasm(
        readFileSync(require.resolve('tmmcore/tmm_kernel.wasm'))
    );

    const { R, T, A } = kernel.rta(stack, 550);
    ```

=== "Browser"

    ```js
    import { initTmmWasmFromUrl, getTmmWasm, rta } from 'tmmcore';

    // Resolves false, without throwing, if the file is missing or fails to load.
    const ok = await initTmmWasmFromUrl('/tmm_kernel.wasm');

    const { R } = ok ? getTmmWasm().rta(stack, 550) : rta(stack, 550);
    ```

    Your bundler must serve `tmm_kernel.wasm` as a static asset.

Instances are not shared between threads, since there is no shared memory. Each worker instantiates its own from the same bytes.

On the instance, a `Spectrum` call crosses from JavaScript into WebAssembly per grid rather than per wavelength. That is where the speed is.

## Next

- [API reference](api.md), every call and what it returns
- [Examples](examples.md), runnable scripts from a single layer to a needle scan
- [Python](python.md), the same calls from Python
- [Validation](validation.md), how to rerun the accuracy checks
