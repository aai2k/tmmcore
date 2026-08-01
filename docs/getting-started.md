# Getting started

## Install

```bash
npm install tmmcore
```

The WebAssembly binary ships prebuilt. Nothing is compiled at install time and
no toolchain is required.

## Conventions

Mismatched conventions are the single most common reason two transfer-matrix
codes disagree, and the disagreement looks like a bug in one of them until
someone checks.

| | |
|---|---|
| Refractive index | $\tilde n = n + ik$, with $k \ge 0$ for absorbing media |
| Time factor | $\exp(-i\omega t)$, so a wave $\exp(i(kz - \omega t))$ decays for $k > 0$ |
| Wavelength | nanometres |
| Thickness | nanometres |
| Angle of incidence | degrees from normal |
| Complex numbers | `[re, im]` pairs |
| Layer order | incident medium first, substrate last |

This is the complex conjugate of the convention in Macleod's *Thin-Film Optical
Filters* ($\tilde n = n - ik$, $\exp(+i\omega t)$). Reflectance, transmittance
and absorptance are identical under conjugation. Phase-sensitive quantities are
not, so mind the sign if you compute phase from these results.

!!! warning "k is positive for loss"

    Some codes and datasets use $\tilde n = n - ik$. Feeding a negative $k$ to
    tmmcore describes a **gain** medium, which will produce transmittance above
    one rather than an error.

## Describing a stack

A layer is `{ n: [re, im], d }`. A stack is an array of them, ordered from the
incident medium toward the substrate.

```js
const layers = [
    { n: [1.38, 0],      d: 92.4 },    // outermost, facing the incident medium
    { n: [2.20, 0],      d: 115.9 },
    { n: [1.70, 0],      d: 75.0 },    // innermost, against the substrate
];
```

The incident medium and substrate are passed separately, as bare `[re, im]`
pairs, since they are semi-infinite and have no thickness.

Dispersion is your responsibility: tmmcore takes the index **at the wavelength
you are evaluating**. Sweeping a spectrum means recomputing each layer's index
per wavelength.

```js
const nSiO2 = lam => [1.46 + 3000 / (lam * lam), 0];

for (const lam of wavelengths) {
    const layers = design.map(({ mat, d }) => ({ n: mat(lam), d }));
    const { R } = tmm(lam, 0, 's', [1, 0], substrate(lam), layers);
}
```

Keeping material data out of the kernel is what makes the comparison against
other packages meaningful. The same numbers go into every implementation, so
only the mathematics differs.

## Unpolarized light

There is no `'unpolarized'` option. Compute both and average:

```js
const s = tmm(lam, theta, 's', n0, ns, layers);
const p = tmm(lam, theta, 'p', n0, ns, layers);
const R = (s.R + p.R) / 2;
```

At normal incidence the two are identical, so this costs nothing to write and
stays correct when you change the angle.

## Turning on WebAssembly

The JavaScript path works immediately. WebAssembly is opt-in: load the bytes,
instantiate, and call methods on the instance.

=== "Node"

    ```js
    import { readFileSync } from 'node:fs';
    import { createRequire } from 'node:module';
    import { instantiateTmmWasm } from 'tmmcore';

    const require = createRequire(import.meta.url);
    const kernel = await instantiateTmmWasm(
        readFileSync(require.resolve('tmmcore/tmm_kernel.wasm'))
    );

    // polarization is 0 for s, 1 for p
    const { R, T, A } = kernel.tmmOne(550, 0, 0, [1, 0], [1.52, 0], layers);
    ```

=== "Browser"

    ```js
    import { initTmmWasmFromUrl, getTmmWasm, tmm } from 'tmmcore';

    // Resolves false if the asset is missing or instantiation fails; it does
    // not throw, so that a deployment problem degrades to the JS path instead
    // of breaking the page.
    const ok = await initTmmWasmFromUrl('/tmm_kernel.wasm');

    const { R } = ok
        ? getTmmWasm().tmmOne(550, 0, 0, [1, 0], [1.52, 0], layers)
        : tmm(550, 0, 's', [1, 0], [1.52, 0], layers);
    ```

    Your bundler must serve `tmm_kernel.wasm` as a static asset.

Instances are not shared between threads, since there is no shared memory. Each
worker instantiates its own from the same bytes.

For a whole spectrum, `tmmSpectrum` evaluates every wavelength and both
polarizations in one call. That is where the speed advantage comes from: it
amortizes the JavaScript-to-WebAssembly boundary crossing.

## Next

- [API reference](api.md), every function and what it returns
- [Examples](examples.md), runnable scripts including optimization with the analytic Jacobian
- [Validation](validation.md), how to reproduce the accuracy claims
