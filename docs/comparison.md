# Comparison with other packages

The timings below were measured with tmmcore 0.1.0 on an AMD Ryzen 7 260 (8 cores, 16 threads), Node 24.18, Python 3.14.6 and Windows 11, and have not been rerun since. The harness is in `benchmarks/`, and [its README](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md) gives the method and its limits. Read the [caveats](#caveats) with the speed table.

## The packages

| Package | Language | Notes |
|---|---|---|
| [`tmm`](https://github.com/sbyrnes321/tmm) (Byrnes) | Python and NumPy | One wavelength per call. |
| [`tmm_faster`](https://github.com/clembron/tmm_faster) | C++ core | float64, multithreaded by default. |
| [`tmm_fast`](https://github.com/MLResearchAtOSRAM/tmm_fast) | PyTorch | Batched, GPU-capable, autograd. |
| [`tmmax`](https://github.com/bahremsd/tmmax) | JAX | JIT and vmap. Bundles its own material database. |
| **tmmcore** | JavaScript, C, WebAssembly | This package. |

## Method

Every implementation reads the same precomputed file: one wavelength grid, one set of complex indices, one set of thicknesses. Dispersion and material lookup sit **outside** the timed region, so only the transfer-matrix calculation is timed.

One spectrum is every wavelength of the grid in both polarizations.

The stacks have 4, 7, 21 and 40 layers; the 7-layer stack contains absorbing silver. The grids are 400-1100 nm at 10 nm (71 points) and at 1 nm (701 points).

## Accuracy

Byrnes' outputs are committed, so this needs only Node:

```bash
npm run compare
```

Largest absolute difference in R and T from Byrnes, over the whole grid, both polarizations:

| Case | `tmm_faster` | tmmcore (JS and WASM) |
|---|--:|--:|
| 4-layer, 71 pts | 2.9e-16 | 2.1e-15 |
| 21-layer, 71 pts | 6.7e-15 | 1.2e-14 |
| 7-layer with Ag, 71 pts | 3.2e-15 | 1.4e-15 |
| 40-layer, 71 pts | 1.8e-14 | 1.7e-14 |
| 4-layer, 701 pts | 5.8e-16 | 2.6e-15 |
| 21-layer, 701 pts | 1.3e-14 | 2.4e-14 |
| 7-layer with Ag, 701 pts | 3.6e-15 | 3.1e-15 |
| 40-layer, 701 pts | 1.8e-14 | 8.6e-14 |

Every case is lit from air. An absorbing incident medium has no cross-library reference here: Byrnes' `tmm` asserts that $n_0 \sin\theta_0$ is real and refuses a complex $n_0$ at a real angle. tmmcore evaluates that case from the real invariant $\mathrm{Re}(n_0)\sin\theta_0$ and checks it against an exact energy identity instead; see [Validation](validation.md#an-absorbing-incident-medium).

tmmcore and `tmm_faster` agree with Byrnes to float64 round-off. `tmm_fast` and `tmmax` agree only to about $10^{-7}$ whatever dtype they are given, the size of single-precision round-off; `tmmax` runs on JAX, which computes in single precision unless `jax_enable_x64` is set. So neither can serve as an accuracy reference.

## Speed

Single-threaded, best of N runs, one spectrum, in milliseconds. Lower is better.

| Case | byrnes | tmm_fast | tmmax | tmm_faster | tmmcore JS | **tmmcore WASM** |
|---|--:|--:|--:|--:|--:|--:|
| 4-layer, 71 pts | 7.760 | 7.418 | 11.56 | 0.083 | 0.521 | **0.039** |
| 21-layer, 71 pts | 25.16 | 8.782 | 11.80 | 0.343 | 2.577 | **0.169** |
| 7-layer Ag, 71 pts | 10.82 | 7.779 | 11.44 | 0.143 | 0.883 | **0.064** |
| 40-layer, 71 pts | 43.56 | 10.26 | 12.21 | 0.840 | 4.974 | **0.313** |
| 4-layer, 701 pts | 78.79 | 30.78 | 11.52 | 0.800 | 5.850 | **0.367** |
| 21-layer, 701 pts | 248.6 | 38.00 | 21.46 | 4.489 | 27.10 | **1.723** |
| 7-layer Ag, 701 pts | 109.5 | 33.06 | 12.85 | 1.466 | 9.591 | **0.619** |
| 40-layer, 701 pts | 434.3 | 47.27 | 32.62 | 8.366 | 52.05 | **3.182** |

## Caveats

!!! warning "Threading"

    `tmm_faster` is multithreaded by default, at a cpu/wall ratio of 13 to 15 on this 16-thread machine, and in that configuration it is faster than tmmcore on large problems: 1.27 ms against 3.18 ms on the 40-layer, 701-point case. The table above pins it to one thread, where tmmcore's WebAssembly is 2.0 to 2.7 times faster.

!!! warning "Batching"

    The vectorized libraries exist to evaluate many stacks at once. tmmcore has no batch call and loops. Single-threaded on 128 different 40-layer stacks it takes 45.9 ms, against 82.9 ms for `tmm_faster` and 209.8 ms for `tmm_fast`. On a GPU with a large enough batch, `tmm_fast` and `tmmax` would pull ahead, and **that case was not measured**.

From the harness notes:

- `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering`, more than R and T, and its API has no cheaper call.
- Byrnes' `tmm` must be given Python `complex`, not NumPy scalars; `coh_tmm` is much slower on `np.complex128` elements, which would have favoured every other package.
- `tmmax` reads materials from its own database rather than taking n and k, so its timing includes an interpolation the others do not do.
- A control run puts the Python driver loop under 1% of the total time, so the Byrnes figures time the library rather than the harness.

## What the others do that tmmcore does not

- `tmm` (Byrnes) handles incoherent layers, ellipsometric parameters and absorption per layer.
- `tmm_fast` and `tmmax` differentiate through the whole calculation, not only in the thicknesses, and run on GPUs.
- `tmmax` ships a 35-material database, so it computes without being given n and k.

The [examples](examples.md) build an incoherent back side and ellipsometric $\Psi$ and $\Delta$ on tmmcore's calls.
