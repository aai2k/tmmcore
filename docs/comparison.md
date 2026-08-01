# Comparison with other packages

Numbers below were measured on an AMD Ryzen 7 260 (8C/16T), Node 24.18, Python
3.14.6, Windows 11. The harness lives in `benchmarks/`; the method and its known
limitations are described there and summarised here.

Two of the [caveats](#caveats) materially change what the figures mean.

## The packages

| Package | Language | Notes |
|---|---|---|
| [`tmm`](https://github.com/sbyrnes321/tmm) (Byrnes) | Python + NumPy | The de-facto standard. Scalar in wavelength. Cited by most machine-learning-for-thin-films work. |
| [`tmm_faster`](https://github.com/clembron/tmm_faster) | C++ core | Genuinely fast, float64. Multithreaded by default. |
| [`tmm_fast`](https://github.com/MLResearchAtOSRAM/tmm_fast) | PyTorch | Batched, GPU-capable, autograd. |
| [`tmmax`](https://github.com/bahremsd/tmmax) | JAX | JIT and vmap. Bundles its own material database. |
| **tmmcore** | JS / C / WASM | This package. |

## Method

Every implementation reads the same precomputed file: one wavelength grid, one
set of complex indices, one set of thicknesses. Dispersion evaluation and
material lookup therefore sit **outside** the timed region, so what is measured
is the transfer-matrix core alone and nothing else.

One "spectrum" means every wavelength on the grid for **both** polarizations,
the unit an optimizer evaluates thousands of times per design.

Stacks are 4-, 7-, 21- and 40-layer; the 7-layer contains absorbing silver.
Grids are 400–1100 nm at 10 nm (71 points) and at 1 nm (701 points).

## Accuracy

Reproducible with Node alone, since Byrnes' outputs are committed and no Python
is needed:

```bash
npm run compare
```

Maximum absolute difference in R and T against Byrnes, over the whole grid, both
polarizations:

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

tmmcore and `tmm_faster` both agree with Byrnes at the level of float64
accumulation noise. The other two cannot: **both `tmm_fast` and `tmmax` are
single precision**. `tmm_fast` casts internally regardless of the dtype you hand
it; `tmmax` runs on JAX with `jax_enable_x64` left at its default of false. Both
are fine for generating machine-learning training data. Neither is suitable as
an accuracy reference, and neither can be checked against one.

## Speed

Single-threaded, best-of-N, one full spectrum, in milliseconds. Lower is better.

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

As speedup factors, single-threaded:

| Against | tmmcore WASM is |
|---|---|
| `tmm` (Byrnes) | 136–215× faster |
| `tmm_fast` (PyTorch) | 15–192× faster |
| `tmmax` (JAX) | 10–300× faster |
| `tmm_faster` (C++) | 2.0–2.7× faster |

## Caveats

!!! warning "Threading: tmmcore is not the fastest TMM library"

    `tmm_faster` is **multithreaded by default**, measured at a cpu/wall ratio
    of 13–15 on a 16-thread machine. In default configuration it *beats*
    tmmcore on large problems: 1.27 ms against 3.18 ms on the 40-layer,
    701-point case.

    The 2.0–2.7× advantage above is **per core**, with `OMP_NUM_THREADS=1`. That
    is the honest comparison of the code itself, and it is the relevant one when
    you parallelize at a level above the kernel. "Fastest TMM library" is not a
    claim this data supports, and it is not one made here.

!!! warning "Batching"

    The vectorized libraries exist to evaluate many stacks at once. tmmcore has
    no batch entry point and simply loops. It still wins single-threaded on 128
    distinct stacks (45.9 ms against 82.9 ms for `tmm_faster` and 209.8 ms for
    `tmm_fast`, on the 40-layer case). On a GPU with a large enough batch,
    `tmm_fast` and `tmmax` would pull ahead, and **that case was not measured**.

Further points of fairness, carried from the harness notes:

- `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering`,
  strictly more than R and T, with no cheaper entry point available in its API.
- Byrnes must be fed Python `complex`, not NumPy scalars; `coh_tmm` is markedly
  slower on `np.complex128` elements, which would have flattered every other
  contender.
- `tmmax` reads materials from its own bundled database rather than accepting
  explicit n,k, so its timing includes an interpolation the others do not
  perform.
- A control run confirms the Python driver loop is under 1% of total time, so
  the Byrnes figures measure the library rather than the harness.

## What the others do that tmmcore does not

- `tmm` (Byrnes) handles incoherent layers, ellipsometric parameters, and
  absorption per layer. It is a far broader package, and it is the reference the
  whole field is built on.
- `tmm_fast` and `tmmax` offer autograd through the whole calculation, not just
  thickness derivatives, and run on GPUs.
- `tmmax` ships a 35-material database, so you can compute something without
  supplying n,k at all.

tmmcore is deliberately narrow: the kernel, exact thickness derivatives, and
nothing else.
