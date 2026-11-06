# Comparison with other packages

This page sets tmmcore 0.5.0 beside four Python packages that do the same calculation: Byrnes' `tmm` 0.2.0, `tmm_faster` 0.1.3, `tmm_fast` 0.2.1 on PyTorch 2.13.0 (CPU build) and `tmmax` 1.1.4 on JAX 0.11.0, under Python 3.14.6 with NumPy 2.5.1. tmmcore runs under Node 24.18.0 and, from its PyPI wheel, under the same Python; the wheel's library was built with MSVC 14.51. Everything was measured on 6 November 2026 on an AMD Ryzen 7 260 (8 cores, 16 threads) under Windows 11, plugged in, on the Performance power plan. The harness is in `benchmarks/`, and [its README](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md#timing-tier) gives the commands that rerun every table on this page.

## What each package does

| Package | Language | One call covers | Derivatives | Incoherent layers | GPU |
|---|---|---|---|---|---|
| [`tmm`](https://github.com/sbyrnes321/tmm) (Byrnes) | Python and NumPy | one wavelength, one angle, one polarization | none | yes | no |
| [`tmm_faster`](https://github.com/clembron/tmm_faster) | C++, called from Python | one stack over many wavelengths and angles, s and p | none | yes | no |
| [`tmm_fast`](https://github.com/MLResearchAtOSRAM/tmm_fast) | PyTorch | many stacks over many wavelengths and angles, one polarization | automatic, through PyTorch | no | yes |
| [`tmmax`](https://github.com/bahremsd/tmmax) | JAX | one stack over many wavelengths and angles, one polarization | automatic, through JAX | yes | yes |
| **tmmcore** | JavaScript, WebAssembly and C, and Python on the same C | one stack over a wavelength grid; R and T in s and p together | analytic: [in the thicknesses](api.md#thickness-derivatives) to second order, the [needle scan](api.md#needle-scan), and [GD, GDD and TOD](api.md#phase-gd-gdd-and-tod) with their thickness derivatives | no | no |

Byrnes' `tmm` also gives ellipsometric $\Psi$ and $\Delta$, the absorption in each layer and the field at any depth. `tmmax` ships a 35-material database, so it computes without being given $n$ and $k$. tmmcore gives the [monitoring curve](api.md#monitoring-and-deposition) of a coating as it is deposited, and its [examples](examples.md) build an incoherent back side and ellipsometric $\Psi$ and $\Delta$ on its calls.

## Accuracy

Byrnes' outputs are committed, so checking tmmcore against them needs only Node:

```bash
npm run compare
```

Largest absolute difference in $R$ and $T$, in s and p, from Byrnes' `tmm` over the whole grid:

| Case | Angle | `tmm_faster` | tmmcore |
|---|--:|--:|--:|
| 4-layer, 71 pts | 0° | 2.2e-15 | 2.1e-15 |
| 21-layer, 71 pts | 0° | 1.1e-14 | 1.2e-14 |
| 7-layer with Ag, 71 pts | 0° | 3.2e-15 | 1.6e-15 |
| 40-layer, 71 pts | 0° | 2.2e-14 | 1.7e-14 |
| 4-layer, 701 pts | 0° | 2.4e-15 | 2.6e-15 |
| 21-layer, 701 pts | 0° | 2.1e-14 | 2.4e-14 |
| 7-layer with Ag, 701 pts | 0° | 3.6e-15 | 3.1e-15 |
| 40-layer, 701 pts | 0° | 9.2e-14 | 8.6e-14 |
| 4-layer, 71 pts | 30° | 2.0e-15 | 2.3e-15 |
| 21-layer, 71 pts | 30° | 1.3e-14 | 1.7e-14 |
| 7-layer with Ag, 71 pts | 30° | 3.0e-15 | 2.0e-15 |
| 40-layer, 71 pts | 30° | 3.6e-14 | 5.4e-14 |
| 4-layer, 71 pts | 60° | 2.9e-15 | 2.4e-15 |
| 21-layer, 71 pts | 60° | 1.7e-14 | 1.5e-14 |
| 7-layer with Ag, 71 pts | 60° | 3.1e-15 | 2.6e-15 |
| 40-layer, 71 pts | 60° | 1.0e-13 | 8.1e-14 |

The tmmcore column is the larger of the JavaScript and WebAssembly figures of `npm run compare`, whose whole output is on the [Validation](validation.md#against-byrnes-tmm) page. The `tmm_faster` column comes from `compare_tmm_faster.py`, which runs `tmm_faster` on the same inputs and checks it against the same committed Byrnes outputs. Both agree with Byrnes to float64 round-off. `tmm_fast` and `tmmax` are timed only; the harness does not compare their output with Byrnes. The Python and C speed figures below time the library in tmmcore's wheel, a native build of the same C, which the package's own tests compare with the WebAssembly kernel; the C timing driver checks that its calls return the wheel's outputs to the bit.

Every case is lit from air. An absorbing incident medium has no cross-library reference here: Byrnes' `tmm` asserts that $n_0 \sin\theta_0$ is real and refuses a complex $n_0$ at a real angle. tmmcore evaluates that case from the real invariant $\mathrm{Re}(n_0)\sin\theta_0$ and checks it against an exact energy identity instead; see [Validation](validation.md#an-absorbing-incident-medium).

## Speed on one stack

One spectrum, every wavelength of the grid in both polarizations, at normal incidence, in milliseconds. Every package but `tmmax` runs on one thread; the harness does not set JAX's threads. Lower is better.

| Case | byrnes | tmm_fast | tmmax | tmm_faster | tmmcore Python | tmmcore JS | tmmcore WASM | **tmmcore C** |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| 4-layer, 71 pts | 8.126 | 5.479 | 8.095 | 0.083 | 0.087 | 0.734 | 0.044 | **0.035** |
| 21-layer, 71 pts | 26.12 | 6.582 | 8.019 | 0.346 | 0.309 | 3.204 | 0.186 | **0.168** |
| 7-layer with Ag, 71 pts | 11.37 | 5.768 | 7.906 | 0.130 | 0.126 | 1.170 | 0.071 | **0.058** |
| 40-layer, 71 pts | 46.13 | 7.528 | 8.120 | 0.636 | 0.561 | 5.998 | 0.347 | **0.318** |
| 4-layer, 701 pts | 82.01 | 22.03 | 8.338 | 0.784 | 0.419 | 7.721 | 0.395 | **0.350** |
| 21-layer, 701 pts | 260.0 | 26.89 | 16.66 | 3.388 | 1.935 | 33.47 | 1.851 | **1.722** |
| 7-layer with Ag, 701 pts | 113.9 | 22.62 | 9.927 | 1.230 | 0.676 | 12.07 | 0.679 | **0.579** |
| 40-layer, 701 pts | 456.9 | 34.84 | 23.10 | 6.356 | 3.644 | 62.55 | 3.445 | **3.201** |

`tmmax` takes its own material data rather than the shared $n$ and $k$, so its column times stacks of the same size, not the same coatings; see [How the timing was done](#how-the-timing-was-done).

!!! warning "Threading"

    `tmm_faster` and `tmm_fast` are multithreaded by default, and the table above pins them to one thread. With its default threads `tmm_faster` keeps about 15 of the 16 hardware threads busy, as `probe_threads.py` measures, and it is faster than tmmcore in every case. `tmm_fast` gains little from its threads on one stack.

    | Case | `tmm_faster`, default threads | `tmm_fast`, default threads | tmmcore C, one thread |
    |---|--:|--:|--:|
    | 4-layer, 71 pts | 0.023 | 5.477 | 0.035 |
    | 21-layer, 71 pts | 0.052 | 6.509 | 0.168 |
    | 7-layer with Ag, 71 pts | 0.028 | 5.772 | 0.058 |
    | 40-layer, 71 pts | 0.083 | 7.760 | 0.318 |
    | 4-layer, 701 pts | 0.139 | 24.52 | 0.350 |
    | 21-layer, 701 pts | 0.386 | 25.10 | 1.722 |
    | 7-layer with Ag, 701 pts | 0.182 | 22.92 | 0.579 |
    | 40-layer, 701 pts | 0.660 | 29.31 | 3.201 |

    tmmcore's calls are single-threaded and leave threads to the caller, which can run one instance per worker.

## Speed on 128 stacks

The 71-point cases again, now as 128 stacks each, every stack's thicknesses scattered at random with a 10% standard deviation: the whole batch, every wavelength in both polarizations, in milliseconds. `tmm_fast` takes the 128 stacks in one call. Byrnes' `tmm`, `tmm_faster` and tmmcore take one stack per call, and the harness loops over the stacks. `tmmax` is not timed on this workload. Every package runs on one thread, except in the last two columns, which give `tmm_faster` and `tmm_fast` their default threads.

| Case | byrnes | tmm_faster | tmm_fast | tmmcore Python | tmmcore JS | tmmcore WASM | **tmmcore C** | tmm_faster, default threads | tmm_fast, default threads |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| 4-layer, 71 pts | 1083.2 | 11.17 | 28.35 | 11.98 | 104.4 | 6.054 | **4.692** | 3.078 | 17.54 |
| 21-layer, 71 pts | 3387.9 | 45.93 | 95.88 | 43.47 | 438.4 | 28.25 | **23.66** | 7.275 | 45.77 |
| 7-layer with Ag, 71 pts | 1499.5 | 17.15 | 38.26 | 17.26 | 155.8 | 10.43 | **7.942** | 3.794 | 20.16 |
| 40-layer, 71 pts | 6013.2 | 84.40 | 186.9 | 78.84 | 798.2 | 52.51 | **44.56** | 13.21 | 87.21 |

On a GPU with a large enough batch, `tmm_fast` and `tmmax` may pull ahead; **that case was not measured**.

## How the timing was done

- Every package but `tmmax` reads one precomputed file per case: the wavelength grid, the complex index of every layer at every wavelength, and the thicknesses. Dispersion and material lookup sit outside the timed region, so only the transfer-matrix calculation is timed.
- `tmmax` takes material names and its own data rather than $n$ and $k$. The harness gives it TiO2 and SiO2 in turn, at each case's layer count, thicknesses and grid, so its 7-layer row has no silver, and its timing includes an interpolation of its data that the others do not do. It runs on JAX, which computes in single precision unless `jax_enable_x64` is set; the harness does not set it.
- A package that takes one polarization per call is called once for each.
- tmmcore is timed through `rtaSpectrum` on the JavaScript module and on a WebAssembly instance, through `rta_spectrum` from its Python wheel, and through `tmm_rta_spectrum` called from C, each with the stack built before timing, as the other drivers build their inputs. The Python and C columns run the same library, the one in the wheel. The C column calls it from a C program with no Python in the process, so the gap between the two is what the Python layer costs per call.
- Each figure is the fastest of repeated timings, with a floor on their number and on the time spent. Byrnes' batches outlast the time floor and are timed once each.
- `TMMBENCH_THREADS=1` pins `tmm_faster` and `tmm_fast` to one thread. The default-thread figures come from the same script run without the pin.
- In the batch, the Python and C drivers time the same stacks. The JavaScript driver scatters the thicknesses with its own seeded generator, so its stacks have the same spread but are not the same stacks.
- Byrnes' `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering`, more than $R$ and $T$, and its API has no cheaper call.
- `overhead_py.py` puts the Python driver loop under 1% of Byrnes' time, so the Byrnes figures time the library rather than the harness.
- Timings move by tens of percent between sessions, so every timing on this page comes from one session on an otherwise idle machine. The whole set was run twice back to back, and the tables give the second run. Between the runs, the figures on one thread moved by 5% at most, except the JavaScript batch, by up to 14%; `tmmax` and the default-thread figures moved by up to 11%.
