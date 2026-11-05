# Comparison with other packages

This page sets tmmcore 0.5.0 beside four Python packages that do the same calculation: Byrnes' `tmm` 0.2.0, `tmm_faster` 0.1.3, `tmm_fast` 0.2.1 on PyTorch 2.13.0 (CPU build) and `tmmax` 1.1.4 on JAX 0.11.0, under Python 3.14.6 with NumPy 2.5.1, with tmmcore under Node 24.18.0. Everything was measured on 5 November 2026 on an AMD Ryzen 7 260 (8 cores, 16 threads) under Windows 11, plugged in, on the Performance power plan. The harness is in `benchmarks/`, and [its README](https://github.com/aai2k/tmmcore/blob/main/benchmarks/README.md#timing-tier) gives the commands that rerun every table on this page.

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

Byrnes' outputs are committed, so this needs only Node:

```bash
npm run compare
```

Largest absolute difference from Byrnes' `tmm` over the whole grid, at normal incidence:

| Case | `tmm_faster`, $R_s$ | tmmcore, $R$ and $T$ in s and p |
|---|--:|--:|
| 4-layer, 71 pts | 2.9e-16 | 2.1e-15 |
| 21-layer, 71 pts | 6.7e-15 | 1.2e-14 |
| 7-layer with Ag, 71 pts | 3.2e-15 | 1.6e-15 |
| 40-layer, 71 pts | 1.8e-14 | 1.7e-14 |
| 4-layer, 701 pts | 5.8e-16 | 2.6e-15 |
| 21-layer, 701 pts | 1.3e-14 | 2.4e-14 |
| 7-layer with Ag, 701 pts | 3.6e-15 | 3.1e-15 |
| 40-layer, 701 pts | 1.8e-14 | 8.6e-14 |

The tmmcore column is the larger of the JavaScript and WebAssembly figures of `npm run compare`; [Validation](validation.md#against-byrnes-tmm) shows its whole output, with the cases at 30° and 60°. The `tmm_faster` column comes from `bench_fast.py`, which runs it beside Byrnes' `tmm` on the same inputs and compares $R_s$ only. Both agree with Byrnes to float64 round-off. `tmm_fast` and `tmmax` are timed only; the harness does not compare their output with Byrnes.

Every case is lit from air. An absorbing incident medium has no cross-library reference here: Byrnes' `tmm` asserts that $n_0 \sin\theta_0$ is real and refuses a complex $n_0$ at a real angle. tmmcore evaluates that case from the real invariant $\mathrm{Re}(n_0)\sin\theta_0$ and checks it against an exact energy identity instead; see [Validation](validation.md#an-absorbing-incident-medium).

## Speed on one stack

One spectrum, every wavelength of the grid in both polarizations, at normal incidence, in milliseconds. Every package but `tmmax` runs on one thread; the harness does not set JAX's threads. Lower is better.

| Case | byrnes | tmm_fast | tmmax | tmm_faster | tmmcore JS | **tmmcore WASM** |
|---|--:|--:|--:|--:|--:|--:|
| 4-layer, 71 pts | 8.136 | 5.579 | 8.230 | 0.085 | 0.749 | **0.044** |
| 21-layer, 71 pts | 26.22 | 6.663 | 8.447 | 0.353 | 3.253 | **0.187** |
| 7-layer with Ag, 71 pts | 11.30 | 5.830 | 8.281 | 0.132 | 1.172 | **0.072** |
| 40-layer, 71 pts | 45.65 | 7.929 | 8.453 | 0.639 | 5.950 | **0.351** |
| 4-layer, 701 pts | 82.13 | 23.25 | 8.472 | 0.800 | 7.752 | **0.406** |
| 21-layer, 701 pts | 266.5 | 29.46 | 17.69 | 3.446 | 34.03 | **1.852** |
| 7-layer with Ag, 701 pts | 119.1 | 23.85 | 10.93 | 1.261 | 12.37 | **0.674** |
| 40-layer, 701 pts | 454.9 | 36.21 | 23.81 | 6.430 | 63.35 | **3.478** |

`tmmax` takes its own material data rather than the shared $n$ and $k$, so its column times stacks of the same size, not the same coatings; see [How the timing was done](#how-the-timing-was-done).

!!! warning "Threading"

    `tmm_faster` and `tmm_fast` are multithreaded by default, and the table above pins them to one thread. With its default threads `tmm_faster` keeps about 15 of the 16 hardware threads busy, as `probe_threads.py` measures, and it is faster than tmmcore in every case. `tmm_fast` runs no faster with its threads on one stack.

    | Case | `tmm_faster`, default threads | `tmm_fast`, default threads | tmmcore WASM, one thread |
    |---|--:|--:|--:|
    | 4-layer, 71 pts | 0.023 | 5.572 | 0.044 |
    | 21-layer, 71 pts | 0.052 | 6.644 | 0.187 |
    | 7-layer with Ag, 71 pts | 0.029 | 5.792 | 0.072 |
    | 40-layer, 71 pts | 0.083 | 7.979 | 0.351 |
    | 4-layer, 701 pts | 0.140 | 24.23 | 0.406 |
    | 21-layer, 701 pts | 0.383 | 29.12 | 1.852 |
    | 7-layer with Ag, 701 pts | 0.181 | 26.53 | 0.674 |
    | 40-layer, 701 pts | 0.657 | 32.87 | 3.478 |

    tmmcore's calls are single-threaded and leave threads to the caller, which can run one instance per worker.

## Speed on 128 stacks

The 71-point cases again, now as 128 stacks each, every stack's thicknesses scattered at random with a 10% standard deviation: the whole batch, every wavelength in both polarizations, in milliseconds. `tmm_fast` takes the 128 stacks in one call. Byrnes' `tmm`, `tmm_faster` and tmmcore take one stack per call, and the harness loops over the stacks. `tmmax` is not timed on this workload. Every package runs on one thread, except in the last two columns, which give `tmm_faster` and `tmm_fast` their default threads.

| Case | byrnes | tmm_faster | tmm_fast | tmmcore JS | **tmmcore WASM** | tmm_faster, default threads | tmm_fast, default threads |
|---|--:|--:|--:|--:|--:|--:|--:|
| 4-layer, 71 pts | 1059.8 | 11.46 | 29.12 | 103.5 | **6.236** | 3.091 | 18.10 |
| 21-layer, 71 pts | 3307.5 | 46.95 | 98.49 | 439.7 | **28.24** | 8.054 | 47.64 |
| 7-layer with Ag, 71 pts | 1464.7 | 17.68 | 39.30 | 159.5 | **10.47** | 3.915 | 22.09 |
| 40-layer, 71 pts | 5889.4 | 86.33 | 191.7 | 805.3 | **52.98** | 14.67 | 91.06 |

On a GPU with a large enough batch, `tmm_fast` and `tmmax` may pull ahead; **that case was not measured**.

## How the timing was done

- Every package but `tmmax` reads one precomputed file per case: the wavelength grid, the complex index of every layer at every wavelength, and the thicknesses. Dispersion and material lookup sit outside the timed region, so only the transfer-matrix calculation is timed.
- `tmmax` takes material names and its own data rather than $n$ and $k$. The harness gives it TiO2 and SiO2 in turn, at each case's layer count, thicknesses and grid, so its 7-layer row has no silver, and its timing includes an interpolation of its data that the others do not do. It runs on JAX, which computes in single precision unless `jax_enable_x64` is set; the harness does not set it.
- A package that takes one polarization per call is called once for each.
- tmmcore is timed through `rtaSpectrum`, on the JavaScript module and on a WebAssembly instance, with the stack built before timing, as the other drivers build their inputs.
- Each figure is the fastest of repeated timings, with a floor on their number and on the time spent. Byrnes' batches outlast the time floor and are timed once each.
- `TMMBENCH_THREADS=1` pins `tmm_faster` and `tmm_fast` to one thread. The default-thread figures come from the same script run without the pin.
- In the batch, the Python and JavaScript drivers scatter the thicknesses with different seeded generators, so the stacks have the same spread but are not the same stacks.
- Byrnes' `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering`, more than $R$ and $T$, and its API has no cheaper call.
- `overhead_py.py` puts the Python driver loop under 1% of Byrnes' time, so the Byrnes figures time the library rather than the harness.
- Timings move by tens of percent between sessions, so every timing on this page comes from one session on an otherwise idle machine. The whole set was run twice back to back, and the two runs agree within 10% on every figure; the tables give the second.
