# Benchmarks and cross-library comparison

Two tiers. They differ in what you have to install, not in which libraries are covered.

| | Command | Needs |
|---|---|---|
| **Accuracy** | `npm run compare` | Node only |
| **Timing** | see [below](#timing-tier) | A ~1.3 GB Python environment |

Results are published in [`../docs/comparison.md`](../docs/comparison.md), which is the only place they live. `summarize.mjs` emits `tables.md` here as an intermediate; its numbers are transcribed into the published page along with the caveats that give them meaning.

## Implementations compared

| Label | What it is |
|---|---|
| `byrnes` | The [`tmm`](https://github.com/sbyrnes321/tmm) PyPI package, `coh_tmm`. Pure Python and NumPy, scalar in wavelength. |
| `tmm_faster` | [`tmm-faster`](https://github.com/clembron/tmm_faster), a C++ core with a Python wrapper. Returns both polarizations per call. Multithreaded by default. |
| `tmm_fast` | [`tmm_fast`](https://github.com/MLResearchAtOSRAM/tmm_fast), on PyTorch: batched, GPU-capable, autograd. |
| `tmmax` | [`tmmax`](https://github.com/bahremsd/tmmax), on JAX with JIT and vmap. |
| `tmmcore js` | `rtaSpectrum` from the JavaScript module: one call for the whole grid, s and p together, which runs the one-wavelength JavaScript kernel at each point. |
| `tmmcore wasm` | `rtaSpectrum` on a WebAssembly instance: one call for the whole grid, s and p together, which runs the C kernel's `tmm_rta_spectrum()` once. |
| `tmmcore python` | `rta_spectrum` from the tmmcore wheel on PyPI: the C kernel compiled with the package, with MSVC on Windows, and called through ctypes. One call for the whole grid, s and p together. |
| `tmmcore c` | The same library, the wheel's `_tmmcore.dll` or `_tmmcore.so`, loaded by `bench_c.c` and called from C through `tmm_rta_spectrum()`, with no Python in the process. |

## The shared input set

`gen_cases.py` precomputes the wavelength grid, per-layer complex indices and thicknesses into `cases.json`, which every driver reads. Dispersion evaluation and material lookup therefore sit **outside** the timed region: what is measured is the transfer-matrix core alone.

Indices use ñ = n + ik with k ≥ 0, the convention tmmcore and Byrnes share, so values transfer verbatim with no conjugation and there is no material-data confound.

Cases are 4-, 7-, 21- and 40-layer stacks (the 7-layer contains absorbing silver) on two grids: 400-1100 nm at 10 nm (71 points) and at 1 nm (701 points), at normal incidence. The coarse grid runs again at 30° and 60°, where s and p part.

Adding a case is a pull request against `gen_cases.py`.

## Accuracy tier

```bash
npm run compare
```

Feeds `cases.json` to tmmcore's JavaScript and WebAssembly paths and diffs against `reference_byrnes.json`, output committed from Byrnes' `tmm`, so no Python is needed. R and T come from `rta()`, the complex reflection coefficient from `phase()`. Exits non-zero if either disagrees by more than 1e-12.

r_s is compared as it stands and r_p with its sign reversed, for the reason the [conventions](../docs/getting-started.md#conventions) give.

To regenerate the reference file:

```bash
python -m venv .venv
.venv/Scripts/python -m pip install tmm numpy
.venv/Scripts/python gen_cases.py       # regenerate the shared inputs
.venv/Scripts/python gen_reference.py   # runs Byrnes, writes reference_byrnes.json
```

`npm run compare` shows that tmmcore agrees with Byrnes to float64 round-off. `compare_tmm_faster.py` checks R and T of `tmm_faster`, in s and p, against the same file; it runs in the Python environment of the [timing tier](#timing-tier). `tmm_fast` and `tmmax` are timed without a comparison.

## Timing tier

```bash
python -m venv .venv
.venv/Scripts/python -m pip install tmm numpy tmmax tmm-faster tmmcore
.venv/Scripts/python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
.venv/Scripts/python -m pip install tmm-fast

.venv/Scripts/python gen_cases.py                        # shared inputs
.venv/Scripts/python bench_py.py                         # Byrnes    -> results_py.json
TMMBENCH_THREADS=1 .venv/Scripts/python bench_fast.py    # fast libs -> results_fast_t1.json
node bench_js.mjs                                        # tmmcore, single stack
node bench_js_batch.mjs                                  # tmmcore, 128-stack batch
.venv/Scripts/python bench_c.py                          # tmmcore's library from C -> results_c.json
node summarize.mjs                                       #          -> tables.md
.venv/Scripts/python bench_fast.py                       # default threads -> results_fast.json
.venv/Scripts/python probe_threads.py                    # cpu/wall of tmm_faster
.venv/Scripts/python overhead_py.py                      # cost of the Python driver loop
.venv/Scripts/python compare_tmm_faster.py               # tmm_faster against Byrnes, no timing
```

On PowerShell set the thread pin with `$env:TMMBENCH_THREADS="1"` beforehand, and remove it with `Remove-Item Env:TMMBENCH_THREADS` before the default-thread run.

`bench_js*.mjs` need `src/wasm/tmm_kernel.wasm`, which ships prebuilt. `bench_c.py` builds `bench_c.c` with `gcc`, or with the compiler `CC` names; on Windows the gcc of MSYS2's ucrt64 does.

Timing is best-of-N with both a repetition floor and a wall-time floor per case. One "spectrum" is every wavelength on the grid for **both** polarizations, the unit an optimizer evaluates thousands of times per design.

Run the whole tier in one session on an otherwise idle machine. Timings drift by tens of percent between sessions, so mixing runs produces a table whose rows are not comparable with one another.

### Threading

`tmm_faster` and `tmm_fast` are multithreaded by default; tmmcore's kernels are single-threaded and expect the caller to parallelize above them. `TMMBENCH_THREADS=1` pins `tmm_faster` and `tmm_fast` to one thread for a per-core comparison, which is what the main published tables report; it sets no thread count for JAX. The run without the pin writes `results_fast.json`, the source of the published default-thread figures. `probe_threads.py` measures the cpu/wall ratio that motivated the pinning.

With its default threads, `tmm_faster` is faster than tmmcore on every case.

## Fairness notes

- `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering`: strictly more than R and T, with no cheaper entry point in its API.
- `overhead_py.py` confirms the Python driver loop is under 1% of total time, so the Byrnes figures measure the library rather than the harness.
- `tmmax` reads materials from its own bundled database rather than accepting explicit n,k, so it is timed at matching layer and grid sizes but not value-checked, and its timing includes an interpolation the others do not perform. It runs on JAX, which computes in single precision unless `jax_enable_x64` is set; the harness does not set it.
- `tmmcore python` and `tmmcore c` run the same machine code, the library in the wheel, on the same stacks, so the difference between them is what the Python layer costs per call. `bench_c.py` checks that the C calls return the wheel's outputs to the bit.
- Byrnes' single-stack timings come from `bench_py.py`, in a process of its own. Its batch timings come from `bench_fast.py`, which shares its process with torch and JAX and times each batch once, since one batch outlasts the time floor.

## Files

| File | Role |
|---|---|
| `compare.mjs` | Accuracy tier, Node only |
| `compare_tmm_faster.py` | `tmm_faster` against `reference_byrnes.json`, in the timing tier's Python |
| `reference_byrnes.json` | Committed reference spectra |
| `gen_reference.py` | Runs Byrnes on `cases.json`, writes `reference_byrnes.json` |
| `gen_cases.py` | Emits `cases.json`, the shared input set |
| `bench_py.py` | Byrnes timing and reference spectra |
| `bench_fast.py` | `tmm_faster` / tmmcore wheel / `tmm_fast` / `tmmax`, single and batch |
| `bench_js.mjs` | tmmcore JS and WASM, single stack, plus agreement check |
| `bench_js_batch.mjs` | tmmcore JS and WASM, 128-stack batch |
| `bench_c.c` | The tmmcore wheel's library timed from C, single stack and batch |
| `bench_c.py` | Builds and runs `bench_c.c`, checks its outputs against the wheel, writes `results_c.json` |
| `summarize.mjs` | Joins every result file into `tables.md` |
| `overhead_py.py` | Control: the Python driver's own cost |
| `probe_apis.py` | Pins down each library's index and thickness conventions |
| `probe_threads.py` | Measures cpu/wall to detect multithreading |
