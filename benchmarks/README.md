# Benchmarks and cross-library comparison

Two tiers. They differ in what you have to install, not in which libraries are
covered.

| | Command | Needs |
|---|---|---|
| **Accuracy** | `npm run compare` | Node only |
| **Timing** | see [below](#timing-tier) | A ~1.3 GB Python environment |

Results are published in [`../docs/comparison.md`](../docs/comparison.md),
which is the only place they live. `summarize.mjs` emits `tables.md` here as an
intermediate; its numbers are transcribed into the published page along with the
caveats that give them meaning.

## Implementations compared

| Label | What it is |
|---|---|
| `byrnes` | The [`tmm`](https://github.com/sbyrnes321/tmm) PyPI package, `coh_tmm`. Pure Python and NumPy, scalar in wavelength. The de-facto standard. |
| `tmm_faster` | [`tmm-faster`](https://github.com/clembron/tmm_faster) — C++ core with a Python wrapper. Returns both polarizations per call. Multithreaded by default. |
| `tmm_fast` | [`tmm_fast`](https://github.com/MLResearchAtOSRAM/tmm_fast) — PyTorch, batched, GPU-capable, autograd. |
| `tmmax` | [`tmmax`](https://github.com/bahremsd/tmmax) — JAX with JIT and vmap. |
| `tmmcore js` | The JavaScript `tmm()`, one call per (λ, polarization). |
| `tmmcore wasm` | The C kernel's `tmm_spectrum()`, one call for the whole grid and both polarizations. |

## The shared input set

`gen_cases.py` precomputes the wavelength grid, per-layer complex indices and
thicknesses into `cases.json`, which every driver reads. Dispersion evaluation
and material lookup therefore sit **outside** the timed region: what is measured
is the transfer-matrix core alone.

Indices use ñ = n + ik with k ≥ 0 — the convention tmmcore and Byrnes share — so
values transfer verbatim with no conjugation and there is no material-data
confound.

Cases are 4-, 7-, 21- and 40-layer stacks (the 7-layer contains absorbing
silver) on two grids: 400–1100 nm at 10 nm (71 points) and at 1 nm (701 points).

Adding a case is a pull request against `cases.json`.

## Accuracy tier

```bash
npm run compare
```

Feeds `cases.json` to tmmcore's JavaScript and WebAssembly paths and diffs
against `reference_byrnes.json` — output committed from Byrnes' `tmm`, so no
Python is needed. Exits non-zero if disagreement exceeds 1e-12.

Distrusting a committed reference file is the correct instinct. To regenerate it
yourself:

```bash
python -m venv .venv
.venv/Scripts/python -m pip install tmm numpy
.venv/Scripts/python gen_cases.py     # regenerate the shared inputs
.venv/Scripts/python bench_py.py      # runs Byrnes, writes results_py.json
```

`results_py.json` carries the spectra; `reference_byrnes.json` is the same
values with the timing fields stripped.

Only Byrnes and `tmm_faster` are usable as accuracy references. `tmm_fast` casts
to single precision internally, and `tmmax` runs on JAX with `jax_enable_x64`
left at its default of false — both agree only to ~1e-7 regardless of the input
dtype.

## Timing tier

```bash
python -m venv .venv
.venv/Scripts/python -m pip install tmm numpy tmmax tmm-faster
.venv/Scripts/python -m pip install torch --index-url https://download.pytorch.org/whl/cpu
.venv/Scripts/python -m pip install tmm-fast

.venv/Scripts/python gen_cases.py                        # shared inputs
.venv/Scripts/python bench_py.py                         # Byrnes    -> results_py.json
TMMBENCH_THREADS=1 .venv/Scripts/python bench_fast.py    # fast libs -> results_fast_t1.json
node bench_js.mjs                                        # tmmcore, single stack
node bench_js_batch.mjs                                  # tmmcore, 128-stack batch
node summarize.mjs                                       #          -> tables.md
```

On PowerShell set the thread pin with `$env:TMMBENCH_THREADS="1"` beforehand.

`bench_js*.mjs` need `src/tmm_kernel.wasm`, which ships prebuilt.

Timing is best-of-N with both a repetition floor and a wall-time floor per case.
One "spectrum" is every wavelength on the grid for **both** polarizations — the
unit an optimizer evaluates thousands of times per design.

Run the whole tier in one session on an otherwise idle machine. Timings drift by
tens of percent between sessions, so mixing runs produces a table whose rows are
not comparable with one another.

### Threading

`tmm_faster` and `tmm_fast` are multithreaded by default; tmmcore's kernels are
single-threaded and expect the caller to parallelize above them.
`TMMBENCH_THREADS=1` pins the Python side to one thread for a per-core
comparison, which is what the published tables report. `probe_threads.py` measures the
cpu/wall ratio that motivated the pinning.

In default configuration `tmm_faster` beats tmmcore on large problems. That is
stated plainly in the published comparison and should stay that way.

## Fairness notes

- `coh_tmm` also returns `vw_list`, `kz_list`, `r`, `t` and `power_entering` —
  strictly more than R and T, with no cheaper entry point in its API.
- `overhead_py.py` confirms the Python driver loop is under 1% of total time, so
  the Byrnes figures measure the library rather than the harness.
- Byrnes must be fed Python `complex`, not NumPy scalars; `coh_tmm` is markedly
  slower on `np.complex128` elements, which would flatter every other contender.
- `tmmax` reads materials from its own bundled database rather than accepting
  explicit n,k, so it is timed at matching layer and grid sizes but not
  value-checked, and its timing includes an interpolation the others do not
  perform.
- Canonical Byrnes timings come from `bench_py.py` in a dedicated process; the
  copy inside `results_fast_t1.json` shares a process with torch and JAX.

## Files

| File | Role |
|---|---|
| `compare.mjs` | Accuracy tier — Node only |
| `reference_byrnes.json` | Committed reference spectra |
| `gen_cases.py` | Emits `cases.json`, the shared input set |
| `bench_py.py` | Byrnes timing and reference spectra |
| `bench_fast.py` | `tmm_faster` / `tmm_fast` / `tmmax`, single and batch |
| `bench_js.mjs` | tmmcore JS and WASM, single stack, plus agreement check |
| `bench_js_batch.mjs` | tmmcore JS and WASM, 128-stack batch |
| `summarize.mjs` | Joins every result file into `tables.md` |
| `overhead_py.py` | Control: the Python driver's own cost |
| `probe_apis.py` | Pins down each library's index and thickness conventions |
| `probe_threads.py` | Measures cpu/wall to detect multithreading |
