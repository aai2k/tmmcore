"""Benchmark the fast Python TMM libraries on the shared cases.

Contenders
----------
  byrnes      tmm            pure Python + NumPy, scalar in wavelength (the de-facto standard)
  tmm_faster  tmm_faster     C++ core, one call per stack, both polarizations
  tmmcore     tmmcore        C kernel in the wheel, one call per stack, both polarizations
  tmm_fast    tmm_fast       PyTorch, batched over stacks/wavelengths/angles
  tmmax       tmmax          JAX + JIT + vmap

Two workloads
-------------
  A (single)  one stack, full wavelength grid, both polarizations.
              This is the unit an optimizer evaluates thousands of times.
  B (batch)   NBATCH distinct stacks at once. This is the vectorized libraries'
              home turf, and matches ML dataset generation / DE populations.

Nothing here is checked against Byrnes: compare_tmm_faster.py checks the values
of tmm_faster, and tmm_fast is timed without a comparison. tmmax reads materials
from its own bundled database, so it is timed at matching layer/grid sizes but
not value-checked, and its timing also includes an nk interpolation the others
do not perform.
"""

import json
import os
import platform
import time
from importlib.metadata import version

# Thread pinning must happen before numpy/torch import their native backends.
THREADS = int(os.environ.get("TMMBENCH_THREADS", "0"))   # 0 = library default
if THREADS:
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS"):
        os.environ[var] = str(THREADS)

import numpy as np
import tmm

NBATCH = 128
DEG = np.pi / 180.0
SUFFIX = f"_t{THREADS}" if THREADS else ""


def best_of(fn, min_reps=3, min_seconds=1.0):
    best, reps = float('inf'), 0
    t0 = time.perf_counter()
    while reps < min_reps or (time.perf_counter() - t0) < min_seconds:
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
        reps += 1
    return best, reps


def index_matrix(case):
    """[nLam x (nLayers+2)] complex, semi-infinite media on both ends."""
    lam = case["lambdas"]
    nL = case["nLayers"]
    m = np.empty((len(lam), nL + 2), dtype=complex)
    for j in range(len(lam)):
        m[j, 0] = complex(*case["n0"][j])
        for k in range(nL):
            m[j, k + 1] = complex(*case["layerNK"][k][j])
        m[j, -1] = complex(*case["ns"][j])
    return m


# ── byrnes ───────────────────────────────────────────────────────────────────

def run_byrnes(case, nmat, thick):
    lam = case["lambdas"]
    d_list = [np.inf] + list(thick) + [np.inf]
    Rs = np.empty(len(lam))
    for i, l in enumerate(lam):
        n_list = [complex(z) for z in nmat[i]]
        Rs[i] = tmm.coh_tmm('s', n_list, d_list, 0.0, l)['R']
        tmm.coh_tmm('p', n_list, d_list, 0.0, l)
    return Rs


# ── tmm_faster (C++ core) ────────────────────────────────────────────────────

import tmm_faster

def run_tmm_faster(case, nmat, thick):
    d = np.concatenate([[0.0], thick, [0.0]])
    r = tmm_faster.calc_coherent(nmat, d, np.array([0.0]), np.array(case["lambdas"]))
    return np.asarray(r['R_s']).ravel()


# ── tmmcore (the wheel's C kernel, through ctypes) ───────────────────────────

import tmmcore

def tmmcore_medium(column):
    """One column of the index matrix as tmmcore's {n, k}, each an array over the grid."""
    return {"n": column.real.copy(), "k": column.imag.copy()}

def tmmcore_stack(nmat, thick):
    """The stack as tmmcore takes it, built before timing as the other inputs are."""
    return {"incident": tmmcore_medium(nmat[:, 0]),
            "substrate": tmmcore_medium(nmat[:, -1]),
            "layers": [{**tmmcore_medium(nmat[:, j + 1]), "d": float(d)} for j, d in enumerate(thick)]}

def run_tmmcore(stack, lam):
    return tmmcore.rta_spectrum(stack, lam)["Rs"]


# ── tmm_fast (PyTorch) ───────────────────────────────────────────────────────

import torch
from tmm_fast import coh_tmm as fast_coh_tmm

if THREADS:
    torch.set_num_threads(THREADS)

def _fast_tensors(case, nmat, thick, nstack):
    nLam = len(case["lambdas"])
    nL = case["nLayers"]
    N = torch.zeros((nstack, nL + 2, nLam), dtype=torch.cdouble)
    for k in range(nL + 2):
        N[:, k, :] = torch.tensor(nmat[:, k], dtype=torch.cdouble)
    T = torch.tensor(np.concatenate([[np.inf], thick, [np.inf]]),
                     dtype=torch.float64).repeat(nstack, 1)
    th = torch.tensor([0.0], dtype=torch.float64)
    wl = torch.tensor(case["lambdas"], dtype=torch.float64)
    return N, T, th, wl

def run_tmm_fast(N, T, th, wl):
    fast_coh_tmm('s', N, T, th, wl, device='cpu')
    r = fast_coh_tmm('p', N, T, th, wl, device='cpu')
    return np.asarray(r['R'])


# ── tmmax (JAX) ──────────────────────────────────────────────────────────────

import jax.numpy as jnp
from tmmax.tmm import tmm_coh

TMMAX_H, TMMAX_L = "TiO2", "SiO2"

def tmmax_materials(nL):
    return ["Air"] + [TMMAX_H if i % 2 == 0 else TMMAX_L for i in range(nL)] + ["SiO2"]

def run_tmmax(mats, thick_j, wl_j, ang_j):
    tmm_coh(mats, thick_j, wl_j, ang_j, 's')
    return tmm_coh(mats, thick_j, wl_j, ang_j, 'p')


# ── Driver ───────────────────────────────────────────────────────────────────

def main():
    # Every contender below is called at normal incidence, so the oblique cases
    # of the accuracy tier are left out rather than timed at the wrong angle.
    cases = [c for c in json.load(open("cases.json"))["cases"] if c["theta_deg"] == 0]
    out = []

    print("=== A: ONE stack, full grid, both polarizations (ms) ===")
    print(f"{'case':13s} {'N':>3s} {'nLam':>5s} {'byrnes':>9s} {'tmm_faster':>11s} "
          f"{'tmmcore':>9s} {'tmm_fast':>9s} {'tmmax':>9s}")
    print("-" * 75)

    for case in cases:
        nmat = index_matrix(case)
        thick = np.array(case["thick"])
        nL, nLam = case["nLayers"], len(case["lambdas"])

        run_byrnes(case, nmat, thick)              # warm-up
        t_by, _ = best_of(lambda: run_byrnes(case, nmat, thick))

        run_tmm_faster(case, nmat, thick)          # warm-up
        t_cpp, _ = best_of(lambda: run_tmm_faster(case, nmat, thick))

        stack, lam = tmmcore_stack(nmat, thick), np.array(case["lambdas"])
        run_tmmcore(stack, lam)                    # warm-up
        t_tc, _ = best_of(lambda: run_tmmcore(stack, lam))

        N, T, th, wl = _fast_tensors(case, nmat, thick, 1)
        r_ft = run_tmm_fast(N, T, th, wl)
        t_ft, _ = best_of(lambda: run_tmm_fast(N, T, th, wl))

        mats = tmmax_materials(nL)
        thick_j = jnp.array(thick * 1e-9)          # tmmax works in metres
        wl_j = jnp.array(np.array(case["lambdas"]) * 1e-9)
        ang_j = jnp.array([0.0])
        run_tmmax(mats, thick_j, wl_j, ang_j)      # warm-up: JIT compile
        t_ax, _ = best_of(lambda: run_tmmax(mats, thick_j, wl_j, ang_j))

        print(f"{case['name']:13s} {nL:3d} {nLam:5d} {t_by*1e3:9.3f} {t_cpp*1e3:11.3f} "
              f"{t_tc*1e3:9.3f} {t_ft*1e3:9.3f} {t_ax*1e3:9.3f}")

        out.append({"name": case["name"], "nLayers": nL, "nLam": nLam,
                    "single": {"byrnes": t_by, "tmm_faster": t_cpp, "tmmcore": t_tc,
                               "tmm_fast": t_ft, "tmmax": t_ax}})

    print(f"\n=== B: BATCH of {NBATCH} stacks, full grid, both pols (ms) ===")
    print(f"{'case':13s} {'N':>3s} {'nLam':>5s} {'byrnes':>9s} {'tmm_faster':>11s} "
          f"{'tmmcore':>9s} {'tmm_fast':>9s}   (per-stack ms in parens)")
    print("-" * 105)

    for case in cases:
        if not case["name"].endswith("g71"):
            continue
        nmat = index_matrix(case)
        thick = np.array(case["thick"])
        nL, nLam = case["nLayers"], len(case["lambdas"])

        # Distinct stacks: jitter thicknesses, as a real population would.
        rng = np.random.default_rng(0)
        thicks = thick[None, :] * (1 + 0.1 * rng.standard_normal((NBATCH, nL)))

        def by_batch():
            for s in range(NBATCH):
                run_byrnes(case, nmat, thicks[s])

        def cpp_batch():
            for s in range(NBATCH):
                run_tmm_faster(case, nmat, thicks[s])

        stacks, lam = [tmmcore_stack(nmat, t) for t in thicks], np.array(case["lambdas"])

        def tc_batch():
            for stack in stacks:
                run_tmmcore(stack, lam)

        N, T, th, wl = _fast_tensors(case, nmat, thick, NBATCH)
        T = torch.tensor(
            np.concatenate([np.full((NBATCH, 1), np.inf), thicks,
                            np.full((NBATCH, 1), np.inf)], axis=1), dtype=torch.float64)

        t_by, _ = best_of(by_batch, min_reps=1, min_seconds=0.5)
        t_cpp, _ = best_of(cpp_batch, min_reps=2, min_seconds=0.5)
        t_tc, _ = best_of(tc_batch, min_reps=2, min_seconds=0.5)
        run_tmm_fast(N, T, th, wl)
        t_ft, _ = best_of(lambda: run_tmm_fast(N, T, th, wl), min_reps=2, min_seconds=0.5)

        print(f"{case['name']:13s} {nL:3d} {nLam:5d} "
              f"{t_by*1e3:9.1f} {t_cpp*1e3:11.1f} {t_tc*1e3:9.1f} {t_ft*1e3:9.1f}   "
              f"({t_by/NBATCH*1e3:.3f} / {t_cpp/NBATCH*1e3:.3f} / {t_tc/NBATCH*1e3:.3f} "
              f"/ {t_ft/NBATCH*1e3:.3f})")

        for rec in out:
            if rec["name"] == case["name"]:
                rec["batch"] = {"n": NBATCH, "byrnes": t_by, "tmm_faster": t_cpp,
                                "tmmcore": t_tc, "tmm_fast": t_ft}

    fn = f"results_fast{SUFFIX}.json"
    versions = {"python": platform.python_version(),
                **{name: version(name) for name in ("numpy", "tmm", "tmm-faster", "tmmcore",
                                                    "tmm-fast", "torch", "tmmax", "jax")}}
    json.dump({"nbatch": NBATCH, "threads": THREADS or "default", "versions": versions,
               "results": out}, open(fn, "w"))
    print(f"\nwrote {fn}")


if __name__ == "__main__":
    main()
