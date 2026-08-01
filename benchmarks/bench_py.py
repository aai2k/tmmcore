"""Benchmark Byrnes' `tmm` package (coh_tmm) on the shared cases.

Emits results_py.json: per-case timings plus the R/T spectra, which the Node
side cross-checks against both tmmcore kernels before reporting any speedups.

Byrnes' coh_tmm is scalar in wavelength, so one "spectrum evaluation" is
nLam x 2 calls (s and p) — the same unit of work the other two kernels are
timed on.
"""

import json
import time

import numpy as np
import tmm

DEG = np.pi / 180.0


def spectrum(case):
    """One full R/T spectrum, both polarizations. Returns (Rs, Ts, Rp, Tp)."""
    lambdas = case["lambdas"]
    thick = case["thick"]
    n0 = case["n0"]
    ns = case["ns"]
    nk = case["layerNK"]
    th0 = case["theta_deg"] * DEG
    nL = len(thick)
    nLam = len(lambdas)

    # d_list is fixed across lambda; inf on the semi-infinite bounding media.
    d_list = [np.inf] + list(thick) + [np.inf]

    Rs = np.empty(nLam); Ts = np.empty(nLam)
    Rp = np.empty(nLam); Tp = np.empty(nLam)

    for i, lam in enumerate(lambdas):
        n_list = [complex(n0[i][0], n0[i][1])]
        for k in range(nL):
            n_list.append(complex(nk[k][i][0], nk[k][i][1]))
        n_list.append(complex(ns[i][0], ns[i][1]))

        s = tmm.coh_tmm('s', n_list, d_list, th0, lam)
        p = tmm.coh_tmm('p', n_list, d_list, th0, lam)
        Rs[i] = s['R']; Ts[i] = s['T']
        Rp[i] = p['R']; Tp[i] = p['T']

    return Rs, Ts, Rp, Tp


def time_it(fn, min_reps, min_seconds):
    """Best-of timing: repeat until both the rep and time floors are met."""
    best = float('inf')
    reps = 0
    t_start = time.perf_counter()
    while reps < min_reps or (time.perf_counter() - t_start) < min_seconds:
        t0 = time.perf_counter()
        fn()
        dt = time.perf_counter() - t0
        best = min(best, dt)
        reps += 1
    return best, reps


def main():
    with open("cases.json") as fh:
        data = json.load(fh)

    out = []
    print(f"{'case':14s} {'N':>4s} {'nLam':>5s} {'spectrum ms':>12s} "
          f"{'us/eval':>9s} {'reps':>5s}")
    print("-" * 56)

    for case in data["cases"]:
        # Warm-up (import-time lazies, numpy dispatch caches).
        spectrum(case)

        best, reps = time_it(lambda: spectrum(case), min_reps=3, min_seconds=1.5)
        Rs, Ts, Rp, Tp = spectrum(case)

        nLam = len(case["lambdas"])
        evals = nLam * 2                     # s and p
        us_per_eval = best / evals * 1e6

        print(f"{case['name']:14s} {case['nLayers']:4d} {nLam:5d} "
              f"{best * 1e3:12.3f} {us_per_eval:9.3f} {reps:5d}")

        out.append({
            "name": case["name"],
            "nLayers": case["nLayers"],
            "nLam": nLam,
            "seconds": best,
            "evals": evals,
            "us_per_eval": us_per_eval,
            "reps": reps,
            "Rs": Rs.tolist(), "Ts": Ts.tolist(),
            "Rp": Rp.tolist(), "Tp": Tp.tolist(),
        })

    with open("results_py.json", "w") as fh:
        json.dump({"impl": "byrnes-tmm", "version": tmm.__version__ if hasattr(tmm, "__version__") else "0.2.0",
                   "results": out}, fh)
    print("\nwrote results_py.json")


if __name__ == "__main__":
    main()
