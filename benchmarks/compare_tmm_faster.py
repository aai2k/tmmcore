"""Accuracy of tmm_faster against Byrnes' tmm, on the shared cases.

Runs tmm_faster on every case of cases.json, at the case's angle, and compares
its R and T in s and p with Byrnes' outputs in reference_byrnes.json, the file
compare.mjs checks tmmcore against. Times nothing.

    python compare_tmm_faster.py
"""

import json
from importlib.metadata import version

import numpy as np
import tmm_faster


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


def main():
    with open("cases.json") as fh:
        cases = json.load(fh)["cases"]
    with open("reference_byrnes.json") as fh:
        ref = json.load(fh)

    src = ref["source"]
    print(f"tmm_faster {version('tmm_faster')}, NumPy {version('numpy')}")
    print(f"Reference: {src['package']} {src['version']}, {src['author']}\n")
    print(f"{'case':18s}{'deg':>6s}{'layers':>7s}{'points':>8s}{'max |diff|':>12s}")
    print("-" * 51)

    worst = 0.0
    for case in cases:
        r = ref["cases"][case["name"]]
        # tmm_faster takes the semi-infinite media as zero-thickness entries,
        # and the angle in degrees.
        d = np.concatenate([[0.0], case["thick"], [0.0]])
        out = tmm_faster.calc_coherent(index_matrix(case), d,
                                       np.array([float(case["theta_deg"])]),
                                       np.array(case["lambdas"], dtype=float))
        diff = max(float(np.max(np.abs(np.ravel(out[ours]) - np.array(r[theirs]))))
                   for ours, theirs in (("R_s", "Rs"), ("T_s", "Ts"),
                                        ("R_p", "Rp"), ("T_p", "Tp")))
        worst = max(worst, diff)
        print(f"{case['name']:18s}{case['theta_deg']:>6g}{case['nLayers']:7d}"
              f"{len(case['lambdas']):8d}{diff:12.1e}")

    print(f"\nLargest difference from Byrnes over {len(cases)} cases: "
          f"{worst:.1e} in R and T, s and p")


if __name__ == "__main__":
    main()
