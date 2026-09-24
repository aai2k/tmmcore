"""Emit reference_byrnes.json: Byrnes' tmm on every case in cases.json.

Values only, no timing, so the accuracy tier (compare.mjs) needs no Python.
For each case and wavelength: R and T in s and p, and the complex reflection
coefficients r_s and r_p as [re, im].

Byrnes and tmmcore share n~ = n + ik and the time factor exp(-i w t), so r_s
agrees as it stands. r_p does not: Byrnes follows the Fresnel convention,
r_p = -r_s at normal incidence, and tmmcore the tilted admittances, r_p = r_s.
compare.mjs accounts for that sign; nothing here converts.

    python gen_reference.py
"""

import json

import numpy as np
import tmm

DEG = np.pi / 180.0


def spectra(case):
    d_list = [np.inf] + list(case["thick"]) + [np.inf]
    th0 = case["theta_deg"] * DEG
    out = {"Rs": [], "Ts": [], "Rp": [], "Tp": [], "rs": [], "rp": []}
    for i, lam in enumerate(case["lambdas"]):
        # Python complex, not NumPy scalars: coh_tmm is markedly slower on the latter.
        n_list = ([complex(*case["n0"][i])]
                  + [complex(*layer[i]) for layer in case["layerNK"]]
                  + [complex(*case["ns"][i])])
        for pol in ("s", "p"):
            res = tmm.coh_tmm(pol, n_list, d_list, th0, lam)
            out["R" + pol].append(float(res["R"]))
            out["T" + pol].append(float(res["T"]))
            out["r" + pol].append([float(res["r"].real), float(res["r"].imag)])
    return out


def main():
    with open("cases.json") as fh:
        cases = json.load(fh)["cases"]
    reference = {
        "_comment": "Reference spectra from an independent implementation. "
                    "Regenerate with gen_reference.py, see benchmarks/README.md. "
                    "Values only; timings live elsewhere.",
        "source": {
            "package": "tmm",
            "author": "Steven Byrnes",
            "url": "https://github.com/sbyrnes321/tmm",
            "paper": "S. J. Byrnes, Multilayer optical calculations, arXiv:1603.02720",
            "version": getattr(tmm, "__version__", "0.2.0"),
            "convention": "n~ = n + ik, k >= 0, identical to tmmcore, so inputs transfer "
                          "verbatim; r_p takes the Fresnel sign, opposite to tmmcore's",
        },
        "inputs": "cases.json",
        "cases": {c["name"]: spectra(c) for c in cases},
    }
    with open("reference_byrnes.json", "w") as fh:
        json.dump(reference, fh, indent=1)
    print(f"wrote reference_byrnes.json: {len(cases)} cases")


if __name__ == "__main__":
    main()
