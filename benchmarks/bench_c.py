"""tmmcore's C kernel timed from C, on the shared cases.

bench_c.c loads the library inside the installed tmmcore wheel, the one
tmmcore's Python calls go through, and times tmm_rta_spectrum with no Python
in its process. This script builds it, writes it the normal-incidence cases of
cases.json with the batch stacks of bench_fast.py, and checks the outputs of
its single-stack calls against tmmcore.rta_spectrum on the same stacks: they
must agree to the bit. Writes results_c.json for summarize.mjs.

The C compiler, gcc unless CC names another, builds only the timing loop; the
timed code is the wheel's.

    python bench_c.py
"""

import json
import os
import platform
import subprocess
import sys
import tempfile
from importlib.metadata import version

import numpy as np
import tmmcore

NBATCH = 128
HERE = os.path.dirname(os.path.abspath(__file__))
LIBRARY = os.path.join(os.path.dirname(tmmcore.__file__),
                       "_tmmcore.dll" if sys.platform == "win32" else "_tmmcore.so")
QUANTITIES = ("Rs", "Ts", "As", "Rp", "Tp", "Ap")


def batch_thicknesses(thick):
    """The batch bench_fast.py times: the thicknesses scattered by 10%, from the same seed."""
    rng = np.random.default_rng(0)
    return thick[None, :] * (1 + 0.1 * rng.standard_normal((NBATCH, len(thick))))


def stack_of(case):
    """The case as tmmcore.rta_spectrum takes it, n and k as arrays over the grid."""
    def medium(pairs):
        pairs = np.asarray(pairs, dtype=np.float64)
        return {"n": pairs[:, 0].copy(), "k": pairs[:, 1].copy()}
    return {"incident": medium(case["n0"]), "substrate": medium(case["ns"]),
            "layers": [{**medium(row), "d": d} for row, d in zip(case["layerNK"], case["thick"])]}


def write_stacks(path, cases):
    with open(path, "wb") as file:
        for case in cases:
            thick = np.array(case["thick"], dtype=np.float64)
            batch = batch_thicknesses(thick) if case["name"].endswith("g71") else np.empty((0, len(thick)))
            np.array([len(thick), len(case["lambdas"]), len(batch)], dtype="<i4").tofile(file)
            for values in (case["theta_deg"], case["lambdas"], case["n0"], case["ns"],
                           case["layerNK"], thick, batch):
                np.asarray(values, dtype="<f8").tofile(file)


def main():
    with open(os.path.join(HERE, "cases.json")) as file:
        cases = [c for c in json.load(file)["cases"] if c["theta_deg"] == 0]

    with tempfile.TemporaryDirectory() as work:
        program = os.path.join(work, "bench_c.exe" if sys.platform == "win32" else "bench_c")
        build = [os.environ.get("CC", "gcc"), "-O2", "-o", program, os.path.join(HERE, "bench_c.c")]
        subprocess.run(build + ([] if sys.platform == "win32" else ["-ldl"]), check=True)
        stacks, outputs = os.path.join(work, "stacks.bin"), os.path.join(work, "outputs.bin")
        write_stacks(stacks, cases)
        run = subprocess.run([program, LIBRARY, stacks, outputs], check=True, capture_output=True, text=True)
        produced = np.fromfile(outputs, dtype="<f8")

    print(f"tmmcore {version('tmmcore')}, its library called from C (ms)\n")
    print(f"{'case':13s} {'N':>3s} {'nLam':>5s} {'one stack':>10s} {f'{NBATCH} stacks':>12s}")
    print("-" * 47)
    out, at = [], 0
    for case, line in zip(cases, run.stdout.splitlines(), strict=True):
        n_lam = len(case["lambdas"])
        ours = produced[at:at + 6 * n_lam].reshape(6, n_lam)
        at += 6 * n_lam
        theirs = tmmcore.rta_spectrum(stack_of(case), case["lambdas"])
        for row, key in zip(ours, QUANTITIES):
            if not np.array_equal(row, theirs[key]):
                sys.exit(f"{case['name']}: {key} from C differs from tmmcore.rta_spectrum "
                         f"by up to {np.abs(row - theirs[key]).max():.1e}")

        t_one, t_batch = map(float, line.split())
        rec = {"name": case["name"], "nLayers": case["nLayers"], "nLam": n_lam, "single_s": t_one}
        if t_batch:
            rec["batch_s"] = t_batch
        out.append(rec)
        batch_ms = f"{t_batch * 1e3:12.2f}" if t_batch else f"{'':12s}"
        print(f"{case['name']:13s} {case['nLayers']:3d} {n_lam:5d} {t_one * 1e3:10.4f} {batch_ms}")

    print("\nThe outputs of every single-stack call match tmmcore.rta_spectrum to the bit.")
    versions = {"python": platform.python_version(), "numpy": version("numpy"), "tmmcore": version("tmmcore")}
    with open(os.path.join(HERE, "results_c.json"), "w") as file:
        json.dump({"impl": "tmmcore-c", "nbatch": NBATCH, "versions": versions, "results": out}, file)
    print("wrote results_c.json")


if __name__ == "__main__":
    main()
