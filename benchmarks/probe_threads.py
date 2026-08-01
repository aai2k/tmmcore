"""Is tmm_faster's C++ core using multiple threads?

Compares process CPU time against wall time for a long run. A ratio near 1
means single-threaded; near the core count means it is parallelizing (which
would make a wall-clock comparison against a single-threaded kernel unfair).
"""

import json
import time

import numpy as np
import tmm_faster

case = [c for c in json.load(open("cases.json"))["cases"] if c["name"] == "BIG40/g701"][0]
lam = np.array(case["lambdas"])
nL = case["nLayers"]
thick = np.array(case["thick"])

nmat = np.empty((len(lam), nL + 2), dtype=complex)
for j in range(len(lam)):
    nmat[j, 0] = complex(*case["n0"][j])
    for k in range(nL):
        nmat[j, k + 1] = complex(*case["layerNK"][k][j])
    nmat[j, -1] = complex(*case["ns"][j])
d = np.concatenate([[0.0], thick, [0.0]])
angles = np.array([0.0])

REPS = 400
tmm_faster.calc_coherent(nmat, d, angles, lam)          # warm-up

w0, c0 = time.perf_counter(), time.process_time()
for _ in range(REPS):
    tmm_faster.calc_coherent(nmat, d, angles, lam)
wall = time.perf_counter() - w0
cpu = time.process_time() - c0

print(f"reps       {REPS}")
print(f"wall       {wall:.3f} s   ({wall / REPS * 1e3:.3f} ms/call)")
print(f"cpu        {cpu:.3f} s")
print(f"cpu/wall   {cpu / wall:.2f}   (1.0 = single-threaded)")

# Same probe with many angles, where parallelism would be easiest to exploit.
angles_many = np.linspace(0.0, 60.0, 32)
tmm_faster.calc_coherent(nmat, d, angles_many, lam)
w0, c0 = time.perf_counter(), time.process_time()
for _ in range(20):
    tmm_faster.calc_coherent(nmat, d, angles_many, lam)
wall2 = time.perf_counter() - w0
cpu2 = time.process_time() - c0
print(f"\n32 angles: cpu/wall {cpu2 / wall2:.2f}  ({wall2 / 20 * 1e3:.3f} ms/call)")
