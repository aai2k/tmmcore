"""How much of the Python time is loop/marshalling overhead vs. actual TMM?

Two probes on the same cases:
  harness  — the per-lambda n_list construction only, coh_tmm never called.
             This is the Python-interpreter floor the library cannot go below.
  full     — the real bench_py.py inner loop.

Also reports what coh_tmm returns, since it computes more than R and T and
there is no cheaper entry point in the API.
"""

import json
import time

import numpy as np
import tmm

DEG = np.pi / 180.0


def build_only(case):
    lambdas = case["lambdas"]; thick = case["thick"]
    n0 = case["n0"]; ns = case["ns"]; nk = case["layerNK"]
    nL = len(thick)
    d_list = [np.inf] + list(thick) + [np.inf]
    sink = 0.0
    for i, lam in enumerate(lambdas):
        n_list = [complex(n0[i][0], n0[i][1])]
        for k in range(nL):
            n_list.append(complex(nk[k][i][0], nk[k][i][1]))
        n_list.append(complex(ns[i][0], ns[i][1]))
        sink += n_list[0].real + d_list[0]
    return sink


def full(case):
    lambdas = case["lambdas"]; thick = case["thick"]
    n0 = case["n0"]; ns = case["ns"]; nk = case["layerNK"]
    th0 = case["theta_deg"] * DEG
    nL = len(thick)
    d_list = [np.inf] + list(thick) + [np.inf]
    acc = 0.0
    for i, lam in enumerate(lambdas):
        n_list = [complex(n0[i][0], n0[i][1])]
        for k in range(nL):
            n_list.append(complex(nk[k][i][0], nk[k][i][1]))
        n_list.append(complex(ns[i][0], ns[i][1]))
        acc += tmm.coh_tmm('s', n_list, d_list, th0, lam)['R']
        acc += tmm.coh_tmm('p', n_list, d_list, th0, lam)['R']
    return acc


def best_of(fn, reps=5):
    b = float('inf')
    for _ in range(reps):
        t = time.perf_counter(); fn(); b = min(b, time.perf_counter() - t)
    return b


with open("cases.json") as fh:
    cases = json.load(fh)["cases"]

print(f"{'case':14s} {'harness ms':>11s} {'full ms':>10s} {'harness %':>10s}")
print("-" * 48)
for c in cases:
    if not c["name"].endswith("g71"):
        continue
    build_only(c); full(c)
    h = best_of(lambda: build_only(c))
    f = best_of(lambda: full(c))
    print(f"{c['name']:14s} {h*1e3:11.3f} {f*1e3:10.3f} {h/f*100:9.1f}%")

# What coh_tmm actually hands back.
c = cases[0]
n_list = [complex(c['n0'][0][0], c['n0'][0][1])] \
    + [complex(c['layerNK'][k][0][0], c['layerNK'][k][0][1]) for k in range(c['nLayers'])] \
    + [complex(c['ns'][0][0], c['ns'][0][1])]
d_list = [np.inf] + list(c['thick']) + [np.inf]
res = tmm.coh_tmm('s', n_list, d_list, 0.0, c['lambdas'][0])
print("\ncoh_tmm returns keys:", sorted(res.keys()))
