"""Probe the layout conventions of tmm_faster and tmm_fast against a known case.

Both take a full complex-index matrix, but neither documents whether the
semi-infinite bounding media are included or how thicknesses line up. Compare
against Byrnes on the same stack to pin it down before benchmarking.
"""

import json

import numpy as np
import tmm

case = json.load(open("cases.json"))["cases"][0]     # AR4/g71
lam = np.array(case["lambdas"])
nL = case["nLayers"]
thick = np.array(case["thick"])

# Reference: Byrnes, first wavelength only.
i = 0
n_list = ([complex(*case["n0"][i])]
          + [complex(*case["layerNK"][k][i]) for k in range(nL)]
          + [complex(*case["ns"][i])])
d_list = [np.inf] + list(thick) + [np.inf]
ref_s = tmm.coh_tmm('s', n_list, d_list, 0.0, lam[i])
print(f"byrnes   Rs={ref_s['R']:.12f}  Ts={ref_s['T']:.12f}")

# Full [nLam x (nL+2)] index matrix including both semi-infinite media.
nmat = np.empty((len(lam), nL + 2), dtype=complex)
for j in range(len(lam)):
    nmat[j, 0] = complex(*case["n0"][j])
    for k in range(nL):
        nmat[j, k + 1] = complex(*case["layerNK"][k][j])
    nmat[j, -1] = complex(*case["ns"][j])

# --- tmm_faster ---------------------------------------------------------
import tmm_faster
for d_variant, label in [
    (np.concatenate([[0.0], thick, [0.0]]), "d padded with 0"),
    (thick, "d = layers only"),
]:
    try:
        r = tmm_faster.calc_coherent(nmat, d_variant, np.array([0.0]), lam)
        print(f"tmm_faster [{label}]  Rs={np.asarray(r['R_s']).ravel()[0]:.12f} "
              f" Ts={np.asarray(r['T_s']).ravel()[0]:.12f}  shape={np.asarray(r['R_s']).shape}")
    except Exception as e:
        print(f"tmm_faster [{label}] -> {type(e).__name__}: {e}")

# --- tmm_fast -----------------------------------------------------------
import torch
from tmm_fast import coh_tmm as fast_coh_tmm

# tmm_fast wants N [stacks, layers, wavelengths] and T [stacks, layers],
# with inf thickness on the bounding media.
N = torch.zeros((1, nL + 2, len(lam)), dtype=torch.cfloat)
for k in range(nL + 2):
    N[0, k, :] = torch.tensor(nmat[:, k], dtype=torch.cfloat)
T = torch.tensor(np.concatenate([[np.inf], thick, [np.inf]]), dtype=torch.float32)[None, :]
th = torch.tensor([0.0], dtype=torch.float32)
wl = torch.tensor(lam, dtype=torch.float32)
try:
    r = fast_coh_tmm('s', N, T, th, wl, device='cpu')
    print(f"tmm_fast   Rs={float(np.asarray(r['R']).ravel()[0]):.12f} "
          f" Ts={float(np.asarray(r['T']).ravel()[0]):.12f}  shape={np.asarray(r['R']).shape}")
except Exception as e:
    print(f"tmm_fast -> {type(e).__name__}: {e}")
