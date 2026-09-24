"""Emit cases.json — the shared, byte-identical input set for all three kernels.

Everything (wavelength grid, per-layer complex indices, thicknesses) is
pre-computed here so the benchmark times the TMM core alone: no dispersion
evaluation, no material lookup, no unit conversion inside the timed region.

Convention: n~ = n + ik, k >= 0 for absorbing media (matches both tmmcore's
JavaScript implementation and Byrnes' tmm package, so indices transfer verbatim).
"""

import json

# Analytic dispersions, chosen to be representative rather than exact.
def n_high(lam):  return (2.35 + 8000.0 / (lam * lam), 0.0005)   # TiO2-like
def n_low(lam):   return (1.46 + 3000.0 / (lam * lam), 0.0)      # SiO2-like
def n_metal(lam): return (0.15 + 0.0006 * lam, 3.2 + 0.004 * lam)  # Ag-like
def n_glass(lam): return (1.52 + 4200.0 / (lam * lam), 0.0)      # BK7-like
def n_air(lam):   return (1.0, 0.0)


def grid(start, end, step):
    n = int(round((end - start) / step)) + 1
    return [start + i * step for i in range(n)]


# (name, [(dispersion, thickness_nm), ...])
STACKS = {
    # 4-layer broadband AR on glass
    "AR4": [(n_high, 12.0), (n_low, 35.0), (n_high, 118.0), (n_low, 95.0)],
    # 21-layer quarter-wave high reflector at 550 nm
    "HR21": [(n_high if i % 2 == 0 else n_low, 58.5 if i % 2 == 0 else 94.2)
             for i in range(21)],
    # metal-dielectric absorber / beamsplitter (absorbing layers)
    "AG7": [(n_high, 45.0), (n_low, 88.0), (n_metal, 18.0), (n_low, 130.0),
            (n_metal, 12.0), (n_high, 62.0), (n_low, 101.0)],
    # 40-layer stack — upper end of routine design sizes
    "BIG40": [(n_high if i % 2 == 0 else n_low, 61.0 if i % 2 == 0 else 99.0)
              for i in range(40)],
}

GRIDS = {
    "g71":  grid(400.0, 1100.0, 10.0),   # OptoGPT's grid
    "g701": grid(400.0, 1100.0, 1.0),    # typical optical design grid
}

# Normal incidence on both grids. Oblique incidence, where s and p part and the
# p admittance is n / cos(theta), on the coarse grid: the angle changes what is
# checked, not how long a spectrum takes.
ANGLES = {"g71": [0.0, 30.0, 60.0], "g701": [0.0]}


def case(sname, stack, gname, lambdas, theta):
    return {
        "name": f"{sname}/{gname}" + (f"/{theta:g}deg" if theta else ""),
        "stack": sname,
        "grid": gname,
        "nLayers": len(stack),
        "theta_deg": theta,
        "lambdas": lambdas,
        "n0": [list(n_air(l)) for l in lambdas],
        "ns": [list(n_glass(l)) for l in lambdas],
        "thick": [d for (_, d) in stack],
        # layerNK[layer][lambda] = [re, im]
        "layerNK": [[list(f(l)) for l in lambdas] for (f, _) in stack],
    }


cases = [case(sname, stack, gname, lambdas, 0.0)
         for gname, lambdas in GRIDS.items()
         for sname, stack in STACKS.items()]
cases += [case(sname, stack, gname, GRIDS[gname], theta)
          for gname, angles in ANGLES.items()
          for theta in angles if theta
          for sname, stack in STACKS.items()]

with open("cases.json", "w") as fh:
    json.dump({"cases": cases}, fh)

print(f"wrote cases.json: {len(cases)} cases")
for c in cases:
    print(f"  {c['name']:18s}  N={c['nLayers']:3d}  nLam={len(c['lambdas']):4d}")
