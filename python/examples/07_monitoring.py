# %% [markdown]
# # The monitoring trace of a quarter-wave mirror
#
# A nine-layer quarter-wave mirror at 600 nm, deposited layer by layer on
# glass, watched in reflection at 600 nm. Each quarter wave ends at a turning
# point of the trace, which is how turning-point monitoring knows when to
# stop a layer (Macleod, *Thin-Film Optical Filters* 5th ed., ch. 13).
#
# At a turning point the reflectance has a closed form. A quarter-wave layer
# of index n turns the admittance Y under it into n²/Y (the characteristic
# matrix of Macleod ch. 2 at a phase thickness of π/2), so the surface
# admittance after each layer follows from the substrate's, and
# R = ((n0 − Y)/(n0 + Y))².
#
#     pip install tmmcore matplotlib
#     python 07_monitoring.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore

n0, ns = 1.0, 1.52                # air, glass
lam0 = 600.0                      # design and monitoring wavelength, nm
tio2 = {"name": "TiO2", "n": 2.35}
sio2 = {"name": "SiO2", "n": 1.46}
design = [tio2, sio2] * 4 + [tio2]   # from the air side

# %%
deposited = []                    # finished layers, from the air side
Y = ns
total = 0.0
trace_d, trace_R = [], []

print("layer   material   d (nm)   R at the end, tmmcore   closed form")
for k, material in enumerate(reversed(design), start=1):
    d = lam0 / (4 * material["n"])
    t = np.linspace(0.0, d, 201)
    out = tmmcore.monitor_curve({"incident": n0, "substrate": ns, "layers": deposited}, lam0,
                                top=material["n"], thicknesses=t)
    trace_d.append(total + t)
    trace_R.append(out["Rs"])

    Y = material["n"] ** 2 / Y
    closed = ((n0 - Y) / (n0 + Y)) ** 2
    print(f"  {k}      {material['name']}     {d:6.2f}   {out['Rs'][-1]:.15f}     {closed:.15f}")
    assert abs(out["Rs"][-1] - closed) < 1e-14

    deposited.insert(0, {"n": material["n"], "d": d})
    total += d

# %%
plt.plot(np.concatenate(trace_d), 100 * np.concatenate(trace_R))
plt.xlabel("deposited thickness (nm)")
plt.ylabel(f"reflectance at {lam0:.0f} nm (%)")
plt.show()
