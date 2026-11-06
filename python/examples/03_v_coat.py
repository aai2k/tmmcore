# %% [markdown]
# # A two-layer V-coat for a 1064 nm laser window
#
# Two layers, SiO2 facing air and Ta2O5 against the glass, give zero
# reflectance at one wavelength. Catalan's closed form (Macleod, *Thin-Film
# Optical Filters* 5th ed., Eq. 4.9) gives the thicknesses directly, and there
# are two designs. tmmcore checks each, and shows how wide its low-reflectance
# band is.
#
#     pip install tmmcore matplotlib
#     python 03_v_coat.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore

lam0 = 1064.0                 # design wavelength, nm
n0, ns = 1.0, 1.507           # air, N-BK7 at 1064 nm
n1, n2 = 1.45, 2.05           # SiO2 facing air, Ta2O5 against the glass


def stack(d1, d2):
    return {"incident": n0, "substrate": ns,
            "layers": [{"n": n1, "d": d1}, {"n": n2, "d": d2}]}


# %% [markdown]
# Normal incidence, no loss. Eq. 4.9 gives tan² of each phase thickness, so
# tan δ2 is either root; tan δ1 · tan δ2, from the same derivation, then fixes
# δ1. Each root is one design.

# %%
tan_d2 = np.sqrt((ns - n0) * (n0 * ns - n1**2) * n2**2 / ((n1**2 * ns - n0 * n2**2) * (n2**2 - n0 * ns)))
product = n1 * n2 * (ns - n0) / (n1**2 * ns - n0 * n2**2)

designs = {}
for name, t2 in (("design 1", tan_d2), ("design 2", -tan_d2)):
    delta1, delta2 = np.mod(np.arctan([product / t2, t2]), np.pi)
    designs[name] = (delta1 * lam0 / (2 * np.pi * n1), delta2 * lam0 / (2 * np.pi * n2))

lam = np.arange(900.0, 1250.0, 0.1)
print("          SiO2 (nm)  Ta2O5 (nm)   R at 1064 nm   R < 0.25 %")
for name, (d1, d2) in designs.items():
    R0 = tmmcore.rta(stack(d1, d2), lam0)["R"]
    R = tmmcore.rta_spectrum(stack(d1, d2), lam)["R"]
    band = lam[R < 0.0025]
    print(f"{name}  {d1:9.2f}  {d2:10.2f}   {R0:12.1e}   {band[0]:.0f}-{band[-1]:.0f} nm")
    # r is zero to rounding, about 1e-16, so R = |r|² falls below 1e-30.
    assert R0 < 1e-30

print(f"\nbare glass reflects {tmmcore.rta(stack(0, 0), lam0)['R']:.2%}")

# %%
for name, (d1, d2) in designs.items():
    R = tmmcore.rta_spectrum(stack(d1, d2), lam)["R"]
    plt.plot(lam, 100 * R, label=f"{name}: SiO2 {d1:.1f} nm, Ta2O5 {d2:.1f} nm")
plt.axhline(0.25, color="gray", linestyle=":", label="0.25 %")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.ylim(0, 2)
plt.legend(loc="upper right")
plt.show()
