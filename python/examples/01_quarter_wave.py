# %% [markdown]
# # A quarter-wave layer on glass
#
# One layer of MgF2, a quarter wave thick at 550 nm, on glass at normal
# incidence. At 550 nm its reflectance has a closed form (Macleod, *Thin-Film
# Optical Filters* 5th ed., Eq. 4.5):
#
#     R = ((n0·ns − n1²) / (n0·ns + n1²))²
#
#     pip install tmmcore matplotlib
#     python 01_quarter_wave.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore

n0, n1, ns = 1.0, 1.38, 1.52      # air, MgF2, glass
lam0 = 550.0                      # nm
d = lam0 / (4 * n1)               # quarter wave at lam0, nm

coated = {"incident": n0, "substrate": ns, "layers": [{"n": n1, "d": d}]}
bare = {"incident": n0, "substrate": ns, "layers": []}

# %%
R = tmmcore.rta(coated, lam0)["R"]
exact = ((n0 * ns - n1**2) / (n0 * ns + n1**2)) ** 2

print(f"MgF2 thickness   {d:.2f} nm")
print(f"R, tmmcore       {R:.17f}")
print(f"R, Eq. 4.5       {exact:.17f}")
print(f"difference       {abs(R - exact):.1e}")
print(f"bare glass R     {tmmcore.rta(bare, lam0)['R']:.2%}, coated {R:.2%}")
assert abs(R - exact) < 1e-15

# %%
lam = np.arange(400.0, 801.0, 1.0)
plt.plot(lam, 100 * tmmcore.rta_spectrum(bare, lam)["R"], "--", label="bare glass")
plt.plot(lam, 100 * tmmcore.rta_spectrum(coated, lam)["R"], label="quarter-wave MgF2")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.ylim(0, 5)
plt.legend()
plt.show()
