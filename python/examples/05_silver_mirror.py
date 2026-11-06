# %% [markdown]
# # A silver fold mirror at 45°
#
# 150 nm of silver on glass, used at 45°. The silver indices are measured
# data, Johnson and Christy, *Phys. Rev. B* **6**, 4370 (1972), from the
# refractiveindex.info database, interpolated linearly onto the wavelength
# grid. Give any material this way: a table of n and k, interpolated onto the
# wavelengths you compute at.
#
#     pip install tmmcore matplotlib
#     python 05_silver_mirror.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore

#                 λ (nm)    n      k
silver = np.array([[381.5, 0.05, 1.864], [397.4, 0.05, 2.070], [413.3, 0.05, 2.275],
                   [430.5, 0.04, 2.462], [450.9, 0.04, 2.657], [471.4, 0.05, 2.869],
                   [495.9, 0.05, 3.093], [520.9, 0.05, 3.324], [548.6, 0.06, 3.586],
                   [582.1, 0.05, 3.858], [616.8, 0.06, 4.152], [659.5, 0.05, 4.483],
                   [704.5, 0.04, 4.838], [756.0, 0.03, 5.242], [821.1, 0.04, 5.727]])

lam = np.arange(400.0, 801.0, 1.0)
n_ag = np.interp(lam, silver[:, 0], silver[:, 1])
k_ag = np.interp(lam, silver[:, 0], silver[:, 2])

mirror = {"incident": 1.0, "substrate": 1.52, "layers": [{"n": n_ag, "k": k_ag, "d": 150.0}]}
at45 = tmmcore.rta_spectrum(mirror, lam, angle=45)
at0 = tmmcore.rta_spectrum(mirror, lam)

# %%
print("  λ (nm)    R 0°    Rs 45°   Rp 45°   A 45°, unpolarized")
for wl in (450, 550, 650, 750):
    i = np.searchsorted(lam, wl)
    print(f"  {wl:6.0f}  {at0['R'][i]:6.2%}   {at45['Rs'][i]:6.2%}   {at45['Rp'][i]:6.2%}   {at45['A'][i]:6.2%}")
print(f"\nlargest T through 150 nm of silver: {at0['T'].max():.1e}")

# Silver as the substrate, with no film, against the reflectance of one
# absorbing face at normal incidence, |(n0 − ñ)/(n0 + ñ)|² with ñ = n + ik
# (Macleod Eqs. 2.75 and 2.79).
bulk = tmmcore.rta_spectrum({"incident": 1.0, "substrate": {"n": n_ag, "k": k_ag}, "layers": []}, lam)["R"]
fresnel = np.abs((1 - (n_ag + 1j * k_ag)) / (1 + (n_ag + 1j * k_ag))) ** 2
print(f"bulk silver against the closed form, largest difference {np.abs(bulk - fresnel).max():.1e}")
assert np.abs(bulk - fresnel).max() < 1e-14

# %%
plt.plot(lam, 100 * at45["Rs"], label="Rs, 45°")
plt.plot(lam, 100 * at45["Rp"], label="Rp, 45°")
plt.plot(lam, 100 * at0["R"], "--", label="R, 0°")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.legend()
plt.show()
