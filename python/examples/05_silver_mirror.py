# %% [markdown]
# # A silver fold mirror at 45°
#
# 150 nm of silver on glass, used at 45°. The silver indices are measured
# data, Yang et al., *Phys. Rev. B* **91**, 235137 (2015), from the
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

#                  λ (nm)     n       k
silver = np.array([[399.9, 0.05817, 2.090], [409.9, 0.05647, 2.203], [419.9, 0.05515, 2.310],
                   [429.9, 0.05411, 2.413], [439.9, 0.05332, 2.513], [449.9, 0.05271, 2.609],
                   [459.9, 0.05227, 2.703], [469.9, 0.05198, 2.795], [479.9, 0.05181, 2.886],
                   [489.9, 0.05175, 2.974], [499.9, 0.05180, 3.062], [509.9, 0.05193, 3.148],
                   [519.9, 0.05215, 3.233], [529.9, 0.05245, 3.318], [539.9, 0.05282, 3.401],
                   [549.9, 0.05326, 3.484], [559.9, 0.05376, 3.566], [569.9, 0.05431, 3.648],
                   [579.9, 0.05493, 3.728], [589.9, 0.05560, 3.809], [599.9, 0.05632, 3.889],
                   [609.9, 0.05709, 3.968], [619.9, 0.05790, 4.048], [629.9, 0.05876, 4.126],
                   [639.9, 0.05966, 4.205], [649.9, 0.06061, 4.283], [659.9, 0.06159, 4.361],
                   [669.9, 0.06262, 4.438], [679.9, 0.06368, 4.516], [689.9, 0.06478, 4.593],
                   [699.9, 0.06591, 4.669], [709.9, 0.06708, 4.746], [719.9, 0.06829, 4.823],
                   [729.9, 0.06953, 4.899], [739.9, 0.07080, 4.975], [749.9, 0.07210, 5.051],
                   [759.9, 0.07344, 5.126], [769.9, 0.07480, 5.202], [779.9, 0.07620, 5.278],
                   [789.9, 0.07762, 5.353], [799.9, 0.07908, 5.428], [809.9, 0.08057, 5.503]])

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
