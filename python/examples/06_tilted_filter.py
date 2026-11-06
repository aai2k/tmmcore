# %% [markdown]
# # Tuning a narrowband filter by tilting it
#
# A first-order single-cavity filter: quarter-wave reflectors (HL)⁴ and (LH)⁴
# either side of a half-wave high-index spacer, 17 layers on glass, the type
# of Macleod's Figure 8.8 (*Thin-Film Optical Filters* 5th ed., §8.2.2).
# Tilting it moves the passband to shorter wavelengths. Macleod puts the peak
# at λ0 √(1 − sin²θ / n*²), with n* = √(nH·nL) for this filter (§8.2.5.1,
# Eq. 8.36).
#
# At 0° every layer is a quarter or half wave at λ0, and the filter transmits
# exactly what bare glass does, 4·n0·ns / (n0 + ns)² (Eq. 2.53).
#
#     pip install tmmcore matplotlib
#     python 06_tilted_filter.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore

nH, nL, ns = 2.35, 1.46, 1.52     # TiO2-like, SiO2-like, glass; no dispersion
lam0 = 550.0                      # nm

H = {"n": nH, "d": lam0 / (4 * nH)}
L = {"n": nL, "d": lam0 / (4 * nL)}
spacer = {"n": nH, "d": lam0 / (2 * nH)}
filt = {"incident": 1.0, "substrate": ns, "layers": [H, L] * 4 + [spacer] + [L, H] * 4}

# %%
T0 = tmmcore.rta(filt, lam0)["T"]
bare = 4 * ns / (1 + ns) ** 2
print(f"T at {lam0:.0f} nm and 0°: {T0:.15f}, bare glass {bare:.15f}, difference {abs(T0 - bare):.1e}")
assert abs(T0 - bare) < 1e-14

lam = np.arange(500.0, 560.0, 0.005)
T = tmmcore.rta_spectrum(filt, lam)["T"]
print(f"halfwidth at 0°: {np.count_nonzero(T > T.max() / 2) * 0.005:.1f} nm")

n_star = np.sqrt(nH * nL)
print("\n  AOI   peak s (nm)   peak p (nm)   formula (nm)")
for aoi in (0, 10, 20, 30):
    out = tmmcore.rta_spectrum(filt, lam, angle=aoi)
    formula = lam0 * np.sqrt(1 - (np.sin(np.radians(aoi)) / n_star) ** 2)
    print(f"  {aoi:2d}°   {lam[out['Ts'].argmax()]:11.2f}   {lam[out['Tp'].argmax()]:11.2f}   {formula:12.2f}")

# %%
for aoi in (0, 15, 30):
    plt.plot(lam, 100 * tmmcore.rta_spectrum(filt, lam, angle=aoi)["T"], label=f"{aoi}°, unpolarized")
plt.xlim(520, 556)
plt.xlabel("wavelength (nm)")
plt.ylabel("transmittance (%)")
plt.legend()
plt.show()
