# %% [markdown]
# # The spectrum of an antireflection coating
#
# A quarter-half-quarter antireflection coating at 510 nm on BK7-like glass,
# with every index given by a Cauchy formula, so each one is a list with one
# value per wavelength (Macleod, *Thin-Film Optical Filters* 5th ed., ch. 4).
# The Cauchy coefficients are typical, not measured data.
#
#     pip install tmmcore matplotlib
#     python 02_ar_spectrum.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore


def cauchy(a, b):
    """n(λ) = a + b / λ², λ in nm."""
    return lambda lam: a + b / lam**2


low = cauchy(1.377, 2000)         # MgF2-like
high = cauchy(2.170, 8000)        # ZrO2-like
mid = cauchy(1.680, 5000)         # Al2O3-like
glass = cauchy(1.504, 4200)       # BK7-like

lam0 = 510.0                      # nm
lam = np.arange(400.0, 801.0, 1.0)

# Layers from the air side: a quarter wave of low, a half wave of high, a
# quarter wave of mid. Thicknesses are set at lam0 and stay fixed.
layers = [
    {"n": low(lam), "d": 0.25 * lam0 / low(lam0)},
    {"n": high(lam), "d": 0.5 * lam0 / high(lam0)},
    {"n": mid(lam), "d": 0.25 * lam0 / mid(lam0)},
]
coated = tmmcore.rta_spectrum({"incident": 1.0, "substrate": {"n": glass(lam)}, "layers": layers}, lam)
bare = tmmcore.rta_spectrum({"incident": 1.0, "substrate": {"n": glass(lam)}, "layers": []}, lam)

# %%
print("thicknesses, nm: " + ", ".join(f"{layer['d']:.2f}" for layer in layers))
print("\n  λ (nm)   R coated   T coated   R bare")
for i in range(0, len(lam), 50):
    print(f"  {lam[i]:6.0f}   {coated['R'][i]:7.3%}   {coated['T'][i]:7.2%}   {bare['R'][i]:6.2%}")

band = lam <= 700
print(f"\nmean R over 400-700 nm: coated {coated['R'][band].mean():.2%}, bare {bare['R'][band].mean():.2%}")

# Bare glass against the Fresnel reflectance of one face, Eq. 2.53.
fresnel = ((1.0 - glass(lam)) / (1.0 + glass(lam))) ** 2
print(f"bare glass against Fresnel, largest difference {np.abs(bare['R'] - fresnel).max():.1e}")
assert np.abs(bare["R"] - fresnel).max() < 1e-15

# %%
plt.plot(lam, 100 * bare["R"], "--", label="bare glass")
plt.plot(lam, 100 * coated["R"], label="coated")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.ylim(0, 5)
plt.legend()
plt.show()
