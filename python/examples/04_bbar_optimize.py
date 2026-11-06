# %% [markdown]
# # Optimizing a broadband antireflection coating
#
# Three layers on glass, with the mean reflectance over 450-650 nm minimized
# by scipy. `jacobian_spectrum` returns R and its derivative in every layer's
# thickness at every wavelength in one call, so each step of the optimizer
# costs one call for the merit and its exact gradient.
#
#     pip install tmmcore scipy matplotlib
#     python 04_bbar_optimize.py

# %%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import tmmcore

# Cauchy indices n = a + b/λ², λ in nm, from the air side: fused silica,
# TiO2, Al2O3. Typical coefficients, not measured data.
materials = [(1.458, 3540.0), (2.20, 24000.0), (1.76, 5500.0)]
lam = np.arange(450.0, 651.0, 1.0)


def stack(d, lam):
    layers = [{"n": a + b / lam**2, "d": t} for (a, b), t in zip(materials, d)]
    return {"incident": 1.0, "substrate": 1.52, "layers": layers}


def merit(d):
    """Mean R over the band, and its gradient in the thicknesses, per nm."""
    out = tmmcore.jacobian_spectrum(stack(d, lam), lam)
    return out["R"].mean(), out["dRdd"].mean(axis=0)


# %%
start = np.array([100.0, 100.0, 100.0])
# scipy's default gtol stops about 0.1 nm short of the minimum here; 1e-8
# brings the thicknesses to well within the 0.01 nm printed.
res = minimize(merit, start, jac=True, method="L-BFGS-B", options={"gtol": 1e-8})
R_refined = tmmcore.rta_spectrum(stack(res.x, lam), lam)["R"]

print("          SiO2 (nm)   TiO2 (nm)  Al2O3 (nm)   mean R")
print(f"start     {start[0]:9.2f}  {start[1]:10.2f}  {start[2]:10.2f}   {merit(start)[0]:.3%}")
print(f"refined   {res.x[0]:9.2f}  {res.x[1]:10.2f}  {res.x[2]:10.2f}   {res.fun:.3%}")
print(f"\nlargest R over 450-650 nm: {R_refined.max():.3%}; {res.nfev} calls to merit()")

# %% [markdown]
# The gradient scipy used is exact. Against a central difference at the start:

# %%
h = 1e-3                                   # nm
F, grad = merit(start)
central = np.array([(merit(start + h * e)[0] - merit(start - h * e)[0]) / (2 * h) for e in np.eye(3)])
rel = np.abs(grad - central) / np.abs(grad)

print("layer    analytic      central diff   rel. difference")
for j in range(3):
    print(f"  {j + 1}    {grad[j]: .6e}   {central[j]: .6e}   {rel[j]:.1e}")
assert rel.max() < 1e-8

# %%
plot_lam = np.arange(400.0, 701.0, 1.0)
plt.plot(plot_lam, 100 * tmmcore.rta_spectrum(stack(start, plot_lam), plot_lam)["R"], "--", label="start")
plt.plot(plot_lam, 100 * tmmcore.rta_spectrum(stack(res.x, plot_lam), plot_lam)["R"], label="refined")
plt.axvspan(450, 650, color="green", alpha=0.1, label="optimized band")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.ylim(0, None)
plt.legend()
plt.show()
