# %% [markdown]
# # What thickness errors do to an antireflection coating
#
# The three-layer coating optimized in example 04, deposited 2000 times with
# random thickness errors: each layer off by a normally distributed fraction
# of its thickness, independently of the others. Macleod finds that
# antireflection coatings usually tolerate errors of around 2% standard
# deviation (*Thin-Film Optical Filters* 5th ed., §13.7). The spread of the
# mean reflectance over 450-650 nm shows what 1% and 2% do to this design.
#
#     pip install tmmcore matplotlib
#     python 08_thickness_errors.py

# %%
import matplotlib.pyplot as plt
import numpy as np

import tmmcore


def cauchy(a, b):
    """n(λ) = a + b / λ², λ in nm."""
    return lambda lam: a + b / lam**2


# The design of example 04. Typical coefficients, not measured data.
sio2 = cauchy(1.458, 3540.0)
tio2 = cauchy(2.20, 24000.0)
al2o3 = cauchy(1.76, 5500.0)
nominal = np.array([94.77, 128.87, 87.22])    # nm, SiO2, TiO2, Al2O3
lam = np.arange(450.0, 651.0, 1.0)


def spectrum(d):
    stack = {
        "incident": 1.0,
        "substrate": 1.52,
        "layers": [
            {"n": sio2(lam), "d": d[0]},     # facing air
            {"n": tio2(lam), "d": d[1]},
            {"n": al2o3(lam), "d": d[2]},    # against the glass
        ],
    }
    return tmmcore.rta_spectrum(stack, lam)["R"]


# %%
rng = np.random.default_rng(1)
runs = 2000
spread = {}

print(f"nominal design: mean R over 450-650 nm {spectrum(nominal).mean():.3%}")
print("\n  σ     median mean R   90% of runs below   99% below")
for sigma in (0.01, 0.02):
    R = np.array([spectrum(nominal * (1 + e)) for e in rng.normal(0.0, sigma, size=(runs, 3))])
    spread[sigma] = R
    p50, p90, p99 = np.percentile(R.mean(axis=1), [50, 90, 99])
    print(f"  {sigma:.0%}   {p50:13.3%}   {p90:17.3%}   {p99:9.3%}")

# %%
plt.plot(lam, 100 * spectrum(nominal), color="black", label="nominal")
for sigma, R in spread.items():
    low, high = np.percentile(R, [10, 90], axis=0)
    plt.fill_between(lam, 100 * low, 100 * high, alpha=0.3, label=f"σ = {sigma:.0%}, 10-90% of runs")
plt.xlabel("wavelength (nm)")
plt.ylabel("reflectance (%)")
plt.ylim(0, None)
plt.legend()
plt.show()
