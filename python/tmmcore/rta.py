"""Reflectance, transmittance and absorptance of a stack, at one wavelength and over a grid."""

import numpy as np

from . import _stack
from ._library import lib
from ._stack import address


def rta(stack, wavelength, *, angle=0.0, pol="avg"):
    """
    R, T and A at one wavelength.

    stack       {"incident", "substrate", "layers": [{"n", "k", "d"}]}, d in nm
    wavelength  in nm
    angle       angle of incidence in degrees from the normal
    pol         "s", "p", or "avg" for the mean of the two

    Returns {"R", "T", "A"}.
    """
    lam = _stack.read_wavelength(wavelength)
    theta = _stack.read_angle(angle)
    pol = _stack.read_rta_pol(pol)
    n0, ns, layers = _stack.point_stack(stack)

    def run(side):
        out = np.empty(3)
        lib.tmm_rta(lam, theta, _stack.POL_CODE[side], n0.real, n0.imag, ns.real, ns.imag,
                    address(layers), len(layers), address(out))
        return {"R": float(out[0]), "T": float(out[1]), "A": float(out[2])}

    return _stack.polarized(pol, run)


def rta_spectrum(stack, wavelengths, *, angle=0.0, pol="avg"):
    """
    R, T and A over a wavelength grid.

    Takes the arguments of rta, with the wavelengths as a list or array, and
    the angle as one number or one per wavelength.

    Returns {"R", "T", "A"} for the polarization asked and {"Rs", "Ts", "As",
    "Rp", "Tp", "Ap"} for s and p, each an array over the wavelengths.
    """
    lambdas = _stack.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles = _stack.read_angle(angle, n_lam)
    pol = _stack.read_rta_pol(pol)
    n0, ns, table, thick = _stack.grid_stack(stack, n_lam)
    outs = [np.empty(n_lam) for _ in range(6)]
    lib.tmm_rta_spectrum(address(lambdas), n_lam, address(angles), len(angles), address(n0), address(ns),
                         address(table), address(thick), len(thick), *map(address, outs))
    both = dict(zip(["Rs", "Ts", "As", "Rp", "Tp", "Ap"], outs))

    def pick(quantity):
        s, p = both[quantity + "s"], both[quantity + "p"]
        return _stack.mean(s, p) if pol == "avg" else (s if pol == "s" else p).copy()

    return {"R": pick("R"), "T": pick("T"), "A": pick("A"), **both}
