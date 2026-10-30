"""
Phase, group delay (GD), group delay dispersion (GDD) and third-order
dispersion (TOD) of the reflected and transmitted amplitudes, and their
thickness derivatives, at one wavelength and over a grid.

The quantities are derivatives in angular frequency ω, carried exactly
through the stack as Taylor jets (jet.py). The phase is in radians with
Macleod's sign convention; GD, GDD and TOD are in fs, fs² and fs³.
"""

import math

import numpy as np

from . import _options, _stack
from ._library import lib
from ._stack import address

_QUANTITIES = ["phase_rad", "gd", "gdd", "tod", "magnitude_squared"]
_DERIVATIVES = ["d_phase_deg", "d_gd", "d_gdd", "d_tod", "d_log_magnitude_squared"]


def omega_from_lambda_nm(lambda_nm):
    """Angular frequency in rad/fs of a vacuum wavelength in nm, ω = 2πc/λ."""
    lambdas = np.array(lambda_nm, dtype=np.float64, order="C")
    out = np.empty(lambdas.shape)
    lib.tmm_omega_from_lambda_nm(address(lambdas), lambdas.size, address(out))
    return float(out) if out.ndim == 0 else out


# One side of a point result; the kernel writes NaN where the coefficient is
# exactly zero and has no phase.
def _side(values):
    if math.isnan(values[4]):
        return None
    side = {key: float(value) for key, value in zip(_QUANTITIES, values)}
    return {"phase_rad": side["phase_rad"], "phase_deg": side["phase_rad"] * 180 / math.pi, **side}


def phase(stack, wavelength, *, pol=None, angle=None, sine=None):
    """
    Phase, GD, GDD and TOD of r and t at one wavelength.

    Takes the arguments of rta, with pol "s" or "p", which may be left out at
    normal incidence, and each index as a number or a jet (jet.py). In place
    of angle it takes sine, the sine of the angle of incidence as a jet, for
    an angle held fixed outside a dispersive incident medium. Returns
    {"r", "t"}, each {"phase_rad", "phase_deg", "gd", "gdd", "tod",
    "magnitude_squared"}, or None where the coefficient is exactly zero.
    """
    lam = _options.read_wavelength(wavelength)
    theta, sine, code = _options.read_phase_options(angle, sine, pol)
    n0, ns, layers, thick = _stack.point_jet_stack(stack)
    out = np.empty(10)
    lib.tmm_phase(lam, theta, code, address(n0), address(ns), address(layers), address(thick), len(thick),
                  omega_from_lambda_nm(lam), _stack.optional_address(sine), address(out))
    return {"r": _side(out[:5]), "t": _side(out[5:])}


def phase_jacobian(stack, wavelength, *, pol=None, angle=None, sine=None):
    """
    The phase quantities of phase() at one wavelength with their thickness
    derivatives: "d_phase_deg", "d_gd", "d_gdd", "d_tod" and
    "d_log_magnitude_squared", each an array of N, one per layer.
    """
    lam = _options.read_wavelength(wavelength)
    theta, sine, code = _options.read_phase_options(angle, sine, pol)
    n0, ns, layers, thick = _stack.point_jet_stack(stack)
    out = np.empty((2, 5))
    deriv = np.empty((2, 5, len(thick)))
    lib.tmm_phase_jacobian(lam, theta, code, address(n0), address(ns), address(layers), address(thick),
                           len(thick), omega_from_lambda_nm(lam), _stack.optional_address(sine), address(out),
                           address(deriv))

    def side(index):
        base = _side(out[index])
        if base is None:
            return None
        return {**base, **{key: deriv[index, q].copy() for q, key in enumerate(_DERIVATIVES)}}

    return {"r": side(0), "t": side(1)}


def _phase_grid(stack, wavelengths, options, with_derivatives):
    kernel = lib.tmm_phase_jacobian_spectrum if with_derivatives else lib.tmm_phase_spectrum
    lambdas = _options.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles, sines, code = _options.read_phase_options(*options, n_lam)
    n0, ns, table, thick = _stack.grid_jet_stack(stack, n_lam)
    omegas = omega_from_lambda_nm(lambdas)
    # Per wavelength the kernel writes the five quantities of r, then of t,
    # and the derivatives as [side][quantity][layer].
    outs = [np.empty((n_lam, 2, 5))] + ([np.empty((n_lam, 2, 5, len(thick)))] if with_derivatives else [])
    kernel(address(lambdas), n_lam, address(angles), len(angles), code, address(n0), address(ns),
           address(table), address(thick), len(thick), address(omegas), _stack.optional_address(sines),
           *map(address, outs))
    values = outs[0]

    def side(index):
        quantities = {key: values[:, index, q].copy() for q, key in enumerate(_QUANTITIES)}
        if with_derivatives:
            quantities.update({key: outs[1][:, index, q].copy() for q, key in enumerate(_DERIVATIVES)})
        return quantities

    return {"r": side(0), "t": side(1)}


def phase_spectrum(stack, wavelengths, *, pol=None, angle=None, sine=None):
    """
    Phase, GD, GDD and TOD over a wavelength grid. Takes the arguments of
    rta_spectrum, with pol as in phase(), each index as a number or one jet
    per wavelength, and sine, if given, as one jet per wavelength. Returns
    {"r", "t"}, each {"phase_rad", "gd", "gdd", "tod", "magnitude_squared"} as
    arrays over the grid, NaN where the coefficient is exactly zero.
    """
    return _phase_grid(stack, wavelengths, (angle, sine, pol), with_derivatives=False)


def phase_jacobian_spectrum(stack, wavelengths, *, pol=None, angle=None, sine=None):
    """
    The phase grid of phase_spectrum with the thickness derivatives of
    phase_jacobian, each derivative a (nLam, N) array.
    """
    return _phase_grid(stack, wavelengths, (angle, sine, pol), with_derivatives=True)
