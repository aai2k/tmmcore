"""
Exact analytic derivatives of R, T and A with respect to the layer
thicknesses: the Jacobian, the Hessian, and the needle scan, at one
wavelength and over a grid. Derivatives are per nm.
"""

import numpy as np

from . import _stack
from ._library import lib
from ._stack import address


def _derivatives(stack, wavelength, angle, pol, second):
    lam = _stack.read_wavelength(wavelength)
    theta = _stack.read_angle(angle)
    pol = _stack.read_rta_pol(pol)
    n0, ns, layers = _stack.point_stack(stack)
    N = len(layers)
    keys = ["dRdd", "dTdd", "dAdd"] + (["d2Rdd", "d2Tdd", "d2Add"] if second else [])

    def run(side):
        outs = [np.empty(N) for _ in range(3)] + [np.empty((N, N)) for _ in keys[3:]]
        base = np.empty(3)
        kernel = lib.tmm_hessian if second else lib.tmm_jacobian
        kernel(lam, theta, _stack.POL_CODE[side], n0.real, n0.imag, ns.real, ns.imag,
               address(layers), N, *map(address, outs), address(base))
        return {"R": float(base[0]), "T": float(base[1]), "A": float(base[2]), **dict(zip(keys, outs)), "N": N}

    return _stack.polarized(pol, run)


def _derivatives_spectrum(stack, wavelengths, angle, pol, second):
    lambdas = _stack.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles = _stack.read_angle(angle, n_lam)
    pol = _stack.read_rta_pol(pol)
    n0, ns, table, thick = _stack.grid_stack(stack, n_lam)
    N = len(thick)
    keys = ["dRdd", "dTdd", "dAdd"] + (["d2Rdd", "d2Tdd", "d2Add"] if second else [])

    def run(side):
        outs = [np.empty((n_lam, N)) for _ in range(3)] + [np.empty((n_lam, N, N)) for _ in keys[3:]]
        base = np.empty((n_lam, 3))
        kernel = lib.tmm_hessian_spectrum if second else lib.tmm_jacobian_spectrum
        kernel(address(lambdas), n_lam, address(angles), len(angles), _stack.POL_CODE[side],
               address(n0), address(ns), address(table), address(thick), N,
               *map(address, outs), address(base))
        rta = {"R": base[:, 0].copy(), "T": base[:, 1].copy(), "A": base[:, 2].copy()}
        return {**rta, **dict(zip(keys, outs)), "N": N}

    return _stack.polarized(pol, run)


def jacobian(stack, wavelength, *, angle=0.0, pol="avg"):
    """
    R, T and A at one wavelength with their thickness derivatives.

    Takes the arguments of rta. Returns {"R", "T", "A", "N"} and "dRdd",
    "dTdd", "dAdd", each an array of N derivatives, one per layer. A layer of
    zero thickness keeps its place and its derivative.
    """
    return _derivatives(stack, wavelength, angle, pol, second=False)


def jacobian_spectrum(stack, wavelengths, *, angle=0.0, pol="avg"):
    """
    The Jacobian over a wavelength grid. Takes the arguments of rta_spectrum.
    Returns R, T and A over the grid and the derivatives as (nLam, N) arrays.
    """
    return _derivatives_spectrum(stack, wavelengths, angle, pol, second=False)


def hessian(stack, wavelength, *, angle=0.0, pol="avg"):
    """
    The Jacobian at one wavelength, and the second derivatives "d2Rdd",
    "d2Tdd", "d2Add", each a symmetric (N, N) array.
    """
    return _derivatives(stack, wavelength, angle, pol, second=True)


def hessian_spectrum(stack, wavelengths, *, angle=0.0, pol="avg"):
    """The Hessian over a wavelength grid: first derivatives (nLam, N), second (nLam, N, N)."""
    return _derivatives_spectrum(stack, wavelengths, angle, pol, second=True)


def _needle_scan_kernel(point, pol_code, candidates, fractions):
    """
    The needle scan at one wavelength on the kernel's own arrays. point is
    (λ in nm, θ in degrees, ñ0, ñs, layers as N rows [n, k, d]), candidates
    are complex indices, fractions split points inside each layer. Returns
    {"base": [R, T, A], "gaps": (N + 1, nCand, 3), "intra": (N, nFrac, nCand, 3)},
    the last axis [dR, dT, dA].
    """
    lambda_nm, theta_deg, n0, ns, layers = point
    layers = np.ascontiguousarray(layers, dtype=np.float64)
    candidates = np.ascontiguousarray(candidates, dtype=np.complex128)
    fractions = np.ascontiguousarray(fractions, dtype=np.float64)
    N, n_cand, n_frac = len(layers), len(candidates), len(fractions)
    base = np.empty(3)
    gaps = np.empty((N + 1, n_cand, 3))
    intra = np.empty((N, n_frac, n_cand, 3))
    lib.tmm_needle_scan(float(lambda_nm), float(theta_deg), int(pol_code), n0.real, n0.imag, ns.real, ns.imag,
                        address(layers), N, address(candidates), n_cand, address(fractions), n_frac,
                        address(base), address(gaps), address(intra))
    return {"base": base, "gaps": gaps, "intra": intra}


def _needle_scan_spectrum_kernel(grid, pol_code, candidates, fractions):
    """
    The needle scan over a grid on the kernel's own arrays. grid is (λs in
    nm, θs in degrees, then the media, the layer table and the thicknesses as
    _stack.grid_stack gives them); candidates are (nCand, nLam) complex.
    Returns {"base": (nLam, 3), "gaps": (nLam, N + 1, nCand, 3),
    "intra": (nLam, N, nFrac, nCand, 3)}.
    """
    lambdas, angles, n0, ns, table, thick = grid
    candidates = np.ascontiguousarray(candidates, dtype=np.complex128)
    fractions = np.ascontiguousarray(fractions, dtype=np.float64)
    n_lam, N, n_cand, n_frac = len(lambdas), len(thick), len(candidates), len(fractions)
    base = np.empty((n_lam, 3))
    gaps = np.empty((n_lam, N + 1, n_cand, 3))
    intra = np.empty((n_lam, N, n_frac, n_cand, 3))
    lib.tmm_needle_scan_spectrum(address(lambdas), n_lam, address(angles), len(angles), int(pol_code),
                                 address(n0), address(ns), address(table), address(thick), N,
                                 address(candidates), n_cand, address(fractions), n_frac,
                                 address(base), address(gaps), address(intra))
    return {"base": base, "gaps": gaps, "intra": intra}
