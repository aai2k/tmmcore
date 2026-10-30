"""
Exact analytic derivatives of R, T and A with respect to the layer
thicknesses: the Jacobian, the Hessian, and the needle scan, at one
wavelength and over a grid. Derivatives are per nm.
"""

import numpy as np

from . import _options, _stack
from ._library import lib
from ._stack import address


def _derivatives(stack, wavelength, angle, pol, second):
    lam = _options.read_wavelength(wavelength)
    theta = _options.read_angle(angle)
    pol = _options.read_rta_pol(pol)
    n0, ns, layers = _stack.point_stack(stack)
    N = len(layers)
    keys = ["dRdd", "dTdd", "dAdd"] + (["d2Rdd", "d2Tdd", "d2Add"] if second else [])

    def run(side):
        outs = [np.empty(N) for _ in range(3)] + [np.empty((N, N)) for _ in keys[3:]]
        base = np.empty(3)
        kernel = lib.tmm_hessian if second else lib.tmm_jacobian
        kernel(lam, theta, _options.POL_CODE[side], n0.real, n0.imag, ns.real, ns.imag,
               address(layers), N, *map(address, outs), address(base))
        return {"R": float(base[0]), "T": float(base[1]), "A": float(base[2]), **dict(zip(keys, outs)), "N": N}

    return _options.polarized(pol, run)


def _derivatives_spectrum(stack, wavelengths, angle, pol, second):
    lambdas = _options.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles = _options.read_angle(angle, n_lam)
    pol = _options.read_rta_pol(pol)
    n0, ns, table, thick = _stack.grid_stack(stack, n_lam)
    N = len(thick)
    keys = ["dRdd", "dTdd", "dAdd"] + (["d2Rdd", "d2Tdd", "d2Add"] if second else [])

    def run(side):
        outs = [np.empty((n_lam, N)) for _ in range(3)] + [np.empty((n_lam, N, N)) for _ in keys[3:]]
        base = np.empty((n_lam, 3))
        kernel = lib.tmm_hessian_spectrum if second else lib.tmm_jacobian_spectrum
        kernel(address(lambdas), n_lam, address(angles), len(angles), _options.POL_CODE[side],
               address(n0), address(ns), address(table), address(thick), N,
               *map(address, outs), address(base))
        rta = {"R": base[:, 0].copy(), "T": base[:, 1].copy(), "A": base[:, 2].copy()}
        return {**rta, **dict(zip(keys, outs)), "N": N}

    return _options.polarized(pol, run)


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


def _needle_triples(values):
    """[dR, dT, dA] on the last axis, as three arrays."""
    return {"dR": values[..., 0].copy(), "dT": values[..., 1].copy(), "dA": values[..., 2].copy()}


def needle_scan(stack, wavelength, *, candidates, fractions=(), angle=0.0, pol="avg"):
    """
    The derivatives of R, T and A, per nm of needle, for inserting a layer of
    zero thickness of each candidate material at every interface of the stack
    and, with fractions, inside its layers.

    Takes the arguments of rta, and candidates, a list of indices, each a
    number or {"n", "k"} as a medium is, and fractions, positions inside each
    layer as fractions of its thickness from the incident side. Returns
    {"R", "T", "A", "gaps", "intra", "N"}. "gaps" holds "dR", "dT" and "dA",
    each (N + 1, nCand), interface 0 in front of the first layer and N against
    the substrate; "intra" holds the same, each (N, nFrac, nCand).
    """
    lam = _options.read_wavelength(wavelength)
    theta = _options.read_angle(angle)
    pol = _options.read_rta_pol(pol)
    candidates = _stack.read_candidates(candidates)
    fractions = _options.read_fractions(fractions)
    n0, ns, layers = _stack.point_stack(stack)
    N, n_cand, n_frac = len(layers), len(candidates), len(fractions)

    def run(side):
        base = np.empty(3)
        gaps = np.empty((N + 1, n_cand, 3))
        intra = np.empty((N, n_frac, n_cand, 3))
        lib.tmm_needle_scan(lam, theta, _options.POL_CODE[side], n0.real, n0.imag, ns.real, ns.imag,
                            address(layers), N, address(candidates), n_cand, address(fractions), n_frac,
                            address(base), address(gaps), address(intra))
        return {"R": float(base[0]), "T": float(base[1]), "A": float(base[2]),
                "gaps": _needle_triples(gaps), "intra": _needle_triples(intra), "N": N}

    return _options.polarized(pol, run)


def needle_scan_spectrum(stack, wavelengths, *, candidates, fractions=(), angle=0.0, pol="avg"):
    """
    The needle scan over a wavelength grid. Takes the arguments of
    rta_spectrum, and candidates and fractions as needle_scan does, where a
    candidate's n and k may also be one value per wavelength. Returns R, T and
    A over the grid, and "gaps" and "intra" with a leading wavelength axis,
    (nLam, N + 1, nCand) and (nLam, N, nFrac, nCand).
    """
    lambdas = _options.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles = _options.read_angle(angle, n_lam)
    pol = _options.read_rta_pol(pol)
    candidates = _stack.read_candidates(candidates, n_lam)
    fractions = _options.read_fractions(fractions)
    n0, ns, table, thick = _stack.grid_stack(stack, n_lam)
    N, n_cand, n_frac = len(thick), len(candidates), len(fractions)

    def run(side):
        base = np.empty((n_lam, 3))
        gaps = np.empty((n_lam, N + 1, n_cand, 3))
        intra = np.empty((n_lam, N, n_frac, n_cand, 3))
        lib.tmm_needle_scan_spectrum(address(lambdas), n_lam, address(angles), len(angles), _options.POL_CODE[side],
                                     address(n0), address(ns), address(table), address(thick), N,
                                     address(candidates), n_cand, address(fractions), n_frac,
                                     address(base), address(gaps), address(intra))
        return {"R": base[:, 0].copy(), "T": base[:, 1].copy(), "A": base[:, 2].copy(),
                "gaps": _needle_triples(gaps), "intra": _needle_triples(intra), "N": N}

    return _options.polarized(pol, run)
