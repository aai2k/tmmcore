"""
A stack whose top layer grows, as in a deposition run: the monitoring curve
at one wavelength, the spectra after each deposited layer, and a growing top
layer sampled over a wavelength grid.

Each returns, per polarization, the forward R and T of the coated surface and
its reflectance from the substrate side, for a monitoring simulation that
treats the witness as an incoherent slab with a bare back face:
{"Rs", "Ts", "Rp", "Tp", "Rrs", "Rrp"}.

The calls here take the kernel's own arrays. At one wavelength, point is
(λ in nm, θ in degrees, ñ0, ñs, layers as N rows [n, k, d]). Over a grid,
grid is (λs in nm, θs in degrees, then the media, the layer table and the
thicknesses as _stack.grid_stack gives them).
"""

import numpy as np

from ._library import lib
from ._stack import address

_KEYS = ["Rs", "Ts", "Rp", "Tp", "Rrs", "Rrp"]


def _monitor_curve_kernel(point, top, thicknesses):
    """
    The signal of a layer of index `top` as it grows to each of `thicknesses`
    (nm) on the finished layers, ordered from the incident medium. Each output
    has one value per thickness; 0 gives the finished layers alone.
    """
    lambda_nm, theta_deg, n0, ns, layers = point
    layers = np.ascontiguousarray(layers, dtype=np.float64)
    thicknesses = np.ascontiguousarray(thicknesses, dtype=np.float64)
    outs = [np.empty(len(thicknesses)) for _ in _KEYS]
    lib.tmm_monitor_curve(float(lambda_nm), float(theta_deg), n0.real, n0.imag, ns.real, ns.imag,
                          address(layers), len(layers), top.real, top.imag,
                          address(thicknesses), len(thicknesses), *map(address, outs))
    return dict(zip(_KEYS, outs))


def _deposition_spectra_kernel(grid):
    """
    The spectrum after each layer is deposited, the layers in deposition
    order, first deposited (next to the substrate) first. Each output is
    (N, nLam), one spectrum per step.
    """
    lambdas, angles, n0, ns, table, thick = grid
    n_lam, N = len(lambdas), len(thick)
    outs = [np.empty((N, n_lam)) for _ in _KEYS]
    done = lib.tmm_deposition_spectra(address(lambdas), n_lam, address(angles), len(angles), address(n0),
                                      address(ns), address(table), address(thick), N, *map(address, outs))
    if not done:
        raise MemoryError("tmmcore: the deposition spectra could not allocate their working memory")
    return dict(zip(_KEYS, outs))


class _GrowingStack:
    """
    A finished stack, folded once over a wavelength grid, under a growing top
    layer: set_top() gives the top layer's index, sample(d) the outputs at
    thickness d (nm) over the whole grid. Frees its kernel memory on close(),
    or at the end of a with block.
    """

    def __init__(self, grid):
        lambdas, angles, n0, ns, table, thick = grid
        self._n_lam = len(lambdas)
        self._handle = lib.tmm_growing_stack_create(address(lambdas), self._n_lam, address(angles), len(angles),
                                                    address(n0), address(ns), address(table), address(thick),
                                                    len(thick))
        if not self._handle:
            raise MemoryError("tmmcore: the growing stack could not allocate its working memory")

    def _live(self):
        if not self._handle:
            raise ValueError("tmmcore: the growing stack is closed")
        return self._handle

    def set_top(self, top):
        """The growing layer's index, one complex value per wavelength."""
        top = np.ascontiguousarray(np.broadcast_to(top, (self._n_lam,)), dtype=np.complex128)
        lib.tmm_growing_stack_set_top(self._live(), address(top))

    def sample(self, d):
        """The outputs with the top layer at thickness d (nm), each over the grid."""
        outs = [np.empty(self._n_lam) for _ in _KEYS]
        if not lib.tmm_growing_stack_sample(self._live(), float(d), *map(address, outs)):
            raise ValueError("tmmcore: set_top() must give the growing layer before a sample with d > 0")
        return dict(zip(_KEYS, outs))

    def close(self):
        """Frees the kernel memory; the stack cannot be sampled after."""
        if getattr(self, "_handle", None):
            lib.tmm_growing_stack_free(self._handle)
            self._handle = None

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __del__(self):
        self.close()
