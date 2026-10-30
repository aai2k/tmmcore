"""
A stack whose top layer grows, as in a deposition run: the monitoring curve
at one wavelength, the spectra after each deposited layer, and a growing top
layer sampled over a wavelength grid.

Each returns, per polarization, the forward R and T of the coated surface and
its reflectance from the substrate side, for a monitoring simulation that
treats the witness as an incoherent slab with a bare back face:
{"Rs", "Ts", "Rp", "Tp", "Rrs", "Rrp"}.
"""

import numpy as np

from . import _options, _stack
from ._library import lib
from ._stack import address

_KEYS = ["Rs", "Ts", "Rp", "Tp", "Rrs", "Rrp"]


def monitor_curve(stack, wavelength, *, top, thicknesses, angle=0.0):
    """
    The signal of a layer as it grows on the stack, at one wavelength.

    stack        the finished layers, {"incident", "substrate", "layers"}
    wavelength   in nm
    top          the growing layer's index, given as a medium is: a number or {"n", "k"}
    thicknesses  of the growing layer to sample, in nm; 0 gives the stack alone
    angle        angle of incidence in degrees from the normal

    Returns {"Rs", "Ts", "Rp", "Tp", "Rrs", "Rrp"}, each an array with one
    value per thickness.
    """
    lam = _options.read_wavelength(wavelength)
    theta = _options.read_angle(angle)
    top = _stack.read_index(top, "top")
    thicknesses = _options.read_thicknesses(thicknesses)
    n0, ns, layers = _stack.point_stack(stack)
    outs = [np.empty(len(thicknesses)) for _ in _KEYS]
    lib.tmm_monitor_curve(lam, theta, n0.real, n0.imag, ns.real, ns.imag, address(layers), len(layers),
                          top.real, top.imag, address(thicknesses), len(thicknesses), *map(address, outs))
    return dict(zip(_KEYS, outs))


def deposition_spectra(stack, wavelengths, *, angle=0.0):
    """
    The spectrum after each layer of the stack is deposited, from the
    substrate outward: the last layer listed first, the first layer last.

    Takes the arguments of rta_spectrum but pol. Returns {"Rs", "Ts", "Rp",
    "Tp", "Rrs", "Rrp"}, each (N, nLam): row k with k + 1 layers deposited,
    the last row the whole stack.
    """
    lambdas = _options.read_wavelengths(wavelengths)
    n_lam = len(lambdas)
    angles = _options.read_angle(angle, n_lam)
    n0, ns, table, thick = _stack.grid_stack(stack, n_lam)
    table, thick = np.ascontiguousarray(table[::-1]), np.ascontiguousarray(thick[::-1])
    N = len(thick)
    outs = [np.empty((N, n_lam)) for _ in _KEYS]
    done = lib.tmm_deposition_spectra(address(lambdas), n_lam, address(angles), len(angles), address(n0),
                                      address(ns), address(table), address(thick), N, *map(address, outs))
    if not done:
        raise MemoryError("tmmcore: the deposition spectra could not allocate their working memory")
    return dict(zip(_KEYS, outs))


def growing_stack(stack, wavelengths, *, angle=0.0):
    """
    A layer growing on the stack, over a wavelength grid. The stack is
    multiplied out once, so each sample of the growing layer costs one layer
    matrix per wavelength and polarization.

    Takes the arguments of rta_spectrum but pol, the stack the finished
    layers. Returns a GrowingStack.
    """
    return GrowingStack(stack, wavelengths, angle=angle)


class GrowingStack:
    """
    A finished stack, multiplied out once over a wavelength grid, under a
    growing top layer: set_top() gives the growing layer's index, sample(d)
    the outputs with it at thickness d (nm) over the whole grid. Frees its
    kernel memory on close(), or at the end of a with block. Takes the
    arguments of growing_stack.
    """

    def __init__(self, stack, wavelengths, *, angle=0.0):
        lambdas = _options.read_wavelengths(wavelengths)
        angles = _options.read_angle(angle, len(lambdas))
        n0, ns, table, thick = _stack.grid_stack(stack, len(lambdas))
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
        """The growing layer's index, given as a medium is: a number or {"n", "k"}, each value or one per wavelength."""
        top = _stack.read_index(top, "top", self._n_lam)
        lib.tmm_growing_stack_set_top(self._live(), address(top))

    def sample(self, d):
        """The outputs with the top layer at thickness d (nm), each over the grid; d <= 0 gives the stack alone."""
        d = _options.read_thickness(d)
        outs = [np.empty(self._n_lam) for _ in _KEYS]
        if not lib.tmm_growing_stack_sample(self._live(), d, *map(address, outs)):
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
