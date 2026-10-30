"""
The stack every call takes, read into the arrays the C kernel uses.

A stack is a dict {"incident": ..., "substrate": ..., "layers": [...]}, the
layers ordered from the incident medium toward the substrate, each
{"n": ..., "k": ..., "d": ...} with d in nm. A medium is a number or
{"n": ..., "k": ...}. Each n and k is a number, or in a Spectrum call one value
per wavelength; k may be left out and is then 0. n may also be complex, n + ik,
with k left out. k ≥ 0 for absorbing media.

The phase calls also take an index as a Taylor jet in ω (jet.py), which
carries the material's dispersion: one jet at one wavelength, one per
wavelength in a Spectrum call, as complex arrays of shape (4,) and (nLam, 4).
A number there is a constant index with no dispersion. Values per wavelength
without derivatives are refused, since group delay computed from them would
leave out the material's dispersion.
"""

import numbers

import numpy as np

NO_DISPERSION = ("has no derivatives, so GD would miss the material's dispersion: give a jet per "
                 "wavelength (see jet_from_derivatives), or one number for a constant index")


def address(array):
    """The address of a C-contiguous array, for a pointer argument."""
    return array.ctypes.data


def is_real(value):
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


def is_number(value):
    return isinstance(value, numbers.Complex) and not isinstance(value, bool)


def as_array(value):
    """value as a numpy array, or None where numpy cannot make one, as from a ragged list."""
    try:
        return np.asarray(value)
    except ValueError:
        return None


def numeric_array(value, shape, complex_allowed=False):
    """value as an array if it has this shape and holds numbers, real ones unless complex_allowed; else None."""
    array = as_array(value)
    if array is None or array.shape != shape or not np.issubdtype(array.dtype, np.number):
        return None
    return array if complex_allowed or not np.iscomplexobj(array) else None


def optional_address(array):
    """The address of an array, or None, a null pointer, for no array."""
    return None if array is None else address(array)


def _read_medium(value, where):
    if is_number(value):
        return value, None
    if isinstance(value, dict) and "n" in value:
        return value["n"], value.get("k")
    raise TypeError(f"{where} must be a number or {{ n, k }}")


def _layer_keys(layer, i):
    if not isinstance(layer, dict) or "n" not in layer or "d" not in layer:
        raise TypeError(f"layer {i} needs {{n, d}}")
    return layer


def _thickness(layer, i):
    if not is_real(layer["d"]):
        raise TypeError(f"layer {i} d must be a number (nm)")
    return float(layer["d"])


def _read_stack(stack):
    """
    [(where, n, k)] for the incident medium, the substrate and each layer, and
    the thicknesses. Checked in the JavaScript's order, so the same stack gets
    the same message: every layer's keys, then every thickness, then the media.
    """
    if not isinstance(stack, dict) or not isinstance(stack.get("layers"), (list, tuple)):
        raise TypeError("stack must be { incident, substrate, layers: [{ n, k, d }] }")
    layers = [_layer_keys(layer, i) for i, layer in enumerate(stack["layers"])]
    thick = np.array([_thickness(layer, i) for i, layer in enumerate(layers)], dtype=np.float64)
    indices = [("incident", *_read_medium(stack.get("incident"), "incident")),
               ("substrate", *_read_medium(stack.get("substrate"), "substrate"))]
    indices += [(f"layer {i}", layer["n"], layer.get("k")) for i, layer in enumerate(layers)]
    return indices, thick


# ── R/T calls: ñ as a complex number ─────────────────────────────────────────

def _values(value, n_lam, where, complex_allowed):
    """A number, or with n_lam one value per wavelength, as complex128."""
    if is_number(value) and (complex_allowed or is_real(value)):
        return np.complex128(value) if n_lam is None else np.full(n_lam, value, dtype=np.complex128)
    values = None if n_lam is None else numeric_array(value, (n_lam,), complex_allowed)
    if values is None:
        raise TypeError(f"{where} must be a number at one wavelength" if n_lam is None
                        else f"{where} must be a number or {n_lam} numbers, one per wavelength; n and k are two lists")
    return values.astype(np.complex128)


def _index(n, k, n_lam, where):
    """ñ = n + ik, one value or one per wavelength."""
    index = _values(n, n_lam, f"{where} n", complex_allowed=True)
    if k is None:
        return index
    if np.iscomplexobj(n):
        raise TypeError(f"{where} k goes inside n when n is complex")
    return index + 1j * _values(k, n_lam, f"{where} k", complex_allowed=False).real


def point_stack(stack):
    """
    A stack at one wavelength, as the one-wavelength kernels take it: ñ of the
    two media, and the layers as N rows [n, k, d].
    """
    indices, thick = _read_stack(stack)
    n0, ns, *layer_indices = [_index(n, k, None, where) for where, n, k in indices]
    layers = np.array([[index.real, index.imag, d] for index, d in zip(layer_indices, thick)],
                      dtype=np.float64).reshape(len(thick), 3)
    return n0, ns, layers


def grid_stack(stack, n_lam):
    """
    A stack over a wavelength grid, as the grid kernels take it: ñ per
    wavelength for the two media, (N, nLam) for the layers, and the thicknesses.
    """
    indices, thick = _read_stack(stack)
    n0, ns, *layer_indices = [_index(n, k, n_lam, where) for where, n, k in indices]
    table = np.array(layer_indices, dtype=np.complex128).reshape(len(thick), n_lam)
    return n0, ns, table, thick


def read_index(value, where, n_lam=None):
    """
    An index given as a medium is, a number or {"n", "k"}: ñ, or with n_lam
    one per wavelength. `where` names it in a message.
    """
    return _index(*_read_medium(value, where), n_lam, where)


def read_candidates(candidates, n_lam=None):
    """
    The materials a needle scan tries: a list of indices, each read by
    read_index. ñ of each, or with n_lam an (nCand, nLam) array.
    """
    if isinstance(candidates, np.ndarray) and candidates.ndim == 1:
        candidates = list(candidates)
    if not isinstance(candidates, (list, tuple)):
        raise TypeError("candidates must be a list of indices, each a number or { n, k }")
    indices = [read_index(value, f"candidate {c}", n_lam) for c, value in enumerate(candidates)]
    shape = (len(indices),) if n_lam is None else (len(indices), n_lam)
    return np.array(indices, dtype=np.complex128).reshape(shape)


# ── Phase calls: ñ as a jet ──────────────────────────────────────────────────

def _jet_index(n, k, n_lam, where):
    """One jet, or with n_lam a jet per wavelength, for an index of a phase call."""
    if is_number(n):
        if k is not None and not is_real(k):
            raise TypeError(f"{where} k per wavelength {NO_DISPERSION}")
        jet = np.array([_index(n, k, None, where), 0, 0, 0], dtype=np.complex128)
        return jet if n_lam is None else np.tile(jet, (n_lam, 1))
    if k is not None:
        raise TypeError(f"{where} k goes inside the jet when n is a jet")
    jets = numeric_array(n, (4,) if n_lam is None else (n_lam, 4), complex_allowed=True)
    if jets is not None:
        return np.ascontiguousarray(jets, dtype=np.complex128)
    if n_lam is not None and numeric_array(n, (n_lam,), complex_allowed=True) is not None:
        raise TypeError(f"{where} n per wavelength {NO_DISPERSION}")
    raise TypeError(f"{where} n must be a number or a jet" if n_lam is None
                    else f"{where} n must be a number or {n_lam} jets, one per wavelength")


def point_jet_stack(stack):
    """A stack at one wavelength for the phase kernels: the jets of the media, (N, 4) for the layers, and d."""
    indices, thick = _read_stack(stack)
    n0, ns, *layers = [_jet_index(n, k, None, where) for where, n, k in indices]
    return n0, ns, np.array(layers, dtype=np.complex128).reshape(len(thick), 4), thick


def grid_jet_stack(stack, n_lam):
    """A stack over a grid for the phase kernels: jets per wavelength, (N, nLam, 4) for the layers."""
    indices, thick = _read_stack(stack)
    n0, ns, *layers = [_jet_index(n, k, n_lam, where) for where, n, k in indices]
    return n0, ns, np.array(layers, dtype=np.complex128).reshape(len(thick), n_lam, 4), thick
