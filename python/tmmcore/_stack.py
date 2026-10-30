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

POL_CODE = {"s": 0, "p": 1}

NO_DISPERSION = ("has no derivatives, so GD would miss the material's dispersion: give a jet per "
                 "wavelength (see jet_from_derivatives), or one number for a constant index")


def address(array):
    """The address of a C-contiguous array, for a pointer argument."""
    return array.ctypes.data


def is_real(value):
    return isinstance(value, numbers.Real) and not isinstance(value, bool)


def is_number(value):
    return isinstance(value, numbers.Complex) and not isinstance(value, bool)


def _array(value):
    """value as a numpy array, or None where numpy cannot make one, as from a ragged list."""
    try:
        return np.asarray(value)
    except ValueError:
        return None


def _numbers(value, shape, complex_allowed=False):
    """value as an array if it has this shape and holds numbers, real ones unless complex_allowed; else None."""
    array = _array(value)
    if array is None or array.shape != shape or not np.issubdtype(array.dtype, np.number):
        return None
    return array if complex_allowed or not np.iscomplexobj(array) else None


def read_wavelength(wavelength):
    """The wavelength of a one-wavelength call, in nm."""
    if not is_real(wavelength):
        raise TypeError("wavelength must be a number (nm)")
    return float(wavelength)


def read_wavelengths(wavelengths):
    """The wavelengths of a Spectrum call, in nm."""
    array = _array(wavelengths)
    lambdas = None if array is None or array.ndim != 1 else _numbers(array, array.shape)
    if lambdas is None:
        raise TypeError("wavelengths must be a list of numbers (nm)")
    return np.ascontiguousarray(lambdas, dtype=np.float64)


def read_angle(angle, n_lam=None):
    """
    The angle of incidence in degrees: a number, or in a Spectrum call (n_lam
    given) also one angle per wavelength. A Spectrum call gets it as an array.
    """
    if is_real(angle):
        return float(angle) if n_lam is None else np.array([angle], dtype=np.float64)
    angles = None if n_lam is None else _numbers(angle, (n_lam,))
    if angles is None:
        raise TypeError("angle must be a number (degrees)" if n_lam is None
                        else f"angle must be a number or {n_lam} numbers, one per wavelength (degrees)")
    return np.ascontiguousarray(angles, dtype=np.float64)


def _read_pol(pol, allowed):
    if pol not in allowed:
        raise TypeError(f"pol must be {', '.join(map(repr, allowed[:-1]))} or {allowed[-1]!r}")
    return pol


def read_rta_pol(pol):
    """Polarization of an R/T call: 's', 'p' or 'avg', the mean of the two."""
    return _read_pol(pol, ("s", "p", "avg"))


def read_sine(sine, n_lam=None):
    """
    The sine of the angle of incidence as a jet, which a phase call takes in
    place of the angle: one jet, or with n_lam a jet per wavelength.
    """
    jets = _numbers(sine, (4,) if n_lam is None else (n_lam, 4), complex_allowed=True)
    if jets is None:
        raise TypeError("sine must be a jet" if n_lam is None else f"sine must be {n_lam} jets, one per wavelength")
    return np.ascontiguousarray(jets, dtype=np.complex128)


def read_phase_options(angle, sine, pol, n_lam=None):
    """
    The angle, or the incident sine as jets and the angle 0 the kernel then
    ignores, and the polarization code of a phase call. pol is "s" or "p"; it
    may be left out at normal incidence, where the two are the same, and is
    then "s".
    """
    if sine is None:
        angle = read_angle(0.0 if angle is None else angle, n_lam)
        normal = not np.any(angle)
    elif angle is not None:
        raise TypeError("a phase call takes angle or sine, not both")
    else:
        sine = read_sine(sine, n_lam)
        angle = read_angle(0.0, n_lam)
        normal = not np.any(sine)
    if pol is None and normal:
        pol = "s"
    if pol not in ("s", "p"):
        raise TypeError("pol must be 's' or 'p'; it may be left out at normal incidence")
    return angle, sine, POL_CODE[pol]


def optional_address(array):
    """The address of an array, or None, a null pointer, for no array."""
    return None if array is None else address(array)


def mean(s, p):
    """The mean of the s and p results of one call, entry by entry, in the shape of either."""
    if isinstance(s, dict):
        return {key: mean(s[key], p[key]) for key in s}
    if isinstance(s, int):
        return s
    return (s + p) / 2


def polarized(pol, run):
    """run(side) for the polarization asked, or the mean of s and p for 'avg'."""
    return mean(run("s"), run("p")) if pol == "avg" else run(pol)


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
    values = None if n_lam is None else _numbers(value, (n_lam,), complex_allowed)
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
    jets = _numbers(n, (4,) if n_lam is None else (n_lam, 4), complex_allowed=True)
    if jets is not None:
        return np.ascontiguousarray(jets, dtype=np.complex128)
    if n_lam is not None and _numbers(n, (n_lam,), complex_allowed=True) is not None:
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
