"""
The wavelengths and keyword options of a call, read and checked: each reader
returns the value as the C kernel takes it, or raises TypeError naming what is
wrong. polarized() runs a call in the polarization asked, and for "avg" in s
and p, returning their mean.
"""

import numpy as np

from ._stack import as_array, is_real, numeric_array

POL_CODE = {"s": 0, "p": 1}


def read_wavelength(wavelength):
    """The wavelength of a one-wavelength call, in nm."""
    if not is_real(wavelength):
        raise TypeError("wavelength must be a number (nm)")
    return float(wavelength)


def read_wavelengths(wavelengths):
    """The wavelengths of a Spectrum call, in nm."""
    array = as_array(wavelengths)
    lambdas = None if array is None or array.ndim != 1 else numeric_array(array, array.shape)
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
    angles = None if n_lam is None else numeric_array(angle, (n_lam,))
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
    jets = numeric_array(sine, (4,) if n_lam is None else (n_lam, 4), complex_allowed=True)
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


def read_fractions(fractions):
    """Where a needle scan also opens each layer, as fractions of its thickness from the incident side."""
    array = as_array(fractions)
    values = None if array is None or array.ndim != 1 else numeric_array(array, array.shape)
    if values is None:
        raise TypeError("fractions must be a list of numbers")
    return np.ascontiguousarray(values, dtype=np.float64)


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
