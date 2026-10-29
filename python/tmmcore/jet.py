"""
Taylor jets in angular frequency ω, to give the phase calls a dispersive index.

A jet is the four complex numbers [f, f', f''/2!, f'''/3!] of a function of ω
at one frequency, a complex array of shape (4,); jets over a grid are an
array of shape (..., 4). Ordinary power-series arithmetic on jets
differentiates a formula exactly, so a dispersion formula written with these
calls gives the index jet the phase calls need, with no finite differences.

Every call takes whole grids of jets and real numbers and broadcasts them over
their leading axes, as numpy does. The arithmetic is the C kernel's.
"""

import numpy as np

from ._library import lib
from ._stack import address


def _jets(value, what="jet"):
    jets = np.asarray(value)
    if jets.shape[-1:] != (4,) or not np.issubdtype(jets.dtype, np.number):
        raise TypeError(f"{what} must be a jet, 4 complex numbers [f, f', f''/2!, f'''/3!], or an array of jets")
    return jets.astype(np.complex128)


def _reals(value, what):
    reals = np.asarray(value)
    if not np.issubdtype(reals.dtype, np.number) or np.iscomplexobj(reals):
        raise TypeError(f"{what} must be a real number or an array of them")
    return reals.astype(np.float64)


def _grid(arguments):
    """The leading shape the arguments broadcast to, and each argument broadcast to it.
    `arguments` is (array, own axes) pairs: 1 for a jet, 0 for a real."""
    lead = np.broadcast_shapes(*(array.shape[:array.ndim - own] for array, own in arguments))
    grids = [np.ascontiguousarray(np.broadcast_to(array, lead + array.shape[array.ndim - own:]))
             for array, own in arguments]
    return lead, grids, int(np.prod(lead, dtype=np.int64))


def _apply(kernel, *arguments):
    """kernel(argument grids..., count, out) over the broadcast grid; returns the jets written to out."""
    lead, grids, count = _grid(arguments)
    out = np.empty(lead + (4,), dtype=np.complex128)
    kernel(*map(address, grids), count, address(out))
    return out


def _numbers(value, what):
    numbers = np.asarray(value)
    if not np.issubdtype(numbers.dtype, np.number):
        raise TypeError(f"{what} must be a number or an array of them")
    return numbers


def jet_constant(value, imaginary=0.0):
    """A constant: value + i·imaginary with every derivative zero. value may be complex."""
    value = _numbers(value, "value")
    imaginary = _reals(imaginary, "imaginary")
    return _apply(lib.tmm_jet_constant, (np.real(value).astype(np.float64), 0),
                  (np.imag(value).astype(np.float64) + imaginary, 0))


def jet_from_derivatives(value, first=0.0, second=0.0, third=0.0):
    """A jet from plain derivatives [f, f', f'', f'''], each a number or an array, real or complex."""
    terms = np.broadcast_arrays(*(_numbers(term, name).astype(np.complex128) for name, term in
                                  (("value", value), ("first", first), ("second", second), ("third", third))))
    return _apply(lib.tmm_jet_from_derivatives, (np.stack(terms, axis=-1), 1))


def jet_derivatives(jet):
    """The plain derivatives [f, f', f'', f'''] of a jet."""
    return _apply(lib.tmm_jet_derivatives, (_jets(jet), 1))


def jet_add(left, right):
    return _apply(lib.tmm_jet_add, (_jets(left, "left"), 1), (_jets(right, "right"), 1))


def jet_subtract(left, right):
    return _apply(lib.tmm_jet_subtract, (_jets(left, "left"), 1), (_jets(right, "right"), 1))


def jet_multiply(left, right):
    return _apply(lib.tmm_jet_multiply, (_jets(left, "left"), 1), (_jets(right, "right"), 1))


def jet_divide(numerator, denominator):
    return _apply(lib.tmm_jet_divide, (_jets(numerator, "numerator"), 1),
                  (_jets(denominator, "denominator"), 1))


def jet_scale(jet, scalar):
    """Every order multiplied by a real scalar."""
    return _apply(lib.tmm_jet_scale, (_jets(jet), 1), (_reals(scalar, "scalar"), 0))


def jet_reciprocal(jet):
    return _apply(lib.tmm_jet_reciprocal, (_jets(jet), 1))


def jet_sqrt(jet):
    """Principal square root."""
    return _apply(lib.tmm_jet_sqrt, (_jets(jet), 1))


def jet_exp(jet):
    return _apply(lib.tmm_jet_exp, (_jets(jet), 1))


def jet_log(jet):
    """Principal logarithm."""
    return _apply(lib.tmm_jet_log, (_jets(jet), 1))


def jet_power(jet, exponent):
    """A real power, as exp(exponent · log jet); the branch cut of the log applies."""
    return _apply(lib.tmm_jet_power, (_jets(jet), 1), (_reals(exponent, "exponent"), 0))


def jet_sin_cos(jet):
    """Sine and cosine together, as {"sine", "cosine"}."""
    lead, (jets,), count = _grid([(_jets(jet), 1)])
    sine = np.empty(lead + (4,), dtype=np.complex128)
    cosine = np.empty(lead + (4,), dtype=np.complex128)
    lib.tmm_jet_sin_cos(address(jets), count, address(sine), address(cosine))
    return {"sine": sine, "cosine": cosine}


def jet_sin(jet):
    return _apply(lib.tmm_jet_sin, (_jets(jet), 1))


def jet_cos(jet):
    return _apply(lib.tmm_jet_cos, (_jets(jet), 1))


def jet_compose(value, derivatives, x):
    """
    f(x) for a real function known through its value f(x0) and its
    derivatives [f'(x0), f''(x0), f'''(x0)], x0 the value of the jet x: the
    chain rule that turns a formula differentiated by hand into a jet.
    """
    derivatives = _reals(derivatives, "derivatives")
    if derivatives.shape[-1:] != (3,):
        raise TypeError("derivatives must be [f', f'', f'''] or an array of them")
    return _apply(lib.tmm_jet_compose, (_reals(value, "value"), 0), (derivatives, 1), (_jets(x, "x"), 1))


def wavelength_omega_jet(wavelength, omega):
    """
    The jet of the wavelength as a function of ω, λ(ω) = 2πc/ω, at a
    wavelength and its angular frequency. The wavelength unit of the result is
    the one the wavelength is given in.
    """
    return _apply(lib.tmm_wavelength_omega_jet, (_reals(wavelength, "wavelength"), 0),
                  (_reals(omega, "omega"), 0))


def jet_with_imaginary_part(real_jet, imaginary_jet):
    """The real parts of one jet as the real parts, those of another as the imaginary parts."""
    return _apply(lib.tmm_jet_with_imaginary_part, (_jets(real_jet, "real_jet"), 1),
                  (_jets(imaginary_jet, "imaginary_jet"), 1))


def jet_real_part(jet):
    return _apply(lib.tmm_jet_real_part, (_jets(jet), 1))


def jet_clamp_real_minimum(jet, minimum):
    """The jet while its value is at least `minimum`, otherwise the constant `minimum`."""
    return _apply(lib.tmm_jet_clamp_real_minimum, (_jets(jet), 1), (_reals(minimum, "minimum"), 0))


def jet_clamp_imaginary(jet, limit):
    """
    The value's imaginary part held to ±limit, the real part, which carries the
    phase, kept. The phase calls hold a layer's phase thickness this way.
    """
    return _apply(lib.tmm_jet_clamp_imaginary, (_jets(jet), 1), (_reals(limit, "limit"), 0))


def jet_is_finite(jet):
    """True where all four terms of a jet are finite."""
    lead, (jets,), count = _grid([(_jets(jet), 1)])
    out = np.empty(lead, dtype=np.intc)
    lib.tmm_jet_is_finite(address(jets), count, address(out))
    return bool(out) if out.ndim == 0 else out.astype(bool)
