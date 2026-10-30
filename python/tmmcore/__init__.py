"""
tmmcore: transfer-matrix optics of thin-film stacks, with exact analytic
derivatives.

Every call takes a stack, a wavelength in nm or an array of them, and
keyword options:

    rta, rta_spectrum                           R, T and A
    jacobian, jacobian_spectrum                 with their thickness derivatives
    hessian, hessian_spectrum                   and the second derivatives
    phase, phase_spectrum                       phase, GD, GDD and TOD of r and t
    phase_jacobian, phase_jacobian_spectrum     with their thickness derivatives

    stack = {"incident": 1.0, "substrate": 1.52,
             "layers": [{"n": 2.35, "d": 58.5}, {"n": 1.46, "d": 94.2}]}
    tmmcore.rta_spectrum(stack, numpy.linspace(400, 800, 401))["R"]

A layer is {"n", "k", "d"} with d in nm, ñ = n + ik, k optional; n may be
complex instead. `angle` is the angle of incidence in degrees from the normal,
0 when left out. `pol` is "s", "p" or "avg" (the default, the mean of the two)
in the R/T and derivative calls, and "s" or "p" in the phase calls, which may
leave it out at normal incidence, where the two are the same. In place of
`angle` the phase calls take `sine`, the sine of the angle of incidence as a
jet. The jet calls (jet.py) build the dispersive index jets the phase calls
take.

The numbers come from the same C kernel as the npm package; they agree with
it to round-off, the last digits differing where the platform's math library
does.
"""

from .derivatives import hessian, hessian_spectrum, jacobian, jacobian_spectrum
from .jet import (
    jet_add, jet_clamp_imaginary, jet_clamp_real_minimum, jet_compose, jet_constant, jet_cos,
    jet_derivatives, jet_divide, jet_exp, jet_from_derivatives, jet_is_finite, jet_log, jet_multiply,
    jet_power, jet_real_part, jet_reciprocal, jet_scale, jet_sin, jet_sin_cos, jet_sqrt, jet_subtract,
    jet_with_imaginary_part, wavelength_omega_jet,
)
from .phase import omega_from_lambda_nm, phase, phase_jacobian, phase_jacobian_spectrum, phase_spectrum
from .rta import rta, rta_spectrum

__all__ = [
    "rta", "rta_spectrum", "jacobian", "jacobian_spectrum", "hessian", "hessian_spectrum",
    "phase", "phase_spectrum", "phase_jacobian", "phase_jacobian_spectrum", "omega_from_lambda_nm",
    "jet_constant", "jet_from_derivatives", "jet_derivatives", "jet_add", "jet_subtract", "jet_scale",
    "jet_multiply", "jet_reciprocal", "jet_divide", "jet_sqrt", "jet_exp", "jet_log", "jet_power",
    "jet_sin_cos", "jet_sin", "jet_cos", "jet_compose", "wavelength_omega_jet", "jet_with_imaginary_part",
    "jet_real_part", "jet_clamp_real_minimum", "jet_clamp_imaginary", "jet_is_finite",
]
