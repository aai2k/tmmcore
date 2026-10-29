"""
The C kernel of tmmcore, built with the package as a shared library, and the
argument and result types of its calls.

Every pointer argument is a C-contiguous numpy array passed by its address:
float64 for real values, complex128 for complex ones, whose memory is the
[re, im] pairs the kernel reads.
"""

import ctypes
import os
import sys

_FILE = "_tmmcore.dll" if sys.platform == "win32" else "_tmmcore.so"

lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), _FILE))

_P = ctypes.c_void_p
_D = ctypes.c_double
_I = ctypes.c_int

# The point calls start with (λ, θ, pol, n0 re, n0 im, ns re, ns im, layers, N),
# the grid calls with (λs, nLam, θs, nθ, pol, n0s, nss, index table, d, N),
# pol left out where a call computes both.
_POINT = [_D, _D, _I, _D, _D, _D, _D, _P, _I]
_GRID = [_P, _I, _P, _I, _I, _P, _P, _P, _P, _I]
_GRID_BOTH = [_P, _I, _P, _I, _P, _P, _P, _P, _I]
_PHASE_POINT = [_D, _D, _I, _P, _P, _P, _P, _I, _D, _P]

_SIGNATURES = {
    "tmm_rta": (None, _POINT + [_P]),
    "tmm_rta_spectrum": (None, _GRID_BOTH + [_P] * 6),
    "tmm_jacobian": (None, _POINT + [_P] * 4),
    "tmm_jacobian_spectrum": (None, _GRID + [_P] * 4),
    "tmm_hessian": (None, _POINT + [_P] * 7),
    "tmm_hessian_spectrum": (None, _GRID + [_P] * 7),
    "tmm_needle_scan": (None, _POINT + [_P, _I, _P, _I] + [_P] * 3),
    "tmm_needle_scan_spectrum": (None, _GRID + [_P, _I, _P, _I] + [_P] * 3),
    "tmm_phase": (None, _PHASE_POINT + [_P]),
    "tmm_phase_spectrum": (None, _GRID + [_P, _P, _P]),
    "tmm_phase_jacobian": (None, _PHASE_POINT + [_P, _P]),
    "tmm_phase_jacobian_spectrum": (None, _GRID + [_P, _P, _P, _P]),
    "tmm_omega_from_lambda_nm": (None, [_P, _I, _P]),
    "tmm_monitor_curve": (None, [_D, _D, _D, _D, _D, _D, _P, _I, _D, _D, _P, _I] + [_P] * 6),
    "tmm_deposition_spectra": (_I, _GRID_BOTH + [_P] * 6),
    "tmm_growing_stack_create": (_P, _GRID_BOTH),
    "tmm_growing_stack_set_top": (None, [_P, _P]),
    "tmm_growing_stack_sample": (_I, [_P, _D] + [_P] * 6),
    "tmm_growing_stack_free": (None, [_P]),
    "tmm_jet_constant": (None, [_P, _P, _I, _P]),
    "tmm_jet_from_derivatives": (None, [_P, _I, _P]),
    "tmm_jet_derivatives": (None, [_P, _I, _P]),
    "tmm_jet_compose": (None, [_P, _P, _P, _I, _P]),
    "tmm_jet_sin_cos": (None, [_P, _I, _P, _P]),
    "tmm_jet_is_finite": (None, [_P, _I, _P]),
}
# Jet calls of one jet grid: (jets, count, out).
for _name in ("reciprocal", "sqrt", "exp", "log", "sin", "cos", "real_part"):
    _SIGNATURES["tmm_jet_" + _name] = (None, [_P, _I, _P])
# Jet calls of two grids, of jets or of one real per jet: (a, b, count, out).
for _name in ("jet_add", "jet_subtract", "jet_multiply", "jet_divide", "jet_with_imaginary_part",
              "jet_scale", "jet_power", "jet_clamp_real_minimum", "jet_clamp_imaginary",
              "wavelength_omega_jet"):
    _SIGNATURES["tmm_" + _name] = (None, [_P, _P, _I, _P])

for _name, (_result, _arguments) in _SIGNATURES.items():
    _function = getattr(lib, _name)
    _function.restype = _result
    _function.argtypes = _arguments
