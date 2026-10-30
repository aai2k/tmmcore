"""
The package against Essential Macleod, on the fixtures of tests/macleod.mjs.

tests/macleod/*.json hold Essential Macleod's results for ten designs of
constant-index materials, written by tests/macleod/export.ps1. Each is run
through the package: R, T and A in s, p and natural light, the phase of r
and t, their GD, GDD and TOD, Psi and Delta, and the reflectance from the
substrate side. The conventions, the tolerances and the two places that take
another reference (the reflection delays of the air gap at 30°, and an
absorbing incident medium) are those of tests/macleod.mjs, which gives the
reasons. Two things that file checks are left to it: the monitoring calls,
and the delays under an absorbing incident medium, which it reads from
frequency jets. The package exposes neither.

    python -m unittest discover -s python/tests
"""

import cmath
import glob
import json
import math
import os
import unittest

import numpy as np

import tmmcore

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURES = os.path.join(HERE, "..", "..", "tests", "macleod")

FRACTION = 1e-12  # R, T, A
DEGREES = 1e-10   # phases, Psi and Delta
DELAY = 1e-9      # GD, GDD and TOD, as a fraction of the scan's largest, or of 1 fs^n
T_PHASE_FLOOR = 1e-28


def medium(entry):
    return {"n": entry["n"][0], "k": entry["n"][1]}


def stack_of(case):
    return {
        "incident": medium(case["incident"]),
        "substrate": medium(case["substrate"]),
        "layers": [{**medium(layer), "d": layer["d"]} for layer in case["layers"]],
    }


def reversed_stack(stack):
    return {"incident": stack["substrate"], "substrate": stack["incident"], "layers": stack["layers"][::-1]}


def wrap(degrees):
    return (np.asarray(degrees) + 180) % 360 - 180


def view(incident, angle, pol, rt, ph):
    """R, T, A and the phases as Essential Macleod defines them (tests/macleod.mjs)."""
    if incident["k"] == 0:
        return {"R": rt["R"], "T": rt["T"], "A": rt["A"],
                "r_phase": ph["r"]["phase_rad"], "t_phase": ph["t"]["phase_rad"], "delays": ph}
    n0 = complex(incident["n"], -incident["k"])
    sine = incident["n"] * math.sin(math.radians(angle))
    cosine = cmath.sqrt(n0 * n0 - sine * sine)
    eta0 = cosine if pol == "s" else n0 * n0 / cosine
    a = eta0.real
    r = np.sqrt(ph["r"]["magnitude_squared"]) * np.exp(1j * ph["r"]["phase_rad"])
    denominator = a * (1 + r) + eta0 * (1 - r)
    r_view = (a * (1 + r) - eta0 * (1 - r)) / denominator
    t_ratio = 2 * a / denominator
    return {"R": np.abs(r_view) ** 2, "T": rt["T"] * np.abs(t_ratio) ** 2, "A": rt["A"],
            "r_phase": np.angle(r_view), "t_phase": ph["t"]["phase_rad"] + np.angle(t_ratio), "delays": None}


def film_delays(case, scan, pol):
    """Reflection GD, GDD and TOD of one film between identical media (tests/macleod.mjs)."""
    n_outer, n_film, d = case["incident"]["n"][0], case["layers"][0]["n"][0], case["layers"][0]["d"]
    angle = math.radians(scan["angle"])
    cos_outer = math.cos(angle)
    sine_film = n_outer * math.sin(angle) / n_film
    cos_film = math.sqrt(1 - sine_film * sine_film)
    tilted = (lambda n, c: n * c) if pol == "s" else (lambda n, c: n / c)
    eta_outer, eta_film = tilted(n_outer, cos_outer), tilted(n_film, cos_film)
    r01 = (eta_outer - eta_film) / (eta_outer + eta_film)
    wavelengths = np.asarray(scan["wavelengths"])
    x = 4 * math.pi * n_film * d * cos_film / wavelengths
    a = x / tmmcore.omega_from_lambda_nm(wavelengths)
    u = r01 * r01 * np.exp(-1j * x)
    w = 1 - u
    h1, h2, h3 = 1j * u / w, u / w ** 2, -1j * u * (1 + u) / w ** 3
    return {"gd": a / 2 + a * h1.imag, "gdd": a ** 2 * h2.imag, "tod": a ** 3 * h3.imag}


class AgainstEssentialMacleod(unittest.TestCase):
    def within(self, actual, expected, tolerance, what, mask=None):
        diff = np.abs(np.asarray(actual, dtype=np.float64) - np.asarray(expected, dtype=np.float64))
        if mask is not None:
            diff = np.where(mask, diff, 0.0)
        at = int(np.argmax(diff))
        self.assertLessEqual(diff[at], tolerance, f"{what}, entry {at}: off by {diff[at]!r}")

    def delays(self, ph, em, film, where):
        for side in ("r", "t"):
            for key, name in (("gd", "Gd"), ("gdd", "Gdd"), ("tod", "Tod")):
                exact = side == "r" and film is not None
                reference = np.asarray(film[key] if exact else em[side + name])
                scale = max(1.0, float(np.max(np.abs(reference))))
                self.within(np.asarray(ph[side][key]) / scale, reference / scale, DELAY,
                            f"{where} {side} {key}{' (closed form)' if exact else ''}")

    def scan(self, case, scan):
        stack, angle, wavelengths = stack_of(case), scan["angle"], scan["wavelengths"]
        views = {}
        for pol in ("s", "p"):
            where = f"{case['case']} {angle}° {pol}"
            em = scan[pol]
            rt = tmmcore.rta_spectrum(stack, wavelengths, angle=angle, pol=pol)
            ph = tmmcore.phase_spectrum(stack, wavelengths, pol=pol, angle=angle)
            v = views[pol] = view(stack["incident"], angle, pol, rt, ph)
            back = reversed_stack(stack)
            back_view = view(back["incident"], angle, pol,
                             tmmcore.rta_spectrum(back, wavelengths, angle=angle, pol=pol),
                             tmmcore.phase_spectrum(back, wavelengths, pol=pol, angle=angle))
            for key in ("R", "T", "A"):
                self.within(v[key], em[key], FRACTION, f"{where} {key}")
            self.within(back_view["R"], em["Rback"], FRACTION, f"{where} Rback")
            self.within(wrap(np.degrees(v["r_phase"]) - em["rPhaseDeg"]), 0, DEGREES, f"{where} r phase")
            transmits = np.asarray(em["T"]) > T_PHASE_FLOOR
            self.within(wrap(np.degrees(v["t_phase"]) - em["tPhaseDeg"]), 0, DEGREES, f"{where} t phase", transmits)
            if v["delays"] is not None:
                film = film_delays(case, scan, pol) if case["case"] == "frustrated-tir" and angle == 30 else None
                self.delays(v["delays"], em, film, where)
        s, p, em = views["s"], views["p"], scan["s"]
        where = f"{case['case']} {angle}°"
        transmits = np.asarray(em["T"]) > T_PHASE_FLOOR
        self.within(np.degrees(np.arctan(np.sqrt(p["R"] / s["R"]))), em["rPsiDeg"], DEGREES, f"{where} Psi of r")
        self.within(wrap(np.degrees(p["r_phase"] - s["r_phase"]) + 180 - em["rDeltaDeg"]), 0, DEGREES, f"{where} Delta of r")
        with np.errstate(divide="ignore", invalid="ignore"):
            psi_t = np.degrees(np.arctan(np.sqrt(p["T"] / s["T"])))
        self.within(psi_t, em["tPsiDeg"], DEGREES, f"{where} Psi of t", transmits)
        self.within(wrap(np.degrees(p["t_phase"] - s["t_phase"]) - em["tDeltaDeg"]), 0, DEGREES, f"{where} Delta of t", transmits)
        if stack["incident"]["k"] == 0:
            natural = tmmcore.rta_spectrum(stack, wavelengths, angle=angle)
            for key in ("R", "T", "A"):
                mean = (np.asarray(scan["s"][key]) + np.asarray(scan["p"][key])) / 2
                self.within(natural[key], mean, FRACTION, f"{where} {key} of natural light")

    def test_every_design(self):
        paths = sorted(glob.glob(os.path.join(FIXTURES, "*.json")))
        self.assertTrue(paths, "tests/macleod/*.json are missing")
        for path in paths:
            with open(path, encoding="utf-8") as file:
                case = json.load(file)
            for scan in case["scans"]:
                with self.subTest(design=case["case"], angle=scan["angle"]):
                    self.scan(case, scan)
        print(f"\n{len(paths)} designs against Essential Macleod")


if __name__ == "__main__":
    unittest.main()
