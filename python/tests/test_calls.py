"""
What test_reference.py cannot show: the arguments the package refuses, the
shapes it returns, complex n against n and k, arrays in another memory order
or type than the C reads, plain numbers read without numpy against the same
numbers read through it, the order the layers are deposited in, and the
growing stack's lifetime.

    python -m unittest discover -s python/tests
"""

import unittest

import numpy as np

import tmmcore
from tmmcore import _stack

# A high and a low index layer, each a quarter wave at 550 nm.
PAIR = {"incident": 1.0, "substrate": 1.52,
        "layers": [{"n": 2.35, "d": 58.5}, {"n": 1.46, "d": 94.2}]}
METAL = {"incident": 1.0, "substrate": 1.52, "layers": [{"n": 0.13, "k": 3.0, "d": 40.0}]}
GRID = np.linspace(450.0, 650.0, 5)


def with_layers(stack, *layers):
    return {**stack, "layers": list(layers)}


class Real(float):
    """A float the stack reader does not know by its type, so it reads it through numpy scalars."""


class Complex(complex):
    """A complex number the stack reader does not know by its type, so it reads it through numpy scalars."""


class Arguments(unittest.TestCase):
    def test_pol(self):
        with self.assertRaisesRegex(TypeError, "pol must be 's', 'p' or 'avg'"):
            tmmcore.rta(PAIR, 550, pol="x")
        with self.assertRaisesRegex(TypeError, "pol must be 's' or 'p'"):
            tmmcore.phase(PAIR, 550, pol="avg")
        with self.assertRaisesRegex(TypeError, "pol must be 's' or 'p'; it may be left out at normal incidence"):
            tmmcore.phase_spectrum(PAIR, GRID, angle=30)
        # s and p are the same at normal incidence, where pol may be left out:
        # no angle, an angle of 0 at every wavelength, or a sine of 0.
        s, p, left_out = (tmmcore.phase(PAIR, 550, **kwargs) for kwargs in ({"pol": "s"}, {"pol": "p"}, {}))
        self.assertEqual(p, s)
        self.assertEqual(left_out, s)
        expected = tmmcore.phase_spectrum(PAIR, GRID, pol="p")
        for kwargs in ({"angle": np.zeros(5)}, {"sine": np.zeros((5, 4))}):
            result = tmmcore.phase_spectrum(PAIR, GRID, **kwargs)
            for side in ("r", "t"):
                for key, values in expected[side].items():
                    np.testing.assert_array_equal(result[side][key], values, f"{kwargs} {side} {key}")
        sine = np.zeros((5, 4))
        sine[2, 0] = 0.5
        with self.assertRaisesRegex(TypeError, "pol must be 's' or 'p'"):
            tmmcore.phase_spectrum(PAIR, GRID, sine=sine)

    def test_sine(self):
        sine = tmmcore.jet_constant(0.5)
        with self.assertRaisesRegex(TypeError, "a phase call takes angle or sine, not both"):
            tmmcore.phase(PAIR, 550, pol="s", angle=0, sine=sine)
        with self.assertRaisesRegex(TypeError, "sine must be a jet"):
            tmmcore.phase(PAIR, 550, pol="s", sine=0.5)
        with self.assertRaisesRegex(TypeError, "sine must be 5 jets, one per wavelength"):
            tmmcore.phase_spectrum(PAIR, GRID, pol="s", sine=sine)
        with self.assertRaisesRegex(TypeError, "pol must be 's' or 'p'"):
            tmmcore.phase(PAIR, 550, sine=sine)

    def test_wavelengths(self):
        for bad in (550, [[450, 550]], [450, [550, 650]], ["450", "550"], [550 + 0j]):
            with self.assertRaisesRegex(TypeError, r"wavelengths must be a list of numbers \(nm\)"):
                tmmcore.rta_spectrum(PAIR, bad)
        with self.assertRaisesRegex(TypeError, r"wavelength must be a number \(nm\)"):
            tmmcore.rta(PAIR, [550])

    def test_angle_per_wavelength(self):
        with self.assertRaisesRegex(TypeError, "angle must be a number or 5 numbers, one per wavelength"):
            tmmcore.rta_spectrum(PAIR, GRID, angle=[0, 10])
        np.testing.assert_array_equal(tmmcore.rta_spectrum(PAIR, GRID, angle=np.full(5, 30.0))["R"],
                                      tmmcore.rta_spectrum(PAIR, GRID, angle=30)["R"])

    def test_layer_keys(self):
        with self.assertRaisesRegex(TypeError, r"layer 1 needs \{n, d\}"):
            tmmcore.rta(with_layers(PAIR, {"n": 2.35, "d": 58.5}, {"n": 1.46, "thickness": 94.2}), 550)
        # Every layer's keys come before any thickness, as in the JavaScript.
        with self.assertRaisesRegex(TypeError, r"layer 1 needs \{n, d\}"):
            tmmcore.rta(with_layers(PAIR, {"n": 2.35, "d": "58.5"}, {"n": 1.46}), 550)
        with self.assertRaisesRegex(TypeError, "layer 0 has 'kappa'; a layer takes n, k and d"):
            tmmcore.rta(with_layers(PAIR, {"n": 2.35, "kappa": 0.1, "d": 58.5}), 550)
        with self.assertRaisesRegex(TypeError, "stack has 'angle'; a stack takes incident, substrate and layers"):
            tmmcore.rta({**PAIR, "angle": 30}, 550)

    def test_medium(self):
        with self.assertRaisesRegex(TypeError, r"incident must be a number or \{ n, k \}"):
            tmmcore.rta({**PAIR, "incident": [1.0, 0.0]}, 550)
        with self.assertRaisesRegex(TypeError, "substrate has 'kappa'; a medium takes n and k"):
            tmmcore.rta_spectrum({**PAIR, "substrate": {"n": 1.52, "kappa": 0.0}}, GRID)

    def test_complex_n(self):
        self.assertEqual(tmmcore.rta(with_layers(METAL, {"n": 0.13 + 3.0j, "d": 40.0}), 550), tmmcore.rta(METAL, 550))
        table = with_layers(METAL, {"n": np.full(5, 0.13 + 3.0j), "d": 40.0})
        np.testing.assert_array_equal(tmmcore.rta_spectrum(table, GRID)["R"], tmmcore.rta_spectrum(METAL, GRID)["R"])
        with self.assertRaisesRegex(TypeError, "layer 0 k goes inside n when n is complex"):
            tmmcore.rta(with_layers(METAL, {"n": 0.13 + 3.0j, "k": 3.0, "d": 40.0}), 550)

    def test_needle_options(self):
        np.testing.assert_equal(tmmcore.needle_scan(PAIR, 550, candidates=[2.1 + 0.01j], fractions=[0.5]),
                                tmmcore.needle_scan(PAIR, 550, candidates=[{"n": 2.1, "k": 0.01}], fractions=[0.5]))
        np.testing.assert_equal(tmmcore.needle_scan(PAIR, 550, candidates=[1.38]),
                                tmmcore.needle_scan(PAIR, 550, candidates=[1.38], fractions=[]))
        with self.assertRaisesRegex(TypeError, "candidates"):
            tmmcore.needle_scan(PAIR, 550)
        with self.assertRaisesRegex(TypeError, r"candidates must be a list of indices, each a number or \{ n, k \}"):
            tmmcore.needle_scan(PAIR, 550, candidates=1.38)
        with self.assertRaisesRegex(TypeError, r"candidate 1 must be a number or \{ n, k \}"):
            tmmcore.needle_scan(PAIR, 550, candidates=[1.38, [1.46, 0]])
        with self.assertRaisesRegex(TypeError, "candidate 1 has 'kappa'; a medium takes n and k"):
            tmmcore.needle_scan(PAIR, 550, candidates=[1.38, {"n": 1.46, "kappa": 0.0}])
        with self.assertRaisesRegex(TypeError, "candidate 0 n must be a number or 5 numbers, one per wavelength"):
            tmmcore.needle_scan_spectrum(PAIR, GRID, candidates=[{"n": [1.38, 1.37]}])
        for bad in (0.5, ["0.5"], [[0.5]]):
            with self.assertRaisesRegex(TypeError, "fractions must be a list of numbers"):
                tmmcore.needle_scan(PAIR, 550, candidates=[1.38], fractions=bad)

    def test_phase_index_per_wavelength(self):
        table = with_layers(PAIR, {"n": np.full(5, 2.35), "d": 58.5})
        with self.assertRaisesRegex(TypeError, "layer 0 n per wavelength has no derivatives"):
            tmmcore.phase_spectrum(table, GRID, pol="s")
        # A number is a constant index: the same as a constant jet at every wavelength.
        jets = [{"n": np.tile(tmmcore.jet_constant(layer["n"]), (5, 1)), "d": layer["d"]} for layer in PAIR["layers"]]
        numbers, constants = (tmmcore.phase_spectrum(stack, GRID, pol="s") for stack in (PAIR, with_layers(PAIR, *jets)))
        for side in ("r", "t"):
            for key, values in numbers[side].items():
                np.testing.assert_array_equal(values, constants[side][key], f"{side} {key}")


class Results(unittest.TestCase):
    def test_avg_is_the_mean_of_s_and_p(self):
        s, p, avg = (tmmcore.rta(PAIR, 550, angle=45, pol=pol) for pol in ("s", "p", "avg"))
        for key in ("R", "T", "A"):
            self.assertEqual(avg[key], (s[key] + p[key]) / 2)

    def test_shapes(self):
        N = len(PAIR["layers"])
        jacobian = tmmcore.jacobian_spectrum(PAIR, GRID)
        self.assertEqual((jacobian["R"].shape, jacobian["dRdd"].shape, jacobian["N"]), ((5,), (5, N), N))
        self.assertEqual(tmmcore.hessian_spectrum(PAIR, GRID, pol="s")["d2Rdd"].shape, (5, N, N))
        self.assertEqual(tmmcore.jacobian(PAIR, 550)["dRdd"].shape, (N,))
        self.assertEqual(tmmcore.hessian(PAIR, 550)["d2Tdd"].shape, (N, N))
        phase = tmmcore.phase_jacobian_spectrum(PAIR, GRID, pol="p")
        self.assertEqual((phase["t"]["gd"].shape, phase["r"]["d_gd"].shape), ((5,), (5, N)))
        self.assertEqual(tmmcore.phase_jacobian(PAIR, 550, pol="s")["r"]["d_gdd"].shape, (N,))
        needles = {"candidates": [1.38, {"n": 2.1, "k": 0.01}], "fractions": [0.25, 0.5, 0.75]}
        one = tmmcore.needle_scan(PAIR, 550, **needles)
        self.assertEqual((list(one), one["gaps"]["dR"].shape, one["intra"]["dA"].shape),
                         (["R", "T", "A", "gaps", "intra", "N"], (N + 1, 2), (N, 3, 2)))
        grid = tmmcore.needle_scan_spectrum(PAIR, GRID, pol="p", **needles)
        self.assertEqual((grid["R"].shape, grid["gaps"]["dT"].shape, grid["intra"]["dR"].shape),
                         ((5,), (5, N + 1, 2), (5, N, 3, 2)))

    def test_jets_broadcast(self):
        jets = tmmcore.jet_from_derivatives(np.array([1.5, 1.6]), 0.01)
        scaled = tmmcore.jet_scale(jets, np.array([[1.0], [2.0], [3.0]]))
        self.assertEqual(scaled.shape, (3, 2, 4))
        np.testing.assert_array_equal(scaled[2, 1], tmmcore.jet_scale(jets[1], 3.0))
        self.assertIs(tmmcore.jet_is_finite(jets[0]), True)
        with self.assertRaisesRegex(TypeError, "must be a jet"):
            tmmcore.jet_sqrt([1.0, 2.0, 3.0])


class Buffers(unittest.TestCase):
    """The C reads C-ordered float64 and complex128; arrays in another order or type give the same results."""

    def test_memory_order(self):
        lam = np.array([[450.0, 550.0, 650.0], [500.0, 600.0, 700.0]])
        np.testing.assert_array_equal(tmmcore.omega_from_lambda_nm(lam.T), tmmcore.omega_from_lambda_nm(lam).T)
        jets = tmmcore.jet_from_derivatives(1.5 + 4200 / GRID**2, 0.01)
        for where in ("incident", "substrate"):
            ordered = tmmcore.phase_spectrum({**PAIR, where: {"n": jets}}, GRID, pol="s")
            fortran = tmmcore.phase_spectrum({**PAIR, where: {"n": np.asfortranarray(jets)}}, GRID, pol="s")
            np.testing.assert_array_equal(fortran["r"]["gd"], ordered["r"]["gd"], where)

    def test_types(self):
        expected = tmmcore.jet_constant(np.array([1.5, 2.0]), 0.25)
        for dtype in (np.float32, np.longdouble):
            np.testing.assert_array_equal(tmmcore.jet_constant(np.array([1.5, 2.0], dtype=dtype), 0.25),
                                          expected, np.dtype(dtype).name)

    def test_plain_numbers(self):
        """
        Floats, ints and complex numbers are read without numpy scalars; the
        kernels get the arrays of the same numbers read through numpy.
        """
        def through_numpy(value):
            return Complex(value) if isinstance(value, complex) else Real(value)

        values = [1.46, 2, 0.0, -0.0, -1.5, 1e-310, 1e308]
        layers = [{"n": n, "d": 50.0} for n in values + [0.13 + 3.0j, complex(1.5, -0.0)]]
        layers += [{"n": n, "k": k, "d": 50.0} for n in values for k in values]
        plain = {"incident": 1.0, "substrate": {"n": 1.52, "k": -0.0}, "layers": layers}
        subclassed = {"incident": Real(1.0), "substrate": {"n": Real(1.52), "k": Real(-0.0)},
                      "layers": [{key: through_numpy(value) for key, value in layer.items()} for layer in layers]}
        readers = {"point": _stack.point_stack, "grid": lambda stack: _stack.grid_stack(stack, 3),
                   "point jet": _stack.point_jet_stack, "grid jet": lambda stack: _stack.grid_jet_stack(stack, 3)}
        for name, read in readers.items():
            for part, (a, b) in enumerate(zip(read(plain), read(subclassed))):
                self.assertEqual(np.asarray(a).tobytes(), np.asarray(b).tobytes(), f"{name}, part {part}")


class Monitoring(unittest.TestCase):
    """The H layer of PAIR grows on the L layer."""

    def setUp(self):
        self.top = PAIR["layers"][0]
        self.finished = with_layers(PAIR, PAIR["layers"][1])

    def test_options(self):
        with self.assertRaisesRegex(TypeError, "thicknesses"):
            tmmcore.monitor_curve(self.finished, 550, top=2.35)
        with self.assertRaisesRegex(TypeError, r"top must be a number or \{ n, k \}"):
            tmmcore.monitor_curve(self.finished, 550, top=[2.35, 0], thicknesses=[0])
        with self.assertRaisesRegex(TypeError, "top has 'kappa'; a medium takes n and k"):
            tmmcore.monitor_curve(self.finished, 550, top={"n": 2.35, "kappa": 0.0}, thicknesses=[0])
        with self.assertRaisesRegex(TypeError, r"thicknesses must be a list of numbers \(nm\)"):
            tmmcore.monitor_curve(self.finished, 550, top=2.35, thicknesses=50)
        with self.assertRaisesRegex(TypeError, r"stack must be \{ incident, substrate, layers"):
            tmmcore.deposition_spectra(PAIR["layers"], GRID)
        with tmmcore.growing_stack(self.finished, GRID) as stack:
            with self.assertRaisesRegex(TypeError, "top n must be a number or 5 numbers, one per wavelength"):
                stack.set_top({"n": [2.35, 2.3]})
            with self.assertRaisesRegex(TypeError, r"d must be a number \(nm\)"):
                stack.sample("10")

    def test_the_grown_stack(self):
        curve = tmmcore.monitor_curve(self.finished, 550, top=self.top["n"], thicknesses=[0, self.top["d"]])
        self.assertEqual(curve["Rs"].shape, (2,))
        self.assertAlmostEqual(curve["Rs"][1], tmmcore.rta(PAIR, 550, pol="s")["R"], delta=1e-12)
        with tmmcore.growing_stack(self.finished, GRID) as stack:
            stack.set_top({"n": np.full(5, self.top["n"])})
            np.testing.assert_allclose(stack.sample(self.top["d"])["Rp"], tmmcore.rta_spectrum(PAIR, GRID, pol="p")["R"],
                                       rtol=0, atol=1e-12)

    def test_deposition_order(self):
        """Row k has the last k + 1 layers of the stack deposited."""
        N = len(PAIR["layers"])
        steps = tmmcore.deposition_spectra(PAIR, GRID, angle=30)
        self.assertEqual(steps["Rrs"].shape, (N, 5))
        for k in range(N):
            deposited = with_layers(PAIR, *PAIR["layers"][N - 1 - k:])
            for pol in ("s", "p"):
                rt = tmmcore.rta_spectrum(deposited, GRID, angle=30, pol=pol)
                for key in ("R", "T"):
                    np.testing.assert_allclose(steps[key + pol][k], rt[key], rtol=0, atol=1e-12, err_msg=f"{k} {key}{pol}")

    def test_top_layer_before_a_thickness(self):
        with tmmcore.growing_stack(PAIR, GRID, angle=30) as stack:
            self.assertEqual(stack.sample(0)["Rs"].shape, (5,))
            with self.assertRaisesRegex(ValueError, "set_top"):
                stack.sample(10)

    def test_closed(self):
        with tmmcore.growing_stack(PAIR, GRID, angle=30) as stack:
            stack.set_top(2.35)
            stack.sample(10)
        with self.assertRaisesRegex(ValueError, "closed"):
            stack.sample(10)
        stack.close()


if __name__ == "__main__":
    unittest.main()
