"""
What test_reference.py cannot show: the arguments the package refuses, the
shapes it returns, complex n against n and k, arrays in another memory order
or type than the C reads, and the growing stack's lifetime.

    python -m unittest discover -s python/tests
"""

import unittest

import numpy as np

import tmmcore
from tmmcore._stack import grid_stack
from tmmcore.growing import _GrowingStack

# A high and a low index layer, each a quarter wave at 550 nm.
PAIR = {"incident": 1.0, "substrate": 1.52,
        "layers": [{"n": 2.35, "d": 58.5}, {"n": 1.46, "d": 94.2}]}
METAL = {"incident": 1.0, "substrate": 1.52, "layers": [{"n": 0.13, "k": 3.0, "d": 40.0}]}
GRID = np.linspace(450.0, 650.0, 5)


def with_layers(stack, *layers):
    return {**stack, "layers": list(layers)}


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

    def test_medium(self):
        with self.assertRaisesRegex(TypeError, r"incident must be a number or \{ n, k \}"):
            tmmcore.rta({**PAIR, "incident": [1.0, 0.0]}, 550)

    def test_complex_n(self):
        self.assertEqual(tmmcore.rta(with_layers(METAL, {"n": 0.13 + 3.0j, "d": 40.0}), 550), tmmcore.rta(METAL, 550))
        table = with_layers(METAL, {"n": np.full(5, 0.13 + 3.0j), "d": 40.0})
        np.testing.assert_array_equal(tmmcore.rta_spectrum(table, GRID)["R"], tmmcore.rta_spectrum(METAL, GRID)["R"])
        with self.assertRaisesRegex(TypeError, "layer 0 k goes inside n when n is complex"):
            tmmcore.rta(with_layers(METAL, {"n": 0.13 + 3.0j, "k": 3.0, "d": 40.0}), 550)

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


class GrowingStack(unittest.TestCase):
    def setUp(self):
        self.grid = (GRID, np.array([30.0]), *grid_stack(PAIR, len(GRID)))

    def test_top_layer_before_a_thickness(self):
        with _GrowingStack(self.grid) as stack:
            self.assertEqual(stack.sample(0)["Rs"].shape, (5,))
            with self.assertRaisesRegex(ValueError, "set_top"):
                stack.sample(10)

    def test_closed(self):
        with _GrowingStack(self.grid) as stack:
            stack.set_top(2.35)
            stack.sample(10)
        with self.assertRaisesRegex(ValueError, "closed"):
            stack.sample(10)
        stack.close()


if __name__ == "__main__":
    unittest.main()
