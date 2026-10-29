"""
The package against the npm package's WebAssembly kernel on the cases of the
JavaScript equivalence suite: every call in reference.json run through the
package, its result compared value by value and its keys in order.

reference.json is written from the repository root by

    node python/tests/reference.mjs

Both sides run the same C, built natively here and to WebAssembly for npm.
Their math libraries differ, so the results agree to round-off, not bit for
bit: R, T and A within ABS_RTA, every other value within ABS or REL of its
size. These are the tolerances that tests/equivalence.mjs holds the
JavaScript kernels to against the same WebAssembly, for the same reason.

    python -m unittest discover -s python/tests
"""

import importlib
import json
import os
import unittest

import numpy as np

import tmmcore

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, "reference.json")

ABS_RTA = 1e-9  # R, T and A live in [0, 1]
REL = 1e-7      # derivatives and phase quantities span many orders of magnitude
ABS = 1e-12

# Keys whose values are reflectances, transmittances or absorptances.
RTA_KEYS = {"R", "T", "A", "Rs", "Ts", "As", "Rp", "Tp", "Ap", "Rrs", "Rrp", "base"}


def is_numeric(value):
    """A number, or a list of them at any depth; None stands for NaN."""
    if isinstance(value, list):
        return all(is_numeric(entry) for entry in value)
    return value is None or (isinstance(value, (int, float)) and not isinstance(value, bool))


def decode(value):
    """An argument as the package takes it: {"complex": pairs} as complex numbers, lists of numbers as arrays."""
    if isinstance(value, dict):
        if list(value) == ["complex"]:
            pairs = np.array(value["complex"], dtype=np.float64)
            numbers = pairs[..., 0] + 1j * pairs[..., 1]
            return numbers[()] if numbers.ndim == 0 else numbers
        return {key: decode(entry) for key, entry in value.items()}
    if isinstance(value, list):
        return np.array(value, dtype=np.float64) if is_numeric(value) else [decode(entry) for entry in value]
    return value


def run(case):
    module, _, name = case["call"].rpartition(".")
    function = getattr(importlib.import_module("tmmcore." + module) if module else tmmcore, name)
    args = [decode(arg) for arg in case["args"]]
    if name == "_GrowingStack":
        with function(*args) as stack:
            stack.set_top(decode(case["top"]))
            return [stack.sample(d) for d in case["samples"]]
    return function(*args, **{key: decode(value) for key, value in case["kwargs"].items()})


class Worst:
    """The largest differences seen, for the summary at the end."""
    rta = 0.0
    relative = 0.0
    values = 0
    differing = 0


def close(test, actual, expected, key, where):
    """Two real arrays of the same size within the tolerance of `key`, NaN where the other is NaN."""
    actual = np.ravel(np.asarray(actual, dtype=np.float64))
    expected = np.ravel(np.asarray(expected, dtype=np.float64))
    test.assertEqual(actual.shape, expected.shape, f"{where}: size")
    test.assertTrue(np.array_equal(np.isnan(actual), np.isnan(expected)), f"{where}: NaN where the other is not")
    finite = ~np.isnan(expected)
    actual, expected = actual[finite], expected[finite]
    diff = np.abs(actual - expected)
    Worst.values += diff.size
    Worst.differing += int(np.count_nonzero(diff))
    if key in RTA_KEYS:
        Worst.rta = max(Worst.rta, float(diff.max(initial=0)))
        bad = diff > ABS_RTA
    else:
        size = np.maximum(np.abs(actual), np.abs(expected))
        above = diff > ABS
        if above.any():
            Worst.relative = max(Worst.relative, float((diff[above] / size[above]).max()))
        bad = above & (diff > REL * size)
    if bad.any():
        at = int(np.argmax(bad))
        test.fail(f"{where}: {actual[at]!r} against {expected[at]!r}")


def compare(test, actual, expected, where, key=""):
    """A result against its reference: dicts key by key in order, lists entry by entry, numbers with close()."""
    if expected is None and actual is None:
        return
    if isinstance(expected, dict) and list(expected) == ["complex"]:
        numbers = decode(expected)
        close(test, np.real(actual), np.real(numbers), key, f"{where} real")
        close(test, np.imag(actual), np.imag(numbers), key, f"{where} imaginary")
    elif isinstance(expected, dict):
        test.assertIsInstance(actual, dict, where)
        test.assertEqual(list(actual), list(expected), f"{where}: keys")
        for name, entry in expected.items():
            compare(test, actual[name], entry, f"{where}.{name}", name)
    elif isinstance(expected, list) and not is_numeric(expected):
        test.assertEqual(len(actual), len(expected), f"{where}: length")
        for i, entry in enumerate(expected):
            compare(test, actual[i], entry, f"{where}[{i}]", key)
    elif isinstance(expected, bool):
        test.assertIsInstance(actual, (bool, np.bool_), where)
        test.assertEqual(bool(actual), expected, where)
    else:
        test.assertNotIsInstance(actual, (dict, type(None)), where)
        close(test, actual, np.nan if expected is None else expected, key, where)


class AgainstWebAssembly(unittest.TestCase):
    def test_every_reference_call(self):
        if not os.path.exists(REFERENCE):
            self.fail("reference.json is missing: run `node python/tests/reference.mjs` from the repository root")
        with open(REFERENCE, encoding="utf-8") as file:
            calls = json.load(file)["calls"]
        for i, case in enumerate(calls):
            with self.subTest(call=i, name=case["call"]):
                compare(self, run(case), case["expected"], f"{case['call']} #{i}")
        print(f"\n{len(calls)} calls, {Worst.values} values, {Worst.differing} differ from the WebAssembly; "
              f"worst R/T/A {Worst.rta:.2e}, worst relative above {ABS} {Worst.relative:.2e}")


if __name__ == "__main__":
    unittest.main()
