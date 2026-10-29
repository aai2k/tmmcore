"""
src/c/tmm_kernel.c, which includes the five sources so the kernel compiles as
one file, against the five compiled together with link-time optimization, as
the wheels and the WebAssembly are. Both are built here by the same compiler
with the wheels' flags, into python/build/one-file, and every call in
reference.json runs on each: the results must be equal, NaN where the other
is NaN. Skipped without a C compiler named cc or gcc on PATH.

    python -m unittest discover -s python/tests
"""

import ctypes
import importlib
import json
import os
import shutil
import subprocess
import sys
import unittest
from contextlib import ExitStack
from pathlib import Path
from unittest import mock

import numpy as np

from tmmcore import _library
from test_reference import REFERENCE, run

PYTHON = Path(__file__).resolve().parents[1]
SOURCES = PYTHON.parent / "src" / "c"
FOLDER = PYTHON / "build" / "one-file"
FIVE = ["rta.c", "growing.c", "derivatives.c", "jet.c", "phase.c"]
COMPILER = shutil.which("cc") or shutil.which("gcc")
FLAGS = ["-std=c99", "-O3", "-ffp-contract=off", "-shared", "-fPIC"]
# The modules that call the kernel. The package's names rta and phase are its
# functions, so the modules are looked up by their full names.
MODULES = [importlib.import_module("tmmcore." + name) for name in ("rta", "derivatives", "phase", "growing", "jet")]


def build(name, sources, flags=()):
    """A shared library of these sources, with the argument types the package sets."""
    FOLDER.mkdir(parents=True, exist_ok=True)
    path = FOLDER / (name + (".dll" if sys.platform == "win32" else ".so"))
    subprocess.run([COMPILER, *FLAGS, *flags, *(str(SOURCES / source) for source in sources),
                    "-o", str(path), "-lm"], check=True)
    library = ctypes.CDLL(str(path))
    for function, (result, arguments) in _library._SIGNATURES.items():
        getattr(library, function).restype = result
        getattr(library, function).argtypes = arguments
    return library


def run_on(library, case):
    """A reference call with every module of the package calling `library`."""
    with ExitStack() as patches:
        for module in MODULES:
            patches.enter_context(mock.patch.object(module, "lib", library))
        return run(case)


def equal(a, b):
    """Dicts key by key in order, lists entry by entry, arrays and numbers equal with NaN equal to NaN."""
    if isinstance(a, dict):
        return isinstance(b, dict) and list(a) == list(b) and all(equal(a[key], b[key]) for key in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(map(equal, a, b))
    if a is None or b is None:
        return a is b
    return np.array_equal(np.asarray(a), np.asarray(b), equal_nan=True)


@unittest.skipIf(COMPILER is None, "no C compiler named cc or gcc on PATH")
class OneFile(unittest.TestCase):
    def test_equal_to_the_link_time_build(self):
        if not os.path.exists(REFERENCE):
            self.fail("reference.json is missing: run `node python/tests/reference.mjs` from the repository root")
        with open(REFERENCE, encoding="utf-8") as file:
            calls = json.load(file)["calls"]
        one_file = build("one_file", ["tmm_kernel.c"])
        link_time = build("link_time", FIVE, ["-flto"])
        for i, case in enumerate(calls):
            with self.subTest(call=i, name=case["call"]):
                self.assertTrue(equal(run_on(one_file, case), run_on(link_time, case)))


if __name__ == "__main__":
    unittest.main()
