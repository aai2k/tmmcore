"""
Builds the wheel: the C kernel in ../src/c compiled as one shared library,
tmmcore/_tmmcore.dll on Windows and tmmcore/_tmmcore.so elsewhere, which the
package loads with ctypes. The library uses no Python API, so the wheel is
tagged py3-none-<platform> and one build serves every Python 3. The version is
the npm package's.
"""

import json
import os
import shutil
import sys

from setuptools import Extension, setup
from setuptools.command.bdist_wheel import bdist_wheel
from setuptools.command.build_ext import build_ext

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SOURCES = ["rta.c", "growing.c", "derivatives.c", "jet.c", "phase.c"]

with open(os.path.join(ROOT, "package.json"), encoding="utf-8") as file:
    VERSION = json.load(file)["version"]

# A wheel takes its license files from beside pyproject.toml only.
shutil.copyfile(os.path.join(ROOT, "LICENSE"), os.path.join(HERE, "LICENSE"))


class SharedLibrary(build_ext):
    """The kernel as a plain shared library: no module init function, and a file name ctypes can find."""

    def build_extension(self, ext):
        # Built as the WebAssembly is, all sources optimized together at link
        # time, and without fused multiply-add, which the WebAssembly lacks, so
        # the two round alike. MSVC's defaults already compile with /GL and
        # link with /LTCG, and contract nothing.
        if self.compiler.compiler_type != "msvc":
            ext.extra_compile_args = ["-std=c99", "-O3", "-ffp-contract=off", "-flto"]
            ext.extra_link_args = ["-O3", "-ffp-contract=off", "-flto"]
        super().build_extension(ext)

    def get_export_symbols(self, ext):
        return []

    def get_ext_filename(self, fullname):
        return os.path.join(*fullname.split(".")) + (".dll" if sys.platform == "win32" else ".so")


class PlatformWheel(bdist_wheel):
    """Tagged for any Python 3 on the platform it is built on."""

    def get_tag(self):
        return "py3", "none", super().get_tag()[2]


setup(
    version=VERSION,
    ext_modules=[Extension(
        "tmmcore._tmmcore",
        sources=[os.path.join("..", "src", "c", name) for name in SOURCES],
        libraries=[] if sys.platform == "win32" else ["m"],
    )],
    cmdclass={"build_ext": SharedLibrary, "bdist_wheel": PlatformWheel},
)
