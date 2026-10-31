# Python

The Python package carries the same C kernel as the npm package, compiled for your platform and called through ctypes. The [API reference](api.md) gives each call in Python beside its JavaScript and C forms; this page covers what is particular to Python.

## Install

```bash
pip install tmmcore
```

Wheels are built for Windows (x86-64), macOS (Intel and Apple silicon) and Linux (x86-64 and ARM64, glibc), for Python 3.9 and later. The only dependency is numpy.

## Arguments and results

Each call takes the stack, then a wavelength in nm or, in a grid call, an array of them, and its options as keyword arguments. Wherever a call takes several values, such as the wavelengths of a grid call, an index or an angle per wavelength, or the thicknesses of `monitor_curve`, a list or a numpy array will do. A result is a dict of numbers and numpy arrays, whose keys and shapes the [API reference](api.md#results) gives for each call.

## The stack

A dict, the layers ordered from the incident medium toward the substrate:

```python
stack = {"incident": 1.0,
         "substrate": {"n": 1.52},
         "layers": [{"n": 2.35, "k": 0.0005, "d": 58.5},
                    {"n": 1.46, "d": 94.2}]}   # d in nm
```

`n` may also be complex, `n + ik`, with `k` left out: `{"n": 2.35 + 0.0005j, "d": 58.5}`.

## Jets

A jet is a complex array of shape (4,), the Taylor terms $[f,\ f',\ f''/2!,\ f'''/3!]$ of a function of angular frequency, and jets over a grid are an array of shape (…, 4). Every jet call broadcasts over the leading axes as numpy does, so one call builds a material's jets for a whole grid, as `cauchy_jet` in [Indices as jets](api.md#indices-as-jets) does.

## Agreement with the npm package

The package's tests run every call on the cases of the npm package's equivalence suite and compare each value with its WebAssembly build: R, T and A to within 1e-9, everything else to within 1e-7 of its size or 1e-12, whichever is larger. Both run the same C, so they differ only in the last digits, where the platform's math library rounds differently from the one compiled into the WebAssembly. These are the tolerances the npm package's tests hold its JavaScript kernels to against the same WebAssembly.

They also run the [comparison with Essential Macleod](validation.md#against-essential-macleod) on its committed results, to the same tolerances as the npm package.

## From source

The wheels build `src/c` with the platform's C compiler. To build and test a checkout, with Node for the reference results:

```bash
git clone https://github.com/aai2k/tmmcore
cd tmmcore
pip install ./python
node python/tests/reference.mjs
python -m unittest discover -s python/tests
```

`reference.mjs` runs the npm package on the equivalence cases and writes the results the tests compare with.
