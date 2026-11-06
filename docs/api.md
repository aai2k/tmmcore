# API reference

tmmcore has the same calls in JavaScript, Python and C. The R/T, derivative, needle, phase and monitoring calls take a stack, a wavelength in nm or, in a grid call, a list of them, and options. The grid calls are `depositionSpectra`, `growingStack` and the calls whose names end in `Spectrum`, `_spectrum` in Python and C. Each call below has a tab for each language; choosing a language in one tab switches every tab on the page.

=== "JavaScript"

    ```js
    import { rta, rtaSpectrum } from 'tmmcore';

    const stack = {
        incident: 1,
        substrate: 1.52,
        layers: [
            { n: 2.35, d: 58.5 },   // d in nm
            { n: 1.46, d: 94.2 },
        ],
    };
    const wavelengths = [450, 500, 550, 600, 650];

    rta(stack, 550).R;                                            // natural light, normal incidence
    rtaSpectrum(stack, wavelengths, { angle: 45, pol: 's' }).R;   // one value per wavelength
    ```

=== "Python"

    ```python
    import numpy as np
    import tmmcore

    stack = {
        "incident": 1.0,
        "substrate": 1.52,
        "layers": [
            {"n": 2.35, "d": 58.5},   # d in nm
            {"n": 1.46, "d": 94.2},
        ],
    }
    wavelengths = np.array([450.0, 500.0, 550.0, 600.0, 650.0])

    tmmcore.rta(stack, 550)["R"]                                       # natural light, normal incidence
    tmmcore.rta_spectrum(stack, wavelengths, angle=45, pol="s")["R"]   # one value per wavelength
    ```

=== "C"

    ```c
    #include "tmmcore.h"

    /* n, k and d in nm of each layer, incident side first */
    const double layers[] = {
        2.35, 0, 58.5,
        1.46, 0, 94.2,
    };
    double s[3], p[3];

    tmm_rta(550, 0, 0, 1, 0, 1.52, 0, layers, 2, s);   /* R, T and A in s */
    tmm_rta(550, 0, 1, 1, 0, 1.52, 0, layers, 2, p);   /* and in p */
    double R = (s[0] + p[0]) / 2;                       /* natural light */
    ```

| JavaScript | Python | C | Returns |
|---|---|---|---|
| [`rta`, `rtaSpectrum`](#r-t-and-a) | `rta`, `rta_spectrum` | `tmm_rta`, `tmm_rta_spectrum` | R, T and A |
| [`jacobian`, `jacobianSpectrum`](#thickness-derivatives) | `jacobian`, `jacobian_spectrum` | `tmm_jacobian`, `tmm_jacobian_spectrum` | R, T and A, and their derivatives in each layer's thickness |
| [`hessian`, `hessianSpectrum`](#thickness-derivatives) | `hessian`, `hessian_spectrum` | `tmm_hessian`, `tmm_hessian_spectrum` | the same, and the second derivatives |
| [`needleScan`, `needleScanSpectrum`](#needle-scan) | `needle_scan`, `needle_scan_spectrum` | `tmm_needle_scan`, `tmm_needle_scan_spectrum` | the derivatives of R, T and A for inserting a thin layer anywhere in the stack |
| [`phase`, `phaseSpectrum`](#phase-gd-gdd-and-tod) | `phase`, `phase_spectrum` | `tmm_phase`, `tmm_phase_spectrum` | phase, GD, GDD and TOD of r and t |
| [`phaseJacobian`, `phaseJacobianSpectrum`](#phase-gd-gdd-and-tod) | `phase_jacobian`, `phase_jacobian_spectrum` | `tmm_phase_jacobian`, `tmm_phase_jacobian_spectrum` | the same, and their thickness derivatives |
| [`jetConstant` … `jetIsFinite`](#jet-arithmetic) | `jet_constant` … `jet_is_finite` | `tmm_jet_constant` … `tmm_jet_is_finite` | Taylor-jet arithmetic, which builds the dispersive indices the phase calls take |
| [`monitorCurve`, `depositionSpectra`, `growingStack`](#monitoring-and-deposition) | `monitor_curve`, `deposition_spectra`, `growing_stack` | `tmm_monitor_curve`, `tmm_deposition_spectra`, `tmm_growing_stack_create` | R and T in s and p while the layers are deposited |

In JavaScript the calls also run on a [WebAssembly instance](#webassembly). [Python](python.md) covers installing the package and what is particular to it, and [C](#c) building the kernel into a program. The sign and phase conventions are in [Getting started](getting-started.md#conventions).

## Arguments

### The stack

The incident medium, the substrate, and the layers between them, ordered from the incident medium toward the substrate.

=== "JavaScript"

    `{ incident, substrate, layers }`:

    | | |
    |---|---|
    | `layers` | `{ n, k, d }[]`, with `d` the thickness in nm |
    | `incident`, `substrate` | a number `n`, or `{ n, k }` |
    | `n`, `k` | a number, or in a grid call a list with one value per wavelength |

    A medium given as a bare list is refused, since an `[re, im]` pair would read as two wavelengths.

=== "Python"

    A dict with the keys of the JavaScript object, `{"incident", "substrate", "layers"}`, and each layer `{"n", "k", "d"}` with `d` the thickness in nm. A medium is a number or `{"n", "k"}`. In a grid call `n` and `k` may be arrays with one value per wavelength, and `n` may be [complex](python.md#the-stack) in place of `k`.

=== "C"

    A call at one wavelength takes the incident medium and the substrate as `n0_re, n0_im, ns_re, ns_im`, and the layers as `N` triples:

    ```c
    const double layers[] = {
        2.35, 0, 58.5,   /* n, k, d in nm */
        1.46, 0, 94.2,
    };
    ```

    A grid call takes every index at every wavelength:

    | | |
    |---|---|
    | `lambdas` | the `nLam` wavelengths, in nm |
    | `n0arr`, `nsarr` | the incident medium and the substrate, `nLam` × 2: `[n, k]` at each wavelength |
    | `matNK` | the layers, `N` × `nLam` × 2: `[layer][wavelength][n, k]` |
    | `thick` | the `N` thicknesses, in nm |

$\tilde n = n + ik$, with $k \ge 0$ for absorbing media. In JavaScript and Python `k` may be left out and is then 0. A layer of zero, negative or NaN thickness adds nothing to the spectrum.

In the phase calls an index is a number, which is a constant index with no dispersion, or a [jet](#indices-as-jets) with `k` inside it: one jet at one wavelength, one jet per wavelength in a grid call. A list of plain numbers per wavelength is refused there, since the group delay computed from it would leave out the material's dispersion. In C every index of a phase call is a jet of 8 doubles, and the thicknesses are in `thick`. At one wavelength the indices are `n0jet`, `nsjet` and `layerJets` (N × 8); over a grid they are `n0jets` and `nsjets` (nLam × 8 each) and `matJets` (N × nLam × 8).

### Options

| JavaScript, Python | | C |
|---|---|---|
| `angle` | Angle of incidence in degrees from the normal, 0 when left out. A grid call also takes one angle per wavelength. | `theta_deg`; in a grid call `nTheta` of them, 1 or `nLam` |
| `sine` | Phase calls only, in place of `angle`: the sine of the angle of incidence as a [jet](#indices-as-jets), one jet per wavelength in a grid call. For [an angle fixed outside the incident medium](#an-angle-fixed-outside-the-incident-medium). | `sinJet`, in a grid call `sinJets`, `nLam` × 8; `NULL` to use the angle |
| `pol` | In the R/T, derivative and needle calls `'s'`, `'p'` or `'avg'`, the mean of the two, which is the default. In the phase calls `'s'` or `'p'`; it may be left out at normal incidence, where the two are the same. The monitoring calls take none and return both. | `pol`, 0 for s and 1 for p. For `'avg'`, call with each and take the mean of the results. |
| `candidates` | Needle calls only, and required there: the materials that could be inserted, a list of indices, each a number or `{ n, k }` as a medium is. | `candNs`, `nCand` × 2, `[n, k]` of each; in a grid call `candNK`, `nCand` × `nLam` × 2 |
| `fractions` | Needle calls only: positions inside each layer, as fractions of its thickness from the incident side, where the needle is also tried. None when left out. | `fracs`, `nFrac` of them |
| `top` | `monitorCurve` only, and required there: the index of the growing layer, a number or `{ n, k }` as a medium is. | `ng_re`, `ng_im` |
| `thicknesses` | `monitorCurve` only, and required there: the thicknesses of the growing layer to sample, in nm. | `dArr`, `nD` of them |

=== "JavaScript"

    An object after the wavelength: `rta(stack, 550, { angle: 45, pol: 's' })`.

=== "Python"

    Keyword arguments: `tmmcore.rta(stack, 550, angle=45, pol="s")`.

=== "C"

    Arguments in one order in every call: the wavelength or the grid, the angle, `pol`, the incident medium, the substrate, the layers, any further inputs, then the outputs.

### Results

Each call below gives the layout of what it returns. `N`, where a result has it, is the number of layers.

=== "JavaScript"

    An object. Over a grid a value such as `R` or `gd` is a `Float64Array` with one entry per wavelength, and a first derivative is flat and wavelength-major: wavelength `i`, layer `k` at `[i * N + k]`.

=== "Python"

    A dict with the keys of the JavaScript object, except that the phase keys are in snake_case. A value such as `R` or `gd` is a number at one wavelength and an array over a grid. Derivatives, needle scans and monitoring results are numpy arrays with an axis for each index: a first derivative over a grid has the shape (nLam, N), wavelength `i` and layer `k` at `[i, k]`.

=== "C"

    The call writes into output arrays the caller allocates, in the order its prototype gives them. Each call below gives their sizes and layouts.

### Errors

=== "JavaScript"

    A wrong argument throws a `TypeError` that names it, such as `layer 0 needs {n, d}`, `pol must be 's', 'p' or 'avg'` or `unknown option 'polarization'; the options are angle and pol`. So does a key the stack, a layer or a medium does not take, so a misspelled `k` cannot pass as a lossless layer: `layer 0 has 'kappa'; a layer takes n, k and d`.

=== "Python"

    A wrong argument raises a `TypeError` that names it, such as `layer 0 needs {n, d}`, and so does a key the stack, a layer or a medium does not take: `layer 0 has 'kappa'; a layer takes n, k and d`. An unknown keyword, or a needle scan without `candidates`, raises Python's own `TypeError`.

=== "C"

    The calls do not check their arguments, so each count has to match its array. A call that cannot allocate its working memory writes NaN to its outputs; `tmm_deposition_spectra` returns 0 instead, and `tmm_growing_stack_create` returns `NULL`.

## R, T and A

### `rta`

=== "JavaScript"

    `rta(stack, wavelength, options?)`

    ```js
    const { R, T, A } = rta(stack, 550, { pol: 's' });
    ```

=== "Python"

    `rta(stack, wavelength, *, angle=0.0, pol="avg")`

    ```python
    result = tmmcore.rta(stack, 550, pol="s")   # result["R"], result["T"], result["A"]
    ```

=== "C"

    ```c
    void tmm_rta(double lambda_nm, double theta_deg, int pol,
                 double n0_re, double n0_im, double ns_re, double ns_im,
                 const double *layers, int N, double *out);
    ```

    `out` takes R, T and A.

**Returns** `R`, `T` and `A`: reflectance, transmittance into the substrate, and the absorptance of the layers, each a fraction of the incident irradiance.

$A$ is the net irradiance entering the front of the stack less the irradiance leaving it into the substrate (Macleod, Eq. 2.120), so it is what the layers absorb. With a transparent incident medium $R + T + A = 1$ by construction, so the sum is no check on the calculation.

With an absorbing incident medium they do not sum to one. The incident and reflected waves interfere in the irradiance of the incident medium, and

$$R + T + A = 1 + 2\,\frac{\mathrm{Im}\,\eta_0}{\mathrm{Re}\,\eta_0}\,\mathrm{Im}(r)$$

in tmmcore's sign convention. $T$ can then exceed one: Macleod finds $T = 1 + k_0^2/n_0^2$ for a surface whose admittance is the complex conjugate of the incident one (Eq. 2.84). See [Validation](validation.md#an-absorbing-incident-medium).

### `rtaSpectrum`

R, T and A over a wavelength grid, s and p in one pass.

=== "JavaScript"

    `rtaSpectrum(stack, wavelengths, options?)`

    ```js
    const { R, Rs, Rp } = rtaSpectrum(stack, wavelengths, { angle: 45 });
    ```

=== "Python"

    `rta_spectrum(stack, wavelengths, *, angle=0.0, pol="avg")`

    ```python
    spectrum = tmmcore.rta_spectrum(stack, wavelengths, angle=45)   # spectrum["Rs"], spectrum["Rp"]
    ```

=== "C"

    ```c
    void tmm_rta_spectrum(const double *lambdas, int nLam,
                          const double *theta_deg, int nTheta,
                          const double *n0arr, const double *nsarr,
                          const double *matNK, const double *thick, int N,
                          double *outRs, double *outTs, double *outAs,
                          double *outRp, double *outTp, double *outAp);
    ```

    It takes no `pol` and writes s and p, `nLam` values to each output.

**Returns** `R`, `T` and `A` for the `pol` asked, and the s and p values they come from: `Rs`, `Ts`, `As`, `Rp`, `Tp` and `Ap`.

## Thickness derivatives

The derivatives come from the same characteristic-matrix product as the spectrum, exactly, not by finite differences. They are the derivatives of what `rta` computes, [opaque layers](#opaque-layers) included, and [Validation](validation.md#derivatives-against-finite-differences) checks them against finite differences.

They line up with the layers as given. A layer of zero thickness keeps its derivative, the rate at which the spectrum moves as that layer starts to grow. A layer of negative or NaN thickness is absent, and its derivatives are zero. With `pol: 'avg'` the derivatives are the mean of the s and p ones.

### `jacobian`

=== "JavaScript"

    `jacobian(stack, wavelength, options?)`

    ```js
    import { jacobian } from 'tmmcore';

    const { R, dRdd } = jacobian(stack, 550, { pol: 's' });   // dRdd[j] = ∂R/∂d_j, in 1/nm
    ```

=== "Python"

    `jacobian(stack, wavelength, *, angle=0.0, pol="avg")`

    ```python
    dRdd = tmmcore.jacobian(stack, 550, pol="s")["dRdd"]   # dRdd[j] = ∂R/∂d_j, in 1/nm
    ```

=== "C"

    ```c
    void tmm_jacobian(double lambda_nm, double theta_deg, int pol,
                      double n0_re, double n0_im, double ns_re, double ns_im,
                      const double *layers, int N,
                      double *dRdd, double *dTdd, double *dAdd, double *base);
    ```

    `base` takes R, T and A.

**Returns** `R`, `T`, `A` and `N`, and `dRdd`, `dTdd` and `dAdd`, the derivatives in each layer's thickness, `N` values each.

The spectrum comes with its derivatives, so a Gauss-Newton or Levenberg-Marquardt step costs one evaluation rather than $N + 1$. See the [refinement example](examples.md#refinement-with-the-analytic-jacobian).

### `jacobianSpectrum`

The Jacobian at every wavelength of a grid.

=== "JavaScript"

    `jacobianSpectrum(stack, wavelengths, options?)`

    The derivatives are flat: wavelength `i`, layer `k` at `[i * N + k]`.

=== "Python"

    `jacobian_spectrum(stack, wavelengths, *, angle=0.0, pol="avg")`

    The derivatives have the shape (nLam, N): wavelength `i`, layer `k` at `[i, k]`.

=== "C"

    ```c
    void tmm_jacobian_spectrum(const double *lambdas, int nLam,
                               const double *theta_deg, int nTheta, int pol,
                               const double *n0arr, const double *nsarr,
                               const double *matNK, const double *thick, int N,
                               double *dRdd, double *dTdd, double *dAdd, double *base);
    ```

    The derivatives are nLam × N, `[wavelength][layer]`, and `base` takes R, T and A at each wavelength, nLam × 3.

**Returns** R, T and A over the grid, `dRdd`, `dTdd` and `dAdd` with nLam × N values each, and `N`.

### `hessian`

=== "JavaScript"

    `hessian(stack, wavelength, options?)`

    ```js
    import { hessian } from 'tmmcore';

    const { d2Rdd } = hessian(stack, 550, { pol: 's' });   // d2Rdd[i][j] = ∂²R/∂d_i∂d_j
    ```

    The second derivatives are `N` rows of `N`.

=== "Python"

    `hessian(stack, wavelength, *, angle=0.0, pol="avg")`

    ```python
    d2Rdd = tmmcore.hessian(stack, 550, pol="s")["d2Rdd"]   # d2Rdd[i, j] = ∂²R/∂d_i∂d_j
    ```

=== "C"

    ```c
    void tmm_hessian(double lambda_nm, double theta_deg, int pol,
                     double n0_re, double n0_im, double ns_re, double ns_im,
                     const double *layers, int N,
                     double *dRdd, double *dTdd, double *dAdd,
                     double *d2Rdd, double *d2Tdd, double *d2Add, double *base);
    ```

    The second derivatives are `N` × `N`, row-major.

**Returns** the fields of `jacobian`, and the second derivatives `d2Rdd`, `d2Tdd` and `d2Add`, N × N each and symmetric. The cost is $O(N^2)$ small-matrix operations. See the [tolerance example](examples.md#what-thickness-errors-cost-from-the-exact-hessian).

### `hessianSpectrum`

The Hessian at every wavelength of a grid.

=== "JavaScript"

    `hessianSpectrum(stack, wavelengths, options?)`

    The second derivatives are flat: wavelength `i` and layers `j`, `k` at `[(i * N + j) * N + k]`.

=== "Python"

    `hessian_spectrum(stack, wavelengths, *, angle=0.0, pol="avg")`

    The second derivatives have the shape (nLam, N, N): wavelength `i` and layers `j`, `k` at `[i, j, k]`.

=== "C"

    ```c
    void tmm_hessian_spectrum(const double *lambdas, int nLam,
                              const double *theta_deg, int nTheta, int pol,
                              const double *n0arr, const double *nsarr,
                              const double *matNK, const double *thick, int N,
                              double *dRdd, double *dTdd, double *dAdd,
                              double *d2Rdd, double *d2Tdd, double *d2Add, double *base);
    ```

    The second derivatives are nLam × N × N, `[wavelength][layer][layer]`.

**Returns** the fields of `jacobianSpectrum`, and `d2Rdd`, `d2Tdd` and `d2Add` with nLam × N × N values each.

## Needle scan

The derivatives of R, T and A for inserting a layer of zero thickness of a candidate material, at every interface of the stack at once and, if asked, inside its layers. This is the analytic P-function of Tikhonravov et al. (*Appl. Opt.* **35**, 5493, 1996), the $d \to 0$ limit of the numerical method of Sullivan and Dobrowolski (*Appl. Opt.* **35**, 5484, 1996), so it needs no trial thickness and no second evaluation.

### `needleScan`

=== "JavaScript"

    `needleScan(stack, wavelength, options)`

    ```js
    import { needleScan } from 'tmmcore';

    const scan = needleScan(stack, 550, { pol: 's', candidates: [2.35, 1.46] });
    // scan.gaps[position][candidate] = { dR, dT, dA }
    ```

    `intra[k][i]` is `{ frac, perCand }` for a needle inside layer `k` at fraction `fractions[i]`, with `perCand` shaped like an entry of `gaps`.

=== "Python"

    `needle_scan(stack, wavelength, *, candidates, fractions=(), angle=0.0, pol="avg")`

    ```python
    scan = tmmcore.needle_scan(stack, 550, pol="s", candidates=[2.35, 1.46])
    # scan["gaps"]["dR"][position, candidate]
    ```

    `gaps` and `intra` are each a dict of `dR`, `dT` and `dA`, arrays of shape (N + 1, nCand) in `gaps` and (N, nFrac, nCand) in `intra`.

=== "C"

    ```c
    void tmm_needle_scan(double lambda_nm, double theta_deg, int pol,
                         double n0_re, double n0_im, double ns_re, double ns_im,
                         const double *layers, int N,
                         const double *candNs, int nCand,
                         const double *fracs, int nFrac,
                         double *base, double *gaps, double *intra);
    ```

    `base` takes R, T and A. `gaps` is (N + 1) × nCand × 3, `[position][candidate][dR, dT, dA]`. `intra` is N × nFrac × nCand × 3, `[layer][fraction][candidate][dR, dT, dA]`, and is not written when `nFrac` is 0.

**Returns** `R`, `T`, `A`, `gaps`, `intra` and `N`.

- `gaps` holds `dR`, `dT` and `dA`, per nm of needle, for each candidate at each of the $N + 1$ interfaces: 0 is in front of the first layer and `N` against the substrate. Insert where the derivative of your merit function is negative.
- `intra` holds the same for a needle inside each layer at each of `fractions`. Fractions 0 and 1 give the interfaces either side of the layer.

The positions count the layers as given, zero-thickness ones included. A needle of a layer's own material at its boundary is that layer made thicker, so `dR` in `gaps` at interface 0 for a candidate with layer 0's index equals `dRdd[0]` from `jacobian`; the two are computed independently. With `pol: 'avg'` the derivatives are the mean of the s and p ones. See the [needle example](examples.md#the-needle-p-function-where-the-next-layer-goes).

### `needleScanSpectrum`

The needle scan at every wavelength of a grid, where a candidate's `n` and `k` may also be one value per wavelength.

=== "JavaScript"

    `needleScanSpectrum(stack, wavelengths, options)`

    `gaps` and `intra` are each `{ dR, dT, dA }` of flat, wavelength-major lists: wavelength `i`, interface `p` and candidate `c` at `[(i * (N + 1) + p) * nCand + c]` in `gaps`, and wavelength `i`, layer `k`, fraction `f` and candidate `c` at `[((i * N + k) * nFrac + f) * nCand + c]` in `intra`.

=== "Python"

    `needle_scan_spectrum(stack, wavelengths, *, candidates, fractions=(), angle=0.0, pol="avg")`

    `gaps` and `intra` hold arrays of shape (nLam, N + 1, nCand) and (nLam, N, nFrac, nCand).

=== "C"

    ```c
    void tmm_needle_scan_spectrum(const double *lambdas, int nLam,
                                  const double *theta_deg, int nTheta, int pol,
                                  const double *n0arr, const double *nsarr,
                                  const double *matNK, const double *thick, int N,
                                  const double *candNK, int nCand,
                                  const double *fracs, int nFrac,
                                  double *base, double *gaps, double *intra);
    ```

    The outputs are those of `tmm_needle_scan` at each wavelength in turn: `base` nLam × 3, `gaps` nLam × (N + 1) × nCand × 3 and `intra` nLam × N × nFrac × nCand × 3.

**Returns** R, T and A over the grid, `gaps`, `intra` and `N`.

## Phase, GD, GDD and TOD

Phase, group delay, group delay dispersion and third-order dispersion of r and t. They depend on how each index changes with frequency, so these calls take an index as a number, which has no dispersion, or as a [jet](#indices-as-jets) that carries it.

GD, GDD and TOD are exact frequency derivatives at the wavelength asked, from the frequency derivatives of the characteristic matrix. No phase unwrapping and no wavelength stencil is involved, so a value does not depend on the grid around it.

### `phase`

=== "JavaScript"

    `phase(stack, wavelength, options?)`

    ```js
    import { phase } from 'tmmcore';

    const { r, t } = phase(stack, 800, { pol: 's' });   // r.gdd: GDD on reflection, in fs²
    ```

    `r` and `t` are each `{ phaseRad, phaseDeg, gd, gdd, tod, magnitudeSquared }`, or `null`.

=== "Python"

    `phase(stack, wavelength, *, pol=None, angle=None, sine=None)`

    ```python
    r = tmmcore.phase(stack, 800, pol="s")["r"]   # r["gdd"]: GDD on reflection, in fs²
    ```

    `r` and `t` are each a dict of `phase_rad`, `phase_deg`, `gd`, `gdd`, `tod` and `magnitude_squared`, or `None`.

=== "C"

    ```c
    void tmm_phase(double lambda, double theta_deg, int pol,
                   const double *n0jet, const double *nsjet,
                   const double *layerJets, const double *thick, int N,
                   double omega, const double *sinJet, double *out);
    ```

    `omega` is the angular frequency at `lambda`, from `tmm_omega_from_lambda_nm`. `out` takes 10 values: the phase in radians, GD, GDD, TOD and the squared magnitude of r, then the same of t, all five NaN for a side that is `null`.

**Returns** the phase, GD, GDD, TOD and squared magnitude of `r` and of `t`, in the units of [Units and signs](#units-and-signs). A side is `null`, `None` in Python, where that coefficient is exactly zero and has no phase.

Both are also `null` at the exact critical angle of a layer, and at the substrate's critical angle in p. At the substrate's critical angle in s the phase and the squared magnitude come back, and GD, GDD and TOD are NaN.

### `phaseSpectrum`

`phase` at every wavelength of a grid.

=== "JavaScript"

    `phaseSpectrum(stack, wavelengths, options?)`

    `r` and `t` are each `{ phaseRad, gd, gdd, tod, magnitudeSquared }`.

=== "Python"

    `phase_spectrum(stack, wavelengths, *, pol=None, angle=None, sine=None)`

    `r` and `t` are each a dict of `phase_rad`, `gd`, `gdd`, `tod` and `magnitude_squared`.

=== "C"

    ```c
    void tmm_phase_spectrum(const double *lambdas, int nLam,
                            const double *theta_deg, int nTheta, int pol,
                            const double *n0jets, const double *nsjets,
                            const double *matJets, const double *thick, int N,
                            const double *omegas, const double *sinJets, double *out);
    ```

    `omegas` is the angular frequency at each wavelength, and `out` takes the 10 values of `tmm_phase` at each, nLam × 10.

**Returns** `r` and `t`, each value with one entry per wavelength, NaN where the coefficient vanishes.

### `phaseJacobian`

The phase quantities and their exact derivatives in every layer's thickness.

=== "JavaScript"

    `phaseJacobian(stack, wavelength, options?)`

    ```js
    import { phaseJacobian } from 'tmmcore';

    const { r } = phaseJacobian(stack, 800, { pol: 's' });   // r.dGdd[j] = ∂GDD/∂d_j, in fs²/nm
    ```

    `r` and `t` each add `dPhaseDeg`, `dGd`, `dGdd`, `dTod` and `dLogMagnitudeSquared` to the fields of `phase`.

=== "Python"

    `phase_jacobian(stack, wavelength, *, pol=None, angle=None, sine=None)`

    ```python
    r = tmmcore.phase_jacobian(stack, 800, pol="s")["r"]   # r["d_gdd"][j] = ∂GDD/∂d_j, in fs²/nm
    ```

    `r` and `t` each add `d_phase_deg`, `d_gd`, `d_gdd`, `d_tod` and `d_log_magnitude_squared` to the fields of `phase`.

=== "C"

    ```c
    void tmm_phase_jacobian(double lambda, double theta_deg, int pol,
                            const double *n0jet, const double *nsjet,
                            const double *layerJets, const double *thick, int N,
                            double omega, const double *sinJet, double *out, double *deriv);
    ```

    `out` is that of `tmm_phase`. `deriv` is 2 × 5 × N, `[side][quantity][layer]`: r, then t, and for each the derivatives of the phase in degrees, GD, GDD, TOD and the log of the squared magnitude.

**Returns** the fields of `phase`, and the derivatives of the phase in degrees, GD, GDD and TOD, and of ln |r|² or ln |t|², the relative change in intensity, `N` values each. The derivatives line up with the layers as in the [thickness derivatives](#thickness-derivatives), and a GDD fit costs one evaluation per step rather than $N + 1$.

### `phaseJacobianSpectrum`

`phaseJacobian` at every wavelength of a grid, for fitting measured phase quantities such as group delay or ellipsometric Ψ and Δ.

=== "JavaScript"

    `phaseJacobianSpectrum(stack, wavelengths, options?)`

    The derivatives are flat: wavelength `i`, layer `k` at `[i * N + k]`.

=== "Python"

    `phase_jacobian_spectrum(stack, wavelengths, *, pol=None, angle=None, sine=None)`

    The derivatives have the shape (nLam, N): wavelength `i`, layer `k` at `[i, k]`.

=== "C"

    ```c
    void tmm_phase_jacobian_spectrum(const double *lambdas, int nLam,
                                     const double *theta_deg, int nTheta, int pol,
                                     const double *n0jets, const double *nsjets,
                                     const double *matJets, const double *thick, int N,
                                     const double *omegas, const double *sinJets,
                                     double *out, double *deriv);
    ```

    `out` is that of `tmm_phase_spectrum`, and `deriv` takes the `deriv` of `tmm_phase_jacobian` at each wavelength, nLam × 2 × 5 × N.

**Returns** `r` and `t`, each the five values of `phaseSpectrum` and the five derivatives with nLam × N values each.

### Indices as jets

A jet is a truncated Taylor series in angular frequency, $[f,\ f',\ f''/2!,\ f'''/3!]$, each entry a complex number. An index jet is $[\tilde n,\ d\tilde n/d\omega,\ (d^2\tilde n/d\omega^2)/2,\ (d^3\tilde n/d\omega^3)/6]$.

tmmcore has no material models. You build the jet from your own dispersion formula with the [jet arithmetic](#jet-arithmetic). A Cauchy index:

=== "JavaScript"

    ```js
    import { phaseSpectrum, wavelengthOmegaJet, jetReciprocal, jetMultiply,
             jetScale, jetAdd, jetConstant, omegaFromLambdaNm } from 'tmmcore';

    // n(λ) = A + B/λ² + C/λ⁴ with a constant k, λ in nm
    function cauchyJet(lambda_nm, A, B, C, k) {
        const lambda = wavelengthOmegaJet(lambda_nm, omegaFromLambdaNm(lambda_nm));   // λ(ω)
        const inverseSquare = jetReciprocal(jetMultiply(lambda, lambda));
        return jetAdd(jetConstant(A, k),
            jetAdd(jetScale(inverseSquare, B),
                   jetScale(jetMultiply(inverseSquare, inverseSquare), C)));
    }

    const grid = [700, 750, 800, 850, 900];
    const glass = { n: grid.map(l => cauchyJet(l, 1.5, 4200, 0, 0)) };
    const layers = [{ n: grid.map(l => cauchyJet(l, 2.21, 2.4e4, 0, 0)), d: 85 }];
    const { r } = phaseSpectrum({ incident: 1, substrate: glass, layers }, grid, { pol: 's' });
    ```

=== "Python"

    ```python
    def cauchy_jet(lam, a, b, c, k):
        """n(λ) = a + b/λ² + c/λ⁴ with a constant k, λ in nm, as a jet at each wavelength."""
        lam_jet = tmmcore.wavelength_omega_jet(lam, tmmcore.omega_from_lambda_nm(lam))   # λ(ω)
        inverse_square = tmmcore.jet_reciprocal(tmmcore.jet_multiply(lam_jet, lam_jet))
        return tmmcore.jet_add(tmmcore.jet_constant(a, k),
                               tmmcore.jet_add(tmmcore.jet_scale(inverse_square, b),
                                               tmmcore.jet_scale(tmmcore.jet_multiply(inverse_square, inverse_square), c)))

    grid = np.array([700.0, 750.0, 800.0, 850.0, 900.0])
    glass = {"n": cauchy_jet(grid, 1.5, 4200, 0, 0)}
    layers = [{"n": cauchy_jet(grid, 2.21, 2.4e4, 0, 0), "d": 85}]
    r = tmmcore.phase_spectrum({"incident": 1.0, "substrate": glass, "layers": layers}, grid, pol="s")["r"]
    ```

=== "C"

    ```c
    #define NLAM 5
    static const double grid[NLAM] = { 700, 750, 800, 850, 900 };

    /* n(λ) = a + b/λ² + c/λ⁴ with a constant k, λ in nm, one jet per wavelength */
    static void cauchy_jet(double a, double b, double c, double k, double *jets)
    {
        double omega[NLAM], av[NLAM], kv[NLAM], bv[NLAM], cv[NLAM];
        double inverse_square[8 * NLAM], term[8 * NLAM];
        for (int i = 0; i < NLAM; i++) { av[i] = a; kv[i] = k; bv[i] = b; cv[i] = c; }
        tmm_omega_from_lambda_nm(grid, NLAM, omega);
        tmm_wavelength_omega_jet(grid, omega, NLAM, inverse_square);              /* λ(ω) */
        tmm_jet_multiply(inverse_square, inverse_square, NLAM, inverse_square);
        tmm_jet_reciprocal(inverse_square, NLAM, inverse_square);                /* 1/λ² */
        tmm_jet_multiply(inverse_square, inverse_square, NLAM, term);
        tmm_jet_scale(term, cv, NLAM, term);                                     /* c/λ⁴ */
        tmm_jet_scale(inverse_square, bv, NLAM, inverse_square);                 /* b/λ² */
        tmm_jet_add(inverse_square, term, NLAM, term);
        tmm_jet_constant(av, kv, NLAM, jets);
        tmm_jet_add(jets, term, NLAM, jets);
    }
    ```

    The phase calls take every index as a jet, so a constant index is a jet with only its value set:

    ```c
    double air[8 * NLAM] = { 0 }, glass[8 * NLAM], layer[8 * NLAM], omegas[NLAM], out[10 * NLAM];
    const double theta[1] = { 0 }, thick[1] = { 85 };
    for (int i = 0; i < NLAM; i++) air[8 * i] = 1;
    cauchy_jet(1.5, 4200, 0, 0, glass);
    cauchy_jet(2.21, 2.4e4, 0, 0, layer);
    tmm_omega_from_lambda_nm(grid, NLAM, omegas);
    tmm_phase_spectrum(grid, NLAM, theta, 1, 0, air, glass, layer, thick, 1, omegas, NULL, out);
    ```

Nothing there is differentiated by hand: the wavelength jet carries the derivatives of λ(ω), and every operation after it carries them on. Sellmeier and Lorentz formulas work the same way.

For tabulated data, interpolate with something at least three times differentiable and pass the value and its first three derivatives to `jetCompose`. A PCHIP interpolant is only $C^1$, and its GDD and TOD jump at every knot.

### An angle fixed outside the incident medium

Light that reaches a coating through its own substrate, at an angle fixed in air, meets the coating at an angle that changes with frequency: by Snell's law $\sin\theta_0 = \sin\theta_{\text{air}}/n_0(\omega)$ in the glass, $n_0$ the real part of its index. An `angle` set to the inside angle at each wavelength gets the phase right, but GD, GDD and TOD then leave out how that angle moves. Give the sine as a jet instead, with the Cauchy jet from above:

=== "JavaScript"

    ```js
    import { jetDivide } from 'tmmcore';

    // The two-layer coating seen from inside the glass, air beyond it, 45° in air.
    const lams = [700, 750, 800, 850, 900];
    const n0 = lams.map(l => cauchyJet(l, 1.5, 4000, 0, 0));
    const sine = n0.map(n => jetDivide(jetConstant(Math.sin(Math.PI / 4)), n));
    const back = {
        incident: { n: n0 },
        substrate: 1,
        layers: [
            { n: 1.46, d: 94.2 },
            { n: 2.35, d: 58.5 },
        ],
    };
    const reflected = phaseSpectrum(back, lams, { pol: 'p', sine }).r;
    ```

=== "Python"

    ```python
    # The two-layer coating seen from inside the glass, air beyond it, 45° in air.
    n0 = cauchy_jet(grid, 1.5, 4000, 0, 0)
    sine = tmmcore.jet_divide(tmmcore.jet_constant(np.sin(np.pi / 4)), n0)
    back = {
        "incident": {"n": n0},
        "substrate": 1.0,
        "layers": [
            {"n": 1.46, "d": 94.2},
            {"n": 2.35, "d": 58.5},
        ],
    }
    reflected = tmmcore.phase_spectrum(back, grid, pol="p", sine=sine)["r"]
    ```

=== "C"

    ```c
    /* The two-layer coating seen from inside the glass, air beyond it, 45° in air. */
    double n0[8 * NLAM], sine[8 * NLAM], air[8 * NLAM] = { 0 }, layers[2 * 8 * NLAM] = { 0 };
    double sin45[NLAM], zero[NLAM] = { 0 }, omegas[NLAM], reflected[10 * NLAM];
    const double theta[1] = { 0 }, index[2] = { 1.46, 2.35 }, thick[2] = { 94.2, 58.5 };
    cauchy_jet(1.5, 4000, 0, 0, n0);
    for (int i = 0; i < NLAM; i++) {
        sin45[i] = 0.7071067811865475;
        air[8 * i] = 1;
        for (int j = 0; j < 2; j++) layers[(j * NLAM + i) * 8] = index[j];
    }
    tmm_jet_constant(sin45, zero, NLAM, sine);
    tmm_jet_divide(sine, n0, NLAM, sine);
    tmm_omega_from_lambda_nm(grid, NLAM, omegas);
    tmm_phase_spectrum(grid, NLAM, theta, 1, 1, n0, air, layers, thick, 2, omegas, sine, reflected);
    ```

At 800 nm this gives a GD of 1.3265 fs and a GDD of −0.132 fs² on reflection. The inside angle held at its 800 nm value gives 1.3247 fs and −0.136 fs².

### Units and signs

The phase calls differentiate with respect to the angular frequency $\omega = 2\pi c/\lambda$ in rad/fs, with wavelengths in nm, so GD is in fs, GDD in fs² and TOD in fs³. In C the phase calls take ω as `omega` or `omegas`, and GD, GDD and TOD come out in its time unit: in fs with ω from `tmm_omega_from_lambda_nm`.

| | |
|---|---|
| `omegaFromLambdaNm(lambda_nm)` | angular frequency in rad/fs; `omega_from_lambda_nm` in Python, which also takes an array, and `tmm_omega_from_lambda_nm(lambda_nm, count, out)` in C |
| `C_NM_PER_FS` | the speed of light in vacuum, 299.792458 nm/fs; JavaScript only |

tmmcore uses $\tilde n = n + ik$ with $e^{-i\omega t}$, the complex conjugate of Macleod's convention. The phase is $\phi$ in Macleod's convention, the negative of the argument of $r$ or $t$ as tmmcore computes them, and GD, GDD and TOD are $-d\phi/d\omega$, $-d^2\phi/d\omega^2$ and $-d^3\phi/d\omega^3$, as he defines them. So all of them equal his. Positive GDD means the red end of a pulse arrives first.

$r_p$ has the sign of thin-film practice, $r_p = r_s$ at normal incidence. [Conventions](getting-started.md#conventions) explains the 180° this puts into an ellipsometric $\Delta$, and the [ellipsometry example](examples.md#ellipsometry-of-an-oxide-on-silicon) computes one.

## Jet arithmetic

A jet is four complex numbers: in JavaScript four `[re, im]` pairs, in Python a complex array of shape (4,), and in C 8 doubles. These calls are lower level than the ones above and may change between versions.

| | JavaScript | Python | C |
|---|---|---|---|
| a constant: the value, with zero derivatives | `jetConstant(value, imaginary?)` | `jet_constant(value, imaginary=0.0)` | `tmm_jet_constant(value, imaginary, count, out)` |
| a jet from plain derivatives, applying the factorials | `jetFromDerivatives(value, first?, second?, third?)` | `jet_from_derivatives(value, first=0.0, second=0.0, third=0.0)` | `tmm_jet_from_derivatives(derivatives, count, out)`, `derivatives` `count` × 8 |
| the plain derivatives `[f, f', f'', f''']` back | `jetDerivatives(jet)` | `jet_derivatives(jet)` | `tmm_jet_derivatives(jets, count, out)` |
| arithmetic | `jetAdd`, `jetSubtract`, `jetMultiply`, `jetDivide`, `jetReciprocal` | `jet_add`, `jet_subtract`, `jet_multiply`, `jet_divide`, `jet_reciprocal` | `tmm_jet_add`, `tmm_jet_subtract`, `tmm_jet_multiply`, `tmm_jet_divide`, `tmm_jet_reciprocal` |
| every order times a real number | `jetScale(jet, scalar)` | `jet_scale(jet, scalar)` | `tmm_jet_scale(jets, scalar, count, out)` |
| elementary functions | `jetSqrt`, `jetExp`, `jetLog`, `jetPower`, `jetSin`, `jetCos`, `jetSinCos` | `jet_sqrt`, `jet_exp`, `jet_log`, `jet_power`, `jet_sin`, `jet_cos`, `jet_sin_cos` | `tmm_jet_sqrt`, `tmm_jet_exp`, `tmm_jet_log`, `tmm_jet_power`, `tmm_jet_sin`, `tmm_jet_cos`, `tmm_jet_sin_cos` |
| the chain rule, for a function known through its derivatives | `jetCompose(value, [f', f'', f'''], xJet)` | `jet_compose(value, derivatives, x)` | `tmm_jet_compose(value, derivatives, x, count, out)`, `derivatives` `count` × 3 |
| λ as a jet in ω, the usual starting point | `wavelengthOmegaJet(wavelength, omega)` | `wavelength_omega_jet(wavelength, omega)` | `tmm_wavelength_omega_jet(wavelength, omega, count, out)` |
| the real part, order by order | `jetRealPart(jet)` | `jet_real_part(jet)` | `tmm_jet_real_part(jets, count, out)` |
| a real jet and an imaginary jet combined into one complex jet | `jetWithImaginaryPart(realJet, imaginaryJet)` | `jet_with_imaginary_part(real_jet, imaginary_jet)` | `tmm_jet_with_imaginary_part(realJets, imaginaryJets, count, out)` |
| guards | `jetClampImaginary`, `jetClampRealMinimum`, `jetIsFinite` | `jet_clamp_imaginary`, `jet_clamp_real_minimum`, `jet_is_finite` | `tmm_jet_clamp_imaginary`, `tmm_jet_clamp_real_minimum`, `tmm_jet_is_finite` |
| the highest order, 3 | `JET_ORDER` | | |

`jetSinCos` returns `{ sine, cosine }`, `jet_sin_cos` a dict of the same, and `tmm_jet_sin_cos` writes two arrays. A Python call takes arrays of jets and of real numbers and [broadcasts](python.md#jets) them. A C call applies its operation to `count` jets laid end to end, and a real argument that goes with each jet, such as the factor of `tmm_jet_scale`, is an array of `count`; an output may be the same array as an input of its shape. `tmm_jet_is_finite` writes an `int` per jet.

`wavelengthOmegaJet` needs no value of $c$: given λ and ω, every derivative follows from $\lambda' = -\lambda/\omega$, in the unit λ was given in.

## Monitoring and deposition

What an optical monitor reads while a coating is deposited. These calls take `angle` and no `pol`, and return s and p together: `Rs`, `Ts`, `Rp` and `Tp`, the reflectance and transmittance of the coated surface, and `Rrs` and `Rrp`, its reflectance seen from the substrate side. A caller that treats the substrate as a thick slab combines all three itself; for a semi-infinite substrate, $A = 1 - R - T$.

### `monitorCurve`

The signal of one layer growing on a finished stack, at one wavelength, for each thickness of it in `thicknesses`. `stack` is what is already deposited, and `top` the index of the growing layer.

=== "JavaScript"

    `monitorCurve(stack, wavelength, options)`

    ```js
    import { monitorCurve } from 'tmmcore';

    const finished = { incident: 1, substrate: 1.52, layers: [{ n: 1.46, d: 94.2 }] };
    const { Rs } = monitorCurve(finished, 550, { top: 2.35, thicknesses: [0, 20, 40, 58.5] });
    ```

=== "Python"

    `monitor_curve(stack, wavelength, *, top, thicknesses, angle=0.0)`

    ```python
    finished = {"incident": 1.0, "substrate": 1.52, "layers": [{"n": 1.46, "d": 94.2}]}
    Rs = tmmcore.monitor_curve(finished, 550, top=2.35, thicknesses=[0, 20, 40, 58.5])["Rs"]
    ```

=== "C"

    ```c
    void tmm_monitor_curve(double lambda_nm, double theta_deg,
                           double n0_re, double n0_im, double ns_re, double ns_im,
                           const double *base, int NB,
                           double ng_re, double ng_im,
                           const double *dArr, int nD,
                           double *outRs, double *outTs,
                           double *outRp, double *outTp,
                           double *outRrs, double *outRrp);
    ```

    `base` holds the `NB` layers of the finished stack as `layers` does in `tmm_rta`.

The finished stack's matrix is built once, and each thickness then costs one 2×2 product, so a whole monitoring curve is one call: the incremental control algorithm of Tikhonravov and Trubetskov, *Appl. Opt.* **44**, 6877 (2005).

**Returns** `Rs`, `Ts`, `Rp`, `Tp`, `Rrs` and `Rrp`, each with one value per thickness.

### `depositionSpectra`

The spectrum after each layer of a deposition run, in one call. The layers are deposited from the substrate outward, the last one listed first. Each is added to the running product and the grid evaluated after it, so all N spectra together cost about what the last one costs alone. `n`, `k` and `angle` may be one per wavelength, as in any grid call.

=== "JavaScript"

    `depositionSpectra(stack, wavelengths, options?)`

    Each output is flat, wavelength `i` after `k + 1` layers at `[k * nLam + i]`.

=== "Python"

    `deposition_spectra(stack, wavelengths, *, angle=0.0)`

    Each output has the shape (N, nLam), wavelength `i` after `k + 1` layers at `[k, i]`. It raises `MemoryError` if the kernel cannot allocate its working memory.

=== "C"

    ```c
    int tmm_deposition_spectra(const double *lambdas, int nLam,
                               const double *theta_deg, int nTheta,
                               const double *n0arr, const double *nsarr,
                               const double *matNK, const double *thick, int N,
                               double *outRs, double *outTs,
                               double *outRp, double *outTp,
                               double *outRrs, double *outRrp);
    ```

    It takes the layers in the order they are deposited, the one against the substrate first: `matNK` and `thick` hold the layers of the stack in reverse. Each output is N × nLam, `[step][wavelength]`. It returns 1, or 0 with the outputs untouched when it could not allocate its working memory.

**Returns** `Rs`, `Ts`, `Rp`, `Tp`, `Rrs` and `Rrp` with N × nLam values each: after `k + 1` layers, the last `k + 1` of the stack's layers. A layer of zero thickness repeats the spectrum before it.

### `growingStack`

A layer growing on a finished stack over a whole wavelength grid, one thickness per call: what a broadband monitor reads at each scan. The finished stack is multiplied out once for every wavelength and polarization and kept, so each sample costs one layer matrix and one 2×2 product per wavelength and polarization. `n`, `k` and `angle` may be one per wavelength, as in any grid call.

It returns a handle with three calls. One sets the growing layer's index, given as a medium is, `n` and `k` each a number or one per wavelength; it comes before the first sample with `d > 0`, and again to change the material. One samples the grid with the growing layer at thickness `d` nm, giving `Rs`, `Ts`, `Rp`, `Tp`, `Rrs` and `Rrp` with one value per wavelength, and `d <= 0` gives the finished stack. One releases what the handle holds.

=== "JavaScript"

    `growingStack(stack, wavelengths, options?)`

    | | |
    |---|---|
    | `setTop(top)` | sets the growing layer's index |
    | `sample(d, out?)` | samples the grid. Pass the previous result as `out` to reuse its arrays. Throws if called with `d > 0` before `setTop` |
    | `free()` | releases what the handle holds. On the WebAssembly instance that is kernel memory: a dropped handle is reclaimed by a finalizer in time, but calling `free()` keeps a long run's memory flat. `setTop` and `sample` throw after it |

=== "Python"

    `growing_stack(stack, wavelengths, *, angle=0.0)` returns a `GrowingStack`, with `set_top(top)`, `sample(d)` and `close()`. It is a context manager that closes it, and `sample` returns new arrays each time:

    ```python
    with tmmcore.growing_stack(finished, np.linspace(400, 800, 401), angle=45) as growing:
        growing.set_top(2.35)
        scans = [growing.sample(d)["Rs"] for d in np.arange(0, 60, 5)]
    ```

    `sample` with `d > 0` before `set_top`, or any call after `close()`, raises `ValueError`, and `growing_stack` raises `MemoryError` if the kernel cannot allocate the handle's memory.

=== "C"

    ```c
    tmm_growing_stack *tmm_growing_stack_create(const double *lambdas, int nLam,
                                                const double *theta_deg, int nTheta,
                                                const double *n0arr, const double *nsarr,
                                                const double *matNK, const double *thick,
                                                int NB);
    void tmm_growing_stack_set_top(tmm_growing_stack *h, const double *ngNK);
    int tmm_growing_stack_sample(tmm_growing_stack *h, double d,
                                 double *outRs, double *outTs,
                                 double *outRp, double *outTp,
                                 double *outRrs, double *outRrp);
    void tmm_growing_stack_free(tmm_growing_stack *h);
    ```

    `tmm_growing_stack_create` takes the finished stack as `tmm_rta_spectrum` takes a stack, with `NB` the number of its layers. It copies what it keeps, and returns `NULL` if it cannot allocate. `ngNK` is the growing layer's index at each wavelength, `nLam` × 2. `tmm_growing_stack_sample` writes `nLam` values to each output and returns 1, or 0 with nothing written when `d > 0` comes before `tmm_growing_stack_set_top`.

## C

The kernel in `src/c` is C99 with no dependency beyond libm; the WebAssembly module and the Python package are built from it. `src/c/tmmcore.h` declares every call on this page.

### Building

Compile `tmm_kernel.c`, which includes the other sources, into a program, or into a shared library:

```bash
cc -std=c99 -O2 -c src/c/tmm_kernel.c
cc -std=c99 -O2 -shared -fPIC -o libtmmcore.so src/c/tmm_kernel.c -lm
```

On Windows the calls are marked `__declspec(dllexport)`, so a DLL built from `tmm_kernel.c` exports them.

### Arrays

- Every value is a `double`, except the counts, `pol` and the results of `tmm_jet_is_finite`, which are `int`.
- The caller allocates every array, the outputs too, and a call keeps none of them; `tmm_growing_stack_create` copies what it needs.
- An array of more than one dimension is flat and row-major: `matNK` is `[layer][wavelength][n, k]`, so layer `j` at wavelength `i` starts at `(j * nLam + i) * 2`.
- A jet is 8 doubles, `[re, im]` of each of $f$, $f'$, $f''/2!$ and $f'''/3!$, and a list of them is laid end to end.

### A complete program

Rs and Rp of the two-layer coating at 45°:

```c
#include <stdio.h>
#include "tmmcore.h"

#define NLAM 5
#define NLAYERS 2

int main(void)
{
    const double lambdas[NLAM] = { 450, 500, 550, 600, 650 };
    const double theta[1] = { 45 };                 /* one angle for the whole grid */
    const double index[NLAYERS] = { 2.35, 1.46 };
    const double thick[NLAYERS] = { 58.5, 94.2 };
    double n0[NLAM * 2], ns[NLAM * 2], matNK[NLAYERS * NLAM * 2];
    double Rs[NLAM], Ts[NLAM], As[NLAM], Rp[NLAM], Tp[NLAM], Ap[NLAM];

    for (int i = 0; i < NLAM; i++) {
        n0[2 * i] = 1.0;   n0[2 * i + 1] = 0.0;
        ns[2 * i] = 1.52;  ns[2 * i + 1] = 0.0;
        for (int j = 0; j < NLAYERS; j++) {
            matNK[(j * NLAM + i) * 2] = index[j];
            matNK[(j * NLAM + i) * 2 + 1] = 0.0;
        }
    }
    tmm_rta_spectrum(lambdas, NLAM, theta, 1, n0, ns, matNK, thick, NLAYERS,
                     Rs, Ts, As, Rp, Tp, Ap);
    for (int i = 0; i < NLAM; i++)
        printf("%.0f nm   Rs %.4f   Rp %.4f\n", lambdas[i], Rs[i], Rp[i]);
    return 0;
}
```

```bash
cc -std=c99 -O2 -I src/c -o example example.c src/c/tmm_kernel.c -lm
```

## WebAssembly

The C kernel compiled to WebAssembly, with the same calls as the JavaScript. [Getting started](getting-started.md#turning-on-webassembly) shows how to load it in Node and in a browser.

```js
import { instantiateTmmWasm } from 'tmmcore';

const wasm = await instantiateTmmWasm(bytes);   // bytes: the contents of tmm_kernel.wasm
const { R } = wasm.rtaSpectrum(stack, wavelengths);
```

The grid calls are where the speed is: they cross from JavaScript into WebAssembly per grid, not per wavelength. The loader is also its own entry point, `tmmcore/wasm`, for code that needs only the WebAssembly side.

An instance is not shared between threads, since there is no shared memory. Each worker instantiates its own from the same bytes.

### The instance

| Method | |
|---|---|
| `rta`, `rtaSpectrum`, `jacobian`, `jacobianSpectrum`, `hessian`, `hessianSpectrum`, `needleScan`, `needleScanSpectrum`, `phase`, `phaseSpectrum`, `phaseJacobian`, `phaseJacobianSpectrum`, `monitorCurve`, `depositionSpectra`, `growingStack` | the calls above, with the same arguments and results |
| `jetConstant` … `jetIsFinite`, `wavelengthOmegaJet`, `omegaFromLambdaNm` | the jet arithmetic |
| `hasHessian()`, `hasPhase()`, `hasPhaseJacobianSpectrum()`, `hasGrowingKernels()`, `hasGrowingEval()` | always `true`: a module missing a kernel is refused when it loads |

If the WebAssembly heap is exhausted, a call that cannot allocate memory for its arguments and results throws `Error('tmmWasm: malloc failed')`; the one-wavelength calls reuse theirs from call to call and need more only for a larger stack. A kernel that cannot allocate its own working state writes `NaN` to every output instead: numbers and arrays come back `NaN`, and `phase` and `phaseJacobian` return `null` for both coefficients. `depositionSpectra` and `growingStack` throw.

### Loading

| | |
|---|---|
| `instantiateTmmWasm(bytes)` | instantiates from an `ArrayBuffer` or `Uint8Array`, and resolves to the instance |
| `initTmmWasmFromUrl(url)` | fetches and instantiates. Resolves `true` or `false` and does not throw; a second call returns the first call's promise |
| `getTmmWasm()` | the last instance loaded, or `null`. A failed `initTmmWasmFromUrl` also leaves it `null` |
| `isTmmWasmReady()` | whether an instance is loaded |
| `setTmmWasmEnabled(on)`, `isTmmWasmEnabled()` | a flag for whether to use the instance, off at first |
| `tmmWasmActive()` | whether the flag is on and an instance is loaded |
| `TmmWasmInstance` | the class of the instance |

For workers that cannot fetch the file themselves, the main thread hands them its bytes:

| | |
|---|---|
| `initTmmWasmMainThread(bytes, enabled)` | main thread: keeps the bytes for workers, then instantiates and turns the flag on if `enabled`, or turns it off if not. Resolves `true` when the flag is on |
| `getTmmWasmBytesForWorker()` | main thread: the bytes to send a worker, or `null` while the flag is off |
| `noteTmmWasmBytes(bytes)` | worker: starts instantiating from the bytes it was sent and turns its flag on; does nothing once an instance is loaded |
| `awaitTmmWasmReady()` | worker: resolves when that instantiation ends, `true` if an instance is loaded |

## Edge cases

### Grazing incidence and critical angles

The angle of incidence runs up to and including 90°. $\cos\theta_0$ is taken directly rather than as $\sqrt{1 - \sin^2\theta_0}$, which keeps nothing of it near 90°, and a layer of the incident index keeps that cosine too. At 90° exactly $R$ is one and $T$ zero to rounding.

At a medium's critical angle to the last bit, $\cos\theta$ in it is exactly zero and the layer matrix as written is $0/0$. Its limit is finite and tmmcore uses it: the layer's matrix becomes linear in its thickness, and a substrate at its critical angle reflects everything. The derivative calls take the same limit. The phase calls do not; see [`phase`](#phase-gd-gdd-and-tod).

### Opaque layers

The imaginary part of a layer's phase thickness grows with its thickness, and $\cosh$ of it overflows, which would turn the whole matrix product into `NaN`. So the imaginary phase is held at 50 when the matrix is built. Past that point the held matrix is the true one divided by a real factor, $e^{|\mathrm{Im}\,\delta| - 50}$, to within $e^{-100}$ relative. The factor cancels from $r$ and from every phase, and tmmcore carries it as a log scale for $t$, so $T$ is right however opaque the layer, until it falls below about $10^{-308}$, the range of a double, and reads 0.

The hold depends on the phase, not on $k$. A thick absorbing layer reaches it, and so does a lossless one beyond the critical angle, where the wave in it is evanescent: an air gap of 10 µm under a glass prism at 70° is well past it.

The derivatives include the factor, so $T$ falls with an opaque layer's thickness at the true rate, and $R$ does not move. In a deep stack the partial products outgrow the range of a double; each carries a power-of-two scale that cancels exactly from every derivative, so no wavelength drops to zero. A needle inside an opaque layer sees the field that is there: near the surface its derivatives are those of an ordinary layer, and deep inside $dR$ is zero while $dT$ keeps its value relative to $T$.

## The names of 0.4

0.5.x keeps the JavaScript names of 0.4 with their arguments and results. Each prints one warning on its first call, and 0.6.0 removes them.

| 0.4, module | 0.4, WebAssembly instance | 0.5 |
|---|---|---|
| `tmm` | `tmmOne` | `rta` |
| | `tmmSpectrum` | `rtaSpectrum` |
| `tmmThicknessJacobian` | `tmmJacobian` | `jacobian` |
| `tmmThicknessHessian` | `tmmHessian` | `hessian` |
| `tmmNeedleScan` | `tmmNeedleScan` | `needleScan` |
| `tmmPhaseDispersion` | `tmmPhaseOne` | `phase` |
| | `tmmPhaseSpectrum` | `phaseSpectrum` |
| `tmmPhaseThicknessJacobian` | `tmmPhaseJacobian` | `phaseJacobian` |
| | `tmmPhaseJacobianSpectrum` | `phaseJacobianSpectrum` |
| | `monitorCurve`, `depositionSpectra`, with the arguments of 0.4 | `monitorCurve`, `depositionSpectra` |
| | `growingEval` | `growingStack` |
| `tmmCoefficientJets`, `tmmCoefficientThicknessJets`, `coefficientPhaseDispersion`, `coefficientPhaseThicknessDerivatives` | | none |

The one-wavelength names take `(lambda_nm, theta_deg, pol, n0, ns, layers)`, each index as `[re, im]` and the layers as `{ n, d }[]`, or in the phase names `{ nJet, d }[]` with the media as jets; the instance takes `pol` as `0` or `1`. `tmmNeedleScan` takes after them the candidates, each as `[re, im]`, and the fractions. The [0.4.2 API reference](https://github.com/aai2k/tmmcore/blob/0a432b7/docs/api.md) has them in full.

The instance's `monitorCurve` and `depositionSpectra` tell the two forms apart by the first argument: a stack in 0.5, a wavelength or a list of them in 0.4. `depositionSpectra` and `growingEval` in 0.4 take the wavelengths, the media and the layers as lists of `[re, im]`; `depositionSpectra` takes the layers in the order they are deposited, the first deposited first, where the 0.5 call takes the stack's order.

The 0.4 phase names take the incident sine as the option `sinTheta0Jet`, and over a grid `sinJets`, where the 0.5 calls take `sine`. The 0.5 calls take no `omega`, which sets the time unit; the 0.4 phase names still do.

## Low-level primitives

JavaScript only. Exported for building transfer-matrix variants on the same conventions. These are lower level than the calls above and may change between versions. Complex numbers are `[re, im]` pairs.

| | |
|---|---|
| `cadd`, `csub`, `cmul`, `cdiv` | complex arithmetic |
| `cabs2`, `cconj`, `csqrt`, `ccos`, `csin` | modulus squared, conjugate, square root, cosine and sine |
| `creal`, `cimag` | the parts |
| `matmul(A, B)` | 2×2 complex matrix product |
| `cmatvec(M, v)` | 2×2 matrix times a 2-vector |
| `snellCosTheta(n0, sinTheta0, nj, cosTheta0?)` | complex $\cos\theta$ in a medium, from the real invariant $\mathrm{Re}(n_0)\sin\theta_0$. Pass `cosTheta0` as `[cos θ0, 0]`, as `rta` does, and a medium with an index close to the incident one keeps its cosine at grazing incidence |
| `incidentCosTheta(n0, sinTheta0, cosTheta0?)` | $\cos\theta_0$ of the incident medium itself: its own cosine when it is transparent, from the same invariant when it absorbs |
| `layerMatrix(nj, dj_nm, lambda_nm, cosTheta_j, pol)` | the characteristic matrix of one layer; its limit where `cosTheta_j` is exactly zero, and past the [opaque-layer hold](#opaque-layers) the held matrix |
| `layerLogScale(nj, dj_nm, lambda_nm, cosTheta_j)` | the log of the factor `layerMatrix` leaves out past the hold, 0 short of it |
| `rescaleMatrix(M)` | rescales in place past an overflow threshold and returns the log scale it took out |

`layerLogScale` and `rescaleMatrix` keep opaque stacks finite and their transmittance right. Both are real factors that cancel from reflectance but not from transmittance, so a caller building its own product adds both to one log scale and multiplies T by $e^{-2\,\text{logScale}}$, as `rta` does. A field computed from a partial product carries its own log scale the same way.
