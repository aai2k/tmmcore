/**
 * phase.js : analytic group delay, GDD and TOD of a multilayer.
 *
 * The characteristic matrix is evaluated in third-order Taylor arithmetic with
 * angular frequency as the differentiation variable, so the reflection and
 * transmission coefficients come out with their first three frequency
 * derivatives already attached. The reported quantities are then imaginary parts
 * of logarithmic derivatives of the coefficient:
 *
 *     GD  = Im(r'/r)
 *     GDD = Im(r''/r − (r'/r)²)
 *     TOD = Im(r'''/r − 3 r'r''/r² + 2 (r'/r)³)
 *
 * No phase unwrapping and no wavelength finite-difference stencil takes part, so
 * a value at one wavelength is independent of every neighbouring sample. Sample
 * wherever you like; the numbers do not move.
 *
 * Conventions are those of tmm.js: ñ = n + ik, exp(−iωt), off-diagonals carrying
 * −i. This is the complex conjugate of Macleod's convention, under which each
 * reported order is minus the derivative of physical phase, so all three equal
 * derivatives of the raw transfer-matrix phase computed here.
 *
 * Units. The differentiation variable is angular frequency, and the caller
 * chooses its unit by choosing `omega`. GD comes back in that unit's reciprocal,
 * GDD in its square, TOD in its cube. The default `omega` uses `C_NM_PER_FS`
 * with wavelengths in nm, giving fs, fs² and fs³.
 *
 * Refractive indices arrive as jets, `[ñ, dñ/dω, (d²ñ/dω²)/2, (d³ñ/dω³)/6]`.
 * The package deliberately owns no material models; build the jets with the
 * arithmetic in taylorJet.js, from whatever dispersion formula or interpolant
 * you use. See the docs for a worked Sellmeier example.
 *
 * The C kernel shipped alongside is a port of this file and agrees with it to
 * float64 round-off (see tests/).
 *
 * Reference:
 *   • Birge & Kärtner, "Analysis of the effects of dispersion on the phase of
 *     ultrashort pulses", Appl. Opt. 45, 1478-1483 (2006).
 *     https://doi.org/10.1364/AO.45.001478
 */

import {
    jetAdd,
    jetClampImaginary,
    jetConstant,
    jetDerivatives,
    jetDivide,
    jetMultiply,
    jetRealPart,
    jetScale,
    jetSinCos,
    jetSqrt,
    jetSubtract,
    wavelengthOmegaJet,
} from './taylorJet.js';

/** Speed of light in vacuum, nm/fs. Sets the default time unit to fs. */
export const C_NM_PER_FS = 299.792458;

const MATRIX_RESCALE_THRESHOLD = 1e100;
const BINARY_RESCALE_LIMIT = 2 ** 64;
const MAX_IMAGINARY_PHASE = 50;

/** Angular frequency in rad/fs for a vacuum wavelength in nm. */
export function omegaFromLambdaNm(lambda_nm) {
    return 2 * Math.PI * C_NM_PER_FS / lambda_nm;
}

// ── Jet-valued 2×2 matrices ──────────────────────────────────────────────────

function identityMatrix() {
    return [
        [jetConstant(1), jetConstant(0)],
        [jetConstant(0), jetConstant(1)],
    ];
}

function zeroMatrix() {
    return [
        [jetConstant(0), jetConstant(0)],
        [jetConstant(0), jetConstant(0)],
    ];
}

function matrixMultiply(left, right) {
    return [
        [
            jetAdd(jetMultiply(left[0][0], right[0][0]), jetMultiply(left[0][1], right[1][0])),
            jetAdd(jetMultiply(left[0][0], right[0][1]), jetMultiply(left[0][1], right[1][1])),
        ],
        [
            jetAdd(jetMultiply(left[1][0], right[0][0]), jetMultiply(left[1][1], right[1][0])),
            jetAdd(jetMultiply(left[1][0], right[0][1]), jetMultiply(left[1][1], right[1][1])),
        ],
    ];
}

function matrixMagnitude(matrix) {
    let magnitude = 0;
    for (const row of matrix) {
        for (const jet of row) {
            for (const coefficient of jet) {
                magnitude = Math.max(magnitude, Math.abs(coefficient[0]), Math.abs(coefficient[1]));
            }
        }
    }
    return magnitude;
}

// Divides every order of a jet matrix by a power of two, which is exact, when
// its largest part at any order exceeds BINARY_RESCALE_LIMIT, and returns that
// power's exponent. The prefix and suffix products of the thickness Jacobian
// carry these exponents, as the derivative kernels in tmm.js do.
function rescaleMatrixBinary(matrix) {
    const magnitude = matrixMagnitude(matrix);
    if (!(magnitude > BINARY_RESCALE_LIMIT) || magnitude === Infinity) return 0;
    const exponent = Math.floor(Math.log2(magnitude));
    scaleMatrixInPlace(matrix, 2 ** -exponent);
    return exponent;
}

function scaleMatrixInPlace(matrix, factor) {
    for (const row of matrix) {
        for (let column = 0; column < row.length; column++) {
            row[column] = jetScale(row[column], factor);
        }
    }
}

// The order-0 matrix controls overflow in the physical coefficient. Once
// selected, one plain scalar rescales every jet order and cancels from r.
function rescaleMatrix(matrix, threshold) {
    let scale = 0;
    for (const row of matrix) {
        for (const value of row) {
            scale = Math.max(scale, Math.abs(value[0][0]), Math.abs(value[0][1]));
        }
    }
    if (scale <= threshold) return 0;
    scaleMatrixInPlace(matrix, 1 / scale);
    return Math.log(scale);
}

// The transverse invariant is Re(n0) sinθ0, as in tmm.js: an absorbing incident
// medium carries a wave whose amplitude falls along the normal only. cos²θ is
// 1 − (a/n)² or ((n − n0)(n + n0) + (n0 cosθ0)²)/n², chosen on the values by
// the rule in snellCosTheta there, so that a layer of the incident index keeps
// its cosine at grazing incidence. `incidentCosineJet` is null when the
// incident sine arrives as a jet, and then the first form is used.
function snellCosine(incidentIndex, incidentSine, layerIndex, incidentCosineJet) {
    const realIndex = jetRealPart(incidentIndex);
    if (incidentCosineJet) {
        const nr = incidentIndex[0][0];
        const nc = nr * incidentCosineJet[0][0];
        const q = nc * nc;
        const room = nr * nr - 2 * q;
        const [lr, li] = layerIndex[0];
        const near = room > 0 && Math.abs(lr * lr + li * li - nr * nr) < room;
        const m = near
            ? jetMultiply(jetSubtract(layerIndex, realIndex), jetAdd(layerIndex, realIndex))
            : null;
        if (m && m[0][0] * m[0][0] + m[0][1] * m[0][1] < room * room) {
            const normal = jetMultiply(realIndex, incidentCosineJet);
            return jetSqrt(jetDivide(jetAdd(m, jetMultiply(normal, normal)),
                jetMultiply(layerIndex, layerIndex)));
        }
    }
    const layerSine = jetDivide(jetMultiply(realIndex, incidentSine), layerIndex);
    return jetSqrt(jetSubtract(jetConstant(1), jetMultiply(layerSine, layerSine)));
}

// The incident sine as a jet, and the incident cosine: the constant cosθ0 at a
// fixed angle of incidence, null when the caller gave the sine as a jet.
function incidentAngle(thetaDeg, incidentSineJet) {
    return incidentSineJet
        ? { sine: incidentSineJet, cosine: null }
        : {
            sine: jetConstant(Math.sin(thetaDeg * Math.PI / 180)),
            cosine: jetConstant(Math.cos(thetaDeg * Math.PI / 180)),
        };
}

// cosθ0 of the incident medium: its own cosine while the index is real at
// every order, otherwise from the same real invariant as the layers, so the
// incident admittance is sqrt(N0² − n0² sin²θ0) and matches tmm.js.
function incidentCosine(incidentIndex, { sine, cosine }) {
    if (incidentIndex.some(coefficient => coefficient[1] !== 0)) {
        return snellCosine(incidentIndex, sine, incidentIndex, cosine);
    }
    return cosine ?? jetSqrt(jetSubtract(jetConstant(1), jetMultiply(sine, sine)));
}

function admittance(index, cosine, polarization) {
    return polarization === 's'
        ? jetMultiply(index, cosine)
        : jetDivide(index, cosine);
}

// |Im δ| − MAX_IMAGINARY_PHASE where the bound held the phase, else 0: the log
// of the real factor the held matrix falls short of the true one by, at the
// wavelength itself. It cancels from r and from the phase of t, and only |t|
// needs it; see layerLogScale in tmm.js.
function heldExcess(rawPhase, phase) {
    return rawPhase[0][1] === phase[0][1] ? 0 : Math.abs(rawPhase[0][1]) - MAX_IMAGINARY_PHASE;
}

function layerMatrix(index, thickness, wavelength, cosine, polarization) {
    const rawPhase = jetScale(
        jetDivide(jetMultiply(index, cosine), wavelength),
        2 * Math.PI * thickness,
    );
    const phase = jetClampImaginary(rawPhase, MAX_IMAGINARY_PHASE);
    const { sine, cosine: cosinePhase } = jetSinCos(phase);
    const eta = admittance(index, cosine, polarization);
    const minusI = jetConstant(0, -1);
    return {
        matrix: [
            [cosinePhase, jetMultiply(minusI, jetDivide(sine, eta))],
            [jetMultiply(minusI, jetMultiply(eta, sine)), cosinePhase],
        ],
        excess: heldExcess(rawPhase, phase),
    };
}

function layerMatrixWithThicknessDerivative(index, thickness, wavelength, cosine, polarization) {
    const phasePerUnit = jetScale(
        jetDivide(jetMultiply(index, cosine), wavelength),
        2 * Math.PI,
    );
    const rawPhase = jetScale(phasePerUnit, thickness);
    const phase = jetClampImaginary(rawPhase, MAX_IMAGINARY_PHASE);
    // Past the bound the layer is the held matrix times a real factor, and
    // the full phase derivative, taken with the held sines and cosines, is
    // their derivative to within e^{−100} relative, as derivativeLayer in
    // tmm.js has it. Its imaginary part is what |t| falls by; it cancels from
    // r and from every phase.
    const phaseDerivative = phasePerUnit;
    const { sine, cosine: cosinePhase } = jetSinCos(phase);
    const sineDerivative = jetMultiply(cosinePhase, phaseDerivative);
    const cosineDerivative = jetScale(jetMultiply(sine, phaseDerivative), -1);
    const eta = admittance(index, cosine, polarization);
    const minusI = jetConstant(0, -1);
    return {
        matrix: [
            [cosinePhase, jetMultiply(minusI, jetDivide(sine, eta))],
            [jetMultiply(minusI, jetMultiply(eta, sine)), cosinePhase],
        ],
        thicknessDerivative: [
            [cosineDerivative, jetMultiply(minusI, jetDivide(sineDerivative, eta))],
            [jetMultiply(minusI, jetMultiply(eta, sineDerivative)), cosineDerivative],
        ],
        excess: heldExcess(rawPhase, phase),
    };
}

// `transmissionScale` restores the factor a rescaled matrix gave up: it cancels
// from r but not from t.
function coefficientJetsFromMatrix(matrix, incidentEta, substrateEta, transmissionScale = 1) {
    const boundaryB = jetAdd(matrix[0][0], jetMultiply(matrix[0][1], substrateEta));
    const boundaryC = jetAdd(matrix[1][0], jetMultiply(matrix[1][1], substrateEta));
    const incidentB = jetMultiply(incidentEta, boundaryB);
    const denominator = jetAdd(incidentB, boundaryC);
    const reflectionNumerator = jetSubtract(incidentB, boundaryC);
    const reflection = jetDivide(reflectionNumerator, denominator);
    let transmission = jetDivide(jetScale(incidentEta, 2), denominator);
    if (transmissionScale !== 1) transmission = jetScale(transmission, transmissionScale);
    return { reflection, transmission, denominator, reflectionNumerator };
}

function coefficientThicknessJets(matrixDerivative, coefficientData, incidentEta, substrateEta) {
    const boundaryBDerivative = jetAdd(
        matrixDerivative[0][0],
        jetMultiply(matrixDerivative[0][1], substrateEta),
    );
    const boundaryCDerivative = jetAdd(
        matrixDerivative[1][0],
        jetMultiply(matrixDerivative[1][1], substrateEta),
    );
    const incidentBDerivative = jetMultiply(incidentEta, boundaryBDerivative);
    const denominatorDerivative = jetAdd(incidentBDerivative, boundaryCDerivative);
    const numeratorDerivative = jetSubtract(incidentBDerivative, boundaryCDerivative);
    return {
        reflection: jetDivide(
            jetSubtract(
                numeratorDerivative,
                jetMultiply(coefficientData.reflection, denominatorDerivative),
            ),
            coefficientData.denominator,
        ),
        transmission: jetScale(
            jetDivide(
                jetMultiply(coefficientData.transmission, denominatorDerivative),
                coefficientData.denominator,
            ),
            -1,
        ),
    };
}

// ── Coefficient jets ─────────────────────────────────────────────────────────

/**
 * Complex Taylor jets for r and t at one wavelength.
 *
 * Lower level than `tmmPhaseDispersion`: use this when you want the coefficient
 * jets themselves rather than the phase quantities read off them.
 *
 * A layer past the opaque-layer bound enters the product short by a real
 * factor (see layerLogScale in tmm.js), which makes `transmission` too large
 * by it. The factor is returned as `transmissionLogScale` rather than applied,
 * so that t keeps its phase and every derivative of it however far below the
 * double range |t| falls: the true |t| at this wavelength is |transmission|
 * times e^{−transmissionLogScale}.
 *
 * @param {object} options
 * @param {number[][]} options.wavelengthJet   jet of λ(ω), from `wavelengthOmegaJet`
 * @param {number} options.thetaDeg            angle of incidence, degrees
 * @param {'s'|'p'} options.polarization
 * @param {number[][]} options.incidentIndexJet
 * @param {number[][]} options.substrateIndexJet
 * @param {number[][]} [options.incidentSineJet]  see `tmmPhaseDispersion`
 * @param {{indexJet:number[][], thicknessNm:number}[]} options.layers
 * @returns {{reflection, transmission, transmissionLogScale, incidentEta, substrateEta}}
 */
export function tmmCoefficientJets({
    wavelengthJet,
    thetaDeg,
    polarization,
    incidentIndexJet,
    substrateIndexJet,
    incidentSineJet = null,
    rescaleThreshold = MATRIX_RESCALE_THRESHOLD,
    layers,
}) {
    const angle = incidentAngle(thetaDeg, incidentSineJet);
    const incidentEta = admittance(incidentIndexJet, incidentCosine(incidentIndexJet, angle), polarization);
    const substrateCosine = snellCosine(incidentIndexJet, angle.sine, substrateIndexJet, angle.cosine);
    const substrateEta = admittance(substrateIndexJet, substrateCosine, polarization);

    let matrix = identityMatrix();
    let logScale = 0;
    let transmissionLogScale = 0;
    for (const layer of layers) {
        if (!(layer.thicknessNm > 0)) continue;
        const cosine = snellCosine(incidentIndexJet, angle.sine, layer.indexJet, angle.cosine);
        const held = layerMatrix(
            layer.indexJet,
            layer.thicknessNm,
            wavelengthJet,
            cosine,
            polarization,
        );
        matrix = matrixMultiply(matrix, held.matrix);
        logScale += rescaleMatrix(matrix, rescaleThreshold);
        transmissionLogScale += held.excess;
    }

    const coefficients = coefficientJetsFromMatrix(matrix, incidentEta, substrateEta,
        Math.exp(-logScale));
    return {
        reflection: coefficients.reflection,
        transmission: coefficients.transmission,
        transmissionLogScale,
        incidentEta,
        substrateEta,
    };
}

/**
 * Coefficient jets together with their exact derivatives with respect to every
 * layer thickness. Frequency remains the Taylor variable, so each thickness
 * derivative is itself a third-order frequency jet.
 *
 * Every prefix and suffix product is stored as a mantissa and a binary
 * exponent, so opaque stacks whose products leave double range still return
 * their derivatives. Each derivative is one prefix times one suffix over the
 * full product, and the exponents cancel from it exactly. `transmission` and
 * `transmissionLogScale` are as in `tmmCoefficientJets`; the thickness jets
 * of t are those of `transmission`, so their ratio to it is the true one.
 */
export function tmmCoefficientThicknessJets(options) {
    const {
        wavelengthJet,
        thetaDeg,
        polarization,
        incidentIndexJet,
        substrateIndexJet,
        incidentSineJet = null,
        layers,
    } = options;
    const angle = incidentAngle(thetaDeg, incidentSineJet);
    const incidentEta = admittance(incidentIndexJet, incidentCosine(incidentIndexJet, angle), polarization);
    const substrateCosine = snellCosine(incidentIndexJet, angle.sine, substrateIndexJet, angle.cosine);
    const substrateEta = admittance(substrateIndexJet, substrateCosine, polarization);
    const layerData = layers.map((layer) => {
        // The point evaluator skips values that are not positive. Keep a slot in
        // the Jacobian for the same invalid layer, but give it no optical effect
        // and no derivative so both public paths report the same base result.
        // Zero is different: its nonzero derivative is useful for candidate
        // layers and is why the Jacobian deliberately retains it.
        if (!(layer.thicknessNm >= 0)) {
            return { matrix: identityMatrix(), thicknessDerivative: zeroMatrix(), excess: 0 };
        }
        const cosine = snellCosine(incidentIndexJet, angle.sine, layer.indexJet, angle.cosine);
        return layerMatrixWithThicknessDerivative(
            layer.indexJet,
            layer.thicknessNm,
            wavelengthJet,
            cosine,
            polarization,
        );
    });

    const count = layerData.length;
    const prefix = new Array(count + 1);
    const prefixExponent = new Array(count + 1).fill(0);
    prefix[0] = identityMatrix();
    for (let index = 0; index < count; index++) {
        prefix[index + 1] = matrixMultiply(prefix[index], layerData[index].matrix);
        prefixExponent[index + 1] = prefixExponent[index] + rescaleMatrixBinary(prefix[index + 1]);
    }
    const suffix = new Array(count + 1);
    const suffixExponent = new Array(count + 1).fill(0);
    suffix[count] = identityMatrix();
    for (let index = count - 1; index >= 0; index--) {
        suffix[index] = matrixMultiply(layerData[index].matrix, suffix[index + 1]);
        suffixExponent[index] = suffixExponent[index + 1] + rescaleMatrixBinary(suffix[index]);
    }

    const totalExponent = prefixExponent[count];
    const coefficientData = coefficientJetsFromMatrix(
        prefix[count],
        incidentEta,
        substrateEta,
        2 ** -totalExponent,
    );
    const thicknessDerivatives = layerData.map((layer, index) => {
        const matrixDerivative = matrixMultiply(
            matrixMultiply(prefix[index], layer.thicknessDerivative),
            suffix[index + 1],
        );
        scaleMatrixInPlace(matrixDerivative,
            2 ** (prefixExponent[index] + suffixExponent[index + 1] - totalExponent));
        return coefficientThicknessJets(
            matrixDerivative,
            coefficientData,
            incidentEta,
            substrateEta,
        );
    });
    return {
        reflection: coefficientData.reflection,
        transmission: coefficientData.transmission,
        transmissionLogScale: layerData.reduce((sum, layer) => sum + layer.excess, 0),
        reflectionThickness: thicknessDerivatives.map(value => value.reflection),
        transmissionThickness: thicknessDerivatives.map(value => value.transmission),
        incidentEta,
        substrateEta,
        thicknessDerivatives,
    };
}

// ── Phase quantities ─────────────────────────────────────────────────────────

/**
 * Read phase, GD, GDD and TOD off a coefficient jet.
 *
 * `logScale` is the log of a real factor the jet is too large by, as
 * `transmissionLogScale` from `tmmCoefficientJets`: it leaves every phase
 * quantity alone and divides `magnitudeSquared` by e^{2·logScale}.
 *
 * Returns `null` when the coefficient is exactly zero, where the phase and every
 * derivative of it are undefined.
 */
export function coefficientPhaseDispersion(coefficientJet, logScale = 0) {
    const [value, first, second, third] = jetDerivatives(coefficientJet);
    const magnitudeSquared = value[0] * value[0] + value[1] * value[1];
    if (magnitudeSquared === 0 || !Number.isFinite(magnitudeSquared)) return null;

    const valueJet = jetConstant(value[0], value[1]);
    const firstRatio = jetDivide(jetConstant(first[0], first[1]), valueJet)[0];
    const secondRatio = jetDivide(jetConstant(second[0], second[1]), valueJet)[0];
    const thirdRatio = jetDivide(jetConstant(third[0], third[1]), valueJet)[0];
    const squareFirst = [
        firstRatio[0] * firstRatio[0] - firstRatio[1] * firstRatio[1],
        2 * firstRatio[0] * firstRatio[1],
    ];
    const firstTimesSecond = [
        firstRatio[0] * secondRatio[0] - firstRatio[1] * secondRatio[1],
        firstRatio[0] * secondRatio[1] + firstRatio[1] * secondRatio[0],
    ];
    const cubeFirst = [
        squareFirst[0] * firstRatio[0] - squareFirst[1] * firstRatio[1],
        squareFirst[0] * firstRatio[1] + squareFirst[1] * firstRatio[0],
    ];

    const phaseRad = -Math.atan2(value[1], value[0]);
    return {
        phaseRad,
        phaseDeg: phaseRad * 180 / Math.PI,
        gd: firstRatio[1],
        gdd: secondRatio[1] - squareFirst[1],
        tod: thirdRatio[1] - 3 * firstTimesSecond[1] + 2 * cubeFirst[1],
        magnitudeSquared: logScale > 0 ? magnitudeSquared * Math.exp(-2 * logScale) : magnitudeSquared,
    };
}

/**
 * Exact derivatives of the reported phase quantities with respect to every layer
 * thickness, from the coefficient jet and its thickness-derivative jets.
 *
 * Returns `null` when the thickness jets are absent.
 */
export function coefficientPhaseThicknessDerivatives(coefficientJet, thicknessJets) {
    if (!thicknessJets) return null;
    const result = { phaseDeg: [], gd: [], gdd: [], tod: [], logMagnitudeSquared: [] };
    for (const thicknessJet of thicknessJets) {
        const logarithmicDerivative = jetDivide(thicknessJet, coefficientJet);
        const derivatives = jetDerivatives(logarithmicDerivative);
        // The imaginary part of d(ln c)/dd is the phase derivative; the real
        // part is d(ln |c|)/dd, so twice it is the relative derivative of |c|².
        result.phaseDeg.push(-derivatives[0][1] * 180 / Math.PI);
        result.gd.push(derivatives[1][1]);
        result.gdd.push(derivatives[2][1]);
        result.tod.push(derivatives[3][1]);
        result.logMagnitudeSquared.push(2 * derivatives[0][0]);
    }
    return result;
}

// ── Public entry points ──────────────────────────────────────────────────────

function prepare(lambda, layers, options) {
    const omega = options.omega ?? omegaFromLambdaNm(lambda);
    return {
        wavelengthJet: wavelengthOmegaJet(lambda, omega),
        incidentSineJet: options.sinTheta0Jet || null,
        layers: layers.map(layer => ({ indexJet: layer.nJet, thicknessNm: layer.d })),
    };
}

/**
 * Phase, group delay, GDD and TOD of a multilayer at one wavelength, for both
 * the reflection and the transmission coefficient.
 *
 * @param {number} lambda_nm
 * @param {number} theta_deg   angle of incidence, degrees from normal
 * @param {'s'|'p'} pol
 * @param {number[][]} n0Jet   incident-medium index jet
 * @param {number[][]} nsJet   substrate index jet
 * @param {{nJet:number[][], d:number}[]} layers   incident medium first, d in nm
 * @param {object} [options]
 * @param {number} [options.omega]  angular frequency; sets the time unit of the
 *   result. Defaults to 2πc/λ in rad/fs, giving fs, fs² and fs³.
 * @param {number[][]} [options.sinTheta0Jet]  sine of the incident angle as a
 *   jet. Give this when the stack is embedded in a dispersive medium and the
 *   angle held fixed is the external one, so that the internal angle disperses.
 *   Defaults to the constant sin(theta_deg).
 * @returns {{r: object|null, t: object|null}} each
 *   `{phaseRad, phaseDeg, gd, gdd, tod, magnitudeSquared}`, or `null` where that
 *   coefficient is exactly zero.
 */
export function tmmPhaseDispersion(lambda_nm, theta_deg, pol, n0Jet, nsJet, layers, options = {}) {
    const prepared = prepare(lambda_nm, layers, options);
    const coefficients = tmmCoefficientJets({
        wavelengthJet: prepared.wavelengthJet,
        thetaDeg: theta_deg,
        polarization: pol,
        incidentIndexJet: n0Jet,
        substrateIndexJet: nsJet,
        incidentSineJet: prepared.incidentSineJet,
        layers: prepared.layers,
    });
    return {
        r: coefficientPhaseDispersion(coefficients.reflection),
        t: coefficientPhaseDispersion(coefficients.transmission, coefficients.transmissionLogScale),
    };
}

/**
 * The same quantities plus their exact derivatives with respect to every layer
 * thickness, for gradient-based dispersion design.
 *
 * Arguments are those of `tmmPhaseDispersion`. Zero-thickness layers are kept so
 * derivative indices line up with the design array. Negative or NaN
 * thicknesses are skipped and receive zero derivative entries.
 *
 * @returns {{r, t}} where each side is the phase quantities plus
 *   `{dPhaseDeg, dGd, dGdd, dTod, dLogMagnitudeSquared}`, arrays of length
 *   `layers.length`. `dLogMagnitudeSquared` is d(ln |coefficient|²)/dd, the
 *   relative intensity derivative.
 */
export function tmmPhaseThicknessJacobian(lambda_nm, theta_deg, pol, n0Jet, nsJet, layers, options = {}) {
    const prepared = prepare(lambda_nm, layers, options);
    const coefficients = tmmCoefficientThicknessJets({
        wavelengthJet: prepared.wavelengthJet,
        thetaDeg: theta_deg,
        polarization: pol,
        incidentIndexJet: n0Jet,
        substrateIndexJet: nsJet,
        incidentSineJet: prepared.incidentSineJet,
        layers: prepared.layers,
    });
    const side = (coefficientJet, thicknessJets, logScale = 0) => {
        const base = coefficientPhaseDispersion(coefficientJet, logScale);
        if (!base) return null;
        const derivatives = coefficientPhaseThicknessDerivatives(coefficientJet, thicknessJets);
        return {
            ...base,
            dPhaseDeg: derivatives.phaseDeg,
            dGd: derivatives.gd,
            dGdd: derivatives.gdd,
            dTod: derivatives.tod,
            dLogMagnitudeSquared: derivatives.logMagnitudeSquared,
        };
    };
    return {
        r: side(coefficients.reflection, coefficients.reflectionThickness),
        t: side(coefficients.transmission, coefficients.transmissionThickness,
            coefficients.transmissionLogScale),
    };
}
