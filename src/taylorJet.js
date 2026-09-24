/**
 * taylorJet.js : third-order truncated Taylor arithmetic for complex functions.
 *
 * A jet stores the coefficients of a truncated power series:
 *
 *     [ f, f', f''/2!, f'''/3! ]
 *
 * each entry a complex `[re, im]` pair. Every operation below is ordinary power-
 * series algebra, so composing them differentiates a function exactly, with no
 * finite differences and no step size to choose. Feed a jet through the same
 * code that computes a value and the derivatives come out alongside it.
 *
 * This is what makes the phase-dispersion kernel possible: carry the
 * characteristic matrix in jets of angular frequency and the reflection
 * coefficient emerges with the three derivatives that group delay, GDD and TOD
 * are built from.
 *
 * The differentiation variable is whatever the caller chose. In `phase.js` it is
 * angular frequency, so a jet of the refractive index means
 * `[ñ, dñ/dω, (d²ñ/dω²)/2, (d³ñ/dω³)/6]`. Use `jetFromDerivatives` to build one
 * from plain derivatives without doing the factorials by hand, and
 * `jetDerivatives` to read them back.
 *
 * All arithmetic is double precision.
 */

/** Highest derivative carried. Jets are arrays of `JET_ORDER + 1` entries. */
export const JET_ORDER = 3;

const zero = () => [0, 0];
const addComplex = (a, b) => [a[0] + b[0], a[1] + b[1]];
const subComplex = (a, b) => [a[0] - b[0], a[1] - b[1]];
const scaleComplex = (a, scalar) => [a[0] * scalar, a[1] * scalar];
const multiplyComplex = (a, b) => [
    a[0] * b[0] - a[1] * b[1],
    a[0] * b[1] + a[1] * b[0],
];
const divideComplex = (a, b) => {
    const denominator = b[0] * b[0] + b[1] * b[1];
    return [
        (a[0] * b[0] + a[1] * b[1]) / denominator,
        (a[1] * b[0] - a[0] * b[1]) / denominator,
    ];
};
const sqrtComplex = (value) => {
    const magnitudeRoot = Math.sqrt(Math.hypot(value[0], value[1]));
    const angle = Math.atan2(value[1], value[0]) / 2;
    return [magnitudeRoot * Math.cos(angle), magnitudeRoot * Math.sin(angle)];
};

/** A constant: value with all derivatives zero. */
export function jetConstant(value, imaginary = 0) {
    return [[value, imaginary], zero(), zero(), zero()];
}

/**
 * Build a jet from plain derivatives, applying the factorials.
 * Each argument is a number or a complex `[re, im]` pair.
 */
export function jetFromDerivatives(value, first = 0, second = 0, third = 0) {
    const asComplex = item => Array.isArray(item) ? [...item] : [item, 0];
    return [
        asComplex(value),
        asComplex(first),
        scaleComplex(asComplex(second), 1 / 2),
        scaleComplex(asComplex(third), 1 / 6),
    ];
}

/** Read a jet back as plain derivatives `[f, f', f'', f''']`. */
export function jetDerivatives(jet) {
    return [jet[0], jet[1], scaleComplex(jet[2], 2), scaleComplex(jet[3], 6)];
}

export function jetAdd(left, right) {
    return left.map((value, index) => addComplex(value, right[index]));
}

export function jetSubtract(left, right) {
    return left.map((value, index) => subComplex(value, right[index]));
}

/** Multiply every order by a real scalar. */
export function jetScale(jet, scalar) {
    return jet.map(value => scaleComplex(value, scalar));
}

export function jetMultiply(left, right) {
    const result = Array.from({ length: JET_ORDER + 1 }, zero);
    for (let order = 0; order <= JET_ORDER; order++) {
        for (let index = 0; index <= order; index++) {
            result[order] = addComplex(
                result[order],
                multiplyComplex(left[index], right[order - index]),
            );
        }
    }
    return result;
}

export function jetReciprocal(jet) {
    const result = Array.from({ length: JET_ORDER + 1 }, zero);
    result[0] = divideComplex([1, 0], jet[0]);
    for (let order = 1; order <= JET_ORDER; order++) {
        let sum = zero();
        for (let index = 1; index <= order; index++) {
            sum = addComplex(sum, multiplyComplex(jet[index], result[order - index]));
        }
        result[order] = scaleComplex(divideComplex(sum, jet[0]), -1);
    }
    return result;
}

export function jetDivide(numerator, denominator) {
    return jetMultiply(numerator, jetReciprocal(denominator));
}

/** Principal square root, branch matching `csqrt` in tmm.js. */
export function jetSqrt(jet) {
    const result = Array.from({ length: JET_ORDER + 1 }, zero);
    result[0] = sqrtComplex(jet[0]);
    const twiceRoot = scaleComplex(result[0], 2);
    for (let order = 1; order <= JET_ORDER; order++) {
        let known = zero();
        for (let index = 1; index < order; index++) {
            known = addComplex(known, multiplyComplex(result[index], result[order - index]));
        }
        result[order] = divideComplex(subComplex(jet[order], known), twiceRoot);
    }
    return result;
}

export function jetExp(jet) {
    const result = Array.from({ length: JET_ORDER + 1 }, zero);
    const magnitude = Math.exp(jet[0][0]);
    result[0] = [magnitude * Math.cos(jet[0][1]), magnitude * Math.sin(jet[0][1])];
    for (let order = 1; order <= JET_ORDER; order++) {
        let sum = zero();
        for (let index = 1; index <= order; index++) {
            sum = addComplex(
                sum,
                scaleComplex(multiplyComplex(jet[index], result[order - index]), index),
            );
        }
        result[order] = scaleComplex(sum, 1 / order);
    }
    return result;
}

export function jetLog(jet) {
    const result = Array.from({ length: JET_ORDER + 1 }, zero);
    result[0] = [Math.log(Math.hypot(jet[0][0], jet[0][1])), Math.atan2(jet[0][1], jet[0][0])];
    const derivative = [jet[1], scaleComplex(jet[2], 2), scaleComplex(jet[3], 3), zero()];
    const quotient = jetMultiply(derivative, jetReciprocal(jet));
    for (let order = 1; order <= JET_ORDER; order++) {
        result[order] = scaleComplex(quotient[order - 1], 1 / order);
    }
    return result;
}

/** Real power via exp(p log z); the log branch cut applies. */
export function jetPower(jet, exponent) {
    return jetExp(jetScale(jetLog(jet), exponent));
}

/** Sine and cosine together, since the recurrences share their terms. */
export function jetSinCos(jet) {
    const sine = Array.from({ length: JET_ORDER + 1 }, zero);
    const cosine = Array.from({ length: JET_ORDER + 1 }, zero);
    const [real, imaginary] = jet[0];
    sine[0] = [Math.sin(real) * Math.cosh(imaginary), Math.cos(real) * Math.sinh(imaginary)];
    cosine[0] = [Math.cos(real) * Math.cosh(imaginary), -Math.sin(real) * Math.sinh(imaginary)];
    for (let order = 1; order <= JET_ORDER; order++) {
        let sineSum = zero();
        let cosineSum = zero();
        for (let index = 1; index <= order; index++) {
            sineSum = addComplex(
                sineSum,
                scaleComplex(multiplyComplex(jet[index], cosine[order - index]), index),
            );
            cosineSum = addComplex(
                cosineSum,
                scaleComplex(multiplyComplex(jet[index], sine[order - index]), index),
            );
        }
        sine[order] = scaleComplex(sineSum, 1 / order);
        cosine[order] = scaleComplex(cosineSum, -1 / order);
    }
    return { sine, cosine };
}

export const jetSin = jet => jetSinCos(jet).sine;
export const jetCos = jet => jetSinCos(jet).cosine;

/**
 * Chain rule for a scalar function known only through its derivatives.
 *
 * Given f(x0) and [f'(x0), f''(x0), f'''(x0)], and a jet for x, return the jet
 * for f(x). This is the bridge from a dispersion formula differentiated by hand
 * or by a symbolic tool into the jet arithmetic here.
 */
export function jetCompose(value, derivatives, xJet) {
    const displacement = xJet.map((coefficient, index) => index === 0 ? zero() : coefficient);
    const first = jetScale(displacement, derivatives[0] ?? 0);
    const second = jetScale(jetMultiply(displacement, displacement), (derivatives[1] ?? 0) / 2);
    const third = jetScale(
        jetMultiply(jetMultiply(displacement, displacement), displacement),
        (derivatives[2] ?? 0) / 6,
    );
    return jetAdd(jetAdd(jetConstant(value), first), jetAdd(second, third));
}

/**
 * The jet of wavelength as a function of angular frequency, λ(ω) = 2πc/ω.
 *
 * Needs no value for c: with λ and ω both given, every derivative follows from
 * λ' = -λ/ω. Units are therefore the caller's own, and the wavelength unit of
 * the result is whatever `lambda` was given in.
 */
export function wavelengthOmegaJet(lambda, omega) {
    // Written as repeated multiplication rather than `omega ** 3` so the C
    // kernel can produce the identical double.
    return [
        [lambda, 0],
        [-lambda / omega, 0],
        [lambda / (omega * omega), 0],
        [-lambda / (omega * omega * omega), 0],
    ];
}

/** Combine a real-valued jet and an imaginary-valued jet into one complex jet. */
export function jetWithImaginaryPart(realJet, imaginaryJet) {
    return realJet.map((coefficient, index) => [coefficient[0], imaginaryJet[index][0]]);
}

/** The real part of a complex jet, order by order. */
export function jetRealPart(jet) {
    return jet.map(coefficient => [coefficient[0], 0]);
}

/** Floor the value at `minimum`, flattening the jet to a constant when it bites. */
export function jetClampRealMinimum(jet, minimum) {
    return jet[0][0] >= minimum ? jet : jetConstant(minimum);
}

/**
 * Bound the imaginary part, holding it at the limit at every order when it
 * bites; the real part, which carries the phase, is kept.
 *
 * Guards the same overflow `layerMatrix` guards in tmm.js: the phase thickness
 * of a strongly absorbing layer grows without bound and cosh of it overflows.
 * Past the limit the layer matrix built from the held phase is the true one
 * divided by a real factor, which cancels from r and from every phase; callers
 * that need |t| carry it as a log scale, as tmm.js does with layerLogScale.
 */
export function jetClampImaginary(jet, limit) {
    if (jet[0][1] > limit) {
        return jet.map((coefficient, index) => [coefficient[0], index === 0 ? limit : 0]);
    }
    if (jet[0][1] < -limit) {
        return jet.map((coefficient, index) => [coefficient[0], index === 0 ? -limit : 0]);
    }
    return jet;
}

export function jetIsFinite(jet) {
    return jet.every(coefficient => coefficient.every(Number.isFinite));
}
