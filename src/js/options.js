/**
 * The wavelengths and options of a call, read and checked: each reader returns
 * the value as the kernels take it, or throws a TypeError that names what is
 * wrong.
 */

import { isList, isNumbers, isJet } from './stack.js';

/** A wavelength of a one-wavelength call, in nm. */
export function readWavelength(wavelength) {
    if (typeof wavelength !== 'number') throw new TypeError('wavelength must be a number (nm)');
    return wavelength;
}

// A list of numbers as it was given, or a TypeError with `message`.
function numberList(value, message) {
    if (!isList(value) || !isNumbers(value, value.length)) throw new TypeError(message);
    return value;
}

/** The wavelengths of a grid call, in nm. */
export const readWavelengths = wavelengths => numberList(wavelengths, 'wavelengths must be a list of numbers (nm)');

/**
 * The angle of incidence in degrees, 0 when left out: a number, or in a grid
 * call (nLam given) also one angle per wavelength.
 */
export function readAngle(angle = 0, nLam = null) {
    if (typeof angle === 'number') return angle;
    if (nLam !== null && isNumbers(angle, nLam)) return angle;
    throw new TypeError(nLam === null
        ? 'angle must be a number (degrees)'
        : `angle must be a number or ${nLam} numbers, one per wavelength (degrees)`);
}

/** Polarization of an R/T call: 's', 'p' or 'avg', the mean of the two. */
export function readRtaPol(pol = 'avg') {
    if (pol !== 's' && pol !== 'p' && pol !== 'avg') throw new TypeError("pol must be 's', 'p' or 'avg'");
    return pol;
}

/**
 * Polarization of a phase call: 's' or 'p'. It may be left out at normal
 * incidence, where the two are the same, and is then 's'.
 */
export function readPhasePol(pol, normal) {
    if (pol === undefined && normal) return 's';
    if (pol !== 's' && pol !== 'p') throw new TypeError("pol must be 's' or 'p'; it may be left out at normal incidence");
    return pol;
}

// A list of `length` jets. Array.from reads a hole in a sparse list as undefined, which is not a jet.
const isJets = (value, length) => Array.isArray(value) && value.length === length && Array.from(value).every(isJet);

/**
 * The sine of the angle of incidence as a jet in ω, which a phase call takes in
 * place of the angle: one jet at one wavelength, a jet per wavelength in a grid
 * call (nLam given).
 */
export function readSine(sine, nLam = null) {
    if (nLam === null ? isJet(sine) : isJets(sine, nLam)) return sine;
    throw new TypeError(nLam === null ? 'sine must be a jet' : `sine must be ${nLam} jets, one per wavelength`);
}

const isZeroJet = jet => jet.every(term => term[0] === 0 && term[1] === 0);

/** Light along the normal at every wavelength: every angle 0, or every term of the sine jets 0. */
export function atNormalIncidence(angle, sine, nLam = null) {
    if (sine) return nLam === null ? isZeroJet(sine) : sine.every(isZeroJet);
    return typeof angle === 'number' ? angle === 0 : angle.every(value => value === 0);
}

/**
 * Where a needle scan also opens each layer, as fractions of its thickness from
 * the incident side: none when left out.
 */
export const readFractions = (fractions = []) => numberList(fractions, 'fractions must be a list of numbers');

/** The thicknesses of the growing layer a monitoring curve samples, in nm. */
export const readThicknesses = thicknesses => numberList(thicknesses, 'thicknesses must be a list of numbers (nm)');

/** The thickness of the growing layer a growing stack samples, in nm. */
export function readThickness(d) {
    if (typeof d !== 'number') throw new TypeError('d must be a number (nm)');
    return d;
}
