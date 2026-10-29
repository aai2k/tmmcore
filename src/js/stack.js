/**
 * The stack every call takes, read into the arrays the kernels use.
 *
 * A stack is { incident, substrate, layers }, the layers ordered from the
 * incident medium toward the substrate, each { n, k, d } with d in nm. A
 * medium is a number n or { n, k }. Each n and k is a number, or in a grid call
 * one value per wavelength; k may be left out and is then 0. ñ = n + ik, k ≥ 0
 * for absorbing media.
 *
 * The phase calls also take an index as a Taylor jet in ω (jet.js), which
 * carries the material's dispersion: one jet at one wavelength, a jet per
 * wavelength in a grid call. A number there is a constant index with no
 * dispersion. A list of plain numbers per wavelength is refused, since group
 * delay computed from it would leave out the material's dispersion.
 *
 * A grid call reads the stack into a grid, the layout of the C kernel as
 * JavaScript arrays, which the grid kernels of rta.js, derivatives.js,
 * phase.js and growing.js take:
 *   lambdas            wavelengths (nm)
 *   theta_deg          one angle for the whole grid, or one per wavelength (degrees)
 *   n0List, nsList     incident and substrate ñ per wavelength, [re, im]
 *   layerNK            [layer][λ] ñ, [re, im]
 *   thick              layer thicknesses (nm)
 * and for the phase kernels n0Jets, nsJets and layerJets, the same as Taylor
 * jets, with optional omegas and sinJets per wavelength (see
 * tmmPhaseDispersion). Arrays over the grid come out wavelength-major: the
 * value for wavelength i and layer k at [i·N + k].
 */

import { JET_ORDER, jetConstant } from './jet.js';

const isList = value => Array.isArray(value) || ArrayBuffer.isView(value);

// A list of `length` numbers.
const isNumbers = (value, length) => isList(value) && value.length === length
    && Array.from(value).every(entry => typeof entry === 'number');

/** A jet: JET_ORDER + 1 terms, each [re, im] as any list of two. */
export const isJet = value => Array.isArray(value) && value.length === JET_ORDER + 1
    && value.every(term => isList(term) && term.length === 2);

// A medium as { n, k }; k undefined when it was left out.
function readMedium(value, where) {
    if (typeof value === 'number') return { n: value, k: undefined };
    if (value !== null && typeof value === 'object' && !isList(value) && value.n !== undefined) {
        return { n: value.n, k: value.k };
    }
    throw new TypeError(`${where} must be a number or { n, k }`);
}

/**
 * Every call that takes a layer list reads one index key and one thickness key
 * from each layer: { n, d } for tmm() and its derivative and needle calls,
 * { nJet, d } for the phase calls, { indexJet, thicknessNm } for the
 * coefficient jets. A layer without either key throws here. Read as it stands
 * it would give an undefined thickness, which fails the d > 0 test and drops
 * the layer, so the result would be that of the stack without it.
 *
 * A thickness that is present but zero, negative or NaN passes, and keeps its
 * meaning in each call: an absent layer, and in the derivative calls a zero
 * column (zero thickness keeps its derivative there).
 */
export function checkLayerKeys(layers, indexKey, thicknessKey) {
    for (let i = 0; i < layers.length; i++) {
        const layer = layers[i];
        if (layer == null || layer[indexKey] === undefined || layer[thicknessKey] === undefined) {
            throw new TypeError(`layer ${i} needs {${indexKey}, ${thicknessKey}}`);
        }
    }
}

function readStack(stack) {
    if (stack === null || typeof stack !== 'object' || !Array.isArray(stack.layers)) {
        throw new TypeError('stack must be { incident, substrate, layers: [{ n, k, d }] }');
    }
    checkLayerKeys(stack.layers, 'n', 'd');
    stack.layers.forEach((layer, i) => {
        if (typeof layer.d !== 'number') throw new TypeError(`layer ${i} d must be a number (nm)`);
    });
    return {
        incident: readMedium(stack.incident, 'incident'),
        substrate: readMedium(stack.substrate, 'substrate'),
        layers: stack.layers.map(layer => ({ n: layer.n, k: layer.k, d: layer.d })),
    };
}

// Names each index of a stack for the messages below.
const describeIndices = stack => [
    ['incident', stack.incident], ['substrate', stack.substrate],
    ...stack.layers.map((layer, i) => [`layer ${i}`, layer]),
];

/** A wavelength of a one-wavelength call, in nm. */
export function readWavelength(wavelength) {
    if (typeof wavelength !== 'number') throw new TypeError('wavelength must be a number (nm)');
    return wavelength;
}

/** The wavelengths of a grid call, in nm. */
export function readWavelengths(wavelengths) {
    if (!isList(wavelengths) || !isNumbers(wavelengths, wavelengths.length)) {
        throw new TypeError('wavelengths must be a list of numbers (nm)');
    }
    return wavelengths;
}

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

/** Polarization of a phase call: 's' or 'p'. */
export function readPhasePol(pol) {
    if (pol !== 's' && pol !== 'p') throw new TypeError("pol must be 's' or 'p'");
    return pol;
}

// ── R/T calls: ñ as [re, im] ─────────────────────────────────────────────────

function checkValue(value, nLam, where) {
    if (typeof value === 'number') return;
    if (nLam !== null && isNumbers(value, nLam)) return;
    throw new TypeError(nLam === null
        ? `${where} must be a number at one wavelength`
        : `${where} must be a number or ${nLam} numbers, one per wavelength; n and k are two lists`);
}

const valueAt = (value, i) => (typeof value === 'number' ? value : value[i]);

function checkIndices(stack, nLam) {
    for (const [where, index] of describeIndices(stack)) {
        checkValue(index.n, nLam, `${where} n`);
        if (index.k !== undefined) checkValue(index.k, nLam, `${where} k`);
    }
}

/**
 * A stack at one wavelength, as the one-wavelength kernels take it:
 * { n0, ns, layers: [{ n, d }] } with each ñ as [re, im].
 */
export function pointStack(stack) {
    const s = readStack(stack);
    checkIndices(s, null);
    const nk = index => [index.n, index.k ?? 0];
    return {
        n0: nk(s.incident),
        ns: nk(s.substrate),
        layers: s.layers.map(layer => ({ n: nk(layer), d: layer.d })),
    };
}

/**
 * A stack over a wavelength grid, as the grid kernels take it: ñ per
 * wavelength for the two media, [layer][λ] for the layers, and the
 * thicknesses.
 */
export function gridStack(stack, nLam) {
    const s = readStack(stack);
    checkIndices(s, nLam);
    const nkList = index => Array.from({ length: nLam },
        (_, i) => [valueAt(index.n, i), index.k === undefined ? 0 : valueAt(index.k, i)]);
    return {
        n0List: nkList(s.incident),
        nsList: nkList(s.substrate),
        layerNK: s.layers.map(nkList),
        thick: s.layers.map(layer => layer.d),
    };
}

// ── Phase calls: ñ as a jet ──────────────────────────────────────────────────

const NO_DISPERSION = 'has no derivatives, so GD would miss the material\'s dispersion: give a jet per '
    + 'wavelength (see jetFromDerivatives), or one number for a constant index';

// One jet, or with nLam a jet per wavelength, for an index of a phase call.
function jetIndex(index, nLam, where) {
    if (typeof index.n === 'number') {
        if (index.k !== undefined && typeof index.k !== 'number') {
            throw new TypeError(`${where} k per wavelength ${NO_DISPERSION}`);
        }
        const jet = jetConstant(index.n, index.k ?? 0);
        return nLam === null ? jet : Array.from({ length: nLam }, () => jet);
    }
    if (index.k !== undefined) throw new TypeError(`${where} k goes inside the jet when n is a jet`);
    if (nLam === null && isJet(index.n)) return index.n;
    if (nLam !== null && Array.isArray(index.n) && index.n.length === nLam && index.n.every(isJet)) {
        return index.n;
    }
    if (nLam !== null && isList(index.n) && Array.from(index.n).every(value => typeof value === 'number')) {
        throw new TypeError(`${where} n per wavelength ${NO_DISPERSION}`);
    }
    throw new TypeError(nLam === null
        ? `${where} n must be a number or a jet`
        : `${where} n must be a number or ${nLam} jets, one per wavelength`);
}

/** A stack at one wavelength for the phase kernels: { n0Jet, nsJet, layers: [{ nJet, d }] }. */
export function pointJetStack(stack) {
    const s = readStack(stack);
    return {
        n0Jet: jetIndex(s.incident, null, 'incident'),
        nsJet: jetIndex(s.substrate, null, 'substrate'),
        layers: s.layers.map((layer, i) => ({ nJet: jetIndex(layer, null, `layer ${i}`), d: layer.d })),
    };
}

/** A stack over a grid for the phase kernels: jets per wavelength, [layer][λ] for the layers. */
export function gridJetStack(stack, nLam) {
    const s = readStack(stack);
    return {
        n0Jets: jetIndex(s.incident, nLam, 'incident'),
        nsJets: jetIndex(s.substrate, nLam, 'substrate'),
        layerJets: s.layers.map((layer, i) => jetIndex(layer, nLam, `layer ${i}`)),
        thick: s.layers.map(layer => layer.d),
    };
}

// ── Reading a grid ───────────────────────────────────────────────────────────

/** The angle of a grid at wavelength i: one for the whole grid, or one per wavelength. */
export const thetaAt = (theta_deg, i) => (typeof theta_deg === 'number' ? theta_deg : theta_deg[i]);

const layersAt = ({ layerNK, thick }, i) => Array.from(thick, (d, k) => ({ n: layerNK[k][i], d }));

/** The arguments every one-wavelength kernel takes first, at wavelength i of a grid. */
export const pointArgs = (grid, pol, i) =>
    [grid.lambdas[i], thetaAt(grid.theta_deg, i), pol, grid.n0List[i], grid.nsList[i], layersAt(grid, i)];
