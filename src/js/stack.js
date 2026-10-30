/**
 * The stack every call takes, read into the arrays the kernels use.
 *
 * A stack is { incident, substrate, layers }, the layers ordered from the
 * incident medium toward the substrate, each { n, k, d } with d in nm. A
 * medium is a number n or { n, k }. Each n and k is a number, or in a grid call
 * one value per wavelength; k may be left out and is then 0. ñ = n + ik, k ≥ 0
 * for absorbing media. A key the stack, a layer or a medium does not take
 * throws, so a misspelled k cannot read as a lossless layer.
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

/** An array or a typed array. */
export const isList = value => Array.isArray(value) || ArrayBuffer.isView(value);

/** A list of `length` numbers. A hole in a sparse array is not a number. */
export function isNumbers(value, length) {
    if (!isList(value) || value.length !== length) return false;
    for (let i = 0; i < length; i++) if (typeof value[i] !== 'number') return false;
    return true;
}

/** A jet: JET_ORDER + 1 terms, each [re, im] as any list of two. */
export const isJet = value => Array.isArray(value) && value.length === JET_ORDER + 1
    && value.every(term => isList(term) && term.length === 2);

/** Keys in a message: 'a and b', 'a, b and c'. */
export const listed = keys => `${keys.slice(0, -1).join(', ')} and ${keys[keys.length - 1]}`;

// An index of a stack in a message: 'incident', 'substrate', or a layer by its
// position. A call builds the name only when it throws.
const indexName = where => (typeof where === 'number' ? `layer ${where}` : where);

const STACK_KEYS = ['incident', 'substrate', 'layers'];
const LAYER_KEYS = ['n', 'k', 'd'];
const MEDIUM_KEYS = ['n', 'k'];

// A key that `what`, a stack, a layer or a medium, does not take throws and
// names the keys it does.
function checkKeys(value, keys, where, what) {
    for (const key in value) {
        if (!keys.includes(key)) throw new TypeError(`${indexName(where)} has '${key}'; ${what} takes ${listed(keys)}`);
    }
}

// A medium as { n, k }; k undefined when it was left out.
function readMedium(value, where) {
    if (typeof value === 'number') return { n: value, k: undefined };
    if (value !== null && typeof value === 'object' && !isList(value) && value.n !== undefined) {
        checkKeys(value, MEDIUM_KEYS, where, 'a medium');
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

// The stack with each medium as { n, k }, and its layers as given, checked in
// this order: the stack's keys, every layer's n and d, then each layer's other
// keys and its thickness, then the media.
function readStack(stack) {
    if (stack === null || typeof stack !== 'object' || !Array.isArray(stack.layers)) {
        throw new TypeError('stack must be { incident, substrate, layers: [{ n, k, d }] }');
    }
    checkKeys(stack, STACK_KEYS, 'stack', 'a stack');
    checkLayerKeys(stack.layers, 'n', 'd');
    stack.layers.forEach((layer, i) => {
        checkKeys(layer, LAYER_KEYS, i, 'a layer');
        if (typeof layer.d !== 'number') throw new TypeError(`layer ${i} d must be a number (nm)`);
    });
    return {
        incident: readMedium(stack.incident, 'incident'),
        substrate: readMedium(stack.substrate, 'substrate'),
        layers: stack.layers,
    };
}

// ── R/T calls: ñ as [re, im] ─────────────────────────────────────────────────

function checkValue(value, nLam, where, key) {
    if (typeof value === 'number') return;
    if (nLam !== null && isNumbers(value, nLam)) return;
    const what = `${indexName(where)} ${key}`;
    throw new TypeError(nLam === null
        ? `${what} must be a number at one wavelength`
        : `${what} must be a number or ${nLam} numbers, one per wavelength; n and k are two lists`);
}

/** n or k of an index at wavelength i of a grid: one value for the whole grid, or one per wavelength. */
export const valueAt = (value, i) => (typeof value === 'number' ? value : value[i]);

function checkIndex(index, nLam, where) {
    checkValue(index.n, nLam, where, 'n');
    if (index.k !== undefined) checkValue(index.k, nLam, where, 'k');
}

/** ñ of a checked index at one wavelength, as [re, im]. */
export const pairOf = index => [index.n, index.k ?? 0];

/** ñ of a checked index at each of nLam wavelengths, as [re, im]. */
export const pairsOver = (index, nLam) => Array.from({ length: nLam },
    (_, i) => [valueAt(index.n, i), index.k === undefined ? 0 : valueAt(index.k, i)]);

/**
 * The stack of an R/T call, checked: { incident, substrate, layers }, each
 * medium as { n, k } with k undefined when it was left out, and the layers as
 * given. nLam is the grid size of a Spectrum call, null at one wavelength. A
 * WebAssembly instance writes it into kernel memory as it is (memory.js).
 */
export function checkedStack(stack, nLam) {
    const s = readStack(stack);
    checkIndex(s.incident, nLam, 'incident');
    checkIndex(s.substrate, nLam, 'substrate');
    s.layers.forEach((layer, i) => checkIndex(layer, nLam, i));
    return s;
}

/**
 * A stack at one wavelength, as the one-wavelength kernels take it:
 * { n0, ns, layers: [{ n, d }] } with each ñ as [re, im].
 */
export function pointStack(stack) {
    const s = checkedStack(stack, null);
    return {
        n0: pairOf(s.incident),
        ns: pairOf(s.substrate),
        layers: s.layers.map(layer => ({ n: pairOf(layer), d: layer.d })),
    };
}

/**
 * A stack over a wavelength grid, as the grid kernels take it: ñ per
 * wavelength for the two media, [layer][λ] for the layers, and the
 * thicknesses.
 */
export function gridStack(stack, nLam) {
    const s = checkedStack(stack, nLam);
    return {
        n0List: pairsOver(s.incident, nLam),
        nsList: pairsOver(s.substrate, nLam),
        layerNK: s.layers.map(layer => pairsOver(layer, nLam)),
        thick: s.layers.map(layer => layer.d),
    };
}

/**
 * An index given as a medium is, a number or { n, k }, checked as one and
 * returned as { n, k } with k undefined when it was left out. nLam is the grid
 * size of a Spectrum call, null at one wavelength; `where` names the index in
 * a message.
 */
export function readIndex(value, nLam, where) {
    const index = readMedium(value, where);
    checkIndex(index, nLam, where);
    return index;
}

/** The materials a needle scan tries: a list of indices, each read by readIndex. */
export function readCandidates(candidates, nLam) {
    if (!isList(candidates)) throw new TypeError('candidates must be a list of indices, each a number or { n, k }');
    return Array.from(candidates, (value, c) => readIndex(value, nLam, `candidate ${c}`));
}

// ── Phase calls: ñ as a jet ──────────────────────────────────────────────────

const NO_DISPERSION = 'has no derivatives, so GD would miss the material\'s dispersion: give a jet per '
    + 'wavelength (see jetFromDerivatives), or one number for a constant index';

// One jet, or with nLam a jet per wavelength, for an index of a phase call;
// `where` names the index as indexName reads it.
function jetIndex(index, nLam, where) {
    if (typeof index.n === 'number') {
        if (index.k !== undefined && typeof index.k !== 'number') {
            throw new TypeError(`${indexName(where)} k per wavelength ${NO_DISPERSION}`);
        }
        const jet = jetConstant(index.n, index.k ?? 0);
        return nLam === null ? jet : Array.from({ length: nLam }, () => jet);
    }
    if (index.k !== undefined) throw new TypeError(`${indexName(where)} k goes inside the jet when n is a jet`);
    if (nLam === null && isJet(index.n)) return index.n;
    if (nLam !== null && Array.isArray(index.n) && index.n.length === nLam && index.n.every(isJet)) {
        return index.n;
    }
    if (nLam !== null && isList(index.n) && isNumbers(index.n, index.n.length)) {
        throw new TypeError(`${indexName(where)} n per wavelength ${NO_DISPERSION}`);
    }
    throw new TypeError(nLam === null
        ? `${indexName(where)} n must be a number or a jet`
        : `${indexName(where)} n must be a number or ${nLam} jets, one per wavelength`);
}

/** A stack at one wavelength for the phase kernels: { n0Jet, nsJet, layers: [{ nJet, d }] }. */
export function pointJetStack(stack) {
    const s = readStack(stack);
    return {
        n0Jet: jetIndex(s.incident, null, 'incident'),
        nsJet: jetIndex(s.substrate, null, 'substrate'),
        layers: s.layers.map((layer, i) => ({ nJet: jetIndex(layer, null, i), d: layer.d })),
    };
}

/** A stack over a grid for the phase kernels: jets per wavelength, [layer][λ] for the layers. */
export function gridJetStack(stack, nLam) {
    const s = readStack(stack);
    return {
        n0Jets: jetIndex(s.incident, nLam, 'incident'),
        nsJets: jetIndex(s.substrate, nLam, 'substrate'),
        layerJets: s.layers.map((layer, i) => jetIndex(layer, nLam, i)),
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
