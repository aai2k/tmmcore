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
