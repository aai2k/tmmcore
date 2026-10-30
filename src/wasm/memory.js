/**
 * memory.js : kernel memory for the bindings of a WebAssembly instance: the
 * methods that place inputs and read outputs, the writers of a call's layers
 * at one wavelength, and onGrid, which places a wavelength grid in the layout
 * the grid exports take.
 *
 * Memory growth detaches every view of the old buffer, so a view is made after
 * the allocations, and after the export call, that could grow it. Inputs are
 * read while they are written into such a view, so a getter in them that calls
 * the same instance, or an input that is itself a view of this memory, gives
 * wrong values.
 */

import { writeJet } from './jet.js';
import { valueAt } from '../js/stack.js';

export const memoryMethods = {
    _alloc(nDoubles) {
        const ptr = this.malloc(nDoubles * 8);
        if (!ptr) throw new Error('tmmWasm: malloc failed');
        return ptr;
    },
    // Fresh view : memory.buffer is detached after any growth, so re-create
    // views AFTER all mallocs for a call are done.
    _view(ptr, nDoubles) {
        return new Float64Array(this.memory.buffer, ptr, nDoubles);
    },

    // Persistent scratch arena for the one-wavelength calls.
    // These are invoked thousands of times per optimization run; malloc/free +
    // typed-array churn per call would dominate the JS↔WASM boundary cost and
    // can make per-call WASM SLOWER than JS. Reusing one buffer (grown on demand,
    // never freed between calls) makes each call just write-args / read-result.
    _scratch(nDoubles) {
        if (this._scratchN < nDoubles) {
            if (this._scratchPtr) this.free(this._scratchPtr);
            this._scratchPtr = this._alloc(nDoubles);
            this._scratchN = nDoubles;
        }
        return this._scratchPtr;
    },
};

/**
 * The layers of a call at one wavelength as the kernels take them, n, k and d
 * per layer from view[0]. The methods of 0.4 carry ñ as [re, im] in n
 * (writePairLayers); a stack's layers carry n and k, a missing k read as 0
 * (writeStackLayers, stack.js).
 */
export function writePairLayers(view, layers) {
    for (let i = 0; i < layers.length; i++) {
        view[3 * i] = layers[i].n[0];
        view[3 * i + 1] = layers[i].n[1];
        view[3 * i + 2] = layers[i].d;
    }
}

export function writeStackLayers(view, layers) {
    for (let i = 0; i < layers.length; i++) {
        view[3 * i] = layers[i].n;
        view[3 * i + 1] = layers[i].k ?? 0;
        view[3 * i + 2] = layers[i].d;
    }
}

const writePair = (view, at, pair) => {
    view[at] = pair[0];
    view[at + 1] = pair[1];
};

// Where a grid keeps its indices, and how row `row` of them is written at
// wavelength i: ñ as [re, im] in 2 doubles, a jet as 4 × [re, im] in 8. A
// PAIRS or JETS row holds one entry per wavelength. A STACK row is a medium or
// layer of a checked stack (stack.js), whose n and k are each one value or one
// per wavelength.
export const PAIRS = { incident: 'n0List', substrate: 'nsList', layers: 'layerNK', width: 2,
    write: (view, at, row, i) => writePair(view, at, row[i]) };
export const JETS = { incident: 'n0Jets', substrate: 'nsJets', layers: 'layerJets', width: 8,
    write: (view, at, row, i) => writeJet(view, at, row[i]) };
export const STACK = { incident: 'incident', substrate: 'substrate', layers: 'layers', width: 2,
    write: (view, at, index, i) => {
        view[at] = valueAt(index.n, i);
        view[at + 1] = index.k === undefined ? 0 : valueAt(index.k, i);
    } };

export const named = (keys, values) => Object.fromEntries(keys.map((key, i) => [key, values[i]]));

/**
 * Runs `body` with the grid's inputs in kernel memory. `lead(pol)` gives the
 * arguments every grid export starts with: wavelengths, angles, pol when the
 * export takes one, the two media, the layer table and the thicknesses. The
 * body places its own blocks with `put`, `putTable` and `out` and reads them
 * with `read`; all are freed when it returns or throws. Inputs are written as
 * they are placed, which a later allocation cannot undo, and `read` makes a
 * fresh view, since memory growth detaches older ones.
 *
 * putTable(table, nRows, kind) places row table[r] for r < nRows at every
 * wavelength i, row-major, as `kind` writes it. A PAIRS or JETS row is read by
 * index, so it may be any list, typed or not, and longer than the grid.
 */
export function onGrid(wasm, grid, indices, body) {
    const blocks = [];
    const out = size => {
        const ptr = wasm._alloc(Math.max(1, size));
        blocks.push(ptr);
        return ptr;
    };
    const put = values => {
        const ptr = out(values.length);
        wasm._view(ptr, values.length).set(values);
        return ptr;
    };
    const read = (ptr, size) => Float64Array.from(wasm._view(ptr, size));
    try {
        const nLam = grid.lambdas.length;
        const N = grid.thick.length;
        const putTable = (table, nRows, kind) => {
            const ptr = out(nRows * nLam * kind.width);
            const view = wasm._view(ptr, nRows * nLam * kind.width);
            for (let r = 0; r < nRows; r++) {
                for (let i = 0; i < nLam; i++) kind.write(view, (r * nLam + i) * kind.width, table[r], i);
            }
            return ptr;
        };
        const angles = typeof grid.theta_deg === 'number' ? [grid.theta_deg] : grid.theta_deg;
        const head = [put(grid.lambdas), nLam, put(angles), angles.length];
        const tail = [
            putTable([grid[indices.incident]], 1, indices),
            putTable([grid[indices.substrate]], 1, indices),
            putTable(grid[indices.layers], N, indices),
            put(grid.thick), N,
        ];
        const lead = pol => (pol === undefined ? [...head, ...tail] : [...head, pol, ...tail]);
        return body({ lead, put, putTable, out, read });
    } finally {
        for (const ptr of blocks) wasm.free(ptr);
    }
}
