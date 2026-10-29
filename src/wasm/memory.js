/**
 * memory.js : kernel memory for the bindings of a WebAssembly instance: the
 * methods that place inputs and read outputs, and onGrid, which places a
 * wavelength grid in the layout the grid exports take.
 *
 * Memory growth detaches every view of the old buffer, so a view is made after
 * the allocations, and after the export call, that could grow it.
 */

import { writeJet } from './jet.js';

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

const writePair = (view, at, pair) => {
    view[at] = pair[0];
    view[at + 1] = pair[1];
};

// Where a grid keeps its indices, and how one is written: ñ as [re, im] in 2
// doubles, a jet as 4 × [re, im] in 8.
export const PAIRS = { incident: 'n0List', substrate: 'nsList', layers: 'layerNK', width: 2, write: writePair };
export const JETS = { incident: 'n0Jets', substrate: 'nsJets', layers: 'layerJets', width: 8, write: writeJet };

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
 * putTable(table, nRows, kind) places table[r][i] for r < nRows and every
 * wavelength i, row-major, each entry an index of `kind`. Entries are read by
 * index, so a row may be any list, typed or not, and longer than the grid.
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
                for (let i = 0; i < nLam; i++) kind.write(view, (r * nLam + i) * kind.width, table[r][i]);
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
