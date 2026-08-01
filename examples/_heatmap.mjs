/**
 * Minimal SVG heat map — no dependencies.
 *
 * Companion to _plot.mjs, for quantities that vary over two axes. Cells are
 * quantised to a fixed number of levels and then emitted as one <path> per
 * level, with horizontal runs merged. That keeps a 121 x 81 two-panel grid to
 * tens of kilobytes rather than the ~1 MB a naive <rect> per cell would cost,
 * because the colour is written once per level instead of once per cell.
 *
 * The colour ramp is turbo, sampled at 17 anchors and interpolated. It is a
 * rainbow ramp, so it reads well against both light and dark page backgrounds
 * (GitHub renders markdown in either and strips <style> blocks from inline
 * SVG), while avoiding the perceptual artefacts of jet, which manufactures
 * banding that is not present in the data.
 * Ramp data: Anton Mikhailov, Google LLC, Apache-2.0.
 */

const AXIS = '#8b949e';

const TURBO = [
    [48, 18, 59], [64, 64, 162], [70, 107, 227], [66, 148, 255],
    [40, 188, 235], [24, 221, 194], [50, 242, 152], [109, 254, 98],
    [164, 252, 60], [203, 237, 52], [236, 209, 58], [253, 174, 53],
    [251, 129, 34], [236, 83, 15], [210, 49, 5], [172, 23, 1],
    [122, 4, 3],
];

// A blue / neutral / red ramp for signed quantities, where the midpoint of the
// scale is zero and the two directions must be told apart at a glance. Chosen
// here rather than taken from a published scheme; the neutral centre is light
// enough to read as "no effect" on either page background.
const DIVERGING = [
    [33, 90, 170], [70, 140, 205], [130, 185, 225], [195, 220, 240],
    [240, 240, 240],
    [248, 214, 196], [235, 160, 120], [210, 100, 70], [165, 30, 40],
];

const lerpRamp = (ramp, t) => {
    const u = Math.min(1, Math.max(0, t)) * (ramp.length - 1);
    const i = Math.min(ramp.length - 2, Math.floor(u));
    const f = u - i;
    const c = ramp[i].map((v, k) => Math.round(v + f * (ramp[i + 1][k] - v)));
    return '#' + c.map(v => v.toString(16).padStart(2, '0')).join('');
};

/** value in [0,1] -> "#rrggbb" */
export const turbo = t => lerpRamp(TURBO, t);
export const diverging = t => lerpRamp(DIVERGING, t);

const RAMPS = { turbo, diverging };

/**
 * panels: [{ title, z }] where z[iy][ix] holds values in [zMin, zMax].
 * x, y:   axis coordinate arrays, length matching z's dimensions.
 */
export function heatMap({
    panels, x, y, xLabel = '', yLabel = '', zLabel = '', title = '',
    zMin = 0, zMax = 1, levels = 64, cellW = 3, cellH = 3, ramp = 'turbo', zFormat = null, zTransform = null, zTickValues = null,
}) {
    const m = { top: title ? 62 : 34, right: 92, bottom: 54, left: 62 };
    const gap = 34;
    const pw = x.length * cellW;
    const ph = y.length * cellH;
    const width = m.left + panels.length * pw + (panels.length - 1) * gap + m.right;
    const height = m.top + ph + m.bottom;

    // Quantising before colour lookup is what makes the run merging effective:
    // neighbouring cells that round to the same level collapse into one run.
    // zTransform maps a value to its [0,1] position on the ramp. A non-linear
    // scale is sometimes the only way to show a quantity whose interesting
    // structure spans orders of magnitude; when one is used the caller must
    // also supply zTickValues so the colour bar still tells the truth.
    const toT = zTransform ?? (v => (v - zMin) / (zMax - zMin));
    const quant = v => Math.round(Math.min(1, Math.max(0, toT(v))) * (levels - 1));
    const paint = RAMPS[ramp] ?? turbo;
    const swatch = Array.from({ length: levels }, (_, i) => paint(i / (levels - 1)));

    const out = [];
    out.push(`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" font-family="system-ui, sans-serif" font-size="12">`);
    if (title) out.push(`<text x="${(width / 2).toFixed(1)}" y="26" fill="${AXIS}" font-size="15" text-anchor="middle">${esc(title)}</text>`);

    const fmt = v => String(Math.round(v * 10) / 10);
    const fmtZ = zFormat ?? fmt;

    panels.forEach((panel, p) => {
        const ox = m.left + p * (pw + gap);

        out.push(`<text x="${(ox + pw / 2).toFixed(1)}" y="${m.top - 10}" fill="${AXIS}" font-size="13" text-anchor="middle">${esc(panel.title)}</text>`);
        out.push(`<g shape-rendering="crispEdges">`);

        // Collect runs per quantised level, then emit one path each. Cells are
        // one pixel wider and taller than their pitch so neighbours overlap and
        // no hairline seams show through.
        const byLevel = new Map();
        for (let iy = 0; iy < y.length; iy++) {
            const row = panel.z[iy];
            // y runs bottom-up in data order, top-down in screen order.
            const py = m.top + ph - (iy + 1) * cellH;
            let ix = 0;
            while (ix < row.length) {
                const q = quant(row[ix]);
                let run = 1;
                while (ix + run < row.length && quant(row[ix + run]) === q) run++;
                if (!byLevel.has(q)) byLevel.set(q, []);
                byLevel.get(q).push(`M${ox + ix * cellW} ${py}h${run * cellW + 1}v${cellH + 1}h${-(run * cellW + 1)}z`);
                ix += run;
            }
        }
        for (const [q, subpaths] of [...byLevel].sort((a, b) => a[0] - b[0])) {
            out.push(`<path fill="${swatch[q]}" d="${subpaths.join('')}"/>`);
        }
        out.push('</g>');
        out.push(`<rect x="${ox.toFixed(1)}" y="${m.top}" width="${pw.toFixed(1)}" height="${ph.toFixed(1)}" fill="none" stroke="${AXIS}" stroke-opacity="0.5"/>`);

        const xt = 5, x0 = x[0], x1 = x[x.length - 1];
        for (let i = 0; i <= xt; i++) {
            const v = x0 + (i * (x1 - x0)) / xt;
            const px = ox + ((v - x0) / (x1 - x0)) * pw;
            out.push(`<line x1="${px.toFixed(1)}" y1="${m.top + ph}" x2="${px.toFixed(1)}" y2="${m.top + ph + 5}" stroke="${AXIS}"/>`);
            out.push(`<text x="${px.toFixed(1)}" y="${m.top + ph + 19}" fill="${AXIS}" text-anchor="middle">${fmt(v)}</text>`);
        }

        if (p === 0) {
            const yt = 4, y0 = y[0], y1 = y[y.length - 1];
            for (let i = 0; i <= yt; i++) {
                const v = y0 + (i * (y1 - y0)) / yt;
                const py = m.top + ph - ((v - y0) / (y1 - y0)) * ph;
                out.push(`<line x1="${ox - 5}" y1="${py.toFixed(1)}" x2="${ox}" y2="${py.toFixed(1)}" stroke="${AXIS}"/>`);
                out.push(`<text x="${ox - 9}" y="${(py + 4).toFixed(1)}" fill="${AXIS}" text-anchor="end">${fmt(v)}</text>`);
            }
            out.push(`<text x="18" y="${(m.top + ph / 2).toFixed(1)}" fill="${AXIS}" text-anchor="middle" transform="rotate(-90 18 ${(m.top + ph / 2).toFixed(1)})">${esc(yLabel)}</text>`);
        }
    });

    out.push(`<text x="${(m.left + (width - m.left - m.right) / 2).toFixed(1)}" y="${height - 14}" fill="${AXIS}" text-anchor="middle">${esc(xLabel)}</text>`);

    // Colour bar. Drawn as discrete bands from the same swatch table, so it
    // shows exactly the levels used in the panels.
    const bx = width - m.right + 22, bw = 15, bh = ph;
    const bands = 96;
    out.push('<g shape-rendering="crispEdges">');
    for (let i = 0; i < bands; i++) {
        const t = i / (bands - 1);
        const by = m.top + bh - (i + 1) * (bh / bands);
        out.push(`<rect x="${bx}" y="${by.toFixed(2)}" width="${bw}" height="${(bh / bands + 0.6).toFixed(2)}" fill="${swatch[Math.round(t * (levels - 1))]}"/>`);
    }
    out.push('</g>');
    out.push(`<rect x="${bx}" y="${m.top}" width="${bw}" height="${bh.toFixed(1)}" fill="none" stroke="${AXIS}" stroke-opacity="0.5"/>`);
    for (let i = 0; i <= 4; i++) {
        const v = zTickValues ? zTickValues[i] : zMin + (i * (zMax - zMin)) / 4;
        const py = m.top + bh - (i / 4) * bh;
        out.push(`<text x="${bx + bw + 6}" y="${(py + 4).toFixed(1)}" fill="${AXIS}">${fmtZ(v)}</text>`);
    }
    out.push(`<text x="${bx + bw / 2}" y="${m.top - 10}" fill="${AXIS}" text-anchor="middle">${esc(zLabel)}</text>`);

    out.push('</svg>');
    return out.join('\n');
}

const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
