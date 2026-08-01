/**
 * Minimal SVG line plot — no dependencies.
 *
 * Every figure in docs/ is written by an example script through this helper, so
 * the images cannot drift from the code that produced them. Colours are chosen
 * to stay legible on both light and dark backgrounds, since GitHub renders
 * markdown in either and strips <style> blocks from inline SVG.
 */

const COLORS = ['#2f81f7', '#d29922', '#3fb950', '#f85149', '#a371f7'];
const AXIS = '#8b949e';

export function linePlot({
    series, xLabel = '', yLabel = '', title = '',
    width = 720, height = 380, yMin = null, yMax = null,
}) {
    const m = { top: title ? 44 : 24, right: 20, bottom: 52, left: 64 };
    const w = width - m.left - m.right;
    const h = height - m.top - m.bottom;

    const xs = series.flatMap(s => s.x);
    const ys = series.flatMap(s => s.y);
    const x0 = Math.min(...xs), x1 = Math.max(...xs);
    const y0 = yMin ?? Math.min(...ys);
    const y1 = yMax ?? Math.max(...ys);
    const pad = (y1 - y0) * 0.05 || 0.05;
    const ylo = yMin ?? y0 - pad, yhi = yMax ?? y1 + pad;

    const px = v => m.left + ((v - x0) / (x1 - x0)) * w;
    const py = v => m.top + h - ((v - ylo) / (yhi - ylo)) * h;

    const ticks = (lo, hi, n = 5) =>
        Array.from({ length: n + 1 }, (_, i) => lo + (i * (hi - lo)) / n);
    const fmt = v => (Math.abs(v) >= 1000 || (Math.abs(v) < 0.01 && v !== 0)
        ? v.toExponential(1)
        : String(Math.round(v * 1000) / 1000));

    const out = [];
    out.push(`<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" font-family="system-ui, sans-serif" font-size="12">`);
    if (title) out.push(`<text x="${width / 2}" y="24" fill="${AXIS}" font-size="14" text-anchor="middle">${esc(title)}</text>`);

    for (const t of ticks(ylo, yhi)) {
        out.push(`<line x1="${m.left}" y1="${py(t).toFixed(1)}" x2="${m.left + w}" y2="${py(t).toFixed(1)}" stroke="${AXIS}" stroke-opacity="0.18"/>`);
        out.push(`<text x="${m.left - 8}" y="${(py(t) + 4).toFixed(1)}" fill="${AXIS}" text-anchor="end">${fmt(t)}</text>`);
    }
    for (const t of ticks(x0, x1)) {
        out.push(`<text x="${px(t).toFixed(1)}" y="${m.top + h + 20}" fill="${AXIS}" text-anchor="middle">${fmt(t)}</text>`);
    }

    out.push(`<line x1="${m.left}" y1="${m.top + h}" x2="${m.left + w}" y2="${m.top + h}" stroke="${AXIS}"/>`);
    out.push(`<line x1="${m.left}" y1="${m.top}" x2="${m.left}" y2="${m.top + h}" stroke="${AXIS}"/>`);
    out.push(`<text x="${m.left + w / 2}" y="${height - 12}" fill="${AXIS}" text-anchor="middle">${esc(xLabel)}</text>`);
    out.push(`<text x="16" y="${m.top + h / 2}" fill="${AXIS}" text-anchor="middle" transform="rotate(-90 16 ${m.top + h / 2})">${esc(yLabel)}</text>`);

    series.forEach((s, i) => {
        const c = s.color ?? COLORS[i % COLORS.length];
        const d = s.x.map((xv, j) => `${j ? 'L' : 'M'}${px(xv).toFixed(2)},${py(s.y[j]).toFixed(2)}`).join('');
        out.push(`<path d="${d}" fill="none" stroke="${c}" stroke-width="2" stroke-linejoin="round"${s.dash ? ` stroke-dasharray="${s.dash}"` : ''}/>`);
        if (s.label) {
            const ly = m.top + 14 + i * 18;
            out.push(`<line x1="${m.left + w - 150}" y1="${ly}" x2="${m.left + w - 126}" y2="${ly}" stroke="${c}" stroke-width="2"${s.dash ? ` stroke-dasharray="${s.dash}"` : ''}/>`);
            out.push(`<text x="${m.left + w - 120}" y="${ly + 4}" fill="${AXIS}">${esc(s.label)}</text>`);
        }
    });

    out.push('</svg>');
    return out.join('\n');
}

const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
