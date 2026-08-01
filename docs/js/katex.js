// Renders the math that pymdownx.arithmatex emits, on first load and again
// after each instant-navigation page change.
//
// The delimiters need doubled backslashes: "\\(" is the two-character string
// \( , whereas "\(" collapses to a bare "(" and would treat every parenthesis
// in the page as the start of an equation.
//
// arithmatex in generic mode always writes \(…\) and \[…\], never $…$, so only
// those are listed. Including $ would also catch literal dollar signs in prose.
document$.subscribe(() => {
    renderMathInElement(document.body, {
        delimiters: [
            { left: "\\[", right: "\\]", display: true },
            { left: "\\(", right: "\\)", display: false },
        ],
        throwOnError: false,
    });
});
