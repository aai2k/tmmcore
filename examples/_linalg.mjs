/**
 * Small dense linear algebra for the examples: a linear solve and the
 * eigenvalues of a symmetric matrix. A matrix is an array of rows.
 */

// The row at or below row c with the largest entry in column c.
function pivotRow(a, c) {
    let p = c;
    for (let i = c + 1; i < a.length; i++) if (Math.abs(a[i][c]) > Math.abs(a[p][c])) p = i;
    return p;
}

// Forward elimination with partial pivoting on the augmented matrix a, in
// place. False when a pivot is below 1e-30, far below the entries of the
// systems the examples solve, so the matrix is singular.
function eliminate(a) {
    const n = a.length;
    for (let c = 0; c < n; c++) {
        const p = pivotRow(a, c);
        [a[c], a[p]] = [a[p], a[c]];
        if (Math.abs(a[c][c]) < 1e-30) return false;
        for (let i = c + 1; i < n; i++) {
            const f = a[i][c] / a[c][c];
            for (let j = c; j <= n; j++) a[i][j] -= f * a[c][j];
        }
    }
    return true;
}

function backSubstitute(a) {
    const n = a.length;
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let s = a[i][n];
        for (let j = i + 1; j < n; j++) s -= a[i][j] * x[j];
        x[i] = s / a[i][i];
    }
    return x;
}

// Solve A x = b by Gaussian elimination with partial pivoting; null if A is
// singular.
export function solve(A, b) {
    const a = A.map((row, i) => [...row, b[i]]);
    return eliminate(a) ? backSubstitute(a) : null;
}

// The squared sums of the diagonal and of the entries above it.
function squares(a) {
    let off = 0, diagonal = 0;
    for (let p = 0; p < a.length; p++) {
        diagonal += a[p][p] * a[p][p];
        for (let q = p + 1; q < a.length; q++) off += a[p][q] * a[p][q];
    }
    return { off, diagonal };
}

// The Jacobi rotation that zeroes a[p][q], applied to the columns and rows of
// a and to the columns of the eigenvector matrix v.
function rotate(a, v, p, q) {
    const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
    const t = (theta >= 0 ? 1 : -1) / (Math.abs(theta) + Math.hypot(theta, 1));
    const c = 1 / Math.hypot(t, 1), s = t * c;
    const turn = (x, y) => [c * x - s * y, s * x + c * y];
    for (const row of a) [row[p], row[q]] = turn(row[p], row[q]);
    for (let k = 0; k < a.length; k++) [a[p][k], a[q][k]] = turn(a[p][k], a[q][k]);
    for (const row of v) [row[p], row[q]] = turn(row[p], row[q]);
}

function sweep(a, v) {
    for (let p = 0; p < a.length; p++) {
        for (let q = p + 1; q < a.length; q++) if (a[p][q] !== 0) rotate(a, v, p, q);
    }
}

// Each eigenvector's largest component is made positive, so the signs do not
// depend on the order of the rotations.
function eigenPairs(a, v) {
    return a.map((row, i) => {
        const vector = v.map(r => r[i]);
        const big = vector.reduce((m, x) => (Math.abs(x) > Math.abs(m) ? x : m), 0);
        return { value: row[i], vector: vector.map(x => (big < 0 ? -x : x)) };
    });
}

// Eigenvalues and unit eigenvectors of a symmetric matrix by cyclic Jacobi
// rotations, largest eigenvalue first. The sweeps stop once the off-diagonal
// part is below rounding of the diagonal; they converge quadratically, so the
// cap of 50 sweeps is only a guard.
export function symmetricEigen(A) {
    const n = A.length;
    const a = A.map(row => [...row]);
    const v = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
    for (let round = 0; round < 50; round++) {
        const { off, diagonal } = squares(a);
        if (off <= Number.EPSILON ** 2 * diagonal) break;
        sweep(a, v);
    }
    return eigenPairs(a, v).sort((x, y) => y.value - x.value);
}
