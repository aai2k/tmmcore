/*
 * derivatives.c : exact thickness derivatives of R, T and A: the Jacobian,
 * the Hessian and the needle P-function, at one wavelength and over a grid.
 * Conventions are in tmm_internal.h.
 */

#include "tmm_internal.h"

/* ── Shared pieces of the derivative kernels ─────────────────────────────────
 * Port of the section of the same name in tmm.js. Every prefix
 * Pre[j] = M_0·…·M_{j−1} and suffix Post[j] = M_j·…·M_{N−1}·[1, ηs] is stored
 * as a mantissa times 2^e with the integer e kept alongside; a power of two
 * scales a double exactly, and every derivative, one prefix times one suffix
 * over den = η0 B + C, takes the exponent e_pre + e_post − e_0. */

static inline cx cscale(cx a, double s) { return cmk(a.re * s, a.im * s); }

/* Divides the values in place by a power of two when the largest part exceeds
 * BINARY_RESCALE_LIMIT, and returns that power's exponent. */
static int rescaleBinary(cx **values, int count) {
    double largest = 0.0;
    for (int i = 0; i < count; i++)
        largest = fmax(largest, fmax(fabs(values[i]->re), fabs(values[i]->im)));
    if (!(largest > BINARY_RESCALE_LIMIT) || isinf(largest)) return 0;
    int exponent = (int)floor(log2(largest));
    double factor = ldexp(1.0, -exponent);
    for (int i = 0; i < count; i++) {
        values[i]->re *= factor;
        values[i]->im *= factor;
    }
    return exponent;
}
static int rescaleMat(mat2 *M) {
    cx *e[4] = { &M->a, &M->b, &M->c, &M->d };
    return rescaleBinary(e, 4);
}
static int rescaleVec(vec2 *v) {
    cx *e[2] = { &v->x, &v->y };
    return rescaleBinary(e, 2);
}

/* The pair (P, S) = (Q/η, Qη) of a medium, at its critical-angle limit where
 * cosθ = 0. Mirrors generator in tmm.js. */
static inline generator_t generatorOf(cx nj, cx cosTheta, cx Q, double k0, int pol) {
    if (atCriticalAngle(cosTheta)) return criticalGenerator(nj, k0, pol);
    cx eta = layerAdmittance(nj, cosTheta, pol);
    generator_t g;
    g.P = cdiv(Q, eta);
    g.S = cmul(Q, eta);
    return g;
}

/* Everything the derivative kernels share at one (λ, θ, pol). `layers` is N
 * triples [n_re, n_im, d] used as given, with no d > 0 filter, for index parity
 * with the caller's design. A thickness that is not a number ≥ 0 leaves the
 * layer out: identity matrix, and Q = P = S = 0 so its derivative is zero. Q is
 * dδ/dd in full: past the bound the held matrix times e^{excess} moves with d
 * as the true one does, and Q with the held δ is its derivative. `logScale`
 * sums the excesses, which t carries. Mirrors derivativeLayer and
 * derivativeStack in tmm.js. */
typedef struct {
    int N;
    double k0;
    incidence_t angle;
    cx eta0;
    vec2 substrate;
    cx *n, *cosTheta, *delta, *eta, *Q;
    generator_t *G;
    int *critical;
    double *thickness, *excess;
    double logScale;
    mat2 *M, *Pre;
    vec2 *Post;
    int *preExp, *postExp;
} deriv_stack;

/* The arguments every derivative kernel shares. */
typedef struct {
    double lambda_nm, theta_deg;
    int pol;
    cx n0, ns;
    const double *layers;
    int N;
} deriv_request;

static void deriv_stack_free(deriv_stack *s) {
    free(s->n); free(s->cosTheta); free(s->delta); free(s->eta); free(s->Q);
    free(s->G); free(s->critical);
    free(s->thickness); free(s->excess); free(s->M); free(s->Pre); free(s->Post);
    free(s->preExp); free(s->postExp);
}

/* The matrix of `thickness` of layer k whose phase over it is `delta`. */
static inline mat2 partMatrix(const deriv_stack *s, int k, cx delta, double thickness) {
    return s->critical[k] ? criticalMatrix(s->G[k], thickness) : phaseMatrix(delta, s->eta[k]);
}

/* Returns 1, or 0 with nothing left allocated when the working state could not
 * be allocated. */
static int deriv_stack_build(const deriv_request *q, deriv_stack *out) {
    deriv_stack s;
    int N = q->N, pol = q->pol;
    double lambda_nm = q->lambda_nm;
    const double *layers = q->layers;
    cx n0 = q->n0;
    size_t count = (size_t)(N > 0 ? N : 1);
    s.N = N;
    s.k0 = (2.0 * PI) / lambda_nm;
    s.angle = incidence(q->theta_deg);
    s.eta0 = layerAdmittance(n0, incidentCosTheta(n0, s.angle), pol);
    s.substrate = substrateVector(q->ns, snellCosTheta(n0, s.angle, q->ns), pol);
    s.n         = (cx *)         malloc(sizeof(cx) * count);
    s.cosTheta  = (cx *)         malloc(sizeof(cx) * count);
    s.delta     = (cx *)         malloc(sizeof(cx) * count);
    s.eta       = (cx *)         malloc(sizeof(cx) * count);
    s.Q         = (cx *)         malloc(sizeof(cx) * count);
    s.G         = (generator_t *)malloc(sizeof(generator_t) * count);
    s.critical  = (int *)        malloc(sizeof(int) * count);
    s.thickness = (double *)     malloc(sizeof(double) * count);
    s.excess    = (double *)     malloc(sizeof(double) * count);
    s.M         = (mat2 *)       malloc(sizeof(mat2) * count);
    s.Pre       = (mat2 *)       malloc(sizeof(mat2) * ((size_t)N + 1));
    s.Post      = (vec2 *)       malloc(sizeof(vec2) * ((size_t)N + 1));
    s.preExp    = (int *)        malloc(sizeof(int) * ((size_t)N + 1));
    s.postExp   = (int *)        malloc(sizeof(int) * ((size_t)N + 1));
    void *const blocks[] = { s.n, s.cosTheta, s.delta, s.eta, s.Q, s.G, s.critical,
                             s.thickness, s.excess, s.M, s.Pre, s.Post, s.preExp, s.postExp };
    if (!all_allocated(blocks, (int)(sizeof blocks / sizeof blocks[0]))) {
        deriv_stack_free(&s);
        return 0;
    }

    s.logScale = 0.0;
    for (int k = 0; k < N; k++) {
        double d = layers[3 * k + 2];
        int present = d >= 0.0;
        if (!present) d = 0.0;
        s.n[k] = cmk(layers[3 * k + 0], layers[3 * k + 1]);
        s.thickness[k] = d;
        s.cosTheta[k] = snellCosTheta(n0, s.angle, s.n[k]);
        s.delta[k] = layerPhase(s.n[k], d, lambda_nm, s.cosTheta[k], &s.excess[k]);
        s.logScale += s.excess[k];
        s.eta[k] = layerAdmittance(s.n[k], s.cosTheta[k], pol);
        cx Q = cmul(cmul(s.n[k], cmk(s.k0, 0.0)), s.cosTheta[k]);   /* (2π/λ) n cosθ */
        s.Q[k] = present ? Q : cmk(0.0, 0.0);
        if (present) {
            s.G[k] = generatorOf(s.n[k], s.cosTheta[k], s.Q[k], s.k0, pol);
        } else {
            s.G[k].P = cmk(0.0, 0.0); s.G[k].S = cmk(0.0, 0.0);
        }
        s.critical[k] = atCriticalAngle(s.cosTheta[k]);
        s.M[k] = partMatrix(&s, k, s.delta[k], d);
    }

    s.Pre[0].a = cmk(1.0, 0.0); s.Pre[0].b = cmk(0.0, 0.0);
    s.Pre[0].c = cmk(0.0, 0.0); s.Pre[0].d = cmk(1.0, 0.0);
    s.preExp[0] = 0;
    for (int j = 0; j < N; j++) {
        s.Pre[j + 1] = matmul(s.Pre[j], s.M[j]);
        s.preExp[j + 1] = s.preExp[j] + rescaleMat(&s.Pre[j + 1]);
    }
    s.Post[N] = s.substrate;
    s.postExp[N] = 0;
    for (int j = N - 1; j >= 0; j--) {
        s.Post[j] = cmatvec(s.M[j], s.Post[j + 1]);
        s.postExp[j] = s.postExp[j + 1] + rescaleVec(&s.Post[j]);
    }
    *out = s;
    return 1;
}

/* R, T, A from [B, C] = Post[0] (Macleod Eqs. 2.123–2.125), with r, the true
 * t, b = B/den and c = C/den, and g = 2 Im η0 / Re η0, the weight of Im r in
 * the absorptance. Only t carries the exponent of [B, C] and the log scale of
 * the held layers. */
typedef struct {
    double R, T, A, Tfac, g;
    cx eta0, r, t, b, c, inverseDen;
    int exponent;
} stack_response;

static stack_response stack_response_of(const deriv_stack *s) {
    stack_response o;
    cx B = s->Post[0].x, C = s->Post[0].y;
    cx den = cadd(cmul(s->eta0, B), C);
    o.inverseDen = cdiv(cmk(1.0, 0.0), den);
    o.r = cdiv(csub(cmul(s->eta0, B), C), den);
    o.t = cscale(cmul(cdiv(cmul(cmk(2.0, 0.0), s->eta0), den), s->substrate.x),
                 ldexp(1.0, -s->postExp[0]));
    if (s->logScale > 0.0) o.t = cscale(o.t, exp(-s->logScale));
    o.Tfac = transmittedFlux(s->substrate) / s->eta0.re;
    o.R = cabs2(o.r);
    o.T = o.Tfac * cabs2(o.t); if (o.T < 0.0) o.T = 0.0;
    o.A = layerAbsorptance(o.R, o.T, s->eta0, o.r); if (o.A < 0.0) o.A = 0.0;
    o.g = 2.0 * (s->eta0.im / s->eta0.re);
    o.eta0 = s->eta0;
    o.b = cmul(B, o.inverseDen);
    o.c = cmul(C, o.inverseDen);
    o.exponent = s->postExp[0];
    return o;
}

/* A derivative of [B, C], a mantissa and its exponent, divided by den. */
static inline vec2 per_den(const stack_response *s, vec2 v, int exponent) {
    double scale = ldexp(1.0, exponent - s->exponent);
    vec2 o;
    o.x = cscale(cmul(v.x, s->inverseDen), scale);
    o.y = cscale(cmul(v.y, s->inverseDen), scale);
    return o;
}

/* {dR,dT,dA} from (u, w) = d[B,C]/den: dr = 2η0 (c·u − b·w), dt = −t (η0·u + w),
 * the chain rule with f = 2η0/den² and den divided out once; dA = −(dR + dT)
 * + g Im dr, as layerAbsorptance. */
static void response_derivative(const stack_response *s, vec2 uw, double *o) {
    cx dr = cmul(cmul(cmk(2.0, 0.0), s->eta0), csub(cmul(s->c, uw.x), cmul(s->b, uw.y)));
    double dR = 2.0 * (cmul(cconj(s->r), dr)).re;
    cx dt = cmul(cmk(-1.0, 0.0), cmul(s->t, cadd(cmul(s->eta0, uw.x), uw.y)));
    double dT = s->Tfac * 2.0 * (cmul(cconj(s->t), dt)).re;
    o[0] = dR; o[1] = dT; o[2] = -(dR + dT) + s->g * dr.im;
}

/* dM/dd = [[ −Q sinδ, −iP cosδ ], [ −iS cosδ, −Q sinδ ]], P = Q/η, S = Qη */
static inline mat2 layerMatrixDerivative(cx delta, cx Q, generator_t G) {
    cx cD = ccos_(delta), sD = csin_(delta);
    cx negI = cmk(0.0, -1.0), neg1 = cmk(-1.0, 0.0);
    mat2 dM;
    dM.a = cmul(neg1, cmul(Q, sD));
    dM.b = cmul(negI, cmul(G.P, cD));
    dM.c = cmul(negI, cmul(G.S, cD));
    dM.d = cmul(neg1, cmul(Q, sD));
    return dM;
}

/* ── Exported: analytic thickness Jacobian for one (λ, θ, pol) ────────────────
 * Port of the thickness Jacobian in tmm.js. Returns the exact analytic
 * dR/dd_k, dT/dd_k, dA/dd_k for every layer at one sample.
 *
 * `layers` is N triples [n_re, n_im, d] used as given (see deriv_stack).
 * Outputs (each length N): dRdd, dTdd, dAdd. Also writes base R,T,A to base[0..2]. */

static void jacobian_at(const deriv_request *q,
                        double *dRdd, double *dTdd, double *dAdd, double *base) {
    int N = q->N;
    deriv_stack s;
    if (!deriv_stack_build(q, &s)) {
        fill_nan(base, 3); fill_nan(dRdd, N); fill_nan(dTdd, N); fill_nan(dAdd, N);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    for (int k = 0; k < N; k++) {
        mat2 dMk = layerMatrixDerivative(s.delta[k], s.Q[k], s.G[k]);
        vec2 dV = cmatvec(s.Pre[k], cmatvec(dMk, s.Post[k + 1]));
        double o[3];
        response_derivative(&resp, per_den(&resp, dV, s.preExp[k] + s.postExp[k + 1]), o);
        dRdd[k] = o[0]; dTdd[k] = o[1]; dAdd[k] = o[2];
    }

    deriv_stack_free(&s);
}

TMM_EXPORT("tmm_jacobian")
void tmm_jacobian(double lambda_nm, double theta_deg, int pol,
                  double n0_re, double n0_im, double ns_re, double ns_im,
                  const double *layers, int N,
                  double *dRdd, double *dTdd, double *dAdd, double *base) {
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    jacobian_at(&q, dRdd, dTdd, dAdd, base);
}

/* The request of a grid derivative call at wavelength index `li`; `layers` is
 * the caller's N-triple scratch, filled here. */
static deriv_request grid_request(const grid_layout *g, int pol, int li, double *layers) {
    layers_at(g, li, layers);
    deriv_request q = { g->lambdas[li], theta_at(g->theta_deg, g->nTheta, li), pol,
                        { g->n0arr[2 * li + 0], g->n0arr[2 * li + 1] },
                        { g->nsarr[2 * li + 0], g->nsarr[2 * li + 1] }, layers, g->N };
    return q;
}

/* ── Exported: thickness Jacobian over a wavelength grid ─────────────────────
 * tmm_jacobian at every wavelength of a grid in one call, for one
 * polarization. Inputs are laid out as for tmm_rta_spectrum.
 *
 *   dRdd/dTdd/dAdd : nLam × N each, [λ][layer]
 *   base           : nLam × 3, [λ][R, T, A] */

TMM_EXPORT("tmm_jacobian_spectrum")
void tmm_jacobian_spectrum(const double *lambdas, int nLam,
                           const double *theta_deg, int nTheta, int pol,
                           const double *n0arr, const double *nsarr,
                           const double *matNK, const double *thick, int N,
                           double *dRdd, double *dTdd, double *dAdd, double *base) {
    const grid_layout g = { lambdas, nLam, theta_deg, nTheta, n0arr, nsarr, matNK, thick, N };
    double *layers = (double *)malloc(sizeof(double) * 3 * (N > 0 ? N : 1));
    if (!layers) {
        fill_nan(base, 3L * nLam); fill_nan(dRdd, (long)N * nLam);
        fill_nan(dTdd, (long)N * nLam); fill_nan(dAdd, (long)N * nLam);
        return;
    }
    for (int li = 0; li < nLam; li++) {
        deriv_request q = grid_request(&g, pol, li, layers);
        long at = (long)li * N;
        jacobian_at(&q, &dRdd[at], &dTdd[at], &dAdd[at], &base[3L * li]);
    }
    free(layers);
}

/* ── Analytic needle P-function scan ─────────────────────────────────────────
 * Port of the needle scan in tmm.js : the d→0 limit of
 * Sullivan's pre/post method (Tikhonravov's analytic P-function). Returns the
 * merit-gradient ingredients {dR,dT,dA} of inserting an infinitesimal needle of
 * each candidate index at every gap position (0..N) and, optionally, at intra-
 * layer split fractions. The {dR,dT,dA} chain rule is identical to tmm_jacobian.
 *
 * `layers` used AS-IS (no d>0 filter), index parity with the caller's stack.
 * Outputs (caller-owned f64):
 *   base  : 3                      [R,T,A]
 *   gaps  : (N+1)*nCand*3          layout [pos][cand][dR,dT,dA]
 *   intra : N*nFrac*nCand*3        layout [layer][frac][cand][dR,dT,dA] (nFrac>0)
 */

/* The needles a scan tries: nCand candidate indices as [re, im] pairs, and the
 * nFrac split fractions at which every layer is opened. */
typedef struct {
    const double *candNs; int nCand;
    const double *fracs; int nFrac;
} needle_set;

static void needle_scan_at(const deriv_request *q, const needle_set *needles,
                           double *base, double *gaps, double *intra) {
    const double *candNs = needles->candNs, *fracs = needles->fracs;
    int nCand = needles->nCand, nFrac = needles->nFrac;
    int N = q->N, pol = q->pol;
    double lambda_nm = q->lambda_nm;
    cx n0 = q->n0;
    deriv_stack s;
    mat2 *Ac = (mat2 *)malloc(sizeof(mat2) * (nCand > 0 ? nCand : 1));
    if (!Ac || !deriv_stack_build(q, &s)) {
        free(Ac);
        fill_nan(base, 3);
        fill_nan(gaps, (long)(N + 1) * nCand * 3);
        if (nFrac > 0) fill_nan(intra, (long)N * nFrac * nCand * 3);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    /* Needle derivative-matrix A = G per candidate: [[0,−iQ/ηₐ],[−iηₐQ,0]]. */
    for (int c = 0; c < nCand; c++) {
        cx nA = cmk(candNs[2 * c + 0], candNs[2 * c + 1]);
        cx cthA = snellCosTheta(n0, s.angle, nA);
        cx Q = cmul(cmul(nA, cmk(s.k0, 0.0)), cthA);
        Ac[c] = generatorMatrix(generatorOf(nA, cthA, Q, s.k0, pol));
    }

    /* Gaps: pos = 0..N, every candidate. */
    for (int pos = 0; pos <= N; pos++) {
        for (int c = 0; c < nCand; c++) {
            vec2 dV = cmatvec(s.Pre[pos], cmatvec(Ac[c], s.Post[pos]));
            response_derivative(&resp, per_den(&resp, dV, s.preExp[pos] + s.postExp[pos]),
                                &gaps[((long)pos * nCand + c) * 3]);
        }
    }

    /* Intra-layer positions. Short of the bound the host splits into halves
     * whose product is its own matrix: the front half keeps its own phase and
     * the back half takes the rest. Past it each half is held on its own, and
     * the pair takes e^{front excess + back excess − host excess}, between
     * e^{−50} and 1, the ratio of what the halves and the host leave out. */
    if (nFrac > 0) {
        for (int k = 0; k < N; k++) {
            for (int fi = 0; fi < nFrac; fi++) {
                double frontExcess, backExcess, factor = 1.0;
                cx front = layerPhase(s.n[k], fracs[fi] * s.thickness[k], lambda_nm,
                                      s.cosTheta[k], &frontExcess);
                cx back = csub(s.delta[k], front);
                if (s.excess[k] > 0.0) {
                    back = layerPhase(s.n[k], (1.0 - fracs[fi]) * s.thickness[k], lambda_nm,
                                      s.cosTheta[k], &backExcess);
                    factor = exp(frontExcess + backExcess - s.excess[k]);
                }
                mat2 preIn  = matmul(s.Pre[k], partMatrix(&s, k, front, fracs[fi] * s.thickness[k]));
                vec2 postIn = cmatvec(partMatrix(&s, k, back, (1.0 - fracs[fi]) * s.thickness[k]),
                                      s.Post[k + 1]);
                for (int c = 0; c < nCand; c++) {
                    vec2 dV = cmatvec(preIn, cmatvec(Ac[c], postIn));
                    if (factor != 1.0) { dV.x = cscale(dV.x, factor); dV.y = cscale(dV.y, factor); }
                    long off = (((long)k * nFrac + fi) * nCand + c) * 3;
                    response_derivative(&resp, per_den(&resp, dV, s.preExp[k] + s.postExp[k + 1]),
                                        &intra[off]);
                }
            }
        }
    }

    free(Ac);
    deriv_stack_free(&s);
}

TMM_EXPORT("tmm_needle_scan")
void tmm_needle_scan(double lambda_nm, double theta_deg, int pol,
                     double n0_re, double n0_im, double ns_re, double ns_im,
                     const double *layers, int N,
                     const double *candNs, int nCand,
                     const double *fracs, int nFrac,
                     double *base, double *gaps, double *intra) {
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    const needle_set needles = { candNs, nCand, fracs, nFrac };
    needle_scan_at(&q, &needles, base, gaps, intra);
}

/* ── Exported: needle scan over a wavelength grid ────────────────────────────
 * tmm_needle_scan at every wavelength of a grid in one call, for one
 * polarization. Inputs are laid out as for tmm_rta_spectrum; the candidates
 * are dispersive like the layers:
 *   candNK : nCand × nLam × 2, [candidate][λ][re, im]
 *   fracs  : nFrac split fractions, the same at every wavelength
 * Outputs (caller-owned f64), each block of tmm_needle_scan once per λ:
 *   base  : nLam × 3
 *   gaps  : nLam × (N+1) × nCand × 3
 *   intra : nLam × N × nFrac × nCand × 3   (nFrac > 0) */

TMM_EXPORT("tmm_needle_scan_spectrum")
void tmm_needle_scan_spectrum(const double *lambdas, int nLam,
                              const double *theta_deg, int nTheta, int pol,
                              const double *n0arr, const double *nsarr,
                              const double *matNK, const double *thick, int N,
                              const double *candNK, int nCand,
                              const double *fracs, int nFrac,
                              double *base, double *gaps, double *intra) {
    const grid_layout g = { lambdas, nLam, theta_deg, nTheta, n0arr, nsarr, matNK, thick, N };
    long gapBlock = (long)(N + 1) * nCand * 3;
    long intraBlock = (long)N * nFrac * nCand * 3;
    double *layers = (double *)malloc(sizeof(double) * 3 * (N > 0 ? N : 1));
    double *cand = (double *)malloc(sizeof(double) * 2 * (nCand > 0 ? nCand : 1));
    if (!layers || !cand) {
        free(layers); free(cand);
        fill_nan(base, 3L * nLam);
        fill_nan(gaps, gapBlock * nLam);
        if (nFrac > 0) fill_nan(intra, intraBlock * nLam);
        return;
    }
    const needle_set needles = { cand, nCand, fracs, nFrac };
    for (int li = 0; li < nLam; li++) {
        deriv_request q = grid_request(&g, pol, li, layers);
        for (int c = 0; c < nCand; c++) {
            long at = ((long)c * nLam + li) * 2;
            cand[2 * c + 0] = candNK[at + 0];
            cand[2 * c + 1] = candNK[at + 1];
        }
        needle_scan_at(&q, &needles, &base[3L * li], &gaps[gapBlock * li],
                       nFrac > 0 ? &intra[intraBlock * li] : intra);
    }
    free(layers); free(cand);
}

/* d²M/dd² = Q²·[[ −cosδ, i sinδ/η ], [ i η sinδ, −cosδ ]], with Q²/η = QP and
 * Q²η = QS */
static inline mat2 layerMatrixSecondDerivative(cx delta, cx Q, generator_t G) {
    cx cD = ccos_(delta), sD = csin_(delta);
    cx Q2 = cmul(Q, Q);
    cx posI = cmk(0.0, 1.0), neg1 = cmk(-1.0, 0.0);
    mat2 d2M;
    d2M.a = cmul(neg1, cmul(Q2, cD));
    d2M.b = cmul(posI, cmul(cmul(Q, G.P), sD));
    d2M.c = cmul(posI, cmul(cmul(Q, G.S), sD));
    d2M.d = cmul(neg1, cmul(Q2, cD));
    return d2M;
}

/* d²R, d²T and d²A from (uᵢ,wᵢ), (uⱼ,wⱼ) and the mixed partial (U,W), each a
 * derivative of [B,C] divided by den. With gₖ = η0 uₖ + wₖ and G = η0 U + W:
 *   d²r = 2η0 (wᵢuⱼ + cU − uᵢwⱼ − bW) − 2 drⱼ gᵢ,   d²t = −t G + 2t gᵢgⱼ
 * Port of responseSecondDerivative in tmm.js. */
typedef struct { double R, T, A; } second_derivative;

static second_derivative response_second_derivative(const stack_response *s,
                                                    vec2 fi, vec2 fj, vec2 second) {
    second_derivative o;
    cx twoEta0 = cmul(cmk(2.0, 0.0), s->eta0);
    cx neg1 = cmk(-1.0, 0.0);
    cx dr_i = cmul(twoEta0, csub(cmul(s->c, fi.x), cmul(s->b, fi.y)));
    cx dr_j = cmul(twoEta0, csub(cmul(s->c, fj.x), cmul(s->b, fj.y)));
    cx g_i = cadd(cmul(s->eta0, fi.x), fi.y);
    cx g_j = cadd(cmul(s->eta0, fj.x), fj.y);
    cx innerR = csub(cadd(cmul(fi.y, fj.x), cmul(s->c, second.x)),
                     cadd(cmul(fi.x, fj.y), cmul(s->b, second.y)));
    cx d2r = csub(cmul(twoEta0, innerR), cmul(cmul(cmk(2.0, 0.0), dr_j), g_i));
    o.R = 2.0 * ((cmul(cconj(dr_i), dr_j)).re + (cmul(cconj(s->r), d2r)).re);
    cx dt_i = cmul(neg1, cmul(s->t, g_i));
    cx dt_j = cmul(neg1, cmul(s->t, g_j));
    cx G = cadd(cmul(s->eta0, second.x), second.y);
    cx d2t = cadd(cmul(neg1, cmul(s->t, G)), cmul(cmul(cmk(2.0, 0.0), s->t), cmul(g_i, g_j)));
    o.T = s->Tfac * 2.0 * ((cmul(cconj(dt_i), dt_j)).re + (cmul(cconj(s->t), d2t)).re);
    o.A = -(o.R + o.T) + s->g * d2r.im;
    return o;
}

/* ── Analytic thickness-Hessian kernel ───────────
 * Port of the thickness Hessian in tmm.js : the EXACT analytic second
 * derivatives ∂²{R,T,A}/∂dᵢ∂dⱼ (full N×N symmetric) plus the first
 * derivatives, at one (λ,θ,pol). The middle product M_{i+1}···M_{j-1} carries
 * its own binary exponent, like the prefixes and suffixes.
 *
 * Outputs (caller-owned f64):
 *   base   : 3        [R,T,A]
 *   dRdd/dTdd/dAdd : N each      first derivatives (verbatim tmm_jacobian)
 *   d2Rdd/d2Tdd/d2Add : N*N each row-major, FULL symmetric (both triangles set)
 */
typedef struct {
    double *dRdd, *dTdd, *dAdd;
    double *d2Rdd, *d2Tdd, *d2Add;
    double *base;
} hessian_out;

static void hessian_at(const deriv_request *q, const hessian_out *out) {
    double *dRdd = out->dRdd, *dTdd = out->dTdd, *dAdd = out->dAdd;
    double *d2Rdd = out->d2Rdd, *d2Tdd = out->d2Tdd, *d2Add = out->d2Add;
    double *base = out->base;
    int N = q->N;
    deriv_stack s;
    /* Per-layer first-derivative pieces: dM[k], v[k] = dMₖ·Post[k+1] (exponent
     * postExp[k+1]) and d[B,C]/den. */
    mat2 *dM    = (mat2 *)malloc(sizeof(mat2) * (N > 0 ? N : 1));
    vec2 *v     = (vec2 *)malloc(sizeof(vec2) * (N > 0 ? N : 1));
    vec2 *first = (vec2 *)malloc(sizeof(vec2) * (N > 0 ? N : 1));
    if (!dM || !v || !first || !deriv_stack_build(q, &s)) {
        free(dM); free(v); free(first);
        fill_nan(base, 3); fill_nan(dRdd, N); fill_nan(dTdd, N); fill_nan(dAdd, N);
        fill_nan(d2Rdd, (long)N * N); fill_nan(d2Tdd, (long)N * N); fill_nan(d2Add, (long)N * N);
        return;
    }
    stack_response resp = stack_response_of(&s);
    base[0] = resp.R; base[1] = resp.T; base[2] = resp.A;

    for (int k = 0; k < N; k++) {
        dM[k] = layerMatrixDerivative(s.delta[k], s.Q[k], s.G[k]);
        v[k] = cmatvec(dM[k], s.Post[k + 1]);
        first[k] = per_den(&resp, cmatvec(s.Pre[k], v[k]), s.preExp[k] + s.postExp[k + 1]);
        double o[3];
        response_derivative(&resp, first[k], o);
        dRdd[k] = o[0]; dTdd[k] = o[1]; dAdd[k] = o[2];
    }

    mat2 I; I.a = cmk(1.0, 0.0); I.b = cmk(0.0, 0.0); I.c = cmk(0.0, 0.0); I.d = cmk(1.0, 0.0);
    for (int i = 0; i < N; i++) {
        mat2 Wmat_i = matmul(s.Pre[i], dM[i]); /* Pre[i]·dMᵢ (used for j>i) */
        mat2 Cmid = I;                          /* M_{i+1}···M_{j-1}, empty at j=i+1 */
        int midExp = 0;
        for (int j = i; j < N; j++) {
            vec2 second;
            if (j == i) {
                mat2 d2M = layerMatrixSecondDerivative(s.delta[i], s.Q[i], s.G[i]);
                second = per_den(&resp, cmatvec(s.Pre[i], cmatvec(d2M, s.Post[i + 1])),
                                 s.preExp[i] + s.postExp[i + 1]);
            } else {
                second = per_den(&resp, cmatvec(Wmat_i, cmatvec(Cmid, v[j])),
                                 s.preExp[i] + midExp + s.postExp[j + 1]);
            }
            second_derivative h = response_second_derivative(&resp, first[i], first[j], second);
            double d2Rij = h.R, d2Tij = h.T, d2Aij = h.A;
            d2Rdd[(long)i * N + j] = d2Rdd[(long)j * N + i] = d2Rij;
            d2Tdd[(long)i * N + j] = d2Tdd[(long)j * N + i] = d2Tij;
            d2Add[(long)i * N + j] = d2Add[(long)j * N + i] = d2Aij;
            if (j >= i + 1) {                   /* advance middle: include M_j */
                Cmid = matmul(Cmid, s.M[j]);
                midExp += rescaleMat(&Cmid);
            }
        }
    }

    free(dM); free(v); free(first);
    deriv_stack_free(&s);
}

TMM_EXPORT("tmm_hessian")
void tmm_hessian(double lambda_nm, double theta_deg, int pol,
                 double n0_re, double n0_im, double ns_re, double ns_im,
                 const double *layers, int N,
                 double *dRdd, double *dTdd, double *dAdd,
                 double *d2Rdd, double *d2Tdd, double *d2Add, double *base) {
    deriv_request q = { lambda_nm, theta_deg, pol, { n0_re, n0_im }, { ns_re, ns_im }, layers, N };
    const hessian_out out = { dRdd, dTdd, dAdd, d2Rdd, d2Tdd, d2Add, base };
    hessian_at(&q, &out);
}

/* ── Exported: thickness Hessian over a wavelength grid ──────────────────────
 * tmm_hessian at every wavelength of a grid in one call, for one
 * polarization. Inputs are laid out as for tmm_rta_spectrum.
 *
 *   dRdd/dTdd/dAdd    : nLam × N each, [λ][layer]
 *   d2Rdd/d2Tdd/d2Add : nLam × N × N each, [λ][row][column], full symmetric
 *   base              : nLam × 3, [λ][R, T, A] */

TMM_EXPORT("tmm_hessian_spectrum")
void tmm_hessian_spectrum(const double *lambdas, int nLam,
                          const double *theta_deg, int nTheta, int pol,
                          const double *n0arr, const double *nsarr,
                          const double *matNK, const double *thick, int N,
                          double *dRdd, double *dTdd, double *dAdd,
                          double *d2Rdd, double *d2Tdd, double *d2Add, double *base) {
    const grid_layout g = { lambdas, nLam, theta_deg, nTheta, n0arr, nsarr, matNK, thick, N };
    long square = (long)N * N;
    double *layers = (double *)malloc(sizeof(double) * 3 * (N > 0 ? N : 1));
    if (!layers) {
        fill_nan(base, 3L * nLam);
        fill_nan(dRdd, (long)N * nLam); fill_nan(dTdd, (long)N * nLam); fill_nan(dAdd, (long)N * nLam);
        fill_nan(d2Rdd, square * nLam); fill_nan(d2Tdd, square * nLam); fill_nan(d2Add, square * nLam);
        return;
    }
    for (int li = 0; li < nLam; li++) {
        deriv_request q = grid_request(&g, pol, li, layers);
        long at = (long)li * N, at2 = square * li;
        const hessian_out out = { &dRdd[at], &dTdd[at], &dAdd[at],
                                  &d2Rdd[at2], &d2Tdd[at2], &d2Add[at2], &base[3L * li] };
        hessian_at(&q, &out);
    }
    free(layers);
}
