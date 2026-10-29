/*
 * jet.c : the Taylor-jet arithmetic as exports, each applied to a grid of
 * jets in one call. Conventions are in tmm_internal.h.
 */

#include "tmm_internal.h"

/* ── Exported: Taylor-jet arithmetic over a grid ─────────────────────────────
 * The arithmetic of src/js/jet.js, one operation applied to `count` jets per
 * call, so a caller in another language builds the index jet of a dispersion
 * formula over a whole wavelength grid in a handful of calls. A jet is 8
 * doubles, [re, im] for each of [f, f', f''/2!, f'''/3!]; a grid of jets is
 * count × 8. A real argument that goes with each jet (a factor, an exponent, a
 * bound) is given once per jet. An output may be the same buffer as an input
 * of its shape. */

static inline void jwrite(double *p, jet a) {
    for (int i = 0; i < JET_N; i++) { p[2 * i] = a.c[i].re; p[2 * i + 1] = a.c[i].im; }
}

/* exp(f): e^{f0} at order 0, then (e^f)' = f'·e^f order by order. */
static jet jexp_j(jet a) {
    jet o;
    double magnitude = exp(a.c[0].re);
    o.c[0] = cmk(magnitude * cos(a.c[0].im), magnitude * sin(a.c[0].im));
    for (int order = 1; order < JET_N; order++) {
        cx sum = cmk(0.0, 0.0);
        for (int i = 1; i <= order; i++) {
            cx term = cmul(a.c[i], o.c[order - i]);
            sum = cadd(sum, cmk(term.re * i, term.im * i));
        }
        double inv = 1.0 / (double)order;
        o.c[order] = cmk(sum.re * inv, sum.im * inv);
    }
    return o;
}

/* Principal log: ln|f0| + i arg f0 at order 0, then (ln f)' = f'/f. */
static jet jlog_j(jet a) {
    jet o, derivative;
    o.c[0] = cmk(log(hypot(a.c[0].re, a.c[0].im)), atan2(a.c[0].im, a.c[0].re));
    derivative.c[0] = a.c[1];
    derivative.c[1] = cmk(a.c[2].re * 2.0, a.c[2].im * 2.0);
    derivative.c[2] = cmk(a.c[3].re * 3.0, a.c[3].im * 3.0);
    derivative.c[3] = cmk(0.0, 0.0);
    jet quotient = jmul(derivative, jrecip(a));
    for (int order = 1; order < JET_N; order++) {
        double inv = 1.0 / (double)order;
        o.c[order] = cmk(quotient.c[order - 1].re * inv, quotient.c[order - 1].im * inv);
    }
    return o;
}

/* f(x) from f(x0) and [f'(x0), f''(x0), f'''(x0)]: the Taylor series of f in
 * the displacement x − x0, which is the jet x with its value set to zero. */
static jet jcompose_j(double value, const double *derivatives, jet x) {
    jet displacement = x;
    displacement.c[0] = cmk(0.0, 0.0);
    jet square = jmul(displacement, displacement);
    jet first = jscale(displacement, derivatives[0]);
    jet second = jscale(square, derivatives[1] / 2.0);
    jet third = jscale(jmul(square, displacement), derivatives[2] / 6.0);
    return jadd(jadd(jconst(value, 0.0), first), jadd(second, third));
}

TMM_EXPORT("tmm_jet_constant")
void tmm_jet_constant(const double *value, const double *imaginary, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jconst(value[i], imaginary[i]));
}

/* `derivatives` is count × 8: [f, f', f'', f'''] as [re, im] pairs. */
TMM_EXPORT("tmm_jet_from_derivatives")
void tmm_jet_from_derivatives(const double *derivatives, int count, double *out) {
    for (int i = 0; i < count; i++) {
        const double *d = &derivatives[8L * i];
        jet o;
        o.c[0] = cmk(d[0], d[1]);
        o.c[1] = cmk(d[2], d[3]);
        o.c[2] = cmk(d[4] * (1.0 / 2.0), d[5] * (1.0 / 2.0));
        o.c[3] = cmk(d[6] * (1.0 / 6.0), d[7] * (1.0 / 6.0));
        jwrite(&out[8L * i], o);
    }
}

/* out is count × 8: [f, f', f'', f'''] as [re, im] pairs. */
TMM_EXPORT("tmm_jet_derivatives")
void tmm_jet_derivatives(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) {
        cx d[JET_N];
        jderivs(jread(&jets[8L * i]), d);
        jet o;
        for (int k = 0; k < JET_N; k++) o.c[k] = d[k];
        jwrite(&out[8L * i], o);
    }
}

TMM_EXPORT("tmm_jet_add")
void tmm_jet_add(const double *left, const double *right, int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jadd(jread(&left[8L * i]), jread(&right[8L * i])));
}

TMM_EXPORT("tmm_jet_subtract")
void tmm_jet_subtract(const double *left, const double *right, int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jsub(jread(&left[8L * i]), jread(&right[8L * i])));
}

TMM_EXPORT("tmm_jet_scale")
void tmm_jet_scale(const double *jets, const double *scalar, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jscale(jread(&jets[8L * i]), scalar[i]));
}

TMM_EXPORT("tmm_jet_multiply")
void tmm_jet_multiply(const double *left, const double *right, int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jmul(jread(&left[8L * i]), jread(&right[8L * i])));
}

TMM_EXPORT("tmm_jet_reciprocal")
void tmm_jet_reciprocal(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jrecip(jread(&jets[8L * i])));
}

TMM_EXPORT("tmm_jet_divide")
void tmm_jet_divide(const double *numerator, const double *denominator, int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jdiv(jread(&numerator[8L * i]), jread(&denominator[8L * i])));
}

/* The order-0 root is csqrt_, csqrt of primitives.js; jet.js forms the modulus
 * with Math.hypot, so the two agree to round-off rather than bit for bit while
 * |f| lies between about 1e-154 and 1e154. Outside that range the re² + im² of
 * csqrt_ underflows or overflows where hypot does not. A refractive index or a
 * permittivity is far inside it. */
TMM_EXPORT("tmm_jet_sqrt")
void tmm_jet_sqrt(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jsqrt_j(jread(&jets[8L * i])));
}

TMM_EXPORT("tmm_jet_exp")
void tmm_jet_exp(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jexp_j(jread(&jets[8L * i])));
}

TMM_EXPORT("tmm_jet_log")
void tmm_jet_log(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jlog_j(jread(&jets[8L * i])));
}

/* f^p for real p, as exp(p log f); the branch cut of the log applies. */
TMM_EXPORT("tmm_jet_power")
void tmm_jet_power(const double *jets, const double *exponent, int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jexp_j(jscale(jlog_j(jread(&jets[8L * i])), exponent[i])));
}

TMM_EXPORT("tmm_jet_sin_cos")
void tmm_jet_sin_cos(const double *jets, int count, double *sine, double *cosine) {
    for (int i = 0; i < count; i++) {
        jet s, c;
        jsincos(jread(&jets[8L * i]), &s, &c);
        jwrite(&sine[8L * i], s);
        jwrite(&cosine[8L * i], c);
    }
}

TMM_EXPORT("tmm_jet_sin")
void tmm_jet_sin(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) {
        jet s, c;
        jsincos(jread(&jets[8L * i]), &s, &c);
        jwrite(&out[8L * i], s);
    }
}

TMM_EXPORT("tmm_jet_cos")
void tmm_jet_cos(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) {
        jet s, c;
        jsincos(jread(&jets[8L * i]), &s, &c);
        jwrite(&out[8L * i], c);
    }
}

/* `value` is count reals, `derivatives` count × 3 reals [f', f'', f''']. */
TMM_EXPORT("tmm_jet_compose")
void tmm_jet_compose(const double *value, const double *derivatives, const double *x,
                     int count, double *out) {
    for (int i = 0; i < count; i++)
        jwrite(&out[8L * i], jcompose_j(value[i], &derivatives[3L * i], jread(&x[8L * i])));
}

/* λ(ω) = 2πc/ω, in whatever unit `wavelength` is given in. */
TMM_EXPORT("tmm_wavelength_omega_jet")
void tmm_wavelength_omega_jet(const double *wavelength, const double *omega, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jwavelength(wavelength[i], omega[i]));
}

/* The real parts of one jet with the real parts of another as imaginary parts. */
TMM_EXPORT("tmm_jet_with_imaginary_part")
void tmm_jet_with_imaginary_part(const double *realJets, const double *imaginaryJets, int count,
                                 double *out) {
    for (int i = 0; i < count; i++) {
        jet re = jread(&realJets[8L * i]), im = jread(&imaginaryJets[8L * i]), o;
        for (int k = 0; k < JET_N; k++) o.c[k] = cmk(re.c[k].re, im.c[k].re);
        jwrite(&out[8L * i], o);
    }
}

TMM_EXPORT("tmm_jet_real_part")
void tmm_jet_real_part(const double *jets, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jreal(jread(&jets[8L * i])));
}

/* The jet as given while its value is at least `minimum`, otherwise the
 * constant `minimum`; NaN takes the constant. */
TMM_EXPORT("tmm_jet_clamp_real_minimum")
void tmm_jet_clamp_real_minimum(const double *jets, const double *minimum, int count, double *out) {
    for (int i = 0; i < count; i++) {
        jet a = jread(&jets[8L * i]);
        jwrite(&out[8L * i], a.c[0].re >= minimum[i] ? a : jconst(minimum[i], 0.0));
    }
}

TMM_EXPORT("tmm_jet_clamp_imaginary")
void tmm_jet_clamp_imaginary(const double *jets, const double *limit, int count, double *out) {
    for (int i = 0; i < count; i++) jwrite(&out[8L * i], jclampim(jread(&jets[8L * i]), limit[i]));
}

/* out[i] = 1 when all eight doubles of jet i are finite, else 0. */
TMM_EXPORT("tmm_jet_is_finite")
void tmm_jet_is_finite(const double *jets, int count, int *out) {
    for (int i = 0; i < count; i++) {
        int finite = 1;
        for (int k = 0; k < 8; k++) if (!isfinite(jets[8L * i + k])) finite = 0;
        out[i] = finite;
    }
}
