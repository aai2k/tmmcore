/*
 * tmmcore.h : the calls of the tmmcore C kernel.
 *
 * Compile tmm_kernel.c, which includes the other sources, with any C99
 * compiler and link libm. The arguments and results of every call are in the
 * API reference, https://aai2k.github.io/tmmcore/api/
 *
 *   wavelengths and thicknesses in nm, angles in degrees from the normal
 *   the index n + ik, with k > 0 for absorbing media
 *   pol 0 = s, 1 = p
 *   every array double precision and allocated by the caller; an array of
 *   more than one dimension is flat and row-major, and a complex number is
 *   two doubles, [re, im]
 *
 * A call that cannot allocate its working memory writes NaN to its outputs,
 * except where its comment here says otherwise.
 */

#ifndef TMMCORE_H
#define TMMCORE_H

#ifdef __cplusplus
extern "C" {
#endif

/* R, T and A */

void tmm_rta(double lambda_nm, double theta_deg, int pol,
             double n0_re, double n0_im, double ns_re, double ns_im,
             const double *layers, int N, double *out);

void tmm_rta_spectrum(const double *lambdas, int nLam,
                      const double *theta_deg, int nTheta,
                      const double *n0arr, const double *nsarr,
                      const double *matNK, const double *thick, int N,
                      double *outRs, double *outTs, double *outAs,
                      double *outRp, double *outTp, double *outAp);

/* Thickness derivatives */

void tmm_jacobian(double lambda_nm, double theta_deg, int pol,
                  double n0_re, double n0_im, double ns_re, double ns_im,
                  const double *layers, int N,
                  double *dRdd, double *dTdd, double *dAdd, double *base);

void tmm_jacobian_spectrum(const double *lambdas, int nLam,
                           const double *theta_deg, int nTheta, int pol,
                           const double *n0arr, const double *nsarr,
                           const double *matNK, const double *thick, int N,
                           double *dRdd, double *dTdd, double *dAdd, double *base);

void tmm_hessian(double lambda_nm, double theta_deg, int pol,
                 double n0_re, double n0_im, double ns_re, double ns_im,
                 const double *layers, int N,
                 double *dRdd, double *dTdd, double *dAdd,
                 double *d2Rdd, double *d2Tdd, double *d2Add, double *base);

void tmm_hessian_spectrum(const double *lambdas, int nLam,
                          const double *theta_deg, int nTheta, int pol,
                          const double *n0arr, const double *nsarr,
                          const double *matNK, const double *thick, int N,
                          double *dRdd, double *dTdd, double *dAdd,
                          double *d2Rdd, double *d2Tdd, double *d2Add, double *base);

/* Needle scan */

void tmm_needle_scan(double lambda_nm, double theta_deg, int pol,
                     double n0_re, double n0_im, double ns_re, double ns_im,
                     const double *layers, int N,
                     const double *candNs, int nCand,
                     const double *fracs, int nFrac,
                     double *base, double *gaps, double *intra);

void tmm_needle_scan_spectrum(const double *lambdas, int nLam,
                              const double *theta_deg, int nTheta, int pol,
                              const double *n0arr, const double *nsarr,
                              const double *matNK, const double *thick, int N,
                              const double *candNK, int nCand,
                              const double *fracs, int nFrac,
                              double *base, double *gaps, double *intra);

/* Phase, GD, GDD and TOD. Every index is a jet of 8 doubles; sinJet and
 * sinJets may be NULL. */

void tmm_phase(double lambda, double theta_deg, int pol,
               const double *n0jet, const double *nsjet,
               const double *layerJets, const double *thick, int N,
               double omega, const double *sinJet, double *out);

void tmm_phase_spectrum(const double *lambdas, int nLam,
                        const double *theta_deg, int nTheta, int pol,
                        const double *n0jets, const double *nsjets,
                        const double *matJets, const double *thick, int N,
                        const double *omegas, const double *sinJets, double *out);

void tmm_phase_jacobian(double lambda, double theta_deg, int pol,
                        const double *n0jet, const double *nsjet,
                        const double *layerJets, const double *thick, int N,
                        double omega, const double *sinJet, double *out, double *deriv);

void tmm_phase_jacobian_spectrum(const double *lambdas, int nLam,
                                 const double *theta_deg, int nTheta, int pol,
                                 const double *n0jets, const double *nsjets,
                                 const double *matJets, const double *thick, int N,
                                 const double *omegas, const double *sinJets,
                                 double *out, double *deriv);

void tmm_omega_from_lambda_nm(const double *lambda_nm, int count, double *out);

/* Taylor-jet arithmetic. Each call applies one operation to `count` jets of
 * 8 doubles each; a real argument that goes with each jet is an array of
 * `count`. An output may be the same array as an input of its shape. */

void tmm_jet_constant(const double *value, const double *imaginary, int count, double *out);
void tmm_jet_from_derivatives(const double *derivatives, int count, double *out);
void tmm_jet_derivatives(const double *jets, int count, double *out);
void tmm_jet_add(const double *left, const double *right, int count, double *out);
void tmm_jet_subtract(const double *left, const double *right, int count, double *out);
void tmm_jet_scale(const double *jets, const double *scalar, int count, double *out);
void tmm_jet_multiply(const double *left, const double *right, int count, double *out);
void tmm_jet_reciprocal(const double *jets, int count, double *out);
void tmm_jet_divide(const double *numerator, const double *denominator, int count, double *out);
void tmm_jet_sqrt(const double *jets, int count, double *out);
void tmm_jet_exp(const double *jets, int count, double *out);
void tmm_jet_log(const double *jets, int count, double *out);
void tmm_jet_power(const double *jets, const double *exponent, int count, double *out);
void tmm_jet_sin_cos(const double *jets, int count, double *sine, double *cosine);
void tmm_jet_sin(const double *jets, int count, double *out);
void tmm_jet_cos(const double *jets, int count, double *out);
void tmm_jet_compose(const double *value, const double *derivatives, const double *x,
                     int count, double *out);
void tmm_wavelength_omega_jet(const double *wavelength, const double *omega, int count, double *out);
void tmm_jet_with_imaginary_part(const double *realJets, const double *imaginaryJets,
                                 int count, double *out);
void tmm_jet_real_part(const double *jets, int count, double *out);
void tmm_jet_clamp_real_minimum(const double *jets, const double *minimum, int count, double *out);
void tmm_jet_clamp_imaginary(const double *jets, const double *limit, int count, double *out);
void tmm_jet_is_finite(const double *jets, int count, int *out);

/* Monitoring and deposition */

void tmm_monitor_curve(double lambda_nm, double theta_deg,
                       double n0_re, double n0_im, double ns_re, double ns_im,
                       const double *base, int NB,
                       double ng_re, double ng_im,
                       const double *dArr, int nD,
                       double *outRs, double *outTs,
                       double *outRp, double *outTp,
                       double *outRrs, double *outRrp);

/* Returns 1, or 0 with the outputs untouched when it could not allocate its
 * working memory. */
int tmm_deposition_spectra(const double *lambdas, int nLam,
                           const double *theta_deg, int nTheta,
                           const double *n0arr, const double *nsarr,
                           const double *matNK, const double *thick, int N,
                           double *outRs, double *outTs,
                           double *outRp, double *outTp,
                           double *outRrs, double *outRrp);

/* A finished stack under a growing top layer. */
typedef struct tmm_growing_stack tmm_growing_stack;

/* Returns NULL when it could not allocate its memory. */
tmm_growing_stack *tmm_growing_stack_create(const double *lambdas, int nLam,
                                            const double *theta_deg, int nTheta,
                                            const double *n0arr, const double *nsarr,
                                            const double *matNK, const double *thick,
                                            int NB);

void tmm_growing_stack_set_top(tmm_growing_stack *h, const double *ngNK);

/* Returns 1, or 0 with the outputs untouched for a NULL handle, or for d > 0
 * before tmm_growing_stack_set_top. */
int tmm_growing_stack_sample(tmm_growing_stack *h, double d,
                             double *outRs, double *outTs,
                             double *outRp, double *outTp,
                             double *outRrs, double *outRrp);

void tmm_growing_stack_free(tmm_growing_stack *h);

#ifdef __cplusplus
}
#endif

#endif
