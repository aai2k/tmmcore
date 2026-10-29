/*
 * tmm_kernel.c : the whole kernel as one translation unit, for a build that
 * compiles a single file:
 *
 *   cc -std=c99 -O2 -c tmm_kernel.c
 *
 * The package's own builds compile the five sources together with link-time
 * optimization instead. Conventions and references are in tmm_internal.h.
 */

#include "rta.c"
#include "growing.c"
#include "derivatives.c"
#include "phase.c"
#include "jet.c"
