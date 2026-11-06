/*
 * tmmcore's C kernel timed from C, with no Python in the process. Loads the
 * library inside the installed tmmcore wheel, the one tmmcore's Python calls
 * go through, and times tmm_rta_spectrum on the stacks bench_c.py writes:
 * each case's stack alone, then its batch of stacks, one call per stack, s
 * and p together.
 *
 *   bench_c <library> <stacks file> <outputs file>
 *
 * bench_c.py builds it, writes the stacks file and checks the outputs file.
 * For each case the stacks file holds int32 N, nLam and nBatch, then float64
 * values: the angle in degrees, the wavelengths [nLam], the incident and
 * substrate indices [nLam][re, im], the layer indices [N][nLam][re, im], the
 * thicknesses [N] and the batch's thicknesses [nBatch][N], all little-endian.
 *
 * Prints a line per case, the seconds for one stack and for the batch, each
 * the fastest of repeated timings, and writes the six outputs of the
 * single-stack call, Rs Ts As Rp Tp Ap, to the outputs file.
 */

#ifndef _WIN32
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdio.h>
#include <stdlib.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#include <time.h>
#endif

/* tmm_rta_spectrum as src/c/tmmcore.h declares it. */
typedef void (*rta_spectrum_fn)(const double *lambdas, int nLam,
                                const double *theta_deg, int nTheta,
                                const double *n0arr, const double *nsarr,
                                const double *matNK, const double *thick, int N,
                                double *outRs, double *outTs, double *outAs,
                                double *outRp, double *outTp, double *outAp);

typedef struct {
    rta_spectrum_fn rta_spectrum;
    int N, nLam, nBatch;
    double theta;
    double *lambdas, *n0, *ns, *nk, *thick, *batch;
    double *out[6];
} job;

static double now(void)
{
#ifdef _WIN32
    LARGE_INTEGER count, frequency;
    QueryPerformanceCounter(&count);
    QueryPerformanceFrequency(&frequency);
    return (double)count.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
#endif
}

static rta_spectrum_fn load(const char *path)
{
#ifdef _WIN32
    HMODULE library = LoadLibraryExA(path, NULL, LOAD_WITH_ALTERED_SEARCH_PATH);
    return library ? (rta_spectrum_fn)GetProcAddress(library, "tmm_rta_spectrum") : NULL;
#else
    void *library = dlopen(path, RTLD_NOW);
    return library ? (rta_spectrum_fn)dlsym(library, "tmm_rta_spectrum") : NULL;
#endif
}

static double *read_doubles(FILE *file, size_t count)
{
    double *values = malloc((count ? count : 1) * sizeof *values);
    if (!values || fread(values, sizeof *values, count, file) != count) {
        fprintf(stderr, "bench_c: the stacks file ends early\n");
        exit(1);
    }
    return values;
}

static void spectrum(const job *j, const double *thick)
{
    j->rta_spectrum(j->lambdas, j->nLam, &j->theta, 1, j->n0, j->ns, j->nk, thick, j->N,
                    j->out[0], j->out[1], j->out[2], j->out[3], j->out[4], j->out[5]);
}

static void one_stack(const job *j)
{
    spectrum(j, j->thick);
}

static void whole_batch(const job *j)
{
    for (int s = 0; s < j->nBatch; s++)
        spectrum(j, j->batch + (size_t)s * j->N);
}

/* The fastest of repeated runs, after at least min_reps runs and min_seconds,
 * the floors bench_fast.py sets for tmm_faster and the tmmcore wheel. */
static double best_of(void (*run)(const job *), const job *j, int min_reps, double min_seconds)
{
    double best = 1e300, start = now();
    for (int reps = 0; reps < min_reps || now() - start < min_seconds; reps++) {
        double t = now();
        run(j);
        t = now() - t;
        if (t < best)
            best = t;
    }
    return best;
}

int main(int argc, char **argv)
{
    if (argc != 4) {
        fprintf(stderr, "usage: bench_c <library> <stacks file> <outputs file>\n");
        return 2;
    }
    job j = { load(argv[1]) };
    if (!j.rta_spectrum) {
        fprintf(stderr, "bench_c: no tmm_rta_spectrum in %s\n", argv[1]);
        return 1;
    }
    FILE *in = fopen(argv[2], "rb"), *out = fopen(argv[3], "wb");
    if (!in || !out) {
        fprintf(stderr, "bench_c: cannot open %s or %s\n", argv[2], argv[3]);
        return 1;
    }

    int header[3];
    while (fread(header, sizeof header[0], 3, in) == 3) {
        j.N = header[0];
        j.nLam = header[1];
        j.nBatch = header[2];
        double *theta = read_doubles(in, 1);
        j.theta = *theta;
        free(theta);
        j.lambdas = read_doubles(in, (size_t)j.nLam);
        j.n0 = read_doubles(in, 2 * (size_t)j.nLam);
        j.ns = read_doubles(in, 2 * (size_t)j.nLam);
        j.nk = read_doubles(in, 2 * (size_t)j.N * (size_t)j.nLam);
        j.thick = read_doubles(in, (size_t)j.N);
        j.batch = read_doubles(in, (size_t)j.nBatch * (size_t)j.N);
        for (int q = 0; q < 6; q++)
            j.out[q] = malloc((size_t)j.nLam * sizeof(double));

        one_stack(&j);   /* warm-up, and the outputs bench_c.py checks */
        for (int q = 0; q < 6; q++)
            fwrite(j.out[q], sizeof(double), (size_t)j.nLam, out);
        double t_one = best_of(one_stack, &j, 3, 1.0);
        double t_batch = j.nBatch ? best_of(whole_batch, &j, 2, 0.5) : 0.0;
        printf("%.6e %.6e\n", t_one, t_batch);

        free(j.lambdas);
        free(j.n0);
        free(j.ns);
        free(j.nk);
        free(j.thick);
        free(j.batch);
        for (int q = 0; q < 6; q++)
            free(j.out[q]);
    }
    fclose(in);
    fclose(out);
    return 0;
}
