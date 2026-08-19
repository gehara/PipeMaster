/*
 * compute_sumstats.c — Standalone summary statistics + SFS from haplotype matrices
 *
 * Computes the same statistics as msABC from R integer matrices
 * (e.g., scrm output), without running the coalescent simulator.
 *
 * Two entry points:
 *
 * 1. compute_sumstats_call(hapmat, config, npop)
 *    Single-locus: returns named numeric vector of stats + SFS
 *
 * 2. compute_sumstats_batch_call(haplist, config, npop)
 *    Multi-locus: takes a list of haplotype matrices (one per locus),
 *    computes per-locus stats, aggregates (mean/var/skew/kurtosis),
 *    and returns a list(sumstats=numeric_vector, sfs=numeric_vector).
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "data_sumstat.h"

/* Forward declarations for seg_sites (name mismatch in header vs implementation) */
int seg_sites(int* almap, int n, int segsites);

/* ================================================================
 * Helper: compute SFS from char** list
 * Returns a flat array of joint SFS counts (folded to minor allele).
 * For 1 pop: length = nsam-1 bins (freq 1..nsam-1)
 * For multi-pop: length = prod(config[p]+1 for p in 0..npop-1)
 * ================================================================ */
static void compute_sfs(char **list, int nsam, int nsegsites,
                         int *config, int npop,
                         double *sfs, long long sfs_len) {

    memset(sfs, 0, (size_t)sfs_len * sizeof(double));

    if (nsegsites <= 0) return;

    if (npop == 1) {
        /* 1D SFS: bins for freq 1..nsam-1, folded */
        for (int s = 0; s < nsegsites; s++) {
            int freq = 0;
            for (int i = 0; i < nsam; i++)
                freq += (list[i][s] == '1');
            if (freq >= 1 && freq <= nsam - 1) {
                int folded = (freq <= nsam - freq) ? freq : nsam - freq;
                sfs[folded - 1] += 1.0;
            }
        }
    } else {
        /* Joint SFS: flattened multi-dimensional array */
        int cum_sizes[npop + 1];
        cum_sizes[0] = 0;
        for (int p = 0; p < npop; p++)
            cum_sizes[p + 1] = cum_sizes[p] + config[p];

        int strides[npop];
        strides[0] = 1;
        for (int p = 1; p < npop; p++)
            strides[p] = strides[p - 1] * (config[p - 1] + 1);

        for (int s = 0; s < nsegsites; s++) {
            int counts[npop];
            int total = 0;
            for (int p = 0; p < npop; p++) {
                counts[p] = 0;
                for (int i = cum_sizes[p]; i < cum_sizes[p + 1]; i++)
                    counts[p] += (list[i][s] == '1');
                total += counts[p];
            }
            /* Fold to minor allele side */
            int do_flip = 0;
            if (2 * total > nsam) {
                do_flip = 1;
            } else if (2 * total == nsam) {
                for (int p = 0; p < npop; p++) {
                    int comp = config[p] - counts[p];
                    if (counts[p] > comp) { do_flip = 1; break; }
                    if (counts[p] < comp) { break; }
                }
            }
            if (do_flip) {
                for (int p = 0; p < npop; p++)
                    counts[p] = config[p] - counts[p];
            }
            int flat_idx = 0;
            for (int p = 0; p < npop; p++)
                flat_idx += counts[p] * strides[p];
            if (flat_idx >= 0 && flat_idx < sfs_len)
                sfs[flat_idx] += 1.0;
        }
    }
}


/* ================================================================
 * Helper: convert R integer matrix to char** list
 * R matrix is column-major: hapdata[row + col*nsam]
 * ================================================================ */
static char **make_list(int *hapdata, int nsam, int nsegsites) {
    char **list = (char **)malloc(nsam * sizeof(char *));
    if (!list) return NULL;
    for (int i = 0; i < nsam; i++) {
        list[i] = (char *)malloc((nsegsites + 1) * sizeof(char));
        if (!list[i]) return NULL;
        for (int j = 0; j < nsegsites; j++)
            list[i][j] = (hapdata[i + j * nsam] == 1) ? '1' : '0';
        list[i][nsegsites] = '\0';
    }
    return list;
}

static void free_list(char **list, int nsam) {
    for (int i = 0; i < nsam; i++) free(list[i]);
    free(list);
}


/* ================================================================
 * Per-locus stat computation into a pre-allocated double array
 *
 * Stat layout (by stat type): for each of 5 main stats, per-pop
 * then global; then pairwise block (Fst global + shared/private/
 * fixed/Fst per pair); then nhap/Hd per-pop interleaved + global.
 *
 * Note: thomson and FayWuH excluded (same as sim.all.stats default)
 * ================================================================ */
#define N_MAIN_STATS 5   /* segs, pi, thetaW, tajd, ZnS */
#define N_HAP_STATS  2   /* nhap, Hd */
#define N_PAIR_STATS 4   /* shared, private, fixed, Fst */

static int nstat_total(int npop) {
    int npairs = npop * (npop - 1) / 2;
    return N_MAIN_STATS * (npop + 1)
         + N_HAP_STATS * (npop + 1)
         + ((npop > 1) ? (1 + N_PAIR_STATS * npairs) : 0);
}

static void compute_locus_stats(char **list, int nsam, int nsegsites,
                                 int *config, int npop, int skip_zns,
                                 double *out) {

    int nstats = nstat_total(npop);
    memset(out, 0, nstats * sizeof(double));

    if (nsegsites == 0) return;

    /* Global maps */
    int *almap = (int *)malloc(nsegsites * sizeof(int));
    int *missing_arr = (int *)calloc(nsegsites, sizeof(int));
    int **ders = (int **)malloc(nsam * sizeof(int *));
    for (int i = 0; i < nsam; i++)
        ders[i] = (int *)malloc(nsegsites * sizeof(int));

    set_almap(almap, list, nsam, nsegsites, 0, nsam);
    set_missing(missing_arr, list, nsam, nsegsites, 0, nsam);
    set_ders(ders, list, nsam, nsegsites, 0, nsam);

    double hn, sqhn, bn;
    denominators(nsam, &hn, &sqhn, &bn);

    int n_segs = seg_sites(almap, nsam, nsegsites);
    double val_pi = 0., val_tw = 0., val_tajd = 0., val_zns = 0.;
    double val_dvk = 0., val_dvh = 0.;

    theta_pi(&val_pi, almap, nsam, nsegsites, missing_arr);
    theta_w(&val_tw, almap, nsam, nsegsites, hn, missing_arr);
    if (n_segs > 0) tajD(&val_tajd, n_segs, nsam, val_tw, val_pi, hn, sqhn);
    if (!skip_zns && n_segs > 1) ZnS(&val_zns, list, nsam, nsegsites, 0, almap, ders, missing_arr);
    dvstat(list, nsam, nsegsites, 0, nsam, &val_dvk, &val_dvh);

    /* Per-population stats into arrays */
    int samples_b[npop], samples_e[npop];
    samples_b[0] = 0;
    samples_e[0] = config[0];
    for (int i = 1; i < npop; i++) {
        samples_b[i] = samples_b[i - 1] + config[i - 1];
        samples_e[i] = samples_b[i] + config[i];
    }

    double pop_segs[npop], pop_pi_a[npop], pop_tw[npop];
    double pop_tajd[npop], pop_zns_a[npop];
    double pop_nhap[npop], pop_hd[npop];
    memset(pop_segs, 0, npop * sizeof(double));
    memset(pop_pi_a, 0, npop * sizeof(double));
    memset(pop_tw, 0, npop * sizeof(double));
    memset(pop_tajd, 0, npop * sizeof(double));
    memset(pop_zns_a, 0, npop * sizeof(double));
    memset(pop_nhap, 0, npop * sizeof(double));
    memset(pop_hd, 0, npop * sizeof(double));

    int *pop_almap = (int *)malloc(nsegsites * sizeof(int));
    int *pop_missing = (int *)calloc(nsegsites, sizeof(int));
    int **pop_ders = (int **)malloc(nsam * sizeof(int *));
    for (int i = 0; i < nsam; i++)
        pop_ders[i] = (int *)malloc(nsegsites * sizeof(int));

    for (int p = 0; p < npop; p++) {
        int n_pop = config[p];
        if (n_pop >= 2) {
            set_almap(pop_almap, list, nsam, nsegsites, samples_b[p], samples_e[p]);
            set_missing(pop_missing, list, nsam, nsegsites, samples_b[p], samples_e[p]);
            set_ders(pop_ders, list, nsam, nsegsites, samples_b[p], samples_e[p]);

            double pop_hn, pop_sqhn, pop_bn;
            denominators(n_pop, &pop_hn, &pop_sqhn, &pop_bn);

            int p_segs = seg_sites(pop_almap, n_pop, nsegsites);
            pop_segs[p] = (double)p_segs;
            theta_pi(&pop_pi_a[p], pop_almap, n_pop, nsegsites, pop_missing);
            theta_w(&pop_tw[p], pop_almap, n_pop, nsegsites, pop_hn, pop_missing);
            if (p_segs > 0) tajD(&pop_tajd[p], p_segs, n_pop, pop_tw[p], pop_pi_a[p], pop_hn, pop_sqhn);
            if (!skip_zns && p_segs > 1) ZnS(&pop_zns_a[p], list, n_pop, nsegsites, 0, pop_almap, pop_ders, pop_missing);
            dvstat(list, nsam, nsegsites, samples_b[p], samples_e[p], &pop_nhap[p], &pop_hd[p]);
        }
    }

    /* === Write stats in by-stat-type order === */
    int si = 0;

    /* Main stats: per-pop then global for each */
    for (int p = 0; p < npop; p++) out[si++] = pop_segs[p];
    out[si++] = (double)n_segs;
    for (int p = 0; p < npop; p++) out[si++] = pop_pi_a[p];
    out[si++] = val_pi;
    for (int p = 0; p < npop; p++) out[si++] = pop_tw[p];
    out[si++] = val_tw;
    for (int p = 0; p < npop; p++) out[si++] = pop_tajd[p];
    out[si++] = val_tajd;
    for (int p = 0; p < npop; p++) out[si++] = pop_zns_a[p];
    out[si++] = val_zns;

    /* Pairwise block: global Fst + shared/private/fixed/Fst per pair */
    if (npop > 1) {
        double *weights = (double *)malloc(npop * sizeof(double));
        for (int i = 0; i < npop; i++) weights[i] = 1.0 / (double)npop;

        double piT = 0., piS = 0., piB = 0., piD = 0.;
        double **shared_mat = (double **)malloc(npop * sizeof(double *));
        double **private_mat = (double **)malloc(npop * sizeof(double *));
        double **fixed_mat = (double **)malloc(npop * sizeof(double *));
        for (int i = 0; i < npop; i++) {
            shared_mat[i] = (double *)calloc(npop, sizeof(double));
            private_mat[i] = (double *)calloc(npop, sizeof(double));
            fixed_mat[i] = (double *)calloc(npop, sizeof(double));
        }

        calculations(weights, list, config, npop, nsegsites, nsam,
                      &piT, &piS, &piB, &piD,
                      shared_mat, private_mat, fixed_mat, 0);

        out[si++] = Fst_HBK(piS, piT);
        for (int i = 0; i < npop - 1; i++) {
            for (int j = i + 1; j < npop; j++) {
                out[si++] = shared_mat[i][j];
                out[si++] = private_mat[i][j];
                out[si++] = fixed_mat[i][j];
                double piT_ij = 0., piS_ij = 0., piB_ij = 0., piD_ij = 0.;
                pairwiseFstcalculations(i, j, weights, list, config, npop,
                                        nsegsites, nsam,
                                        &piT_ij, &piS_ij, &piB_ij, &piD_ij);
                double fst_ij = Fst_HBK(piS_ij, piT_ij);
                out[si++] = isnan(fst_ij) ? 0. : fst_ij;
            }
        }

        for (int i = 0; i < npop; i++) {
            free(shared_mat[i]); free(private_mat[i]); free(fixed_mat[i]);
        }
        free(shared_mat); free(private_mat); free(fixed_mat);
        free(weights);
    }

    /* nhap/Hd: per-pop interleaved, then global */
    for (int p = 0; p < npop; p++) {
        out[si++] = pop_nhap[p];
        out[si++] = pop_hd[p];
    }
    out[si++] = val_dvk;
    out[si++] = val_dvh;

    /* Cleanup */
    for (int i = 0; i < nsam; i++) { free(ders[i]); free(pop_ders[i]); }
    free(almap); free(missing_arr); free(ders);
    free(pop_almap); free(pop_missing); free(pop_ders);
}


/* ================================================================
 * Build stat names vector (used by both single and batch)
 * ================================================================ */
static SEXP make_stat_names(int npop) {
    int nstats = nstat_total(npop);
    SEXP names = PROTECT(allocVector(STRSXP, nstats));
    int si = 0;
    char buf[64];

    /* Main stats: per-pop then global for each */
    const char *main_stats[] = {"segs", "pi", "thetaW", "tajd", "ZnS"};
    for (int s = 0; s < N_MAIN_STATS; s++) {
        for (int p = 0; p < npop; p++) {
            snprintf(buf, 64, "%s_%d", main_stats[s], p+1);
            SET_STRING_ELT(names, si++, mkChar(buf));
        }
        SET_STRING_ELT(names, si++, mkChar(main_stats[s]));
    }

    /* Pairwise block: global Fst + shared/private/fixed/Fst per pair */
    if (npop > 1) {
        SET_STRING_ELT(names, si++, mkChar("Fst"));
        for (int i = 0; i < npop - 1; i++) {
            for (int j = i + 1; j < npop; j++) {
                snprintf(buf, 64, "shared_%d_%d", i+1, j+1); SET_STRING_ELT(names, si++, mkChar(buf));
                snprintf(buf, 64, "private_%d_%d", i+1, j+1); SET_STRING_ELT(names, si++, mkChar(buf));
                snprintf(buf, 64, "fixed_%d_%d", i+1, j+1); SET_STRING_ELT(names, si++, mkChar(buf));
                snprintf(buf, 64, "Fst_%d_%d", i+1, j+1); SET_STRING_ELT(names, si++, mkChar(buf));
            }
        }
    }

    /* nhap/Hd: per-pop interleaved, then global */
    for (int p = 0; p < npop; p++) {
        snprintf(buf, 64, "nhap_%d", p+1); SET_STRING_ELT(names, si++, mkChar(buf));
        snprintf(buf, 64, "Hd_%d", p+1); SET_STRING_ELT(names, si++, mkChar(buf));
    }
    SET_STRING_ELT(names, si++, mkChar("nhap"));
    SET_STRING_ELT(names, si++, mkChar("Hd"));

    UNPROTECT(1);
    return names;
}


/* ================================================================
 * .Call("compute_sumstats_call", hapmat, config, npop)
 *
 * Single-locus entry point.
 * Returns: list(stats = named_numeric, sfs = numeric)
 * ================================================================ */
/* PM-SFSLEN-20260819: see the matching note in scrm_stats.cpp. prod(config+1)
 * was computed as a 32-bit int and overflowed silently for high-npop designs
 * (10 pops x 8 haps wraps NEGATIVE; 10 pops x 10 haps wraps to a wrong
 * positive length). Accumulate in 64 bits and refuse rather than return
 * garbage. */
#define PM_SFS_MAX_BINS 50000000LL   /* 50M bins = 400 MB per array */

static long long pm_sfs_len(int nsam, int *config, int npop, const char *who) {
    if (npop == 1) return (long long)nsam - 1;
    long long n = 1;
    int p;
    for (p = 0; p < npop; p++) {
        n *= (long long)(config[p] + 1);
        if (n > PM_SFS_MAX_BINS)
            Rf_error("%s: joint SFS would need > %lld bins across %d populations "
                     "(%lld and still multiplying). Exceeds the %lld-bin limit and "
                     "would overflow. Use skip.sfs = TRUE, or fewer populations.",
                     who, PM_SFS_MAX_BINS, npop, n, PM_SFS_MAX_BINS);
    }
    return n;
}

SEXP compute_sumstats_call(SEXP hapmat_sexp, SEXP config_sexp, SEXP npop_sexp) {

    if (!isInteger(hapmat_sexp) || !isMatrix(hapmat_sexp))
        Rf_error("compute_sumstats: 'hapmat' must be an integer matrix");

    SEXP dim = getAttrib(hapmat_sexp, R_DimSymbol);
    int nsam = INTEGER(dim)[0];
    int nsegsites = INTEGER(dim)[1];
    int *hapdata = INTEGER(hapmat_sexp);

    int npop = asInteger(npop_sexp);
    if (!isInteger(config_sexp) || length(config_sexp) != npop)
        Rf_error("compute_sumstats: 'config' must be integer vector of length npop");
    int *config = INTEGER(config_sexp);

    int config_sum = 0;
    for (int i = 0; i < npop; i++) config_sum += config[i];
    if (config_sum != nsam)
        Rf_error("compute_sumstats: sum(config)=%d != nsam=%d", config_sum, nsam);

    char **list = make_list(hapdata, nsam, nsegsites);
    if (!list) Rf_error("compute_sumstats: malloc failed");

    /* Compute stats */
    int nstats = nstat_total(npop);
    double *stat_vals = (double *)calloc(nstats, sizeof(double));
    compute_locus_stats(list, nsam, nsegsites, config, npop, 0, stat_vals);

    /* Compute SFS */
    long long sfs_len = pm_sfs_len(nsam, config, npop, "compute_sumstats");
    double *sfs_vals = (double *)calloc((size_t)sfs_len, sizeof(double));
    if (sfs_vals == NULL)
        Rf_error("compute_sumstats: could not allocate %lld SFS bins", sfs_len);
    compute_sfs(list, nsam, nsegsites, config, npop, sfs_vals, sfs_len);

    free_list(list, nsam);

    /* Build result list */
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("stats"));
    SET_STRING_ELT(result_names, 1, mkChar("sfs"));
    setAttrib(result, R_NamesSymbol, result_names);

    /* Stats vector */
    SEXP stats_sexp = PROTECT(allocVector(REALSXP, nstats));
    memcpy(REAL(stats_sexp), stat_vals, nstats * sizeof(double));
    SEXP stat_names = PROTECT(make_stat_names(npop));
    setAttrib(stats_sexp, R_NamesSymbol, stat_names);
    SET_VECTOR_ELT(result, 0, stats_sexp);

    /* SFS vector */
    SEXP sfs_sexp = PROTECT(allocVector(REALSXP, (R_xlen_t)sfs_len));
    memcpy(REAL(sfs_sexp), sfs_vals, (size_t)sfs_len * sizeof(double));
    SET_VECTOR_ELT(result, 1, sfs_sexp);

    free(stat_vals);
    free(sfs_vals);

    UNPROTECT(5);
    return result;
}


/* ================================================================
 * .Call("compute_sumstats_batch_call", haplist, config, npop)
 *
 * Multi-locus entry point. haplist is an R list of integer matrices.
 * Returns: list(
 *   stats = named_numeric (mean/var/skew/kurtosis of each stat),
 *   sfs   = numeric (accumulated SFS across all loci)
 * )
 * ================================================================ */
SEXP compute_sumstats_batch_call(SEXP haplist_sexp, SEXP config_sexp,
                                  SEXP npop_sexp, SEXP skip_zns_sexp,
                                  SEXP skip_sfs_sexp) {

    if (!isNewList(haplist_sexp))
        Rf_error("compute_sumstats_batch: 'haplist' must be a list");

    int nloci = length(haplist_sexp);
    int npop = asInteger(npop_sexp);
    int skip_zns = asLogical(skip_zns_sexp);
    int skip_sfs = asLogical(skip_sfs_sexp);
    if (skip_sfs == NA_LOGICAL) skip_sfs = 0;
    if (!isInteger(config_sexp) || length(config_sexp) != npop)
        Rf_error("compute_sumstats_batch: 'config' must be integer vector of length npop");
    int *config = INTEGER(config_sexp);

    int config_sum = 0;
    for (int i = 0; i < npop; i++) config_sum += config[i];
    int nsam = config_sum;

    int nstats = nstat_total(npop);

    /* Accumulators for moments */
    double *sum1 = (double *)calloc(nstats, sizeof(double));  /* sum of x */
    double *sum2 = (double *)calloc(nstats, sizeof(double));  /* sum of x^2 */
    double *sum3 = (double *)calloc(nstats, sizeof(double));  /* sum of x^3 */
    double *sum4 = (double *)calloc(nstats, sizeof(double));  /* sum of x^4 */
    double *locus_stats = (double *)calloc(nstats, sizeof(double));

    /* SFS accumulator */
    long long sfs_len = skip_sfs ? 0 : pm_sfs_len(nsam, config, npop, "compute_sumstats_batch");
    double *sfs_accum = NULL, *locus_sfs = NULL;
    if (!skip_sfs) {
        sfs_accum = (double *)calloc((size_t)sfs_len, sizeof(double));
        locus_sfs = (double *)calloc((size_t)sfs_len, sizeof(double));
        if (sfs_accum == NULL || locus_sfs == NULL)
            Rf_error("compute_sumstats_batch: could not allocate %lld SFS bins "
                     "(%.2f GB each). Use skip.sfs = TRUE.",
                     sfs_len, sfs_len * 8.0 / (1024*1024*1024));
    }

    /* Process each locus */
    for (int loc = 0; loc < nloci; loc++) {
        SEXP hapmat = VECTOR_ELT(haplist_sexp, loc);

        if (isNull(hapmat) || length(hapmat) == 0) {
            /* Empty locus (no segregating sites) — contributes zeros */
            continue;
        }

        if (!isInteger(hapmat) || !isMatrix(hapmat))
            Rf_error("compute_sumstats_batch: element %d must be an integer matrix", loc + 1);

        SEXP dim = getAttrib(hapmat, R_DimSymbol);
        int loc_nsam = INTEGER(dim)[0];
        int loc_segs = INTEGER(dim)[1];

        if (loc_nsam != nsam)
            Rf_error("compute_sumstats_batch: locus %d has %d samples, expected %d",
                     loc + 1, loc_nsam, nsam);

        char **list = make_list(INTEGER(hapmat), nsam, loc_segs);
        if (!list) Rf_error("compute_sumstats_batch: malloc failed at locus %d", loc + 1);

        /* Stats */
        compute_locus_stats(list, nsam, loc_segs, config, npop, skip_zns, locus_stats);
        for (int k = 0; k < nstats; k++) {
            double v = locus_stats[k];
            sum1[k] += v;
            sum2[k] += v * v;
            sum3[k] += v * v * v;
            sum4[k] += v * v * v * v;
        }

        /* SFS */
        if (!skip_sfs) {
            compute_sfs(list, nsam, loc_segs, config, npop, locus_sfs, sfs_len);
            for (long long k = 0; k < sfs_len; k++)
                sfs_accum[k] += locus_sfs[k];
        }

        free_list(list, nsam);

        /* Check for user interrupt every 100 loci */
        if ((loc + 1) % 100 == 0) R_CheckUserInterrupt();
    }

    /* Compute mean, var, skew, kurtosis */
    int n_agg = 4;  /* mean, var, skew, kurtosis */
    int n_out = nstats * n_agg;
    double *out_stats = (double *)calloc(n_out, sizeof(double));

    double nl = (double)nloci;
    for (int k = 0; k < nstats; k++) {
        double mean = sum1[k] / nl;
        double var = (sum2[k] / nl) - mean * mean;
        double sd = (var > 0) ? sqrt(var) : 0.;

        double skew = 0., kurt = 0.;
        if (sd > 0 && nl > 2) {
            double m3 = (sum3[k] / nl) - 3 * mean * (sum2[k] / nl) + 2 * mean * mean * mean;
            skew = m3 / (sd * sd * sd);
        }
        if (sd > 0 && nl > 3) {
            double m4 = (sum4[k] / nl) - 4 * mean * (sum3[k] / nl)
                         + 6 * mean * mean * (sum2[k] / nl) - 3 * mean * mean * mean * mean;
            kurt = m4 / (sd * sd * sd * sd) - 3.0;  /* excess kurtosis */
        }

        out_stats[k] = mean;
        out_stats[nstats + k] = var;
        out_stats[2 * nstats + k] = skew;
        out_stats[3 * nstats + k] = kurt;
    }

    /* Build output */
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("stats"));
    SET_STRING_ELT(result_names, 1, mkChar("sfs"));
    setAttrib(result, R_NamesSymbol, result_names);

    /* Stats: named vector with s_mean_*, s_var_*, s_skew_*, s_kurt_* */
    SEXP stats_sexp = PROTECT(allocVector(REALSXP, n_out));
    memcpy(REAL(stats_sexp), out_stats, n_out * sizeof(double));

    SEXP stat_base_names = PROTECT(make_stat_names(npop));
    SEXP agg_names = PROTECT(allocVector(STRSXP, n_out));
    const char *prefixes[] = {"s_mean_", "s_var_", "s_skew_", "s_kurt_"};
    char nbuf[128];
    for (int a = 0; a < n_agg; a++) {
        for (int k = 0; k < nstats; k++) {
            snprintf(nbuf, 128, "%s%s", prefixes[a], CHAR(STRING_ELT(stat_base_names, k)));
            SET_STRING_ELT(agg_names, a * nstats + k, mkChar(nbuf));
        }
    }
    setAttrib(stats_sexp, R_NamesSymbol, agg_names);
    SET_VECTOR_ELT(result, 0, stats_sexp);

    /* SFS */
    SEXP sfs_sexp = PROTECT(allocVector(REALSXP, (R_xlen_t)sfs_len));
    if (!skip_sfs)
        memcpy(REAL(sfs_sexp), sfs_accum, (size_t)sfs_len * sizeof(double));
    SET_VECTOR_ELT(result, 1, sfs_sexp);

    /* Cleanup */
    free(sum1); free(sum2); free(sum3); free(sum4);
    free(locus_stats); free(out_stats);
    free(sfs_accum); free(locus_sfs);

    UNPROTECT(6);
    return result;
}
