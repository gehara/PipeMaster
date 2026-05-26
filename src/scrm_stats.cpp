/*
 * scrm_stats.cpp — Combined coalescent simulation + summary statistics
 *
 * Uses scrm's SMC' coalescent engine for fast simulation with recombination,
 * then computes summary statistics directly from the in-memory haplotype data
 * using PipeMaster's C stat functions (data_sumstat.c). No R matrix overhead.
 *
 * scrm source: Staab et al. (2015) Genetics 199:1-6, GPL-3 license.
 * Copyright (C) 2013, 2014 Paul R. Staab, Sha (Joe) Zhu, Dirk Metzler, Gerton Lunter
 *
 * .Call("scrm_stats_call", args, config, npop, nloci, skip_zns)
 *   args:     character string with scrm command-line arguments
 *   config:   integer vector of per-population sample sizes
 *   npop:     integer, number of populations
 *   nloci:    integer, number of loci (must match scrm args)
 *   skip_zns: logical, whether to skip ZnS computation
 *
 * Returns: list(stats = named_numeric_vector, sfs = numeric_vector)
 */

/* Include scrm headers BEFORE R headers to avoid macro conflicts */
#include <cstring>
#include <memory>
#include <cmath>

#include "scrm/param.h"
#include "scrm/forest.h"
#include "scrm/model.h"
#include "scrm/random/fastfunc.h"
#include "scrm/summary_statistics/seg_sites.h"

/* Now include R headers, but undef problematic macros first */
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* R defines 'length' as a macro which breaks C++ std library headers.
 * We already have Rinternals loaded, so undef to avoid conflicts
 * with any remaining C++ standard library usage. */
#undef length
#undef Realloc
#undef Free

#include "r_random_generator.h"

/* C stat functions from data_sumstat.c */
extern "C" {
    int seg_sites(int* almap, int n, int segsites);
    void set_almap(int* almap, char** list, int n, int segsites, int start, int end);
    void set_missing(int* missing, char** list, int n, int segsites, int start, int end);
    void set_ders(int** ders, char** list, int n, int segsites, int start, int end);
    int denominators(int n, double* hn, double* sqhn, double* bn);
    int theta_pi(double* theta, int* almap, int n, int segsites, int* missing);
    int theta_w(double* theta, int* almap, int n, int segsites, double denom, int* missing);
    int tajD(double* tajd, int segsites, int n, double thetaw, double thetap, double hn, double sqhn);
    int ZnS(double* ld, char** list, int n, int segsites, int filter, int* almap, int** ders, int* missing);
    int dvstat(char** list, int n, int segs, int start, int end, double* dvk, double* dvh);
    int calculations(double* weights, char** list, int* config, int npop,
                     int segsites, int n,
                     double* piT, double* piS, double* piB, double* piD,
                     double** shared, double** private_alleles, double** fixed_dif,
                     int derived);
    int pairwiseFstcalculations(int popi, int popj, double* weights,
                                char** list, int* config, int npop,
                                int segsites, int n,
                                double* piT, double* piS, double* piB, double* piD);
    double Fst_HBK(double piS, double piT);
}

/* Stat layout constants
 * Layout (by stat type): for each of 5 main stats, per-pop then global;
 * then pairwise block (Fst global + shared/private/fixed/Fst per pair);
 * then nhap/Hd per-pop interleaved + global.
 */
#define N_MAIN_STATS 5   /* segs, pi, thetaW, tajd, ZnS */
#define N_HAP_STATS  2   /* nhap, Hd */
#define N_PAIR_STATS 4   /* shared, private, fixed, Fst */

static int nstat_total(int npop) {
    int npairs = npop * (npop - 1) / 2;
    return N_MAIN_STATS * (npop + 1)
         + N_HAP_STATS * (npop + 1)
         + ((npop > 1) ? (1 + N_PAIR_STATS * npairs) : 0);
}

/*
 * Convert scrm SegSites haplotype data to char** list format
 * and compute stats + SFS for one locus.
 */
static void process_locus(SegSites *ss, int nsam, int nsegsites,
                           int *config, int npop, int skip_zns,
                           double *stat_out, int nstats,
                           double *sfs_out, int sfs_len) {

    memset(stat_out, 0, nstats * sizeof(double));
    memset(sfs_out, 0, sfs_len * sizeof(double));

    if (nsegsites == 0) return;

    /* Build char** list from scrm's valarray<bool> haplotypes */
    char **list = (char **)malloc(nsam * sizeof(char *));
    for (int i = 0; i < nsam; i++) {
        list[i] = (char *)malloc((nsegsites + 1) * sizeof(char));
        for (int j = 0; j < nsegsites; j++) {
            list[i][j] = (*(ss->getHaplotype(j)))[i] ? '1' : '0';
        }
        list[i][nsegsites] = '\0';
    }

    /* === Compute all statistics into temporary storage === */

    int *almap = (int *)malloc(nsegsites * sizeof(int));
    int *missing_arr = (int *)calloc(nsegsites, sizeof(int));
    int **ders = (int **)malloc(nsam * sizeof(int *));
    for (int i = 0; i < nsam; i++)
        ders[i] = (int *)malloc(nsegsites * sizeof(int));

    set_almap(almap, list, nsam, nsegsites, 0, nsam);
    set_missing(missing_arr, list, nsam, nsegsites, 0, nsam);
    if (!skip_zns)
        set_ders(ders, list, nsam, nsegsites, 0, nsam);

    double hn, sqhn, bn;
    denominators(nsam, &hn, &sqhn, &bn);

    /* Global stats */
    int n_segs = seg_sites(almap, nsam, nsegsites);
    double val_pi = 0., val_tw = 0., val_tajd = 0., val_zns = 0.;
    double val_nhap = 0., val_hd = 0.;

    theta_pi(&val_pi, almap, nsam, nsegsites, missing_arr);
    theta_w(&val_tw, almap, nsam, nsegsites, hn, missing_arr);
    if (n_segs > 0) tajD(&val_tajd, n_segs, nsam, val_tw, val_pi, hn, sqhn);
    if (!skip_zns && n_segs > 1) ZnS(&val_zns, list, nsam, nsegsites, 0, almap, ders, missing_arr);
    dvstat(list, nsam, nsegsites, 0, nsam, &val_nhap, &val_hd);

    /* Per-population stats into arrays */
    int *samples_b = (int *)malloc(npop * sizeof(int));
    int *samples_e = (int *)malloc(npop * sizeof(int));
    samples_b[0] = 0;
    samples_e[0] = config[0];
    for (int i = 1; i < npop; i++) {
        samples_b[i] = samples_b[i - 1] + config[i - 1];
        samples_e[i] = samples_b[i] + config[i];
    }

    double *pop_segs = (double *)calloc(npop, sizeof(double));
    double *pop_pi   = (double *)calloc(npop, sizeof(double));
    double *pop_tw   = (double *)calloc(npop, sizeof(double));
    double *pop_tajd = (double *)calloc(npop, sizeof(double));
    double *pop_zns  = (double *)calloc(npop, sizeof(double));
    double *pop_nhap = (double *)calloc(npop, sizeof(double));
    double *pop_hd   = (double *)calloc(npop, sizeof(double));

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

            double pop_hn, pop_sqhn, pop_bn;
            denominators(n_pop, &pop_hn, &pop_sqhn, &pop_bn);

            int p_segs = seg_sites(pop_almap, n_pop, nsegsites);
            pop_segs[p] = (double)p_segs;
            theta_pi(&pop_pi[p], pop_almap, n_pop, nsegsites, pop_missing);
            theta_w(&pop_tw[p], pop_almap, n_pop, nsegsites, pop_hn, pop_missing);
            if (p_segs > 0) tajD(&pop_tajd[p], p_segs, n_pop, pop_tw[p], pop_pi[p], pop_hn, pop_sqhn);
            if (!skip_zns && p_segs > 1) {
                set_ders(pop_ders, list, nsam, nsegsites, samples_b[p], samples_e[p]);
                ZnS(&pop_zns[p], list, nsam, nsegsites, 0, pop_almap, pop_ders, pop_missing);
            }
            dvstat(list, nsam, nsegsites, samples_b[p], samples_e[p], &pop_nhap[p], &pop_hd[p]);
        }
    }

    /* === Write stats in by-stat-type order === */
    int si = 0;

    /* Main stats: per-pop then global for each */
    for (int p = 0; p < npop; p++) stat_out[si++] = pop_segs[p];
    stat_out[si++] = (double)n_segs;
    for (int p = 0; p < npop; p++) stat_out[si++] = pop_pi[p];
    stat_out[si++] = val_pi;
    for (int p = 0; p < npop; p++) stat_out[si++] = pop_tw[p];
    stat_out[si++] = val_tw;
    for (int p = 0; p < npop; p++) stat_out[si++] = pop_tajd[p];
    stat_out[si++] = val_tajd;
    for (int p = 0; p < npop; p++) stat_out[si++] = pop_zns[p];
    stat_out[si++] = val_zns;

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

        stat_out[si++] = Fst_HBK(piS, piT);
        for (int i = 0; i < npop - 1; i++) {
            for (int j = i + 1; j < npop; j++) {
                stat_out[si++] = shared_mat[i][j];
                stat_out[si++] = private_mat[i][j];
                stat_out[si++] = fixed_mat[i][j];
                double piT_ij = 0., piS_ij = 0., piB_ij = 0., piD_ij = 0.;
                pairwiseFstcalculations(i, j, weights, list, config, npop,
                                        nsegsites, nsam,
                                        &piT_ij, &piS_ij, &piB_ij, &piD_ij);
                double fst_ij = Fst_HBK(piS_ij, piT_ij);
                stat_out[si++] = std::isnan(fst_ij) ? 0. : fst_ij;
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
        stat_out[si++] = pop_nhap[p];
        stat_out[si++] = pop_hd[p];
    }
    stat_out[si++] = val_nhap;
    stat_out[si++] = val_hd;

    free(pop_segs); free(pop_pi); free(pop_tw);
    free(pop_tajd); free(pop_zns);
    free(pop_nhap); free(pop_hd);

    /* === SFS === */
    if (npop == 1) {
        /* 1D folded SFS */
        for (int s = 0; s < nsegsites; s++) {
            int freq = almap[s];
            if (freq >= 1 && freq <= nsam - 1) {
                int folded = (freq <= nsam - freq) ? freq : nsam - freq;
                sfs_out[folded - 1] += 1.0;
            }
        }
    } else {
        /* Joint SFS */
        int cum_sizes[npop + 1];
        cum_sizes[0] = 0;
        for (int p = 0; p < npop; p++)
            cum_sizes[p + 1] = cum_sizes[p] + config[p];

        int strides[npop];
        strides[0] = 1;
        for (int p = 1; p < npop; p++)
            strides[p] = strides[p - 1] * (config[p - 1] + 1);

        /* Need per-pop allele counts — recompute from list */
        for (int s = 0; s < nsegsites; s++) {
            int counts[npop];
            int total = 0;
            for (int p = 0; p < npop; p++) {
                counts[p] = 0;
                for (int i = cum_sizes[p]; i < cum_sizes[p + 1]; i++)
                    counts[p] += (list[i][s] == '1');
                total += counts[p];
            }
            /* Fold to minor allele */
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
                sfs_out[flat_idx] += 1.0;
        }
    }

    /* Cleanup */
    for (int i = 0; i < nsam; i++) {
        free(list[i]); free(ders[i]); free(pop_ders[i]);
    }
    free(list); free(almap); free(missing_arr); free(ders);
    free(pop_almap); free(pop_missing); free(pop_ders);
    free(samples_b); free(samples_e);
}


/* Build stat names — by-stat-type layout matching msABC convention */
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
 * Main .Call entry point
 * ================================================================ */
extern "C" SEXP scrm_stats_call(SEXP args_sexp, SEXP config_sexp,
                                  SEXP npop_sexp, SEXP skip_zns_sexp,
                                  SEXP rec_rates_sexp,
                                  SEXP mu_rates_sexp) {

    if (!isString(args_sexp) || Rf_length(args_sexp) != 1)
        Rf_error("scrm_stats: 'args' must be a single character string");

    const char *args_cstr = CHAR(STRING_ELT(args_sexp, 0));
    int npop = asInteger(npop_sexp);
    int skip_zns = asLogical(skip_zns_sexp);

    if (!isInteger(config_sexp) || Rf_length(config_sexp) != npop)
        Rf_error("scrm_stats: 'config' must be integer vector of length npop");
    int *config = INTEGER(config_sexp);

    int nsam = 0;
    for (int i = 0; i < npop; i++) nsam += config[i];

    /* Optional per-locus recombination rates (per bp per generation).
     * NULL -> use rate baked into args.
     * Numeric vector of length nloci -> override per locus. */
    bool use_per_locus_rates = !isNull(rec_rates_sexp);
    double *per_locus_rates = NULL;
    if (use_per_locus_rates) {
        if (!isReal(rec_rates_sexp))
            Rf_error("scrm_stats: 'rec_rates' must be numeric or NULL");
        per_locus_rates = REAL(rec_rates_sexp);
    }

    /* Optional per-locus mutation rates (per bp per generation).
     * NULL -> use rate baked into args (via -t).
     * Numeric vector of length nloci -> override per locus. mu is post-hoc
     * to the ARG (only used by SegSites to place mutations on the built
     * genealogy), so per-locus override is O(1) and does not affect cost. */
    bool use_per_locus_mu = !isNull(mu_rates_sexp);
    double *per_locus_mu = NULL;
    if (use_per_locus_mu) {
        if (!isReal(mu_rates_sexp))
            Rf_error("scrm_stats: 'mu_rates' must be numeric or NULL");
        per_locus_mu = REAL(mu_rates_sexp);
    }

    /* Parse scrm arguments */
    std::string args_str(args_cstr);
    Param param(args_str);
    Model model = param.parse();

    int nloci = (int)model.loci_number();

    if (use_per_locus_rates && Rf_length(rec_rates_sexp) != nloci)
        Rf_error("scrm_stats: 'rec_rates' length (%d) != nloci (%d)",
                 Rf_length(rec_rates_sexp), nloci);
    if (use_per_locus_mu && Rf_length(mu_rates_sexp) != nloci)
        Rf_error("scrm_stats: 'mu_rates' length (%d) != nloci (%d)",
                 Rf_length(mu_rates_sexp), nloci);

    /* Set up RNG using R's random number generator */
    std::shared_ptr<FastFunc> ff = std::make_shared<FastFunc>();
    RRandomGenerator rrg(ff);

    /* Initialize forest */
    Forest forest(&model, &rrg);

    /* Allocate accumulators */
    int nstats = nstat_total(npop);
    double *sum1 = (double *)calloc(nstats, sizeof(double));
    double *sum2 = (double *)calloc(nstats, sizeof(double));
    double *sum3 = (double *)calloc(nstats, sizeof(double));
    double *sum4 = (double *)calloc(nstats, sizeof(double));
    double *locus_stats = (double *)calloc(nstats, sizeof(double));

    int sfs_len;
    if (npop == 1) {
        sfs_len = nsam - 1;
    } else {
        sfs_len = 1;
        for (int p = 0; p < npop; p++) sfs_len *= (config[p] + 1);
    }
    double *sfs_accum = (double *)calloc(sfs_len, sizeof(double));
    double *locus_sfs = (double *)calloc(sfs_len, sizeof(double));

    /* Get RNG state from R */
    GetRNGstate();

    /* Main loop over loci */
    for (int loc = 0; loc < nloci; loc++) {

        /* Per-locus rate override (if supplied): change the model's rec
         * rate before building this locus's tree. setRecombinationRate
         * with per_locus=false, scaled=false treats the value as per bp
         * per generation. */
        if (use_per_locus_rates) {
            model.setRecombinationRate(per_locus_rates[loc],
                                        false, false, 0.0);
        }
        /* Per-locus mutation rate override (if supplied): mu does not
         * affect the ARG; SegSites reads model.mutation_rate() to place
         * mutations on the already-built genealogy. So this is O(1) and
         * adds negligible cost. setMutationRate(per_locus=false,
         * scaled=false) stores the value as per bp per generation. */
        if (use_per_locus_mu) {
            model.setMutationRate(per_locus_mu[loc],
                                   false, false, 0.0);
        }

        /* Build tree using scrm's SMC' algorithm */
        forest.buildInitialTree();
        while (forest.next_base() < model.loci_length()) {
            forest.sampleNextGenealogy();
        }

        /* Get segregating sites from the forest */
        SegSites *ss = NULL;
        for (size_t i = 0; i < model.countSummaryStatistics(); i++) {
            SummaryStatistic *s = model.getSummaryStatistic(i);
            if (typeid(*s) == typeid(SegSites)) {
                ss = dynamic_cast<SegSites *>(s);
                break;
            }
        }

        if (ss != NULL) {
            int nsegsites = (int)ss->countMutations();

            /* Compute stats + SFS directly from scrm's in-memory data */
            process_locus(ss, nsam, nsegsites, config, npop, skip_zns,
                          locus_stats, nstats, locus_sfs, sfs_len);

            /* Accumulate moments */
            for (int k = 0; k < nstats; k++) {
                double v = locus_stats[k];
                sum1[k] += v;
                sum2[k] += v * v;
                sum3[k] += v * v * v;
                sum4[k] += v * v * v * v;
            }

            /* Accumulate SFS */
            for (int k = 0; k < sfs_len; k++)
                sfs_accum[k] += locus_sfs[k];
        }

        forest.clear();

        /* Check for user interrupt every 100 loci */
        if ((loc + 1) % 100 == 0) R_CheckUserInterrupt();
    }

    /* Return RNG state to R */
    PutRNGstate();

    /* Compute mean, var, skew, kurtosis */
    int n_agg = 4;
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
            kurt = m4 / (sd * sd * sd * sd) - 3.0;
        }

        out_stats[k] = mean;
        out_stats[nstats + k] = var;
        out_stats[2 * nstats + k] = skew;
        out_stats[3 * nstats + k] = kurt;
    }

    /* Build R output */
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("stats"));
    SET_STRING_ELT(result_names, 1, mkChar("sfs"));
    setAttrib(result, R_NamesSymbol, result_names);

    /* Stats vector with names */
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

    /* SFS vector */
    SEXP sfs_sexp = PROTECT(allocVector(REALSXP, sfs_len));
    memcpy(REAL(sfs_sexp), sfs_accum, sfs_len * sizeof(double));
    SET_VECTOR_ELT(result, 1, sfs_sexp);

    /* Cleanup */
    free(sum1); free(sum2); free(sum3); free(sum4);
    free(locus_stats); free(out_stats);
    free(sfs_accum); free(locus_sfs);

    UNPROTECT(6);
    return result;
}


/* ================================================================
 * Multi-length entry point: accepts vector of scrm command strings
 * (one per length group) and accumulates stats across all groups.
 *
 * .Call("scrm_stats_multi_call", args_vec, config, npop, skip_zns, total_nloci)
 *   args_vec:     character vector, one scrm command per length group
 *   config:       integer vector of per-population sample sizes
 *   npop:         integer, number of populations
 *   skip_zns:     logical, whether to skip ZnS computation
 *   total_nloci:  integer, total loci across all groups (for moment computation)
 * ================================================================ */
extern "C" SEXP scrm_stats_multi_call(SEXP args_vec_sexp, SEXP config_sexp,
                                       SEXP npop_sexp, SEXP skip_zns_sexp,
                                       SEXP total_nloci_sexp,
                                       SEXP rec_rates_sexp,
                                       SEXP mu_rates_sexp) {

    if (!isString(args_vec_sexp))
        Rf_error("scrm_stats_multi: 'args_vec' must be a character vector");

    int n_groups = Rf_length(args_vec_sexp);
    int npop = asInteger(npop_sexp);
    int skip_zns = asLogical(skip_zns_sexp);
    int total_nloci = asInteger(total_nloci_sexp);

    if (!isInteger(config_sexp) || Rf_length(config_sexp) != npop)
        Rf_error("scrm_stats_multi: 'config' must be integer vector of length npop");
    int *config = INTEGER(config_sexp);

    int nsam = 0;
    for (int i = 0; i < npop; i++) nsam += config[i];

    /* Optional per-locus recombination rates (per bp per generation).
     * NULL -> use rate baked into args.
     * Numeric vector of length total_nloci -> override per locus.
     * Order must align with the order of groups in args_vec then loci
     * within each group. */
    bool use_per_locus_rates = !isNull(rec_rates_sexp);
    double *per_locus_rates = NULL;
    if (use_per_locus_rates) {
        if (!isReal(rec_rates_sexp))
            Rf_error("scrm_stats_multi: 'rec_rates' must be numeric or NULL");
        if (Rf_length(rec_rates_sexp) != total_nloci)
            Rf_error("scrm_stats_multi: 'rec_rates' length (%d) != total_nloci (%d)",
                     Rf_length(rec_rates_sexp), total_nloci);
        per_locus_rates = REAL(rec_rates_sexp);
    }

    /* Optional per-locus mutation rates (per bp per generation).
     * NULL -> use rate baked into args (via -t).
     * Numeric vector of length total_nloci -> override per locus. Same
     * group/locus order as rec_rates. mu is post-hoc to the ARG (only
     * affects SegSites mutation placement), so per-locus override is
     * O(1) and does not affect simulation cost. */
    bool use_per_locus_mu = !isNull(mu_rates_sexp);
    double *per_locus_mu = NULL;
    if (use_per_locus_mu) {
        if (!isReal(mu_rates_sexp))
            Rf_error("scrm_stats_multi: 'mu_rates' must be numeric or NULL");
        if (Rf_length(mu_rates_sexp) != total_nloci)
            Rf_error("scrm_stats_multi: 'mu_rates' length (%d) != total_nloci (%d)",
                     Rf_length(mu_rates_sexp), total_nloci);
        per_locus_mu = REAL(mu_rates_sexp);
    }

    /* Allocate shared accumulators */
    int nstats = nstat_total(npop);
    double *sum1 = (double *)calloc(nstats, sizeof(double));
    double *sum2 = (double *)calloc(nstats, sizeof(double));
    double *sum3 = (double *)calloc(nstats, sizeof(double));
    double *sum4 = (double *)calloc(nstats, sizeof(double));
    double *locus_stats = (double *)calloc(nstats, sizeof(double));

    int sfs_len;
    if (npop == 1) {
        sfs_len = nsam - 1;
    } else {
        sfs_len = 1;
        for (int p = 0; p < npop; p++) sfs_len *= (config[p] + 1);
    }
    double *sfs_accum = (double *)calloc(sfs_len, sizeof(double));
    double *locus_sfs = (double *)calloc(sfs_len, sizeof(double));

    /* Get RNG state from R */
    GetRNGstate();

    /* Running index across all loci in all groups, used to look up
     * the per-locus rec rate in the supplied vector. */
    int global_loc_idx = 0;

    /* Loop over length groups */
    for (int g = 0; g < n_groups; g++) {

        const char *args_cstr = CHAR(STRING_ELT(args_vec_sexp, g));
        std::string args_str(args_cstr);
        Param param(args_str);
        Model model = param.parse();

        int group_nloci = (int)model.loci_number();

        std::shared_ptr<FastFunc> ff = std::make_shared<FastFunc>();
        RRandomGenerator rrg(ff);
        Forest forest(&model, &rrg);

        for (int loc = 0; loc < group_nloci; loc++) {

            /* Per-locus rate override (if supplied) */
            if (use_per_locus_rates) {
                model.setRecombinationRate(
                    per_locus_rates[global_loc_idx],
                    false, false, 0.0);
            }
            /* Per-locus mutation rate override (if supplied). Cheap:
             * mu only affects SegSites mutation placement, not the ARG. */
            if (use_per_locus_mu) {
                model.setMutationRate(
                    per_locus_mu[global_loc_idx],
                    false, false, 0.0);
            }
            global_loc_idx++;

            forest.buildInitialTree();
            while (forest.next_base() < model.loci_length()) {
                forest.sampleNextGenealogy();
            }

            SegSites *ss = NULL;
            for (size_t i = 0; i < model.countSummaryStatistics(); i++) {
                SummaryStatistic *s = model.getSummaryStatistic(i);
                if (typeid(*s) == typeid(SegSites)) {
                    ss = dynamic_cast<SegSites *>(s);
                    break;
                }
            }

            if (ss != NULL) {
                int nsegsites = (int)ss->countMutations();

                process_locus(ss, nsam, nsegsites, config, npop, skip_zns,
                              locus_stats, nstats, locus_sfs, sfs_len);

                for (int k = 0; k < nstats; k++) {
                    double v = locus_stats[k];
                    sum1[k] += v;
                    sum2[k] += v * v;
                    sum3[k] += v * v * v;
                    sum4[k] += v * v * v * v;
                }
                for (int k = 0; k < sfs_len; k++)
                    sfs_accum[k] += locus_sfs[k];
            }

            forest.clear();

            if ((loc + 1) % 100 == 0) R_CheckUserInterrupt();
        }

        R_CheckUserInterrupt();
    }

    /* Return RNG state to R */
    PutRNGstate();

    /* Compute mean, var, skew, kurtosis using total_nloci */
    int n_agg = 4;
    int n_out = nstats * n_agg;
    double *out_stats = (double *)calloc(n_out, sizeof(double));
    double nl = (double)total_nloci;

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
            kurt = m4 / (sd * sd * sd * sd) - 3.0;
        }

        out_stats[k] = mean;
        out_stats[nstats + k] = var;
        out_stats[2 * nstats + k] = skew;
        out_stats[3 * nstats + k] = kurt;
    }

    /* Build R output (same format as scrm_stats_call) */
    SEXP result = PROTECT(allocVector(VECSXP, 2));
    SEXP result_names = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(result_names, 0, mkChar("stats"));
    SET_STRING_ELT(result_names, 1, mkChar("sfs"));
    setAttrib(result, R_NamesSymbol, result_names);

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

    SEXP sfs_sexp = PROTECT(allocVector(REALSXP, sfs_len));
    memcpy(REAL(sfs_sexp), sfs_accum, sfs_len * sizeof(double));
    SET_VECTOR_ELT(result, 1, sfs_sexp);

    /* Cleanup */
    free(sum1); free(sum2); free(sum3); free(sum4);
    free(locus_stats); free(out_stats);
    free(sfs_accum); free(locus_sfs);

    UNPROTECT(6);
    return result;
}
