/*
 * sfs_sim.c - Native C SFS simulator for PipeMaster
 *
 * Provides a .Call() entry point that runs coalescent simulations via
 * ms/msABC and accumulates the Site Frequency Spectrum in C, avoiding
 * R-level per-locus overhead. Supports:
 *   - 1D SFS (single population)
 *   - Joint SFS (multi-population)
 *   - one.snp mode (sample one segregating site per locus)
 *   - Per-locus theta values (each locus gets its own ms command)
 *
 * The key optimization: all loci are processed in a single .Call() from R.
 * After each locus, the SFS is computed directly from the haplotype matrix
 * (list[]) in C without any string I/O or R round-trips.
 */

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <setjmp.h>
#include "sfs_sim.h"
#include "msABC_capture.h"

/* Global SFS accumulator state */
int sfs_mode_active = 0;
int sfs_one_snp = 0;
int sfs_npop = 0;
int *sfs_pop_sizes = NULL;
double *sfs_accumulator = NULL;
int sfs_accum_len = 0;
double *sfs_theta_array = NULL;
int sfs_nloci = 0;
int sfs_expected_mode = 0;

/* Jump buffer from ms.c */
extern jmp_buf msABC_jmpbuf;

/* From msABC_wrapper.c */
extern void msABC_init_output_stream(void);
extern void msABC_close_output_stream(void);

/* From ms.c */
extern void msABC_reset_streec_statics(void);
extern void msABC_set_jmpbuf_active(int active);
extern int msABC_main(int argc, char **argv);
extern void seedit_r(unsigned short *seedv);

/* RNG for one.snp sampling - uses R's RNG */
static int sfs_rand_int(int n) {
    double u = unif_rand();
    int r = (int)(u * n);
    if (r >= n) r = n - 1;
    return r;
}

/*
 * sfs_accumulate - called from ms.c after each replicate when sfs_mode_active > 0
 *
 * list[i] is a null-terminated string of '0' and '1' for individual i
 * nsam = total sample size
 * segsites = number of segregating sites
 * config[p] = sample size for population p (from pars.cp.config)
 * npop = number of populations
 */
void sfs_accumulate(char **list, int nsam, int segsites, int *config, int npop) {
    if (!sfs_mode_active || sfs_accumulator == NULL) return;
    if (segsites <= 0) return;

    if (sfs_mode_active == 1) {
        /* ---- 1D SFS (single population) ---- */
        /* SFS has nsam-1 bins for derived allele frequencies 1..nsam-1 */

        if (sfs_one_snp) {
            int col = sfs_rand_int(segsites);
            int freq = 0;
            for (int i = 0; i < nsam; i++) {
                freq += (list[i][col] == '1');
            }
            if (freq >= 1 && freq <= nsam - 1) {
                sfs_accumulator[freq - 1] += 1.0;
            }
        } else {
            for (int s = 0; s < segsites; s++) {
                int freq = 0;
                for (int i = 0; i < nsam; i++) {
                    freq += (list[i][s] == '1');
                }
                if (freq >= 1 && freq <= nsam - 1) {
                    sfs_accumulator[freq - 1] += 1.0;
                }
            }
        }
    } else {
        /* ---- Joint SFS (multi-population) ---- */
        /* Accumulator is a flattened multi-dimensional array:
         * dims = (config[0]+1) x (config[1]+1) x ... x (config[npop-1]+1)
         * First population varies fastest, matching R's expand.grid order.
         */
        int cum_sizes[npop + 1];
        cum_sizes[0] = 0;
        for (int p = 0; p < npop; p++) {
            cum_sizes[p + 1] = cum_sizes[p] + config[p];
        }

        int strides[npop];
        strides[0] = 1;
        for (int p = 1; p < npop; p++) {
            strides[p] = strides[p - 1] * (config[p - 1] + 1);
        }

        if (sfs_one_snp) {
            int col = sfs_rand_int(segsites);
            int flat_idx = 0;
            for (int p = 0; p < npop; p++) {
                int count = 0;
                for (int i = cum_sizes[p]; i < cum_sizes[p + 1]; i++) {
                    count += (list[i][col] == '1');
                }
                flat_idx += count * strides[p];
            }
            sfs_accumulator[flat_idx] += 1.0;
        } else {
            for (int s = 0; s < segsites; s++) {
                int flat_idx = 0;
                for (int p = 0; p < npop; p++) {
                    int count = 0;
                    for (int i = cum_sizes[p]; i < cum_sizes[p + 1]; i++) {
                        count += (list[i][s] == '1');
                    }
                    flat_idx += count * strides[p];
                }
                sfs_accumulator[flat_idx] += 1.0;
            }
        }
    }
}


/* ---- Parse command string into argc/argv ---- */

static int parse_command(const char *cmd, char ***argv_out) {
    char *tmp = strdup(cmd);
    if (tmp == NULL) return 0;

    int argc = 0;
    char *token = strtok(tmp, " \t");
    while (token != NULL) {
        argc++;
        token = strtok(NULL, " \t");
    }
    free(tmp);

    int total = argc + 1;
    char **argv = (char **)malloc(total * sizeof(char *));
    if (argv == NULL) return 0;

    argv[0] = strdup("msABC");

    tmp = strdup(cmd);
    token = strtok(tmp, " \t");
    for (int i = 1; i <= argc; i++) {
        argv[i] = strdup(token);
        token = strtok(NULL, " \t");
    }
    free(tmp);

    *argv_out = argv;
    return total;
}

static void free_argv(int argc, char **argv) {
    for (int i = 0; i < argc; i++) {
        free(argv[i]);
    }
    free(argv);
}


/*
 * msABC_sfs_call - .Call() entry point for native SFS simulation
 *
 * Processes all loci in a single C call. Each locus command is a full ms
 * command string (nsam 1 [flags]). The SFS is accumulated across loci in C.
 *
 * Arguments:
 *   commands_sexp: character vector of ms commands (one per locus)
 *                  Each command: "nsam 1 -t theta [demographic flags]"
 *   pop_sizes_sexp: integer vector of per-population sample sizes
 *   one_snp_sexp: logical, TRUE to sample one SNP per locus
 *   seed_sexp: integer vector of length 3 (or NULL for auto-seed)
 *
 * Returns: numeric vector of SFS counts
 *   1-pop: length nsam-1 (frequencies 1 to nsam-1)
 *   multi-pop: length prod(pop_sizes + 1) (flattened joint SFS)
 */
SEXP msABC_sfs_call(SEXP commands_sexp, SEXP pop_sizes_sexp,
                     SEXP one_snp_sexp, SEXP seed_sexp,
                     SEXP method_sexp) {

    /* Validate inputs */
    if (!isString(commands_sexp) || length(commands_sexp) < 1) {
        Rf_error("msABC_sfs_call: 'commands' must be a character vector");
    }
    if (!isInteger(pop_sizes_sexp) || length(pop_sizes_sexp) < 1) {
        Rf_error("msABC_sfs_call: 'pop_sizes' must be an integer vector");
    }
    if (!isLogical(one_snp_sexp) || length(one_snp_sexp) != 1) {
        Rf_error("msABC_sfs_call: 'one_snp' must be a single logical value");
    }

    /* Parse method: "stochastic" (default) or "expected" */
    int use_expected = 0;
    if (!isNull(method_sexp) && isString(method_sexp) && length(method_sexp) == 1) {
        const char *method_str = CHAR(STRING_ELT(method_sexp, 0));
        if (strcmp(method_str, "expected") == 0) use_expected = 1;
    }

    int nloci = length(commands_sexp);
    int npop = length(pop_sizes_sexp);
    int *pop_sizes = INTEGER(pop_sizes_sexp);
    int one_snp = LOGICAL(one_snp_sexp)[0];

    /* Set seed if provided */
    if (!isNull(seed_sexp)) {
        if (!isInteger(seed_sexp) || length(seed_sexp) != 3) {
            Rf_error("msABC_sfs_call: 'seed' must be an integer vector of length 3, or NULL");
        }
        int *sv = INTEGER(seed_sexp);
        unsigned short seedv[3];
        seedv[0] = (unsigned short)sv[0];
        seedv[1] = (unsigned short)sv[1];
        seedv[2] = (unsigned short)sv[2];
        seedit_r(seedv);
    }

    /* Compute total sample size and SFS dimensions */
    int nsam = 0;
    for (int p = 0; p < npop; p++) {
        nsam += pop_sizes[p];
    }

    int accum_len;
    if (npop == 1) {
        accum_len = nsam - 1;  /* freq bins 1..nsam-1 */
    } else {
        accum_len = 1;
        for (int p = 0; p < npop; p++) {
            accum_len *= (pop_sizes[p] + 1);
        }
    }

    /* Allocate SFS accumulator */
    double *accum = (double *)calloc(accum_len, sizeof(double));
    if (accum == NULL) {
        Rf_error("msABC_sfs_call: failed to allocate SFS accumulator (%d entries)", accum_len);
    }

    /* Set up global SFS mode */
    sfs_mode_active = (npop == 1) ? 1 : 2;
    sfs_one_snp = one_snp;
    sfs_npop = npop;
    sfs_pop_sizes = pop_sizes;
    sfs_accumulator = accum;
    sfs_accum_len = accum_len;
    sfs_expected_mode = use_expected;

    /* Initialize R's RNG (for one.snp sampling) */
    GetRNGstate();

    /* Process each locus: parse command, run ms, accumulate SFS */
    for (int loc = 0; loc < nloci; loc++) {
        const char *cmd = CHAR(STRING_ELT(commands_sexp, loc));

        char **argv = NULL;
        int argc = parse_command(cmd, &argv);
        if (argc == 0) {
            PutRNGstate();
            free(accum);
            sfs_mode_active = 0;
            sfs_accumulator = NULL;
            sfs_pop_sizes = NULL;
            Rf_error("msABC_sfs_call: failed to parse command for locus %d", loc + 1);
        }

        /* Reset global state for clean run */
        msABC_reset_streec_statics();

        /* Initialize output capture (ms.c needs a valid stream) */
        msABC_init_output_stream();

        /* Set up error recovery */
        msABC_set_jmpbuf_active(1);
        int jmpval = setjmp(msABC_jmpbuf);

        if (jmpval != 0) {
            msABC_set_jmpbuf_active(0);
            msABC_close_output_stream();
            PutRNGstate();
            free(accum);
            free_argv(argc, argv);
            sfs_mode_active = 0;
            sfs_accumulator = NULL;
            sfs_pop_sizes = NULL;
            Rf_error("msABC_sfs_call: simulation error at locus %d (exit code %d)",
                     loc + 1, jmpval);
            return R_NilValue;
        }

        /* Run one locus (nreps=1). SFS accumulates via sfs_accumulate() hook in ms.c */
        msABC_main(argc, argv);

        msABC_set_jmpbuf_active(0);
        msABC_close_output_stream();
        free_argv(argc, argv);

        /* Check for user interrupt every 100 loci */
        if ((loc + 1) % 100 == 0) {
            R_CheckUserInterrupt();
        }
    }

    PutRNGstate();

    /* Create R result vector */
    SEXP result = PROTECT(allocVector(REALSXP, accum_len));
    memcpy(REAL(result), accum, accum_len * sizeof(double));

    /* Cleanup */
    free(accum);
    sfs_mode_active = 0;
    sfs_accumulator = NULL;
    sfs_pop_sizes = NULL;
    sfs_npop = 0;
    sfs_one_snp = 0;
    sfs_expected_mode = 0;

    UNPROTECT(1);
    return result;
}


/*
 * msABC_sfs_batch_call - .Call() entry point for batch SFS simulation
 *
 * Processes an entire block of simulations in a single C call, avoiding
 * repeated R↔C transitions and disk I/O for the locfile. Each simulation
 * uses a different ms command (with different demographic parameters) and
 * different per-locus mutation rates. The locfile on disk provides static
 * data (locus IDs, sample sizes, lengths); mu values are overridden in
 * memory via frag_mu_override.
 *
 * Arguments:
 *   commands_sexp:  character vector of length nsims (one ms command per sim)
 *   mu_rates_sexp:  numeric matrix, nrow = locfile rows, ncol = nsims
 *   pop_sizes_sexp: integer vector of per-population sample sizes
 *   one_snp_sexp:   logical, TRUE to sample one SNP per locus
 *
 * Returns: numeric matrix, nsims rows × accum_len columns (SFS per simulation)
 */

/* mu override globals in ms.c */
extern double *frag_mu_override;
extern int frag_mu_override_len;

SEXP msABC_sfs_batch_call(SEXP commands_sexp, SEXP mu_rates_sexp,
                           SEXP pop_sizes_sexp, SEXP one_snp_sexp,
                           SEXP method_sexp) {

    /* Validate inputs */
    if (!isString(commands_sexp) || length(commands_sexp) < 1) {
        Rf_error("msABC_sfs_batch_call: 'commands' must be a character vector");
    }
    if (!isReal(mu_rates_sexp) || !isMatrix(mu_rates_sexp)) {
        Rf_error("msABC_sfs_batch_call: 'mu_rates' must be a numeric matrix");
    }
    if (!isInteger(pop_sizes_sexp) || length(pop_sizes_sexp) < 1) {
        Rf_error("msABC_sfs_batch_call: 'pop_sizes' must be an integer vector");
    }
    if (!isLogical(one_snp_sexp) || length(one_snp_sexp) != 1) {
        Rf_error("msABC_sfs_batch_call: 'one_snp' must be a single logical value");
    }

    /* Parse method: "stochastic" (default) or "expected" */
    int use_expected = 0;
    if (!isNull(method_sexp) && isString(method_sexp) && length(method_sexp) == 1) {
        const char *method_str = CHAR(STRING_ELT(method_sexp, 0));
        if (strcmp(method_str, "expected") == 0) use_expected = 1;
    }

    int nsims = length(commands_sexp);
    int npop = length(pop_sizes_sexp);
    int *pop_sizes = INTEGER(pop_sizes_sexp);
    int one_snp = LOGICAL(one_snp_sexp)[0];

    /* mu_rates matrix dimensions */
    SEXP mu_dim = getAttrib(mu_rates_sexp, R_DimSymbol);
    int mu_nrow = INTEGER(mu_dim)[0];
    int mu_ncol = INTEGER(mu_dim)[1];
    double *mu_data = REAL(mu_rates_sexp);

    if (mu_ncol != nsims) {
        Rf_error("msABC_sfs_batch_call: mu_rates ncol (%d) != length(commands) (%d)",
                 mu_ncol, nsims);
    }

    /* Compute total sample size and SFS dimensions */
    int nsam = 0;
    for (int p = 0; p < npop; p++) {
        nsam += pop_sizes[p];
    }

    int accum_len;
    if (npop == 1) {
        accum_len = nsam - 1;  /* freq bins 1..nsam-1 */
    } else {
        accum_len = 1;
        for (int p = 0; p < npop; p++) {
            accum_len *= (pop_sizes[p] + 1);
        }
    }

    /* Allocate SFS accumulator */
    double *accum = (double *)calloc(accum_len, sizeof(double));
    if (accum == NULL) {
        Rf_error("msABC_sfs_batch_call: failed to allocate SFS accumulator");
    }

    /* Allocate result matrix: nsims rows × accum_len columns (column-major for R) */
    SEXP result = PROTECT(allocMatrix(REALSXP, nsims, accum_len));
    double *result_data = REAL(result);

    /* Set up global SFS mode */
    sfs_mode_active = (npop == 1) ? 1 : 2;
    sfs_one_snp = one_snp;
    sfs_npop = npop;
    sfs_pop_sizes = pop_sizes;
    sfs_accumulator = accum;
    sfs_accum_len = accum_len;
    sfs_expected_mode = use_expected;

    /* Initialize R's RNG */
    GetRNGstate();

    /* Process each simulation */
    for (int sim = 0; sim < nsims; sim++) {

        /* Set mu override for this simulation (column sim of mu_rates matrix) */
        frag_mu_override = mu_data + (long)sim * mu_nrow;
        frag_mu_override_len = mu_nrow;

        /* Reset SFS accumulator */
        memset(accum, 0, accum_len * sizeof(double));

        /* Parse command string */
        const char *cmd = CHAR(STRING_ELT(commands_sexp, sim));
        char **argv = NULL;
        int argc = parse_command(cmd, &argv);
        if (argc == 0) {
            frag_mu_override = NULL;
            frag_mu_override_len = 0;
            PutRNGstate();
            free(accum);
            sfs_mode_active = 0;
            sfs_accumulator = NULL;
            sfs_pop_sizes = NULL;
            UNPROTECT(1);
            Rf_error("msABC_sfs_batch_call: failed to parse command for sim %d", sim + 1);
        }

        /* Reset global state for clean run */
        msABC_reset_streec_statics();
        msABC_init_output_stream();

        /* Set up error recovery */
        msABC_set_jmpbuf_active(1);
        int jmpval = setjmp(msABC_jmpbuf);

        if (jmpval != 0) {
            msABC_set_jmpbuf_active(0);
            msABC_close_output_stream();
            frag_mu_override = NULL;
            frag_mu_override_len = 0;
            PutRNGstate();
            free(accum);
            free_argv(argc, argv);
            sfs_mode_active = 0;
            sfs_accumulator = NULL;
            sfs_pop_sizes = NULL;
            UNPROTECT(1);
            Rf_error("msABC_sfs_batch_call: simulation error at sim %d (exit code %d)",
                     sim + 1, jmpval);
            return R_NilValue;
        }

        /* Run simulation (fragment mode processes all loci internally) */
        msABC_main(argc, argv);

        msABC_set_jmpbuf_active(0);
        msABC_close_output_stream();
        free_argv(argc, argv);

        /* Copy SFS to result matrix (column-major: result[sim, col] = result_data[col * nsims + sim]) */
        for (int j = 0; j < accum_len; j++) {
            result_data[j * nsims + sim] = accum[j];
        }

        /* Check for user interrupt every 10 sims */
        if ((sim + 1) % 10 == 0) {
            R_CheckUserInterrupt();
        }
    }

    PutRNGstate();

    /* Cleanup */
    frag_mu_override = NULL;
    frag_mu_override_len = 0;
    free(accum);
    sfs_mode_active = 0;
    sfs_accumulator = NULL;
    sfs_pop_sizes = NULL;
    sfs_npop = 0;
    sfs_one_snp = 0;
    sfs_expected_mode = 0;

    UNPROTECT(1);
    return result;
}


/*
 * obs_sfs_call - .Call() entry point for observed SFS from alignment data
 *
 * Computes the SFS directly from nucleotide alignment matrices (character
 * matrices from read.phylip.loci()). Processes all loci in a single C call,
 * replacing the slow R per-locus pipeline (mat.snp.2ms -> parse -> tabulate).
 *
 * For each locus: scans alignment columns for variable sites (skipping any
 * column with gaps or N), determines the major allele, encodes minor alleles
 * as derived, and accumulates into the SFS.
 *
 * Arguments:
 *   loci_list:       R list of character matrices (nsam x nsites per locus)
 *   sample_idx_sexp: integer vector (0-based) mapping pop-ordered position
 *                    to row index in the alignment matrix
 *   pop_sizes_sexp:  integer vector of per-population sample sizes
 *   one_snp_sexp:    logical, TRUE to sample one SNP per locus
 *
 * Returns: numeric vector of SFS counts
 *   1-pop:     length nsam+1 (frequency bins 0..nsam)
 *   multi-pop: length prod(pop_sizes + 1) (flattened joint SFS)
 */
SEXP obs_sfs_call(SEXP loci_list, SEXP sample_idx_sexp,
                  SEXP pop_sizes_sexp, SEXP one_snp_sexp) {

    /* Validate inputs */
    if (!isNewList(loci_list)) {
        Rf_error("obs_sfs_call: 'loci_list' must be a list");
    }
    if (!isInteger(sample_idx_sexp)) {
        Rf_error("obs_sfs_call: 'sample_indices' must be an integer vector");
    }
    if (!isInteger(pop_sizes_sexp) || length(pop_sizes_sexp) < 1) {
        Rf_error("obs_sfs_call: 'pop_sizes' must be an integer vector");
    }
    if (!isLogical(one_snp_sexp) || length(one_snp_sexp) != 1) {
        Rf_error("obs_sfs_call: 'one_snp' must be a single logical value");
    }

    int nloci = length(loci_list);
    int npop = length(pop_sizes_sexp);
    int *pop_sizes = INTEGER(pop_sizes_sexp);
    int one_snp = LOGICAL(one_snp_sexp)[0];
    int *sample_idx = INTEGER(sample_idx_sexp);
    int nsam = length(sample_idx_sexp);

    /* Population boundaries (cumulative sizes) */
    int cum_sizes[npop + 1];
    cum_sizes[0] = 0;
    for (int p = 0; p < npop; p++) {
        cum_sizes[p + 1] = cum_sizes[p] + pop_sizes[p];
    }

    /* SFS dimensions */
    int accum_len;
    if (npop == 1) {
        accum_len = nsam + 1;  /* freq bins 0..nsam */
    } else {
        accum_len = 1;
        for (int p = 0; p < npop; p++) {
            accum_len *= (pop_sizes[p] + 1);
        }
    }

    /* Strides for multi-pop flat indexing (pop1 varies fastest) */
    int strides[npop];
    strides[0] = 1;
    for (int p = 1; p < npop; p++) {
        strides[p] = strides[p - 1] * (pop_sizes[p - 1] + 1);
    }

    double *accum = (double *)calloc(accum_len, sizeof(double));
    if (accum == NULL) {
        Rf_error("obs_sfs_call: failed to allocate SFS accumulator");
    }

    /* Temporary arrays for SNP column indices and major alleles */
    /* Allocated once, reused across loci */
    int snp_cap = 1024;
    int *snp_cols = (int *)malloc(snp_cap * sizeof(int));
    int *snp_major = (int *)malloc(snp_cap * sizeof(int));
    if (snp_cols == NULL || snp_major == NULL) {
        free(accum); free(snp_cols); free(snp_major);
        Rf_error("obs_sfs_call: memory allocation failed");
    }

    GetRNGstate();

    for (int loc = 0; loc < nloci; loc++) {
        SEXP mat = VECTOR_ELT(loci_list, loc);
        SEXP dim = getAttrib(mat, R_DimSymbol);
        if (isNull(dim) || length(dim) != 2) continue;

        int nrow = INTEGER(dim)[0];
        int ncol = INTEGER(dim)[1];
        if (ncol == 0) continue;

        /* ---- Pass 1: find segregating sites ---- */
        int n_snps = 0;

        for (int j = 0; j < ncol; j++) {
            /* Count nucleotide frequencies, skip column if gap or N found */
            int counts[4] = {0, 0, 0, 0};  /* a, c, g, t */
            int skip = 0;

            for (int ii = 0; ii < nsam; ii++) {
                int row = sample_idx[ii];
                const char *s = CHAR(STRING_ELT(mat, row + j * nrow));
                char c = s[0];
                switch (c) {
                    case 'a': case 'A': counts[0]++; break;
                    case 'c': case 'C': counts[1]++; break;
                    case 'g': case 'G': counts[2]++; break;
                    case 't': case 'T': counts[3]++; break;
                    default: skip = 1; break;  /* N, -, or other */
                }
                if (skip) break;
            }
            if (skip) continue;

            /* Check if variable: more than one non-zero count */
            int n_alleles = (counts[0] > 0) + (counts[1] > 0) +
                            (counts[2] > 0) + (counts[3] > 0);
            if (n_alleles < 2) continue;

            /* Major allele = most frequent base */
            int max_count = 0, max_idx = 0;
            for (int k = 0; k < 4; k++) {
                if (counts[k] > max_count) {
                    max_count = counts[k];
                    max_idx = k;
                }
            }

            /* Grow arrays if needed */
            if (n_snps >= snp_cap) {
                snp_cap *= 2;
                snp_cols = (int *)realloc(snp_cols, snp_cap * sizeof(int));
                snp_major = (int *)realloc(snp_major, snp_cap * sizeof(int));
                if (snp_cols == NULL || snp_major == NULL) {
                    free(accum);
                    PutRNGstate();
                    Rf_error("obs_sfs_call: realloc failed at locus %d", loc + 1);
                }
            }

            snp_cols[n_snps] = j;
            snp_major[n_snps] = max_idx;  /* 0=a, 1=c, 2=g, 3=t */
            n_snps++;
        }

        if (n_snps == 0) continue;

        /* ---- Pass 2: accumulate SFS ---- */
        /* Map base char to index: a=0, c=1, g=2, t=3 */
        int start_snp = 0, end_snp = n_snps;
        if (one_snp && n_snps > 1) {
            start_snp = sfs_rand_int(n_snps);
            end_snp = start_snp + 1;
        }

        if (npop == 1) {
            for (int si = start_snp; si < end_snp; si++) {
                int j = snp_cols[si];
                int major = snp_major[si];
                int freq = 0;
                for (int ii = 0; ii < nsam; ii++) {
                    int row = sample_idx[ii];
                    const char *s = CHAR(STRING_ELT(mat, row + j * nrow));
                    char c = s[0];
                    int base_idx;
                    switch (c) {
                        case 'a': case 'A': base_idx = 0; break;
                        case 'c': case 'C': base_idx = 1; break;
                        case 'g': case 'G': base_idx = 2; break;
                        default:            base_idx = 3; break;
                    }
                    if (base_idx != major) freq++;
                }
                if (freq >= 0 && freq <= nsam) {
                    accum[freq] += 1.0;
                }
            }
        } else {
            /* Joint SFS for multi-pop */
            for (int si = start_snp; si < end_snp; si++) {
                int j = snp_cols[si];
                int major = snp_major[si];
                int flat_idx = 0;
                for (int p = 0; p < npop; p++) {
                    int count = 0;
                    for (int ii = cum_sizes[p]; ii < cum_sizes[p + 1]; ii++) {
                        int row = sample_idx[ii];
                        const char *s = CHAR(STRING_ELT(mat, row + j * nrow));
                        char c = s[0];
                        int base_idx;
                        switch (c) {
                            case 'a': case 'A': base_idx = 0; break;
                            case 'c': case 'C': base_idx = 1; break;
                            case 'g': case 'G': base_idx = 2; break;
                            default:            base_idx = 3; break;
                        }
                        if (base_idx != major) count++;
                    }
                    flat_idx += count * strides[p];
                }
                accum[flat_idx] += 1.0;
            }
        }

        /* Check for user interrupt every 100 loci */
        if ((loc + 1) % 100 == 0) {
            R_CheckUserInterrupt();
        }
    }

    PutRNGstate();

    /* Create R result vector */
    SEXP result = PROTECT(allocVector(REALSXP, accum_len));
    memcpy(REAL(result), accum, accum_len * sizeof(double));

    /* Cleanup */
    free(accum);
    free(snp_cols);
    free(snp_major);

    UNPROTECT(1);
    return result;
}


/*
 * obs_sfs_vcf_call - .Call() entry point for observed SFS from VCF data
 *
 * Reads a VCF file and computes the SFS directly in C. Needed for
 * whole-genome VCFs with millions of SNPs. Biallelic only, complete
 * data only (skip sites with missing genotypes), major allele polarization.
 *
 * Arguments:
 *   vcf_path:         STRSXP - path to VCF file (plain text)
 *   sample_pop_map:   INTSXP - population number (1-based) for each target sample
 *   sample_col_indices: INTSXP - 0-based column index from first sample column
 *   pop_sizes:        INTSXP - haploid sample size per population
 *   npop_sexp:        INTSXP - number of populations
 *
 * Returns: numeric vector of SFS counts
 *   1-pop:     length nsam+1 (frequency bins 0..nsam, bin 0 starts at 0)
 *   multi-pop: length prod(pop_sizes + 1)
 */
SEXP obs_sfs_vcf_call(SEXP vcf_path, SEXP sample_pop_map, SEXP sample_col_indices,
                       SEXP pop_sizes_sexp, SEXP npop_sexp) {

    /* Validate inputs */
    if (!isString(vcf_path) || length(vcf_path) != 1) {
        Rf_error("obs_sfs_vcf_call: 'vcf_path' must be a single string");
    }
    if (!isInteger(sample_pop_map)) {
        Rf_error("obs_sfs_vcf_call: 'sample_pop_map' must be an integer vector");
    }
    if (!isInteger(sample_col_indices)) {
        Rf_error("obs_sfs_vcf_call: 'sample_col_indices' must be an integer vector");
    }
    if (!isInteger(pop_sizes_sexp) || length(pop_sizes_sexp) < 1) {
        Rf_error("obs_sfs_vcf_call: 'pop_sizes' must be an integer vector");
    }
    if (!isInteger(npop_sexp) || length(npop_sexp) != 1) {
        Rf_error("obs_sfs_vcf_call: 'npop' must be a single integer");
    }

    const char *filepath = CHAR(STRING_ELT(vcf_path, 0));
    int n_samples = length(sample_pop_map);
    int *pop_map = INTEGER(sample_pop_map);       /* 1-based pop number per sample */
    int *col_indices = INTEGER(sample_col_indices); /* 0-based from first sample col */
    int npop = INTEGER(npop_sexp)[0];
    int *pop_sizes = INTEGER(pop_sizes_sexp);

    /* Compute total haploid sample size */
    int nsam = 0;
    for (int p = 0; p < npop; p++) {
        nsam += pop_sizes[p];
    }

    /* SFS dimensions */
    int accum_len;
    if (npop == 1) {
        accum_len = nsam + 1;  /* freq bins 0..nsam */
    } else {
        accum_len = 1;
        for (int p = 0; p < npop; p++) {
            accum_len *= (pop_sizes[p] + 1);
        }
    }

    /* Strides for multi-pop flat indexing (pop1 varies fastest) */
    int *strides = (int *)malloc(npop * sizeof(int));
    if (strides == NULL) Rf_error("obs_sfs_vcf_call: memory allocation failed");
    strides[0] = 1;
    for (int p = 1; p < npop; p++) {
        strides[p] = strides[p - 1] * (pop_sizes[p - 1] + 1);
    }

    double *accum = (double *)calloc(accum_len, sizeof(double));
    if (accum == NULL) {
        free(strides);
        Rf_error("obs_sfs_vcf_call: failed to allocate SFS accumulator");
    }

    /* Per-population allele count workspace */
    int *pop_alt_counts = (int *)calloc(npop, sizeof(int));
    if (pop_alt_counts == NULL) {
        free(strides); free(accum);
        Rf_error("obs_sfs_vcf_call: memory allocation failed");
    }

    /* Find the maximum column index to know how many fields to track */
    int max_col_idx = 0;
    for (int s = 0; s < n_samples; s++) {
        if (col_indices[s] > max_col_idx) max_col_idx = col_indices[s];
    }

    /* Open VCF file */
    FILE *fp = fopen(filepath, "r");
    if (fp == NULL) {
        free(strides); free(accum); free(pop_alt_counts);
        Rf_error("obs_sfs_vcf_call: cannot open VCF file '%s'", filepath);
    }

    /* Read line buffer */
    int buf_cap = 65536;
    char *buf = (char *)malloc(buf_cap);
    if (buf == NULL) {
        fclose(fp); free(strides); free(accum); free(pop_alt_counts);
        Rf_error("obs_sfs_vcf_call: memory allocation failed");
    }

    long lines_processed = 0;

    /* Process file line by line */
    while (1) {
        /* Read a line, handling lines longer than buffer */
        int line_len = 0;
        int c;
        while ((c = fgetc(fp)) != EOF && c != '\n') {
            if (line_len + 1 >= buf_cap) {
                buf_cap *= 2;
                char *new_buf = (char *)realloc(buf, buf_cap);
                if (new_buf == NULL) {
                    fclose(fp); free(buf); free(strides); free(accum); free(pop_alt_counts);
                    Rf_error("obs_sfs_vcf_call: realloc failed");
                }
                buf = new_buf;
            }
            buf[line_len++] = (char)c;
        }
        if (line_len == 0 && c == EOF) break;
        buf[line_len] = '\0';

        /* Skip meta-information lines (##) */
        if (line_len >= 2 && buf[0] == '#' && buf[1] == '#') continue;

        /* Skip header line (#CHROM) */
        if (line_len >= 1 && buf[0] == '#') continue;

        /* ---- Data line ---- */
        /* Fields: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample1 sample2 ... */
        /* Tab-delimited. We need fields 3(REF), 4(ALT), and 8(FORMAT) + sample columns */

        /* Parse by walking through tab-separated fields */
        int field = 0;       /* 0-based field index */
        int pos = 0;
        int skip_site = 0;

        /* Field pointers */
        char *ref_start = NULL;
        int ref_len = 0;
        char *alt_start = NULL;
        int alt_len = 0;
        int format_gt_idx = -1;  /* which sub-field of FORMAT is GT */
        int first_sample_field = 9;  /* VCF standard: samples start at field 9 */

        /* Reset per-pop counts */
        memset(pop_alt_counts, 0, npop * sizeof(int));

        int total_alt = 0;
        int total_ref = 0;

        while (pos <= line_len && !skip_site) {
            /* Find next tab or end of line */
            int field_start = pos;
            while (pos < line_len && buf[pos] != '\t') pos++;
            int field_end = pos;
            pos++;  /* skip tab */

            if (field == 3) {
                /* REF field */
                ref_start = buf + field_start;
                ref_len = field_end - field_start;
                /* Single base REF only */
                if (ref_len != 1) { skip_site = 1; break; }
            } else if (field == 4) {
                /* ALT field */
                alt_start = buf + field_start;
                alt_len = field_end - field_start;
                /* Skip multi-allelic (comma in ALT) or missing ALT (.) */
                if (alt_len != 1 || alt_start[0] == '.') { skip_site = 1; break; }
                for (int k = 0; k < alt_len; k++) {
                    if (alt_start[k] == ',') { skip_site = 1; break; }
                }
                if (skip_site) break;
            } else if (field == 8) {
                /* FORMAT field - find GT position */
                /* GT is usually the first sub-field, but check */
                int sub_field = 0;
                int k = field_start;
                int gt_found = 0;
                while (k <= field_end) {
                    int sf_start = k;
                    while (k < field_end && buf[k] != ':') k++;
                    int sf_len = k - sf_start;
                    if (sf_len == 2 && buf[sf_start] == 'G' && buf[sf_start + 1] == 'T') {
                        format_gt_idx = sub_field;
                        gt_found = 1;
                        break;
                    }
                    k++;  /* skip ':' */
                    sub_field++;
                }
                if (!gt_found) { skip_site = 1; break; }
            } else if (field >= first_sample_field) {
                /* Sample column - check if this is a target sample */
                int sample_col = field - first_sample_field;  /* 0-based */

                /* Check if this sample_col is in our target list */
                /* For efficiency, check all target samples that match this column */
                int is_target = 0;
                for (int s = 0; s < n_samples; s++) {
                    if (col_indices[s] == sample_col) {
                        is_target = 1;
                        break;
                    }
                }

                if (is_target) {
                    /* Extract GT sub-field */
                    int sub_field = 0;
                    int k = field_start;
                    char *gt_start = NULL;
                    int gt_len = 0;

                    while (k <= field_end) {
                        int sf_start = k;
                        while (k < field_end && buf[k] != ':') k++;
                        if (sub_field == format_gt_idx) {
                            gt_start = buf + sf_start;
                            gt_len = k - sf_start;
                            break;
                        }
                        k++;
                        sub_field++;
                    }

                    if (gt_start == NULL || gt_len < 3) { skip_site = 1; break; }

                    /* Parse GT: expect X/Y or X|Y where X,Y are single digits */
                    char a1 = gt_start[0];
                    char sep = gt_start[1];
                    char a2 = gt_start[2];

                    /* Check for missing genotype */
                    if (a1 == '.' || a2 == '.') { skip_site = 1; break; }

                    /* Validate separator */
                    if (sep != '/' && sep != '|') { skip_site = 1; break; }

                    int allele1 = a1 - '0';
                    int allele2 = a2 - '0';

                    /* Count ALT alleles for each target sample mapping to this column */
                    for (int s = 0; s < n_samples; s++) {
                        if (col_indices[s] == sample_col) {
                            int pop_idx = pop_map[s] - 1;  /* convert 1-based to 0-based */
                            pop_alt_counts[pop_idx] += allele1 + allele2;
                            total_alt += allele1 + allele2;
                            total_ref += (1 - allele1) + (1 - allele2);
                        }
                    }
                }
            }

            field++;

            /* Early exit once we've passed all needed columns */
            if (field > first_sample_field + max_col_idx) break;
        }

        if (skip_site) continue;

        /* Make sure we processed enough fields */
        if (field <= first_sample_field + max_col_idx) continue;

        /* Polarize by major allele: if ALT count > REF count, flip */
        if (total_alt > total_ref) {
            for (int p = 0; p < npop; p++) {
                pop_alt_counts[p] = pop_sizes[p] - pop_alt_counts[p];
            }
        }

        /* Accumulate into SFS */
        if (npop == 1) {
            int freq = pop_alt_counts[0];
            if (freq >= 0 && freq <= nsam) {
                accum[freq] += 1.0;
            }
        } else {
            int flat_idx = 0;
            for (int p = 0; p < npop; p++) {
                flat_idx += pop_alt_counts[p] * strides[p];
            }
            if (flat_idx >= 0 && flat_idx < accum_len) {
                accum[flat_idx] += 1.0;
            }
        }

        lines_processed++;
        if (lines_processed % 100000 == 0) {
            R_CheckUserInterrupt();
        }
    }

    fclose(fp);

    /* Create R result vector */
    SEXP result = PROTECT(allocVector(REALSXP, accum_len));
    memcpy(REAL(result), accum, accum_len * sizeof(double));

    /* Cleanup */
    free(buf);
    free(strides);
    free(accum);
    free(pop_alt_counts);

    UNPROTECT(1);
    return result;
}
