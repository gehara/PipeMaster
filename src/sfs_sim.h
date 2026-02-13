#ifndef SFS_SIM_H
#define SFS_SIM_H

/* SFS accumulator mode for native SFS simulation.
 * When active, ms.c accumulates SFS counts after each replicate
 * instead of (in addition to) printing text output.
 */

/* Global state for SFS accumulation */
extern int sfs_mode_active;      /* 0 = off (normal ms mode), 1 = 1D SFS, 2 = joint SFS */
extern int sfs_one_snp;          /* 1 = sample one SNP per locus */
extern int sfs_npop;             /* number of populations */
extern int *sfs_pop_sizes;       /* per-population sample sizes (config) */
extern double *sfs_accumulator;  /* accumulated SFS counts */
extern int sfs_accum_len;        /* length of accumulator array */
extern double *sfs_theta_array;  /* per-locus theta values (NULL = use command theta) */
extern int sfs_nloci;            /* total number of loci (= nreps) */
extern int sfs_expected_mode;    /* 0 = stochastic (default), 1 = expected SFS from tree */

/* Called from ms.c after each replicate when sfs_mode_active > 0 */
void sfs_accumulate(char **list, int nsam, int segsites, int *config, int npop);

#ifdef R_INTERNALS_H_
/* .Call() entry point for batch SFS simulation (eliminates per-sim disk I/O) */
SEXP msABC_sfs_batch_call(SEXP commands_sexp, SEXP mu_rates_sexp,
                           SEXP pop_sizes_sexp, SEXP one_snp_sexp,
                           SEXP method_sexp);

/* .Call() entry point for observed SFS from alignment character matrices */
SEXP obs_sfs_call(SEXP loci_list, SEXP sample_idx_sexp,
                  SEXP pop_sizes_sexp, SEXP one_snp_sexp);
#endif

#endif /* SFS_SIM_H */
