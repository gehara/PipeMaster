#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/* Declared in msABC_wrapper.c */
extern SEXP msABC_call(SEXP command_sexp, SEXP seed_sexp);
extern SEXP msABC_batch_call(SEXP commands_sexp, SEXP mu_rates_sexp,
                              SEXP rec_rates_sexp);

/* Declared in sfs_sim.c */
extern SEXP msABC_sfs_call(SEXP command_sexp, SEXP pop_sizes_sexp,
                            SEXP one_snp_sexp, SEXP seed_sexp,
                            SEXP method_sexp);
extern SEXP msABC_sfs_batch_call(SEXP commands_sexp, SEXP mu_rates_sexp,
                                  SEXP pop_sizes_sexp, SEXP one_snp_sexp,
                                  SEXP method_sexp);
extern SEXP obs_sfs_call(SEXP loci_list, SEXP sample_idx_sexp,
                          SEXP pop_sizes_sexp, SEXP one_snp_sexp);
extern SEXP obs_sfs_vcf_call(SEXP vcf_path, SEXP sample_pop_map,
                              SEXP sample_col_indices, SEXP pop_sizes_sexp,
                              SEXP npop_sexp);

static const R_CallMethodDef CallEntries[] = {
    {"msABC_call",           (DL_FUNC) &msABC_call,           2},
    {"msABC_batch_call",     (DL_FUNC) &msABC_batch_call,     3},
    {"msABC_sfs_call",       (DL_FUNC) &msABC_sfs_call,       5},
    {"msABC_sfs_batch_call", (DL_FUNC) &msABC_sfs_batch_call, 5},
    {"obs_sfs_call",         (DL_FUNC) &obs_sfs_call,         4},
    {"obs_sfs_vcf_call",     (DL_FUNC) &obs_sfs_vcf_call,     5},
    {NULL, NULL, 0}
};

void R_init_PipeMaster(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
