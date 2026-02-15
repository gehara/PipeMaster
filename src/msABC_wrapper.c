#include <R.h>
#include <Rinternals.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "msABC_capture.h"

/* Jump buffer defined in ms.c */
extern jmp_buf msABC_jmpbuf;

/* ---- Output capture using open_memstream ---- */

static FILE *capture_stream = NULL;
static char *capture_buf = NULL;
static size_t capture_len = 0;

void msABC_init_output_stream(void) {
    if (capture_stream != NULL) {
        fclose(capture_stream);
        capture_stream = NULL;
    }
    if (capture_buf != NULL) {
        free(capture_buf);
        capture_buf = NULL;
    }
    capture_len = 0;
    capture_stream = open_memstream(&capture_buf, &capture_len);
    if (capture_stream == NULL) {
        Rf_error("msABC: failed to open output capture stream");
    }
}

FILE *msABC_get_output_stream(void) {
    if (capture_stream == NULL) {
        msABC_init_output_stream();
    }
    return capture_stream;
}

char *msABC_get_output_buffer(size_t *len) {
    if (capture_stream != NULL) {
        fflush(capture_stream);
    }
    if (len != NULL) *len = capture_len;
    return capture_buf;
}

void msABC_close_output_stream(void) {
    if (capture_stream != NULL) {
        fclose(capture_stream);
        capture_stream = NULL;
    }
    /* capture_buf is freed by the caller after use */
}

/* ---- Parse command string into argc/argv ---- */

static int parse_command(const char *cmd, char ***argv_out) {
    /* Count tokens first */
    char *tmp = strdup(cmd);
    if (tmp == NULL) return 0;

    int argc = 0;
    char *token = strtok(tmp, " \t");
    while (token != NULL) {
        argc++;
        token = strtok(NULL, " \t");
    }
    free(tmp);

    /* Build argv: prepend "msABC" as argv[0] */
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

/* ---- .Call entry point ---- */

SEXP msABC_call(SEXP command_sexp, SEXP seed_sexp) {
    /* Validate command argument */
    if (!isString(command_sexp) || length(command_sexp) != 1) {
        Rf_error("msABC_call: 'command' must be a single character string");
    }

    const char *cmd = CHAR(STRING_ELT(command_sexp, 0));

    /* Set seed if provided */
    if (!isNull(seed_sexp)) {
        if (!isInteger(seed_sexp) || length(seed_sexp) != 3) {
            Rf_error("msABC_call: 'seed' must be an integer vector of length 3, or NULL");
        }
        int *sv = INTEGER(seed_sexp);
        unsigned short seedv[3];
        seedv[0] = (unsigned short)sv[0];
        seedv[1] = (unsigned short)sv[1];
        seedv[2] = (unsigned short)sv[2];
        seedit_r(seedv);
    }

    /* Parse command into argc/argv */
    char **argv = NULL;
    int argc = parse_command(cmd, &argv);
    if (argc == 0) {
        Rf_error("msABC_call: failed to parse command string");
    }

    /* Reset global state for clean run */
    msABC_reset_streec_statics();

    /* Initialize output capture */
    msABC_init_output_stream();

    /* Set up error recovery via longjmp */
    msABC_set_jmpbuf_active(1);
    int jmpval = setjmp(msABC_jmpbuf);

    if (jmpval != 0) {
        /* We got here via longjmp from EXIT_MSABC */
        msABC_set_jmpbuf_active(0);
        msABC_close_output_stream();
        if (capture_buf != NULL) {
            free(capture_buf);
            capture_buf = NULL;
        }
        capture_len = 0;
        free_argv(argc, argv);
        Rf_error("msABC: simulation error (exit code %d)", jmpval);
        return R_NilValue; /* not reached */
    }

    /* Run the simulation */
    msABC_main(argc, argv);

    /* Deactivate longjmp */
    msABC_set_jmpbuf_active(0);

    /* Retrieve output */
    size_t outlen = 0;
    char *output = msABC_get_output_buffer(&outlen);

    /* Create R string result */
    SEXP result = PROTECT(allocVector(STRSXP, 1));
    if (output != NULL && outlen > 0) {
        SET_STRING_ELT(result, 0, mkCharLen(output, (int)outlen));
    } else {
        SET_STRING_ELT(result, 0, mkChar(""));
    }

    /* Cleanup */
    msABC_close_output_stream();
    if (capture_buf != NULL) {
        free(capture_buf);
        capture_buf = NULL;
    }
    capture_len = 0;
    free_argv(argc, argv);

    UNPROTECT(1);
    return result;
}


/*
 * msABC_batch_call - .Call() entry point for batch msABC simulation
 *
 * Processes a block of simulations in a single C call. Each simulation
 * uses a different ms command (with different demographic parameters).
 * Per-locus mutation and recombination rates are passed directly as
 * matrices, avoiding per-iteration locfile disk writes.
 *
 * Arguments:
 *   commands_sexp:   character vector of length nsims
 *   mu_rates_sexp:   numeric matrix (locfile_rows x nsims), or NULL
 *   rec_rates_sexp:  numeric matrix (locfile_rows x nsims), or NULL
 *
 * Returns: character vector of length nsims (text output per simulation)
 */

/* Override globals in ms.c */
extern double *frag_mu_override;
extern int frag_mu_override_len;
extern double *frag_rec_override;
extern int frag_rec_override_len;

SEXP msABC_batch_call(SEXP commands_sexp, SEXP mu_rates_sexp,
                       SEXP rec_rates_sexp) {

    /* Validate inputs */
    if (!isString(commands_sexp) || length(commands_sexp) < 1) {
        Rf_error("msABC_batch_call: 'commands' must be a character vector");
    }

    int nsims = length(commands_sexp);

    /* mu_rates matrix (optional) */
    int has_mu = !isNull(mu_rates_sexp);
    double *mu_data = NULL;
    int mu_nrow = 0;
    if (has_mu) {
        if (!isReal(mu_rates_sexp) || !isMatrix(mu_rates_sexp)) {
            Rf_error("msABC_batch_call: 'mu_rates' must be a numeric matrix or NULL");
        }
        SEXP mu_dim = getAttrib(mu_rates_sexp, R_DimSymbol);
        mu_nrow = INTEGER(mu_dim)[0];
        int mu_ncol = INTEGER(mu_dim)[1];
        if (mu_ncol != nsims) {
            Rf_error("msABC_batch_call: mu_rates ncol (%d) != nsims (%d)", mu_ncol, nsims);
        }
        mu_data = REAL(mu_rates_sexp);
    }

    /* rec_rates matrix (optional) */
    int has_rec = !isNull(rec_rates_sexp);
    double *rec_data = NULL;
    int rec_nrow = 0;
    if (has_rec) {
        if (!isReal(rec_rates_sexp) || !isMatrix(rec_rates_sexp)) {
            Rf_error("msABC_batch_call: 'rec_rates' must be a numeric matrix or NULL");
        }
        SEXP rec_dim = getAttrib(rec_rates_sexp, R_DimSymbol);
        rec_nrow = INTEGER(rec_dim)[0];
        int rec_ncol = INTEGER(rec_dim)[1];
        if (rec_ncol != nsims) {
            Rf_error("msABC_batch_call: rec_rates ncol (%d) != nsims (%d)", rec_ncol, nsims);
        }
        rec_data = REAL(rec_rates_sexp);
    }

    /* Allocate result: character vector of text outputs */
    SEXP result = PROTECT(allocVector(STRSXP, nsims));

    for (int sim = 0; sim < nsims; sim++) {

        /* Set mu/rec overrides for this simulation */
        if (has_mu) {
            frag_mu_override = mu_data + (long)sim * mu_nrow;
            frag_mu_override_len = mu_nrow;
        }
        if (has_rec) {
            frag_rec_override = rec_data + (long)sim * rec_nrow;
            frag_rec_override_len = rec_nrow;
        }

        /* Parse command */
        const char *cmd = CHAR(STRING_ELT(commands_sexp, sim));
        char **argv = NULL;
        int argc = parse_command(cmd, &argv);
        if (argc == 0) {
            frag_mu_override = NULL;
            frag_mu_override_len = 0;
            frag_rec_override = NULL;
            frag_rec_override_len = 0;
            UNPROTECT(1);
            Rf_error("msABC_batch_call: failed to parse command for sim %d", sim + 1);
        }

        /* Reset and run */
        msABC_reset_streec_statics();
        msABC_init_output_stream();
        msABC_set_jmpbuf_active(1);
        int jmpval = setjmp(msABC_jmpbuf);

        if (jmpval != 0) {
            msABC_set_jmpbuf_active(0);
            msABC_close_output_stream();
            if (capture_buf != NULL) {
                free(capture_buf);
                capture_buf = NULL;
            }
            capture_len = 0;
            frag_mu_override = NULL;
            frag_mu_override_len = 0;
            frag_rec_override = NULL;
            frag_rec_override_len = 0;
            free_argv(argc, argv);
            UNPROTECT(1);
            Rf_error("msABC_batch_call: simulation error at sim %d (exit code %d)",
                     sim + 1, jmpval);
            return R_NilValue;
        }

        msABC_main(argc, argv);
        msABC_set_jmpbuf_active(0);

        /* Capture output */
        size_t outlen = 0;
        char *output = msABC_get_output_buffer(&outlen);
        if (output != NULL && outlen > 0) {
            SET_STRING_ELT(result, sim, mkCharLen(output, (int)outlen));
        } else {
            SET_STRING_ELT(result, sim, mkChar(""));
        }

        msABC_close_output_stream();
        if (capture_buf != NULL) {
            free(capture_buf);
            capture_buf = NULL;
        }
        capture_len = 0;
        free_argv(argc, argv);

        /* Check for user interrupt every 10 sims */
        if ((sim + 1) % 10 == 0) {
            R_CheckUserInterrupt();
        }
    }

    /* Clear overrides */
    frag_mu_override = NULL;
    frag_mu_override_len = 0;
    frag_rec_override = NULL;
    frag_rec_override_len = 0;

    UNPROTECT(1);
    return result;
}
