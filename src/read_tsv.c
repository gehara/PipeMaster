/*
 * read_tsv.c — Fast tab-separated numeric file reader for PipeMaster
 *
 * .Call("read_tsv_call", filename, col_indices, skip, nrows)
 *
 * Reads a tab-separated file and returns a numeric matrix containing only
 * the requested columns. Much faster than R's read.table() for wide files
 * (e.g., 69K-column SFS reftables) because it:
 *   - Parses doubles directly with strtod (no intermediate strings)
 *   - Skips unrequested columns without allocating memory
 *   - Reads in a single pass with no R-level string operations
 *
 * Parameters:
 *   filename   - STRSXP: path to the input file
 *   col_indices - INTSXP: 1-based column indices to extract
 *   skip       - INTSXP scalar: number of header lines to skip (usually 1)
 *   nrows      - INTSXP scalar: max rows to read (-1 = all)
 *
 * Returns: REALSXP matrix (nrows x length(col_indices))
 */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Initial line buffer size — grown as needed */
#define INIT_BUF_SIZE (1 << 20)  /* 1 MB */

SEXP read_tsv_call(SEXP filename_sexp, SEXP col_indices_sexp,
                    SEXP skip_sexp, SEXP nrows_sexp) {

    const char *filename = CHAR(STRING_ELT(filename_sexp, 0));
    int *col_idx = INTEGER(col_indices_sexp);
    int n_cols = length(col_indices_sexp);
    int skip = asInteger(skip_sexp);
    int max_rows = asInteger(nrows_sexp);

    FILE *fp = fopen(filename, "r");
    if (!fp)
        Rf_error("Cannot open file: %s", filename);

    /* Allocate line buffer */
    size_t buf_size = INIT_BUF_SIZE;
    char *buf = (char *)malloc(buf_size);
    if (!buf) { fclose(fp); Rf_error("malloc failed"); }

    /* Skip header lines */
    for (int i = 0; i < skip; i++) {
        if (getline(&buf, &buf_size, fp) == -1) {
            free(buf); fclose(fp);
            Rf_error("File has fewer than %d lines", skip);
        }
    }

    /* Find max column index to know when to stop scanning each line */
    int max_col = 0;
    for (int i = 0; i < n_cols; i++) {
        if (col_idx[i] > max_col) max_col = col_idx[i];
    }

    /* Build a lookup: for each 1-based column, which output column does it map to?
     * -1 means skip this column */
    int *col_map = (int *)calloc(max_col + 1, sizeof(int));
    if (!col_map) { free(buf); fclose(fp); Rf_error("calloc failed"); }
    for (int i = 0; i <= max_col; i++) col_map[i] = -1;
    for (int i = 0; i < n_cols; i++) col_map[col_idx[i]] = i;

    /* First pass: count rows (if max_rows == -1) */
    int n_rows;
    if (max_rows >= 0) {
        n_rows = max_rows;
    } else {
        long start_pos = ftell(fp);
        n_rows = 0;
        while (getline(&buf, &buf_size, fp) != -1) {
            if (buf[0] != '\0' && buf[0] != '\n') n_rows++;
        }
        fseek(fp, start_pos, SEEK_SET);
    }

    /* Allocate output matrix */
    SEXP result = PROTECT(allocMatrix(REALSXP, n_rows, n_cols));
    double *out = REAL(result);
    /* Initialize to NA */
    for (R_xlen_t i = 0; i < (R_xlen_t)n_rows * n_cols; i++)
        out[i] = NA_REAL;

    /* Second pass: read and parse */
    int row = 0;
    while (row < n_rows && getline(&buf, &buf_size, fp) != -1) {
        if (buf[0] == '\0' || buf[0] == '\n') continue;

        char *p = buf;
        int col = 1;  /* 1-based */

        while (*p && col <= max_col) {
            /* Find end of field */
            char *field_start = p;
            while (*p && *p != '\t' && *p != '\n' && *p != '\r') p++;

            if (col_map[col] >= 0) {
                /* Parse this field */
                char saved = *p;
                *p = '\0';
                char *endptr;
                double val = strtod(field_start, &endptr);
                if (endptr == field_start) val = NA_REAL;
                /* Column-major storage: out[row + col_out * n_rows] */
                out[row + (R_xlen_t)col_map[col] * n_rows] = val;
                *p = saved;
            }

            /* Advance past delimiter */
            if (*p == '\t') p++;
            else if (*p == '\n' || *p == '\r') break;
            col++;
        }
        row++;

        /* Progress every 10000 rows */
        if (row % 10000 == 0) R_CheckUserInterrupt();
    }

    /* If we read fewer rows than allocated, trim */
    if (row < n_rows) {
        SEXP trimmed = PROTECT(allocMatrix(REALSXP, row, n_cols));
        double *trim_out = REAL(trimmed);
        for (int j = 0; j < n_cols; j++)
            memcpy(trim_out + (R_xlen_t)j * row,
                   out + (R_xlen_t)j * n_rows,
                   row * sizeof(double));
        free(col_map);
        free(buf);
        fclose(fp);
        UNPROTECT(2);
        return trimmed;
    }

    free(col_map);
    free(buf);
    fclose(fp);
    UNPROTECT(1);
    return result;
}
