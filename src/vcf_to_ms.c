/* vcf_to_ms.c — Fast VCF-to-ms format converter for PipeMaster
 *
 * Converts a multi-locus VCF file (one CHROM per locus) to concatenated
 * ms format.  Called from R via .Call("vcf_to_ms_call", ...).
 *
 * Expects CHROM names like "locus_N" (1-based).  Loci with no variants
 * in the VCF are emitted as monomorphic ms blocks ("segsites: 0").
 *
 * The haplotype matrix is polarized so the major allele is 0.
 */

#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Helpers                                                             */
/* ------------------------------------------------------------------ */

/* Read one line from fp into *buf (null-terminated).
 * Grows *buf / *cap as needed.  Returns line length or -1 on EOF. */
static int vcf_read_line(FILE *fp, char **buf, int *cap) {
    int len = 0;
    int c;
    while ((c = fgetc(fp)) != EOF && c != '\n') {
        if (len + 1 >= *cap) {
            *cap *= 2;
            char *nb = (char *)realloc(*buf, *cap);
            if (!nb) Rf_error("vcf_to_ms_call: realloc failed");
            *buf = nb;
        }
        (*buf)[len++] = (char)c;
    }
    if (len == 0 && c == EOF) return -1;
    (*buf)[len] = '\0';
    return len;
}

/* Write one ms locus block.
 * var_count == 0  →  monomorphic block (positions/geno ignored). */
static void write_ms_locus(FILE *fp, int var_count, const int *positions,
                           signed char *geno, int n_hap, int contig_len) {
    int v, h;

    if (var_count == 0) {
        fputs("segsites: 0\npositions:\n//\n", fp);
        return;
    }

    /* Polarize: major allele → 0 */
    for (v = 0; v < var_count; v++) {
        int sum = 0;
        signed char *row = geno + (long)v * n_hap;
        for (h = 0; h < n_hap; h++) sum += row[h];
        if (sum * 2 > n_hap) {
            for (h = 0; h < n_hap; h++) row[h] = 1 - row[h];
        }
    }

    /* segsites + positions */
    fprintf(fp, "segsites: %d\npositions:", var_count);
    for (v = 0; v < var_count; v++)
        fprintf(fp, " %.10g", (double)positions[v] / contig_len);
    fputc('\n', fp);

    /* Haplotype strings */
    {
        char *hs = (char *)malloc(var_count + 2);
        if (!hs) Rf_error("vcf_to_ms_call: malloc failed");
        for (h = 0; h < n_hap; h++) {
            for (v = 0; v < var_count; v++)
                hs[v] = '0' + geno[(long)v * n_hap + h];
            hs[var_count] = '\0';
            fputs(hs, fp);
            fputc('\n', fp);
        }
        free(hs);
    }
    fputs("//\n", fp);
}

/* Extract the trailing integer from a CHROM name like "locus_123".
 * Returns 0 if no underscore + digits found. */
static int chrom_to_locus_num(const char *buf, int len) {
    int k;
    for (k = len - 1; k >= 0; k--) {
        if (buf[k] == '_') return atoi(buf + k + 1);
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Main entry point                                                    */
/* ------------------------------------------------------------------ */

/*
 * vcf_to_ms_call(vcf_path, sample_names, contig_length, n_loci, output_path)
 *
 *   vcf_path      : STRSXP(1) — path to input VCF
 *   sample_names  : STRSXP(n) — sample names in desired output order
 *                               (sorted by population in R before calling)
 *   contig_length : INTSXP(1) — bp length for position scaling
 *   n_loci        : INTSXP(1) — total expected loci (for trailing mono blocks)
 *   output_path   : STRSXP(1) — path to output ms file
 *
 * Returns: INTSXP(1) — number of loci written.
 */
SEXP vcf_to_ms_call(SEXP vcf_path_sexp, SEXP sample_names_sexp,
                     SEXP contig_len_sexp, SEXP n_loci_sexp,
                     SEXP output_path_sexp) {

    int s, v;

    /* ---- Validate ---- */
    if (!isString(vcf_path_sexp) || length(vcf_path_sexp) != 1)
        Rf_error("vcf_to_ms_call: 'vcf_path' must be a single string");
    if (!isString(sample_names_sexp))
        Rf_error("vcf_to_ms_call: 'sample_names' must be a character vector");
    if (!isInteger(contig_len_sexp) || length(contig_len_sexp) != 1)
        Rf_error("vcf_to_ms_call: 'contig_length' must be a single integer");
    if (!isInteger(n_loci_sexp) || length(n_loci_sexp) != 1)
        Rf_error("vcf_to_ms_call: 'n_loci' must be a single integer");
    if (!isString(output_path_sexp) || length(output_path_sexp) != 1)
        Rf_error("vcf_to_ms_call: 'output_path' must be a single string");

    const char *vcf_path    = CHAR(STRING_ELT(vcf_path_sexp, 0));
    int         n_diploid   = length(sample_names_sexp);
    int         n_hap       = 2 * n_diploid;
    int         contig_len  = INTEGER(contig_len_sexp)[0];
    int         n_loci      = INTEGER(n_loci_sexp)[0];
    const char *output_path = CHAR(STRING_ELT(output_path_sexp, 0));

    /* ---- Open VCF ---- */
    FILE *fp_in = fopen(vcf_path, "r");
    if (!fp_in)
        Rf_error("vcf_to_ms_call: cannot open VCF '%s'", vcf_path);
    setvbuf(fp_in, NULL, _IOFBF, 1 << 20);  /* 1 MB read buffer */

    int buf_cap = 65536;
    char *buf = (char *)malloc(buf_cap);
    if (!buf) { fclose(fp_in); Rf_error("malloc failed"); }

    /* ---- Parse #CHROM header → find sample column indices ---- */
    /* sample_col[s] = absolute field index in VCF (0-based) for sample s */
    int *sample_col = (int *)calloc(n_diploid, sizeof(int));
    if (!sample_col) { free(buf); fclose(fp_in); Rf_error("malloc failed"); }
    for (s = 0; s < n_diploid; s++) sample_col[s] = -1;

    int header_found = 0;
    while (1) {
        int len = vcf_read_line(fp_in, &buf, &buf_cap);
        if (len < 0) break;
        if (len >= 2 && buf[0] == '#' && buf[1] == '#') continue;  /* meta */
        if (buf[0] == '#') {
            /* #CHROM header — parse tab-separated fields */
            int field = 0, pos = 0;
            while (pos <= len) {
                int fstart = pos;
                while (pos < len && buf[pos] != '\t') pos++;
                if (field >= 9) {
                    int flen = pos - fstart;
                    for (s = 0; s < n_diploid; s++) {
                        const char *sn = CHAR(STRING_ELT(sample_names_sexp, s));
                        if ((int)strlen(sn) == flen &&
                            strncmp(buf + fstart, sn, flen) == 0) {
                            sample_col[s] = field;
                        }
                    }
                }
                pos++;
                field++;
            }
            header_found = 1;
            break;
        }
    }
    if (!header_found) {
        free(sample_col); free(buf); fclose(fp_in);
        Rf_error("vcf_to_ms_call: #CHROM header not found");
    }
    for (s = 0; s < n_diploid; s++) {
        if (sample_col[s] < 0) {
            const char *sn = CHAR(STRING_ELT(sample_names_sexp, s));
            free(sample_col); free(buf); fclose(fp_in);
            Rf_error("vcf_to_ms_call: sample '%s' not found in VCF header", sn);
        }
    }

    /* Build fast field→sample lookup */
    int max_col = 0;
    for (s = 0; s < n_diploid; s++)
        if (sample_col[s] > max_col) max_col = sample_col[s];
    int lookup_len = max_col - 9 + 1;
    int *col_to_sample = (int *)malloc(lookup_len * sizeof(int));
    if (!col_to_sample) {
        free(sample_col); free(buf); fclose(fp_in);
        Rf_error("malloc failed");
    }
    for (v = 0; v < lookup_len; v++) col_to_sample[v] = -1;
    for (s = 0; s < n_diploid; s++)
        col_to_sample[sample_col[s] - 9] = s;

    /* ---- Open output ---- */
    FILE *fp_out = fopen(output_path, "w");
    if (!fp_out) {
        free(col_to_sample); free(sample_col); free(buf); fclose(fp_in);
        Rf_error("vcf_to_ms_call: cannot open output '%s'", output_path);
    }
    setvbuf(fp_out, NULL, _IOFBF, 1 << 20);  /* 1 MB write buffer */

    /* ---- Per-locus variant buffers ---- */
    int var_cap = 128;
    int var_count = 0;
    int *var_pos = (int *)malloc(var_cap * sizeof(int));
    signed char *var_geno = (signed char *)malloc((long)var_cap * n_hap *
                                                   sizeof(signed char));
    if (!var_pos || !var_geno) {
        fclose(fp_in); fclose(fp_out);
        free(col_to_sample); free(sample_col); free(buf);
        if (var_pos)  free(var_pos);
        if (var_geno) free(var_geno);
        Rf_error("malloc failed");
    }

    /* ---- State ---- */
    char current_chrom[256] = "";
    int  current_locus_num  = 0;
    int  loci_written       = 0;
    long n_variants_total   = 0;

    /* ---- Stream through VCF data lines ---- */
    while (1) {
        int len = vcf_read_line(fp_in, &buf, &buf_cap);
        if (len < 0) break;
        if (len == 0 || buf[0] == '#') continue;

        /* --- Field 0: CHROM --- */
        int pos = 0;
        int chrom_start = 0;
        while (pos < len && buf[pos] != '\t') pos++;
        int chrom_len = pos - chrom_start;
        pos++;  /* skip tab */

        /* Detect locus change */
        if (chrom_len != (int)strlen(current_chrom) ||
            strncmp(buf + chrom_start, current_chrom, chrom_len) != 0) {

            /* Flush previous locus */
            if (current_chrom[0] != '\0') {
                write_ms_locus(fp_out, var_count, var_pos, var_geno,
                               n_hap, contig_len);
                loci_written++;
            }

            /* Extract locus number for gap detection */
            int new_num = chrom_to_locus_num(buf + chrom_start, chrom_len);
            int expected = (current_chrom[0] == '\0') ? 1
                                                      : current_locus_num + 1;
            while (expected < new_num) {
                write_ms_locus(fp_out, 0, NULL, NULL, n_hap, contig_len);
                loci_written++;
                expected++;
            }

            /* Update current CHROM */
            int clen = chrom_len < 255 ? chrom_len : 255;
            memcpy(current_chrom, buf + chrom_start, clen);
            current_chrom[clen] = '\0';
            current_locus_num = new_num;
            var_count = 0;
        }

        /* --- Field 1: POS --- */
        int site_pos = 0;
        while (pos < len && buf[pos] != '\t') {
            site_pos = site_pos * 10 + (buf[pos] - '0');
            pos++;
        }
        pos++;  /* skip tab */

        /* --- Skip fields 2..7 (ID REF ALT QUAL FILTER INFO) --- */
        {
            int f;
            for (f = 2; f <= 7 && pos < len; f++) {
                while (pos < len && buf[pos] != '\t') pos++;
                pos++;
            }
        }

        /* --- Field 8: FORMAT — locate GT sub-field index --- */
        int gt_idx = -1;
        {
            int sf = 0, k = pos;
            int fmt_end = pos;
            while (fmt_end < len && buf[fmt_end] != '\t') fmt_end++;
            while (k < fmt_end) {
                int sf_start = k;
                while (k < fmt_end && buf[k] != ':') k++;
                if (k - sf_start == 2 && buf[sf_start] == 'G' &&
                    buf[sf_start + 1] == 'T') {
                    gt_idx = sf;
                    break;
                }
                k++;
                sf++;
            }
            pos = fmt_end + 1;  /* advance past FORMAT field */
        }
        if (gt_idx < 0) continue;  /* no GT — skip variant */

        /* --- Ensure variant buffer capacity --- */
        if (var_count >= var_cap) {
            var_cap *= 2;
            var_pos = (int *)realloc(var_pos, var_cap * sizeof(int));
            var_geno = (signed char *)realloc(var_geno,
                           (long)var_cap * n_hap * sizeof(signed char));
            if (!var_pos || !var_geno)
                Rf_error("vcf_to_ms_call: realloc failed");
        }
        var_pos[var_count] = site_pos;
        signed char *grow = var_geno + (long)var_count * n_hap;

        /* --- Parse sample columns --- */
        int field = 9;
        int skip_variant = 0;
        int samples_filled = 0;

        while (pos <= len && field <= max_col) {
            int fstart = pos;
            while (pos < len && buf[pos] != '\t') pos++;
            /* pos now points to the tab (or end of line) */

            int vcf_col = field - 9;
            if (vcf_col >= 0 && vcf_col < lookup_len) {
                int si = col_to_sample[vcf_col];
                if (si >= 0) {
                    /* Navigate to GT sub-field */
                    int k = fstart;
                    int sf = 0;
                    while (sf < gt_idx && k < pos) {
                        while (k < pos && buf[k] != ':') k++;
                        k++;   /* skip ':' */
                        sf++;
                    }
                    /* k now points to start of GT value */
                    if (k + 2 > pos) { skip_variant = 1; break; }

                    char a1c = buf[k];
                    char a2c = buf[k + 2];
                    if (a1c == '.' || a2c == '.') { skip_variant = 1; break; }

                    grow[2 * si]     = (signed char)(a1c - '0');
                    grow[2 * si + 1] = (signed char)(a2c - '0');
                    samples_filled++;
                }
            }

            pos++;   /* skip tab */
            field++;
        }

        if (skip_variant || samples_filled < n_diploid) continue;

        var_count++;
        n_variants_total++;
        if (n_variants_total % 100000 == 0) R_CheckUserInterrupt();
    }

    /* Flush last locus */
    if (current_chrom[0] != '\0') {
        write_ms_locus(fp_out, var_count, var_pos, var_geno,
                       n_hap, contig_len);
        loci_written++;
    }

    /* Trailing monomorphic loci */
    while (loci_written < n_loci) {
        write_ms_locus(fp_out, 0, NULL, NULL, n_hap, contig_len);
        loci_written++;
    }

    /* ---- Cleanup ---- */
    fclose(fp_in);
    fclose(fp_out);
    free(buf);
    free(sample_col);
    free(col_to_sample);
    free(var_pos);
    free(var_geno);

    SEXP result = PROTECT(ScalarInteger(loci_written));
    UNPROTECT(1);
    return result;
}


/* ================================================================== */
/* vcf_one_snp_call — Subsample one random SNP per locus from VCF     */
/* ================================================================== */

/*
 * vcf_one_snp_call(vcf_path, output_path)
 *
 *   vcf_path    : STRSXP(1) — path to input VCF (CHROM = locus)
 *   output_path : STRSXP(1) — path to output VCF (one SNP per locus)
 *
 * Uses reservoir sampling (probability 1/k for k-th variant in locus)
 * so only one line is buffered at a time.  Uses R's RNG (set.seed in R).
 *
 * Returns: INTSXP(1) — number of loci with at least one variant.
 */
SEXP vcf_one_snp_call(SEXP vcf_path_sexp, SEXP output_path_sexp) {

    if (!isString(vcf_path_sexp) || length(vcf_path_sexp) != 1)
        Rf_error("vcf_one_snp_call: 'vcf_path' must be a single string");
    if (!isString(output_path_sexp) || length(output_path_sexp) != 1)
        Rf_error("vcf_one_snp_call: 'output_path' must be a single string");

    const char *vcf_path    = CHAR(STRING_ELT(vcf_path_sexp, 0));
    const char *output_path = CHAR(STRING_ELT(output_path_sexp, 0));

    FILE *fp_in = fopen(vcf_path, "r");
    if (!fp_in)
        Rf_error("vcf_one_snp_call: cannot open VCF '%s'", vcf_path);
    setvbuf(fp_in, NULL, _IOFBF, 1 << 20);

    FILE *fp_out = fopen(output_path, "w");
    if (!fp_out) {
        fclose(fp_in);
        Rf_error("vcf_one_snp_call: cannot open output '%s'", output_path);
    }
    setvbuf(fp_out, NULL, _IOFBF, 1 << 20);

    /* Line buffers: current line + selected line for this locus */
    int buf_cap = 65536;
    char *buf = (char *)malloc(buf_cap);
    int sel_cap = 65536;
    char *selected = (char *)malloc(sel_cap);
    if (!buf || !selected) {
        if (buf) free(buf);
        if (selected) free(selected);
        fclose(fp_in); fclose(fp_out);
        Rf_error("vcf_one_snp_call: malloc failed");
    }
    int sel_len = 0;

    /* CHROM tracking */
    char current_chrom[256] = "";
    int var_in_locus = 0;
    int loci_written = 0;
    long lines_total = 0;

    /* Use R's RNG */
    GetRNGstate();

    while (1) {
        int len = vcf_read_line(fp_in, &buf, &buf_cap);
        if (len < 0) break;

        /* Copy header lines straight through */
        if (len > 0 && buf[0] == '#') {
            fputs(buf, fp_out);
            fputc('\n', fp_out);
            continue;
        }
        if (len == 0) continue;

        /* Extract CHROM (field 0) */
        int chrom_len = 0;
        while (chrom_len < len && buf[chrom_len] != '\t') chrom_len++;

        /* Detect locus change */
        if (chrom_len != (int)strlen(current_chrom) ||
            strncmp(buf, current_chrom, chrom_len) != 0) {

            /* Write selected line from previous locus */
            if (var_in_locus > 0) {
                selected[sel_len] = '\0';
                fputs(selected, fp_out);
                fputc('\n', fp_out);
                loci_written++;
            }

            /* Update current CHROM */
            int clen = chrom_len < 255 ? chrom_len : 255;
            memcpy(current_chrom, buf, clen);
            current_chrom[clen] = '\0';
            var_in_locus = 0;
        }

        /* Reservoir sampling: keep this line with probability 1/k */
        var_in_locus++;
        if (var_in_locus == 1 || unif_rand() < 1.0 / var_in_locus) {
            if (len >= sel_cap) {
                sel_cap = len + 1;
                selected = (char *)realloc(selected, sel_cap);
                if (!selected) Rf_error("vcf_one_snp_call: realloc failed");
            }
            memcpy(selected, buf, len);
            sel_len = len;
        }

        lines_total++;
        if (lines_total % 100000 == 0) R_CheckUserInterrupt();
    }

    /* Flush last locus */
    if (var_in_locus > 0) {
        selected[sel_len] = '\0';
        fputs(selected, fp_out);
        fputc('\n', fp_out);
        loci_written++;
    }

    PutRNGstate();

    fclose(fp_in);
    fclose(fp_out);
    free(buf);
    free(selected);

    SEXP result_vcf = PROTECT(ScalarInteger(loci_written));
    UNPROTECT(1);
    return result_vcf;
}


/* ================================================================== */
/* phylip_one_snp_call — Subsample one random SNP per locus (PHYLIP)  */
/* ================================================================== */

/*
 * phylip_one_snp_call(phylip_path, output_path)
 *
 *   phylip_path : STRSXP(1) — path to multi-locus sequential PHYLIP
 *   output_path : STRSXP(1) — path to output PHYLIP (one site per locus)
 *
 * For each locus, finds segregating sites (no gaps/N, >1 allele),
 * picks one at random via reservoir sampling, and writes a 1-column
 * PHYLIP block.  Loci with zero seg sites are skipped.
 * Uses R's RNG (set.seed in R).
 *
 * Returns: INTSXP(1) — number of loci written.
 */
SEXP phylip_one_snp_call(SEXP phylip_path_sexp, SEXP output_path_sexp) {

    int t, c, k;

    if (!isString(phylip_path_sexp) || length(phylip_path_sexp) != 1)
        Rf_error("phylip_one_snp_call: 'phylip_path' must be a single string");
    if (!isString(output_path_sexp) || length(output_path_sexp) != 1)
        Rf_error("phylip_one_snp_call: 'output_path' must be a single string");

    const char *phylip_path = CHAR(STRING_ELT(phylip_path_sexp, 0));
    const char *output_path = CHAR(STRING_ELT(output_path_sexp, 0));

    FILE *fp_in = fopen(phylip_path, "r");
    if (!fp_in)
        Rf_error("phylip_one_snp_call: cannot open '%s'", phylip_path);
    setvbuf(fp_in, NULL, _IOFBF, 1 << 20);

    FILE *fp_out = fopen(output_path, "w");
    if (!fp_out) {
        fclose(fp_in);
        Rf_error("phylip_one_snp_call: cannot open output '%s'", output_path);
    }
    setvbuf(fp_out, NULL, _IOFBF, 1 << 20);

    int buf_cap = 65536;
    char *buf = (char *)malloc(buf_cap);
    if (!buf) { fclose(fp_in); fclose(fp_out); Rf_error("malloc failed"); }

    /* Per-locus buffers — grown as needed */
    int max_ntax  = 0;
    int max_nchar = 0;
    char *seq_buf  = NULL;   /* flat [ntax * max_nchar], row-major */
    char *name_buf = NULL;   /* flat [ntax * 11], each name 10 chars + '\0' */

    int loci_written = 0;
    int loci_total   = 0;

    GetRNGstate();

    while (1) {
        /* Read next non-blank line → should be locus header */
        int len;
        int got_line = 0;
        while (1) {
            len = vcf_read_line(fp_in, &buf, &buf_cap);
            if (len < 0) goto phy_done;
            /* skip blank lines */
            int blank = 1;
            for (k = 0; k < len; k++) {
                if (buf[k] != ' ' && buf[k] != '\t' && buf[k] != '\r') {
                    blank = 0; break;
                }
            }
            if (!blank) { got_line = 1; break; }
        }
        if (!got_line) break;

        /* Parse header: "ntax nchar" */
        int ntax = 0, nchar_seq = 0;
        if (sscanf(buf, "%d %d", &ntax, &nchar_seq) != 2 ||
            ntax <= 0 || nchar_seq <= 0)
            continue;

        loci_total++;

        /* Grow buffers if needed */
        if (ntax > max_ntax || nchar_seq > max_nchar) {
            if (ntax > max_ntax) max_ntax = ntax;
            if (nchar_seq > max_nchar) max_nchar = nchar_seq;
            seq_buf  = (char *)realloc(seq_buf,
                           (long)max_ntax * max_nchar);
            name_buf = (char *)realloc(name_buf,
                           (long)max_ntax * 11);
            if (!seq_buf || !name_buf)
                Rf_error("phylip_one_snp_call: realloc failed");
        }

        /* Read ntax sequence lines */
        int valid = 1;
        for (t = 0; t < ntax; t++) {
            len = vcf_read_line(fp_in, &buf, &buf_cap);
            if (len < 0) { valid = 0; break; }

            /* Name: first 10 chars (strict PHYLIP) */
            char *nm = name_buf + t * 11;
            int nlen = (len < 10) ? len : 10;
            memcpy(nm, buf, nlen);
            int e;
            for (e = nlen - 1; e >= 0 && nm[e] == ' '; e--) ;
            nm[e + 1] = '\0';

            /* Sequence: chars from position 10 onward, skip whitespace */
            char *sq = seq_buf + (long)t * max_nchar;
            int si = 0;
            for (k = 10; k < len && si < nchar_seq; k++) {
                if (buf[k] != ' ' && buf[k] != '\t')
                    sq[si++] = buf[k];
            }
            while (si < nchar_seq) sq[si++] = 'n';
        }
        if (!valid) break;

        /* Find segregating sites via reservoir sampling */
        int n_seg = 0;
        int selected_col = -1;

        for (c = 0; c < nchar_seq; c++) {
            int has_bad = 0;
            char first = 0;
            int is_var = 0;

            for (t = 0; t < ntax; t++) {
                char base = seq_buf[(long)t * max_nchar + c];
                if (base >= 'A' && base <= 'Z') base += 32;  /* tolower */
                if (base == '-' || base == 'n') { has_bad = 1; break; }
                if (first == 0) first = base;
                else if (base != first) is_var = 1;
            }
            if (has_bad || !is_var) continue;

            n_seg++;
            if (n_seg == 1 || unif_rand() < 1.0 / n_seg)
                selected_col = c;
        }

        if (selected_col < 0) continue;  /* monomorphic — skip */

        /* Write output block: "ntax 1\n" + name + base */
        fprintf(fp_out, "%d 1\n", ntax);
        for (t = 0; t < ntax; t++) {
            char *nm = name_buf + t * 11;
            char base = seq_buf[(long)t * max_nchar + selected_col];
            fprintf(fp_out, "%-10s%c\n", nm, base);
        }

        loci_written++;
        if (loci_total % 10000 == 0) R_CheckUserInterrupt();
    }

phy_done:
    PutRNGstate();

    if (seq_buf)  free(seq_buf);
    if (name_buf) free(name_buf);
    free(buf);
    fclose(fp_in);
    fclose(fp_out);

    SEXP result_phy = PROTECT(ScalarInteger(loci_written));
    UNPROTECT(1);
    return result_phy;
}
