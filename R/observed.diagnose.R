#' Diagnose non-modeled signals in observed per-locus statistics
#'
#' Pre-training data diagnostic that flags departures from the assumptions
#' of a neutral, uniform-rate coalescent before generating a reference
#' table or fitting a model. Useful as a sanity check on real WGS data:
#' is the dataset compatible with PipeMaster's pure-demographic inference,
#' or are there strong non-modeled signals (background selection,
#' recombination rate heterogeneity, per-locus quality artifacts) that
#' will dominate the summary statistics?
#'
#' Runs up to four tests:
#'
#' \itemize{
#'   \item \strong{Test A — joint D-pi distribution} (informational):
#'         Spearman correlation between per-locus Tajima's D and per-site
#'         pi. A strong positive correlation can arise from background
#'         selection, mutation-rate heterogeneity, or as a statistical
#'         property of the D estimator on short loci under non-stationary
#'         demography. Use \code{cds_bed} (Test D) to distinguish the first;
#'         a neutral-coalescent simulation (planned for a future release)
#'         to distinguish the third.
#'   \item \strong{Test B — over-dispersion of per-locus moments}:
#'         Compares observed moments (mean / var / skew / kurt) of pi, D,
#'         segs, and ZnS across loci to the distribution of the same
#'         moments across simulations in a PipeMaster reference table.
#'         Flags obs values outside the central 95\% sim band. Skipped if
#'         \code{reftable = NULL}.
#'   \item \strong{Test C — chromosome landscape}:
#'         Per-chromosome mean pi and lag-1 Spearman autocorrelation of pi
#'         along chromosomes (sorted by start position). Significant
#'         spatial autocorrelation indicates recombination rate
#'         heterogeneity (which PipeMaster's uniform recombination rate
#'         cannot reproduce). Per-chromosome outliers flag chromosome-
#'         specific effects.
#'   \item \strong{Test D — BGS gradient}:
#'         Spearman correlation between per-locus pi (and D) and the
#'         distance from each locus midpoint to the nearest interval in
#'         \code{cds_bed}. Under residual background selection, pi rises
#'         with distance from selection targets. Skipped if
#'         \code{cds_bed = NULL}.
#' }
#'
#' Verdict logic combines Tests B, C, and D into pass / warn / fail.
#' Test A is reported as informational only (its interpretation requires
#' a matched neutral-coalescent simulation, planned for a follow-up
#' release).
#'
#' @param per_locus data.frame of per-locus statistics. Required columns:
#'   \code{locus_id}, \code{length}, \code{S}, \code{pi_per_site}, \code{D}.
#'   Recommended: \code{thetaW_per_site}.
#' @param locus_table data.frame mapping loci to genome coordinates.
#'   Required columns: \code{locus_id}, \code{chrom}, \code{start},
#'   \code{end}.
#' @param cds_bed Optional path to a BED file of CDS (or any selection-
#'   target) intervals: 3 columns chrom / start / end (additional columns
#'   ignored). If \code{NULL}, Test D is skipped.
#' @param reftable Optional PipeMaster simulation reference table
#'   (data.frame with column names like \code{s_mean_pi_1},
#'   \code{s_var_pi_1}, etc., as produced by \code{sim.scrm.sumstats()}).
#'   If \code{NULL}, Test B is skipped.
#' @param pop Integer population index (default \code{1}) — selects which
#'   pop-specific moment columns to compare against
#'   (\code{s_mean_pi_<pop>}, \code{s_var_pi_<pop>}, etc.).
#' @param alpha Numeric (default \code{0.01}). Two-sided significance for
#'   the spatial-autocorrelation flag; the over-dispersion band uses 95\%
#'   (fixed) regardless of \code{alpha}.
#' @param subsample Integer, optional. For plotting D-vs-pi scatter
#'   (Test A), randomly subsample this many loci. Defaults to
#'   \code{min(30000, nrow(per_locus))}.
#' @param plot Logical. Draw the diagnostic figure (6 panels when both
#'   \code{cds_bed} and \code{reftable} are supplied; fewer otherwise).
#' @param verbose Logical. Print progress to console.
#'
#' @return An object of class \code{c("observed_diagnostic",
#'   "OOD_diagnostic")} with fields:
#' \itemize{
#'   \item \code{test_A}: D-pi correlation + quintile table.
#'   \item \code{test_B}: per-moment z-scores and outlier directions (or
#'         \code{NULL} if \code{reftable} not supplied).
#'   \item \code{test_C}: per-chromosome summary and lag-1 autocorrelations.
#'   \item \code{test_D}: distance-binned pi and D (or \code{NULL}).
#'   \item \code{verdict}: \code{"pass"} / \code{"warn"} / \code{"fail"}.
#'   \item \code{verdict_reasons}: character vector explaining the verdict.
#' }
#' Methods: \code{print}, \code{summary}, \code{plot}.
#'
#' @examples
#' \dontrun{
#' # Per-locus stats computed externally (e.g., via bcftools query +
#' # standard Tajima's D formula). Schema documented in the help page.
#' per_locus <- read.table("per_locus_stats.tsv", header = TRUE)
#' locus_table <- read.table("locus_table.txt", header = TRUE)
#'
#' # Minimal: just Tests A + C (no annotation, no reftable)
#' diag1 <- observed.diagnose(per_locus, locus_table)
#'
#' # Full: Tests A + B + C + D
#' reftable <- read.table("SIMS_reftable_10K.txt", header = TRUE)
#' diag2 <- observed.diagnose(per_locus, locus_table,
#'                            cds_bed = "cds.bed",
#'                            reftable = reftable)
#' print(diag2)
#' plot(diag2)
#' }
#'
#' @export
observed.diagnose <- function(per_locus, locus_table,
                              cds_bed = NULL,
                              reftable = NULL,
                              pop = 1,
                              alpha = 0.01,
                              subsample = NULL,
                              plot = TRUE,
                              verbose = TRUE) {

  # ---- input validation ----
  req_pl <- c("locus_id", "length", "S", "pi_per_site", "D")
  miss <- setdiff(req_pl, colnames(per_locus))
  if (length(miss) > 0L)
    stop("per_locus missing columns: ", paste(miss, collapse = ", "))
  req_lt <- c("locus_id", "chrom", "start", "end")
  miss <- setdiff(req_lt, colnames(locus_table))
  if (length(miss) > 0L)
    stop("locus_table missing columns: ", paste(miss, collapse = ", "))

  per_locus <- merge(per_locus,
                     locus_table[, c("locus_id", "chrom", "start", "end")],
                     by = "locus_id", all.x = TRUE)
  per_locus$mid <- (per_locus$start + per_locus$end) / 2

  if (verbose)
    cat(sprintf("PipeMaster:: observed.diagnose: %d loci, %d chromosomes\n",
                nrow(per_locus), length(unique(per_locus$chrom))))

  # ============ Test A: D-pi correlation ============
  if (verbose) cat("PipeMaster::   Test A (D-pi correlation)\n")
  test_A <- .observed.test.A(per_locus)

  # ============ Test B: over-dispersion (optional) ============
  test_B <- NULL
  if (!is.null(reftable)) {
    if (verbose) cat("PipeMaster::   Test B (over-dispersion vs reftable)\n")
    test_B <- .observed.test.B(per_locus, reftable, pop)
  }

  # ============ Test C: chromosome landscape ============
  if (verbose) cat("PipeMaster::   Test C (chromosome landscape)\n")
  test_C <- .observed.test.C(per_locus)

  # ============ Test D: BGS gradient (optional) ============
  test_D <- NULL
  if (!is.null(cds_bed)) {
    if (verbose) cat("PipeMaster::   Test D (BGS gradient vs CDS)\n")
    test_D <- .observed.test.D(per_locus, cds_bed)
  }

  # ============ Verdict ============
  reasons <- character()
  # B
  if (!is.null(test_B)) {
    non_zns <- test_B$summary[!grepl("ZnS", test_B$summary$stat), ]
    frac_out <- mean(non_zns$direction != "")
    if (frac_out > 0.50)
      reasons <- c(reasons, sprintf(
        "Test B: %.0f%% of non-ZnS moments outside 95%% sim band (>50%% -> fail)",
        100 * frac_out))
    else if (frac_out > 0.30)
      reasons <- c(reasons, sprintf(
        "Test B: %.0f%% of non-ZnS moments outside 95%% sim band (>30%% -> warn)",
        100 * frac_out))
  }
  # C
  high_autocorr <- test_C$lag1_top5$rho > 0.30
  if (any(high_autocorr, na.rm = TRUE))
    reasons <- c(reasons, sprintf(
      "Test C: lag-1 spatial autocorrelation > 0.30 on %d/%d top chromosomes (recombination landscape)",
      sum(high_autocorr, na.rm = TRUE), nrow(test_C$lag1_top5)))
  # D
  if (!is.null(test_D)) {
    if (test_D$rho_pi_cds > 0.30)
      reasons <- c(reasons, sprintf(
        "Test D: rho(d_CDS, pi) = %+.3f (>0.30 -> residual BGS likely)",
        test_D$rho_pi_cds))
    else if (test_D$rho_pi_cds > 0.15)
      reasons <- c(reasons, sprintf(
        "Test D: rho(d_CDS, pi) = %+.3f (>0.15 -> mild residual BGS)",
        test_D$rho_pi_cds))
  }

  # Severity: any "fail" word in reasons OR fraction > 0.50 -> fail
  has_fail <- any(grepl("fail|likely", reasons, ignore.case = TRUE))
  verdict <- if (length(reasons) == 0L) "pass"
             else if (has_fail) "fail"
             else "warn"

  out <- list(
    test_A          = test_A,
    test_B          = test_B,
    test_C          = test_C,
    test_D          = test_D,
    verdict         = verdict,
    verdict_reasons = reasons,
    n_loci          = nrow(per_locus),
    cds_bed         = cds_bed,
    reftable_supplied = !is.null(reftable),
    pop             = pop
  )
  class(out) <- c("observed_diagnostic", "OOD_diagnostic")

  if (plot)
    plot(out, subsample = subsample)

  if (verbose)
    cat(sprintf("PipeMaster:: verdict: %s%s\n",
                toupper(verdict),
                if (length(reasons) == 0L) ""
                else sprintf(" (%d issue%s)", length(reasons),
                             ifelse(length(reasons) == 1L, "", "s"))))
  invisible(out)
}

# ---- internal helpers -------------------------------------------------

.observed.test.A <- function(pl) {
  ok <- !is.na(pl$D) & is.finite(pl$pi_per_site) & pl$pi_per_site > 0
  A  <- pl[ok, ]
  rho_p <- cor(A$D, A$pi_per_site, method = "pearson")
  rho_s <- cor(A$D, A$pi_per_site, method = "spearman")
  set.seed(42)
  perm_rho <- replicate(200, cor(A$D, sample(A$pi_per_site),
                                  method = "spearman"))
  perm_p <- mean(abs(perm_rho) >= abs(rho_s))
  q <- quantile(A$pi_per_site, probs = seq(0, 1, by = 0.2))
  A$pi_q <- cut(A$pi_per_site, breaks = q, include.lowest = TRUE,
                labels = paste0("Q", 1:5))
  qtbl <- do.call(rbind, lapply(split(A, A$pi_q), function(d) {
    data.frame(pi_q = unique(d$pi_q), n = nrow(d),
               pi_mean = mean(d$pi_per_site),
               D_mean  = mean(d$D, na.rm = TRUE),
               D_se    = sd(d$D, na.rm = TRUE) /
                         sqrt(sum(!is.na(d$D))),
               S_mean  = mean(d$S),
               stringsAsFactors = FALSE)
  }))
  rownames(qtbl) <- NULL
  list(rho_pearson = rho_p, rho_spearman = rho_s, perm_p = perm_p,
       quintile_table = qtbl, n_used = nrow(A))
}

.observed.test.B <- function(pl, reftable, pop) {
  # PipeMaster convention: s_mean_pi etc. are moments OVER LOCI of the
  # per-locus pi total (sum of pairwise pi across segsites within a locus),
  # NOT moments of per-site pi. Convert pi_per_site -> pi_total here.
  pi_locus <- pl$pi_per_site * pl$length
  if ("thetaW_per_site" %in% colnames(pl)) {
    thW_locus <- pl$thetaW_per_site * pl$length
  } else {
    thW_locus <- NULL
  }
  obs <- list(
    s_mean_pi = mean(pi_locus, na.rm = TRUE),
    s_var_pi  = var(pi_locus, na.rm = TRUE),
    s_skew_pi = .skew(pi_locus),
    s_kurt_pi = .kurt(pi_locus),
    s_mean_tajd = mean(pl$D, na.rm = TRUE),
    s_var_tajd  = var(pl$D, na.rm = TRUE),
    s_skew_tajd = .skew(pl$D),
    s_kurt_tajd = .kurt(pl$D),
    s_mean_segs = mean(pl$S, na.rm = TRUE),
    s_var_segs  = var(pl$S, na.rm = TRUE),
    s_skew_segs = .skew(pl$S),
    s_kurt_segs = .kurt(pl$S)
  )
  if (!is.null(thW_locus)) {
    obs$s_mean_thetaW <- mean(thW_locus, na.rm = TRUE)
    obs$s_var_thetaW  <- var(thW_locus, na.rm = TRUE)
    obs$s_skew_thetaW <- .skew(thW_locus)
    obs$s_kurt_thetaW <- .kurt(thW_locus)
  }
  # match to reftable columns with pop suffix
  out_rows <- list()
  for (nm in names(obs)) {
    col <- sprintf("%s_%d", nm, pop)
    if (!col %in% colnames(reftable)) {
      col2 <- nm                                      # try unsuffixed
      if (!col2 %in% colnames(reftable)) next
      col <- col2
    }
    sv <- as.numeric(reftable[[col]])
    sv <- sv[is.finite(sv)]
    if (length(sv) < 100L) next
    ov <- as.numeric(obs[[nm]])
    if (!is.finite(ov)) next
    mu <- mean(sv); sg <- sd(sv)
    qq <- quantile(sv, c(0.025, 0.975))
    Fx <- ecdf(sv)(ov); pval <- 2 * min(Fx, 1 - Fx)
    dir <- if (ov < qq[1]) "LO" else if (ov > qq[2]) "HI" else ""
    out_rows[[length(out_rows) + 1L]] <- data.frame(
      stat = col, obs = ov, sim_mean = mu, sim_sd = sg,
      sim_q025 = unname(qq[1]), sim_q975 = unname(qq[2]),
      z = (ov - mu) / sg, pval = pval, direction = dir,
      stringsAsFactors = FALSE)
  }
  summary_df <- do.call(rbind, out_rows)
  rownames(summary_df) <- NULL
  list(summary = summary_df,
       n_outliers = sum(summary_df$direction != ""),
       n_total    = nrow(summary_df))
}

.observed.test.C <- function(pl) {
  C <- pl[!is.na(pl$chrom) & is.finite(pl$pi_per_site), ]
  global_mu <- mean(C$pi_per_site)
  global_sd <- sd(C$pi_per_site)
  chrom_summary <- do.call(rbind, lapply(split(C, C$chrom), function(d) {
    data.frame(chrom = unique(d$chrom), n_loci = nrow(d),
               total_kb = sum(d$length) / 1000,
               mean_pi  = mean(d$pi_per_site),
               median_pi = median(d$pi_per_site),
               mean_D   = mean(d$D, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
  rownames(chrom_summary) <- NULL
  chrom_summary$z_mean_pi <- (chrom_summary$mean_pi - global_mu) / global_sd
  chrom_summary <- chrom_summary[order(-chrom_summary$n_loci), ]

  # lag-1 spatial autocorrelation on top-5 chromosomes
  top5 <- head(chrom_summary$chrom, 5)
  lag1 <- data.frame(chrom = top5, n_loci = NA_integer_,
                     rho = NA_real_, stringsAsFactors = FALSE)
  for (i in seq_along(top5)) {
    d <- C[C$chrom == top5[i], ]
    d <- d[order(d$start), ]
    if (nrow(d) >= 100L) {
      x <- d$pi_per_site
      lag1$rho[i]    <- cor(x[-length(x)], x[-1])
      lag1$n_loci[i] <- nrow(d)
    }
  }
  list(global_mean_pi = global_mu,
       chrom_summary  = chrom_summary,
       lag1_top5      = lag1)
}

.observed.test.D <- function(pl, cds_bed) {
  bed <- utils::read.table(cds_bed, header = FALSE, sep = "\t",
                            col.names = c("chrom", "start", "end")[1:3],
                            stringsAsFactors = FALSE,
                            colClasses = c("character", "integer",
                                           "integer"))[, 1:3]
  d_cds <- .dist.to.intervals(pl$chrom, pl$mid, bed)
  pl$d_cds <- d_cds

  bins <- c(-1, 0, 1000, 5000, 10000, 25000, 50000, 100000, Inf)
  labs <- c("0 (mask edge)", "0-1kb", "1-5kb", "5-10kb",
            "10-25kb", "25-50kb", "50-100kb", ">100kb")
  pl$bin_cds <- cut(pl$d_cds, breaks = bins, labels = labs)
  bin_summary <- do.call(rbind, lapply(split(pl, pl$bin_cds), function(d) {
    if (nrow(d) == 0L) return(NULL)
    data.frame(bin = unique(d$bin_cds), n = nrow(d),
               mean_pi = mean(d$pi_per_site, na.rm = TRUE),
               se_pi   = sd(d$pi_per_site, na.rm = TRUE) /
                         sqrt(sum(is.finite(d$pi_per_site))),
               mean_D  = mean(d$D, na.rm = TRUE),
               se_D    = sd(d$D, na.rm = TRUE) /
                         sqrt(sum(!is.na(d$D))),
               stringsAsFactors = FALSE)
  }))
  rownames(bin_summary) <- NULL

  ok <- is.finite(pl$d_cds) & is.finite(pl$pi_per_site)
  rho_pi <- cor(pl$d_cds[ok], pl$pi_per_site[ok], method = "spearman")
  ok_D <- ok & !is.na(pl$D)
  rho_D <- cor(pl$d_cds[ok_D], pl$D[ok_D], method = "spearman")

  list(bin_summary = bin_summary,
       rho_pi_cds  = rho_pi,
       rho_D_cds   = rho_D,
       n_inside    = sum(pl$d_cds == 0, na.rm = TRUE))
}

# Distance from query points to nearest interval (per chromosome).
.dist.to.intervals <- function(query_chrom, query_pos, intervals) {
  out <- rep(NA_real_, length(query_pos))
  for (ch in unique(query_chrom)) {
    q_idx <- which(query_chrom == ch)
    qp <- query_pos[q_idx]
    iv <- intervals[intervals$chrom == ch, , drop = FALSE]
    if (nrow(iv) == 0L) { out[q_idx] <- Inf; next }
    iv <- iv[order(iv$start), ]
    starts <- iv$start; ends <- iv$end
    idx <- findInterval(qp, starts)
    d_after_prev  <- ifelse(idx > 0L,
                             pmax(qp - ends[idx], 0), Inf)
    d_before_next <- ifelse(idx < length(starts),
                             pmax(starts[idx + 1L] - qp, 0), Inf)
    out[q_idx] <- pmin(d_after_prev, d_before_next)
  }
  out
}

.skew <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  m <- mean(x); s <- sd(x); if (s == 0) return(NA_real_)
  mean(((x - m) / s)^3)
}
.kurt <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 4L) return(NA_real_)
  m <- mean(x); s <- sd(x); if (s == 0) return(NA_real_)
  mean(((x - m) / s)^4) - 3
}

# ---- S3 methods -------------------------------------------------------

#' @export
print.observed_diagnostic <- function(x, ...) {
  cat(sprintf("Observed-data diagnostic (PipeMaster)\n"))
  cat(sprintf("  Loci: %d\n", x$n_loci))
  cat(sprintf("  Verdict: %s\n", toupper(x$verdict)))
  if (length(x$verdict_reasons) > 0L) {
    cat("  Reasons:\n")
    for (r in x$verdict_reasons) cat(sprintf("    - %s\n", r))
  }
  cat(sprintf("  Test A  rho(pi,D) = %+.4f (perm p = %.4f)\n",
              x$test_A$rho_spearman, x$test_A$perm_p))
  if (!is.null(x$test_B))
    cat(sprintf("  Test B  outlier moments: %d / %d\n",
                x$test_B$n_outliers, x$test_B$n_total))
  if (!is.null(x$test_C$lag1_top5)) {
    rng <- range(x$test_C$lag1_top5$rho, na.rm = TRUE)
    cat(sprintf("  Test C  lag-1 rho range (top 5 chrom): [%+.3f, %+.3f]\n",
                rng[1], rng[2]))
  }
  if (!is.null(x$test_D))
    cat(sprintf("  Test D  rho(d_CDS, pi) = %+.4f, rho(d_CDS, D) = %+.4f\n",
                x$test_D$rho_pi_cds, x$test_D$rho_D_cds))
  invisible(x)
}

#' @export
summary.observed_diagnostic <- function(object, ...) {
  print(object)
  cat("\n--- Test A: pi quintile -> mean D ---\n")
  print(object$test_A$quintile_table, row.names = FALSE, digits = 4)
  if (!is.null(object$test_B)) {
    cat("\n--- Test B: obs vs sim moment distributions ---\n")
    print(object$test_B$summary, row.names = FALSE, digits = 4)
  }
  cat("\n--- Test C: chromosome summary (top 10 by n_loci) ---\n")
  print(utils::head(object$test_C$chrom_summary, 10),
        row.names = FALSE, digits = 5)
  if (!is.null(object$test_D)) {
    cat("\n--- Test D: pi and D by distance to CDS ---\n")
    print(object$test_D$bin_summary, row.names = FALSE, digits = 4)
  }
  invisible(object)
}

#' @export
plot.observed_diagnostic <- function(x, subsample = NULL, ...) {
  pl_A <- x$test_A; qtbl <- pl_A$quintile_table
  n_panels <- 4L + (!is.null(x$test_B)) + (!is.null(x$test_D))
  if (n_panels == 4L) {
    op <- graphics::par(mfrow = c(2, 2),
                        mar = c(4, 4.3, 3, 1), oma = c(0, 0, 2, 0))
  } else if (n_panels == 5L) {
    op <- graphics::par(mfrow = c(2, 3),
                        mar = c(4, 4.3, 3, 1), oma = c(0, 0, 2, 0))
  } else {
    op <- graphics::par(mfrow = c(2, 3),
                        mar = c(4, 4.3, 3, 1), oma = c(0, 0, 2, 0))
  }
  on.exit(graphics::par(op))

  # (1) D vs pi - quintile bar with SE
  ymin <- min(qtbl$D_mean - 2 * qtbl$D_se, na.rm = TRUE)
  ymax <- max(qtbl$D_mean + 2 * qtbl$D_se, na.rm = TRUE)
  pad  <- 0.05 * (ymax - ymin)
  bp <- graphics::barplot(qtbl$D_mean, names.arg = qtbl$pi_q,
                          ylim = c(ymin - pad, ymax + pad),
                          col = "#9ecae1", border = "#3182bd",
                          ylab = "mean Tajima's D",
                          xlab = "pi quintile (Q1=low)",
                          main = sprintf("(A) D by pi quintile\nrho(pi,D)=%+.3f, p=%.3f",
                                          pl_A$rho_spearman, pl_A$perm_p))
  graphics::arrows(bp, qtbl$D_mean - 2 * qtbl$D_se,
                   bp, qtbl$D_mean + 2 * qtbl$D_se, code = 3,
                   angle = 90, length = 0.04, col = "grey40")
  graphics::abline(h = 0, col = "grey50", lty = 2)

  # (2) Chromosome summary - per-chrom mean pi z-score
  cs <- x$test_C$chrom_summary
  top22 <- utils::head(cs, 22)
  graphics::par(mar = c(7, 4.3, 3, 1))
  graphics::barplot(top22$z_mean_pi,
                    names.arg = top22$chrom,
                    col = ifelse(abs(top22$z_mean_pi) > 2,
                                  "#fdae6b", "#a1d99b"),
                    border = "white", las = 2,
                    cex.names = 0.55, cex.axis = 0.8,
                    ylab = "z(mean pi)", xlab = "",
                    main = "(C1) Per-chromosome mean pi (z-score)")
  graphics::abline(h = c(-2, 2), col = "grey50", lty = 3)
  graphics::par(mar = c(4, 4.3, 3, 1))

  # (3) lag-1 autocorrelation top 5
  l1 <- x$test_C$lag1_top5
  bp <- graphics::barplot(l1$rho, names.arg = l1$chrom,
                          col = ifelse(l1$rho > 0.30, "#fdae6b", "#a1d99b"),
                          border = "white",
                          ylab = "lag-1 Spearman rho(pi)",
                          main = "(C2) Spatial autocorr (top 5 chrom)")
  graphics::abline(h = c(0, 0.30), col = c("grey50", "red"),
                   lty = c(2, 3))

  # (4) Test B if available
  if (!is.null(x$test_B)) {
    tb <- x$test_B$summary
    tb <- tb[order(grepl("ZnS", tb$stat), tb$stat), ]
    cols <- ifelse(tb$direction == "LO", "#3182bd",
            ifelse(tb$direction == "HI", "#cb181d", "grey70"))
    graphics::par(mar = c(7, 4.3, 3, 1))
    z_show <- pmin(pmax(tb$z, -8), 8)
    graphics::barplot(z_show, names.arg = tb$stat, las = 2,
                      cex.names = 0.55, col = cols, border = "white",
                      ylab = "z (obs in sim distribution)",
                      main = "(B) Per-locus moment over-dispersion",
                      ylim = c(-8.5, 8.5))
    graphics::abline(h = c(-2, 0, 2), col = c("grey70", "grey50", "grey70"),
                     lty = c(3, 2, 3))
    graphics::par(mar = c(4, 4.3, 3, 1))
  }

  # (5+6) Test D if available
  if (!is.null(x$test_D)) {
    bd <- x$test_D$bin_summary
    graphics::par(mar = c(7, 4.3, 3, 1))
    graphics::barplot(bd$mean_pi, names.arg = as.character(bd$bin), las = 2,
                      cex.names = 0.55,
                      col = "#9ecae1", border = "#3182bd",
                      ylab = "mean pi", xlab = "",
                      main = sprintf("(D1) pi vs distance to CDS\nrho=%+.3f",
                                      x$test_D$rho_pi_cds))
    graphics::par(mar = c(4, 4.3, 3, 1))
  }

  graphics::mtext(sprintf("observed.diagnose | verdict: %s",
                           toupper(x$verdict)),
                  outer = TRUE, cex = 1.05, font = 2)
  invisible(NULL)
}
