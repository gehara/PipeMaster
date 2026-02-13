#!/usr/bin/env Rscript
#
# Simulation-based calibration (coverage test) for ABC parameter estimation
#
# Tests whether the ABC pipeline is well-calibrated by:
#   1. Sampling rows from the reference table as pseudo-observed data
#   2. Running ABC on the remaining rows
#   3. Checking if the 95% CI contains the true parameter values
#   4. Repeating 100 times and computing coverage frequencies
#
# A well-calibrated pipeline should have ~95% coverage for all parameters.
# If migration parameters have significantly lower coverage than Ne/time
# parameters, it indicates a statistical power issue, not a pipeline bug.

library(abc)

# ============================================================
# Configuration
# ============================================================

n_reps     <- 100     # number of calibration replicates
tol        <- 0.005   # ABC tolerance (0.5% -> 500 accepted from 100K)
ci_level   <- 0.95    # credible interval level
seed       <- 42      # reproducibility

set.seed(seed)

# ============================================================
# Load data
# ============================================================

cat("=== Migration Pipeline Calibration Test ===\n\n")
cat("Loading reference table and model...\n")

load("data_to_test/sfs_reftable_results.RData")  # sim_sfs_ooa
load("data_to_test/test_models.RData")            # OutOfAfrica_3G09

# Parameter and SFS column names
ooa_params <- c("Ne0.pop1", "Ne0.pop2", "Ne1.pop2", "Ne1.pop1",
                "join1", "t.Ne1.pop2", "t.Ne1.pop1",
                "mig0.1_2", "mig0.2_1")
sfs_cols <- grep("^sfs_", colnames(sim_sfs_ooa), value = TRUE)

cat(sprintf("  Reference table: %d simulations, %d SFS bins, %d parameters\n",
            nrow(sim_sfs_ooa), length(sfs_cols), length(ooa_params)))

# Pre-extract full matrices
full_params <- as.matrix(sim_sfs_ooa[, ooa_params])
full_stats  <- as.matrix(sim_sfs_ooa[, sfs_cols])

# Remove zero-variance SFS columns (same as in test_abc_sfs.R)
col_sd <- apply(full_stats, 2, sd)
keep <- col_sd > 1e-10
n_removed <- sum(!keep)
if (n_removed > 0) {
  cat(sprintf("  Removed %d zero-variance SFS bins (%d -> %d)\n",
              n_removed, ncol(full_stats), sum(keep)))
  full_stats <- full_stats[, keep]
}

# Remove rows with NA/Inf
bad <- apply(full_stats, 1, function(x) any(!is.finite(x)))
if (sum(bad) > 0) {
  cat(sprintf("  Removed %d rows with NA/Inf\n", sum(bad)))
  full_stats  <- full_stats[!bad, ]
  full_params <- full_params[!bad, ]
}

n_total <- nrow(full_stats)
cat(sprintf("  Usable simulations: %d\n", n_total))

# ============================================================
# Calibration loop
# ============================================================

cat(sprintf("\nRunning %d calibration replicates (tol=%.3f, CI=%.0f%%)...\n",
            n_reps, tol, ci_level * 100))

alpha <- 1 - ci_level
q_lo  <- alpha / 2
q_hi  <- 1 - alpha / 2

# Storage: coverage (TRUE/FALSE), relative error, CI width
coverage_mat <- matrix(FALSE, nrow = n_reps, ncol = length(ooa_params),
                       dimnames = list(NULL, ooa_params))
error_mat    <- matrix(NA, nrow = n_reps, ncol = length(ooa_params),
                       dimnames = list(NULL, ooa_params))
ci_width_mat <- matrix(NA, nrow = n_reps, ncol = length(ooa_params),
                       dimnames = list(NULL, ooa_params))

# Sample pseudo-observed indices (without replacement for independence)
pod_indices <- sample(n_total, n_reps, replace = FALSE)

t0 <- proc.time()

for (r in 1:n_reps) {
  idx <- pod_indices[r]

  # Pseudo-observed data
  obs_stats  <- matrix(full_stats[idx, ], nrow = 1)
  true_vals  <- full_params[idx, ]

  # Reference table (leave one out)
  ref_stats  <- full_stats[-idx, ]
  ref_params <- full_params[-idx, ]

  # ABC rejection
  abc_res <- abc(target  = obs_stats,
                 param   = ref_params,
                 sumstat = ref_stats,
                 tol     = tol,
                 method  = "rejection")

  # Check coverage and error for each parameter
  for (j in seq_along(ooa_params)) {
    posterior <- abc_res$unadj.values[, j]
    ci <- quantile(posterior, c(q_lo, q_hi))
    tv <- true_vals[j]

    coverage_mat[r, j] <- (tv >= ci[1]) && (tv <= ci[2])
    error_mat[r, j]    <- abs(median(posterior) - tv) / abs(tv) * 100
    ci_width_mat[r, j] <- (ci[2] - ci[1]) / abs(tv) * 100
  }

  if (r %% 10 == 0) {
    elapsed <- (proc.time() - t0)[3]
    rate <- r / elapsed
    remaining <- (n_reps - r) / rate
    cat(sprintf("  rep %3d/%d  (%.0f sec elapsed, ~%.0f sec remaining)\n",
                r, n_reps, elapsed, remaining))
  }
}

elapsed_total <- (proc.time() - t0)[3]
cat(sprintf("\nCompleted %d replicates in %.1f seconds (%.2f sec/rep)\n",
            n_reps, elapsed_total, elapsed_total / n_reps))

# ============================================================
# Results summary
# ============================================================

cat("\n===================================================\n")
cat("  CALIBRATION RESULTS: OutOfAfrica_3G09\n")
cat("===================================================\n\n")

# Coverage table
coverage_pct <- colMeans(coverage_mat) * 100
mean_error   <- colMeans(error_mat)
mean_ci_width <- colMeans(ci_width_mat)

cat(sprintf("  %-15s %10s %12s %12s %8s\n",
            "Parameter", "Coverage%", "Mean.Err%", "Mean.CI.W%", "Status"))
cat(sprintf("  %s\n", paste(rep("-", 62), collapse = "")))

# Classify parameters
ne_params  <- c("Ne0.pop1", "Ne0.pop2", "Ne1.pop2", "Ne1.pop1")
mig_params <- c("mig0.1_2", "mig0.2_1")
time_params <- c("join1", "t.Ne1.pop2", "t.Ne1.pop1")

for (p in ooa_params) {
  # Expected coverage with 100 reps: 95% +/- ~4.4% (binomial SE)
  status <- if (coverage_pct[p] >= 85) "OK" else "LOW"
  cat(sprintf("  %-15s %9.1f%% %11.1f%% %11.1f%% %8s\n",
              p, coverage_pct[p], mean_error[p], mean_ci_width[p], status))
}

# Group summaries
cat(sprintf("\n  --- Group Summaries ---\n"))
cat(sprintf("  Ne parameters:   mean coverage = %.1f%%  mean error = %.1f%%\n",
            mean(coverage_pct[ne_params]), mean(mean_error[ne_params])))
cat(sprintf("  Time parameters: mean coverage = %.1f%%  mean error = %.1f%%\n",
            mean(coverage_pct[time_params]), mean(mean_error[time_params])))
cat(sprintf("  Mig parameters:  mean coverage = %.1f%%  mean error = %.1f%%\n",
            mean(coverage_pct[mig_params]), mean(mean_error[mig_params])))

# Binomial confidence interval for coverage estimate
binom_se <- sqrt(0.95 * 0.05 / n_reps) * 100
cat(sprintf("\n  Expected coverage: 95%% +/- %.1f%% (binomial SE with %d reps)\n",
            binom_se, n_reps))

# Interpretation
cat("\n  --- Interpretation ---\n")

all_ok <- all(coverage_pct >= 85)
mig_ok <- all(coverage_pct[mig_params] >= 85)
ne_ok  <- all(coverage_pct[ne_params] >= 85)

if (all_ok) {
  cat("  All parameters have adequate coverage (>=85%).\n")
  cat("  The migration pipeline is correctly calibrated.\n")
  if (mean(mean_error[mig_params]) > 2 * mean(mean_error[ne_params])) {
    cat("  Note: Migration has higher estimation error than Ne, consistent\n")
    cat("  with a statistical power issue (not a bug).\n")
  }
} else if (ne_ok && !mig_ok) {
  cat("  WARNING: Migration parameters have low coverage while Ne is OK.\n")
  cat("  This could indicate a pipeline issue OR extreme power limitation.\n")
  cat("  Consider: narrower priors, more loci, or longer sequences.\n")
} else {
  cat("  WARNING: Multiple parameter groups have low coverage.\n")
  cat("  This may indicate an overall calibration issue or insufficient\n")
  cat("  reference table size.\n")
}

cat("\nDone!\n")
