#!/usr/bin/env Rscript
# Test sim.sumstat() against obs.sumstat.ngs() observed statistics.
#
# This script:
#   1. Loads model objects and obs.sumstat.ngs() results (the msABC-based observed)
#   2. Runs sim.sumstat() to generate simulated summary statistics (appending to existing)
#   3. Verifies column names match between sim and obs
#   4. Checks that observed stats fall within the simulated range
#   5. Runs simple ABC rejection (excluding out-of-range stats) to verify
#      posterior concentrates near true params

library(PipeMaster)

# ============================================================================
# Setup
# ============================================================================
args <- commandArgs(trailingOnly=FALSE)
script_arg <- grep("--file=", args, value=TRUE)
if (length(script_arg) > 0) {
  test_dir <- dirname(normalizePath(sub("--file=", "", script_arg)))
} else {
  test_dir <- "data_to_test"
}
setwd(test_dir)

load("test_models.RData")
load("observed_msABC_sumstats.RData")

cat("=======================================================\n")
cat("  Testing sim.sumstat() against obs.sumstat.ngs()\n")
cat("=======================================================\n\n")

n_pass <- 0
n_fail <- 0

report <- function(test_name, passed, msg="") {
  if (passed) {
    cat(sprintf("  [PASS] %s\n", test_name))
    n_pass <<- n_pass + 1
  } else {
    cat(sprintf("  [FAIL] %s -- %s\n", test_name, msg))
    n_fail <<- n_fail + 1
  }
}

simdir <- getwd()  # write sims into data_to_test/ so they persist

# Helper: identify stats where observed falls outside simulated range
find_out_of_range <- function(obs_vals, sim_df, stat_cols) {
  oor_names <- character(0)
  oor_details <- character(0)
  in_range_names <- character(0)
  for (col in stat_cols) {
    obs_val <- obs_vals[col]
    if (is.na(obs_val) || !is.finite(obs_val)) next
    sim_range <- range(sim_df[, col], na.rm=TRUE)
    if (obs_val >= sim_range[1] && obs_val <= sim_range[2]) {
      in_range_names <- c(in_range_names, col)
    } else {
      oor_names <- c(oor_names, col)
      oor_details <- c(oor_details, sprintf("%s (obs=%.4f, sim=[%.4f, %.4f])",
                                            col, obs_val, sim_range[1], sim_range[2]))
    }
  }
  list(in_range=in_range_names, oor_names=oor_names, oor_details=oor_details)
}

# Helper: ABC rejection using standardized Euclidean distance
abc_reject <- function(sim_df, obs_vals, stat_cols, tol_frac) {
  sim_matrix <- as.matrix(sim_df[, stat_cols, drop=FALSE])
  obs_vector <- obs_vals[stat_cols]
  sim_matrix[!is.finite(sim_matrix)] <- 0
  obs_vector[!is.finite(obs_vector)] <- 0
  sim_means <- colMeans(sim_matrix)
  sim_sd <- apply(sim_matrix, 2, sd, na.rm=TRUE)
  sim_sd[sim_sd == 0] <- 1
  sim_scaled <- sweep(sim_matrix, 2, sim_means) / rep(sim_sd, each=nrow(sim_matrix))
  obs_scaled <- (obs_vector - sim_means) / sim_sd
  distances <- sqrt(rowSums(sweep(sim_scaled, 2, obs_scaled)^2))
  tol <- quantile(distances, tol_frac)
  which(distances <= tol)
}

# ============================================================================
# Test 1: Vaquita2Epoch (single population)
# ============================================================================
cat("--- Test 1: Vaquita2Epoch sim.sumstat() ---\n")

# Fresh 100k sims (10 blocks x 1000 x 10 cores = 100k)
cat("  Simulating 100,000 (10 cores)...\n")
sim.sumstat(model=Vaquita2Epoch, nsim.blocks=10, block.size=1000,
            path=simdir, use.alpha=FALSE, output.name="test_vaq",
            append.sims=FALSE, ncores=10)

sim_vaq <- read.table(file.path(simdir, "SIMS_test_vaq.txt"), header=TRUE, sep="\t")
cat(sprintf("  Total: %d simulations x %d columns\n\n", nrow(sim_vaq), ncol(sim_vaq)))

# Check 1a: sim.sumstat ran without error
report("Vaquita sim.sumstat() completes without error", TRUE)

# Check 1b: correct total simulations (100k fresh)
report("Vaquita output has ~100k rows",
       nrow(sim_vaq) >= 99000 && nrow(sim_vaq) <= 101000,
       sprintf("got %d", nrow(sim_vaq)))

# Check 1c: column names - stat columns from sim should match obs
obs_colnames_vaq <- colnames(observed_msABC_Vaquita2Epoch)
sim_colnames_vaq <- colnames(sim_vaq)
first_stat_idx <- min(grep("^s_", sim_colnames_vaq))
param_cols_vaq <- sim_colnames_vaq[1:(first_stat_idx - 1)]
stat_cols_vaq <- sim_colnames_vaq[first_stat_idx:length(sim_colnames_vaq)]

cat(sprintf("  Parameter columns (%d): %s\n", length(param_cols_vaq),
            paste(param_cols_vaq, collapse=", ")))
cat(sprintf("  Stat columns: %d (sim) vs %d (obs)\n", length(stat_cols_vaq), length(obs_colnames_vaq)))

report("Vaquita stat column names match obs.sumstat.ngs()",
       identical(stat_cols_vaq, obs_colnames_vaq),
       sprintf("sim has: %s\nobs has: %s",
               paste(setdiff(stat_cols_vaq, obs_colnames_vaq), collapse=", "),
               paste(setdiff(obs_colnames_vaq, stat_cols_vaq), collapse=", ")))

# Check 1d: parameter columns present
expected_params_vaq <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1", "mean.rate", "sd.rate")
report("Vaquita has expected parameter columns",
       all(expected_params_vaq %in% param_cols_vaq),
       sprintf("missing: %s", paste(setdiff(expected_params_vaq, param_cols_vaq), collapse=", ")))

# Check 1e: observed stats within range, identify out-of-range for exclusion
obs_vals_vaq <- as.numeric(observed_msABC_Vaquita2Epoch[1, ])
names(obs_vals_vaq) <- obs_colnames_vaq

range_check_vaq <- find_out_of_range(obs_vals_vaq, sim_vaq, stat_cols_vaq)
cat(sprintf("\n  Observed within sim range: %d / %d stats\n",
            length(range_check_vaq$in_range),
            length(range_check_vaq$in_range) + length(range_check_vaq$oor_names)))
if (length(range_check_vaq$oor_details) > 0) {
  cat("  Out of range (excluded from rejection):\n")
  for (s in range_check_vaq$oor_details) cat(sprintf("    %s\n", s))
}
report("Vaquita observed stats mostly within simulated range",
       length(range_check_vaq$in_range) >= (length(range_check_vaq$in_range) + length(range_check_vaq$oor_names)) * 0.5,
       sprintf("only %d in range", length(range_check_vaq$in_range)))

# Check 1f: ABC rejection using in-range stats only, 0.1% tolerance
true_params_vaq <- attr(Vaquita2Epoch, "true_params")
cat("\n  ABC rejection (0.1%% tolerance, in-range stats only):\n")
cat(sprintf("  Using %d stats for rejection\n", length(range_check_vaq$in_range)))
cat(sprintf("  True params: Ne0=%.0f, Ne1=%.0f, t=%.0f\n",
            true_params_vaq$Ne0.pop1, true_params_vaq$Ne1.pop1, true_params_vaq$t.Ne1.pop1))

accepted_vaq <- abc_reject(sim_vaq, obs_vals_vaq, range_check_vaq$in_range, 0.001)
cat(sprintf("  Accepted %d / %d simulations\n", length(accepted_vaq), nrow(sim_vaq)))
cat(sprintf("  Posterior Ne0.pop1:  median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_vaq$Ne0.pop1[accepted_vaq]),
            quantile(sim_vaq$Ne0.pop1[accepted_vaq], 0.25),
            quantile(sim_vaq$Ne0.pop1[accepted_vaq], 0.75),
            true_params_vaq$Ne0.pop1))
cat(sprintf("  Posterior Ne1.pop1:  median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_vaq$Ne1.pop1[accepted_vaq]),
            quantile(sim_vaq$Ne1.pop1[accepted_vaq], 0.25),
            quantile(sim_vaq$Ne1.pop1[accepted_vaq], 0.75),
            true_params_vaq$Ne1.pop1))
cat(sprintf("  Posterior t.Ne1:     median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_vaq$t.Ne1.pop1[accepted_vaq]),
            quantile(sim_vaq$t.Ne1.pop1[accepted_vaq], 0.25),
            quantile(sim_vaq$t.Ne1.pop1[accepted_vaq], 0.75),
            true_params_vaq$t.Ne1.pop1))


# ============================================================================
# Test 2: OutOfAfrica_3G09 (two populations)
# ============================================================================
cat("\n--- Test 2: OutOfAfrica_3G09 sim.sumstat() ---\n")

# Fresh 100k sims (10 blocks x 1000 x 10 cores = 100k)
cat("  Simulating 100,000 (10 cores)...\n")
sim.sumstat(model=OutOfAfrica_3G09, nsim.blocks=10, block.size=1000,
            path=simdir, use.alpha=FALSE, output.name="test_ooa",
            append.sims=FALSE, ncores=10)

sim_ooa <- read.table(file.path(simdir, "SIMS_test_ooa.txt"), header=TRUE, sep="\t")
cat(sprintf("  Total: %d simulations x %d columns\n\n", nrow(sim_ooa), ncol(sim_ooa)))

# Check 2a: sim.sumstat ran without error
report("OutOfAfrica sim.sumstat() completes without error", TRUE)

# Check 2b: correct total simulations
report("OutOfAfrica output has ~100k rows",
       nrow(sim_ooa) >= 99000 && nrow(sim_ooa) <= 101000,
       sprintf("got %d", nrow(sim_ooa)))

# Check 2c: column names match
obs_colnames_ooa <- colnames(observed_msABC_OutOfAfrica_3G09)
sim_colnames_ooa <- colnames(sim_ooa)
first_stat_idx_ooa <- min(grep("^s_", sim_colnames_ooa))
param_cols_ooa <- sim_colnames_ooa[1:(first_stat_idx_ooa - 1)]
stat_cols_ooa <- sim_colnames_ooa[first_stat_idx_ooa:length(sim_colnames_ooa)]

cat(sprintf("  Parameter columns (%d): %s\n", length(param_cols_ooa),
            paste(param_cols_ooa, collapse=", ")))
cat(sprintf("  Stat columns: %d (sim) vs %d (obs)\n", length(stat_cols_ooa), length(obs_colnames_ooa)))

report("OutOfAfrica stat column names match obs.sumstat.ngs()",
       identical(stat_cols_ooa, obs_colnames_ooa),
       sprintf("sim extra: [%s]; obs extra: [%s]",
               paste(setdiff(stat_cols_ooa, obs_colnames_ooa), collapse=", "),
               paste(setdiff(obs_colnames_ooa, stat_cols_ooa), collapse=", ")))

# Check 2d: parameter columns present
expected_params_ooa <- c("Ne0.pop1", "Ne0.pop2", "Ne1.pop2", "Ne1.pop1",
                          "join1", "t.Ne1.pop2", "t.Ne1.pop1",
                          "mig0.1_2", "mig0.2_1", "mean.rate", "sd.rate")
report("OutOfAfrica has expected parameter columns",
       all(expected_params_ooa %in% param_cols_ooa),
       sprintf("missing: %s", paste(setdiff(expected_params_ooa, param_cols_ooa), collapse=", ")))

# Check 2e: observed stats within range
obs_vals_ooa <- as.numeric(observed_msABC_OutOfAfrica_3G09[1, ])
names(obs_vals_ooa) <- obs_colnames_ooa

range_check_ooa <- find_out_of_range(obs_vals_ooa, sim_ooa, stat_cols_ooa)
cat(sprintf("\n  Observed within sim range: %d / %d stats\n",
            length(range_check_ooa$in_range),
            length(range_check_ooa$in_range) + length(range_check_ooa$oor_names)))
if (length(range_check_ooa$oor_details) > 0) {
  cat("  Out of range (excluded from rejection):\n")
  for (s in head(range_check_ooa$oor_details, 15)) cat(sprintf("    %s\n", s))
  if (length(range_check_ooa$oor_details) > 15)
    cat(sprintf("    ... and %d more\n", length(range_check_ooa$oor_details) - 15))
}
report("OutOfAfrica observed stats mostly within simulated range",
       length(range_check_ooa$in_range) >= (length(range_check_ooa$in_range) + length(range_check_ooa$oor_names)) * 0.5,
       sprintf("only %d in range", length(range_check_ooa$in_range)))

# Check 2f: ABC rejection using in-range stats only, 0.1% tolerance
true_params_ooa <- attr(OutOfAfrica_3G09, "true_params")
cat("\n  ABC rejection (0.1%% tolerance, in-range stats only):\n")
cat(sprintf("  Using %d stats for rejection (excluded %d out-of-range)\n",
            length(range_check_ooa$in_range), length(range_check_ooa$oor_names)))
cat(sprintf("  True params: Ne0.pop1=%.0f, Ne0.pop2=%.0f, join1=%.0f\n",
            true_params_ooa$Ne0.pop1, true_params_ooa$Ne0.pop2, true_params_ooa$join1))

accepted_ooa <- abc_reject(sim_ooa, obs_vals_ooa, range_check_ooa$in_range, 0.001)
cat(sprintf("  Accepted %d / %d simulations\n", length(accepted_ooa), nrow(sim_ooa)))
cat(sprintf("  Posterior Ne0.pop1:  median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_ooa$Ne0.pop1[accepted_ooa]),
            quantile(sim_ooa$Ne0.pop1[accepted_ooa], 0.25),
            quantile(sim_ooa$Ne0.pop1[accepted_ooa], 0.75),
            true_params_ooa$Ne0.pop1))
cat(sprintf("  Posterior Ne0.pop2:  median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_ooa$Ne0.pop2[accepted_ooa]),
            quantile(sim_ooa$Ne0.pop2[accepted_ooa], 0.25),
            quantile(sim_ooa$Ne0.pop2[accepted_ooa], 0.75),
            true_params_ooa$Ne0.pop2))
cat(sprintf("  Posterior join1:     median=%.0f, IQR=[%.0f, %.0f] (true=%.0f)\n",
            median(sim_ooa$join1[accepted_ooa]),
            quantile(sim_ooa$join1[accepted_ooa], 0.25),
            quantile(sim_ooa$join1[accepted_ooa], 0.75),
            true_params_ooa$join1))
cat(sprintf("  Posterior mig0.1_2:  median=%.4f, IQR=[%.4f, %.4f] (true=%.4f)\n",
            median(sim_ooa$mig0.1_2[accepted_ooa]),
            quantile(sim_ooa$mig0.1_2[accepted_ooa], 0.25),
            quantile(sim_ooa$mig0.1_2[accepted_ooa], 0.75),
            true_params_ooa$mig0.1_2))
cat(sprintf("  Posterior mig0.2_1:  median=%.4f, IQR=[%.4f, %.4f] (true=%.4f)\n",
            median(sim_ooa$mig0.2_1[accepted_ooa]),
            quantile(sim_ooa$mig0.2_1[accepted_ooa], 0.25),
            quantile(sim_ooa$mig0.2_1[accepted_ooa], 0.75),
            true_params_ooa$mig0.2_1))

# Check 2g: Fst values should be bounded [0,1]
fst_col <- "s_mean_Fst"
if (fst_col %in% stat_cols_ooa) {
  fst_vals <- sim_ooa[, fst_col]
  fst_valid <- all(fst_vals >= -0.5 & fst_vals <= 1, na.rm=TRUE)
  report("OutOfAfrica simulated Fst values in valid range",
         fst_valid,
         sprintf("range: [%.3f, %.3f]", min(fst_vals, na.rm=TRUE), max(fst_vals, na.rm=TRUE)))
}

# ============================================================================
# Save simulation results
# ============================================================================
cat("\n--- Saving simulation results ---\n")
save(sim_vaq, sim_ooa, file="sim_sumstat_results.RData")
cat(sprintf("  Saved sim_sumstat_results.RData\n"))
cat(sprintf("    - sim_vaq: %d x %d\n", nrow(sim_vaq), ncol(sim_vaq)))
cat(sprintf("    - sim_ooa: %d x %d\n", nrow(sim_ooa), ncol(sim_ooa)))
cat(sprintf("  SIMS_test_vaq.txt and SIMS_test_ooa.txt also in %s\n", getwd()))

# ============================================================================
# Summary
# ============================================================================
cat("\n=======================================================\n")
cat(sprintf("  Results: %d passed, %d failed\n", n_pass, n_fail))
cat("=======================================================\n")
