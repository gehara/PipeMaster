#!/usr/bin/env Rscript
#
# Test emulator.posterior() MCMC on the Vaquita dataset.
# Uses saved emulator from emulator_vaquita_results/ and observed stats.
#
# The Vaquita2Epoch model is a 1-population bottleneck with 3 parameters:
#   Ne0.pop1     = 2807  (current Ne)
#   Ne1.pop1     = 4485  (ancestral Ne)
#   t.Ne1.pop1   = 2162  (bottleneck time in generations)
#
# Usage (run from tests/):
#   cd tests && Rscript test_emulator_posterior_vaquita.R
#   cd tests && Rscript test_emulator_posterior_vaquita.R 2>&1 | tee test_emulator_posterior_vaquita.log
#

devtools::load_all("..", quiet = TRUE)

cat("========================================\n")
cat("  Emulator MCMC Posterior — Vaquita2Epoch\n")
cat("========================================\n\n")

# =====================================================================
# Load saved emulator + observed stats
# =====================================================================
cat("=== 1. Loading saved emulator ===\n")
emulator <- load.emulator.result("emulator_vaquita_results")

cat("\n=== 2. Loading observed stats ===\n")
load("obs_vaquita_100Kloci.RData")  # obs_raw, stat_cols, true_vals, param_cols

cat(sprintf("Observed: %d summary statistics\n", length(obs_raw)))
cat(sprintf("Parameters: %s\n", paste(param_cols, collapse = ", ")))
cat(sprintf("True values: %s\n",
            paste(sprintf("%s=%.0f", names(true_vals), true_vals), collapse = " ")))

# =====================================================================
# Load model for prior bounds
# =====================================================================
data_dir <- file.path("data", "Vaquita2Epoch", "vaquita_100K")
load(file.path(data_dir, "model.RData"))
model <- Vaquita2Epoch_100K

ptab <- get.prior.table(model)
cat("\nPrior bounds:\n")
print(ptab[ptab$Parameter %in% param_cols, ])

# =====================================================================
# 3. Quick run (5K samples, 1K burn-in)
# =====================================================================
cat("\n=== 3. Quick MCMC run (5K samples) ===\n")

t0 <- proc.time()
mcmc <- emulator.posterior(
  emulator,
  observed  = obs_raw,
  model     = model,
  n_samples = 5000L,
  burnin    = 1000L,
  thin      = 1L,
  n_chains  = 10L,
  adaptive  = TRUE,
  verbose   = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nQuick run completed in %.1f sec (%.1f min)\n",
            elapsed, elapsed / 60))

# =====================================================================
# 4. Summary
# =====================================================================
cat("\n=== 4. Posterior summary ===\n")
print(summary(mcmc))

# =====================================================================
# 5. Compare to true values and other methods
# =====================================================================
cat("\n=== 5. Comparison ===\n")

# Gradient-based estimate
opt <- emulator.optimize(emulator, obs_raw, model = model,
                         n_starts = 50L, n_steps = 1000L,
                         verbose = FALSE)

# ABC estimate (if available)
abc_est <- NULL
if (file.exists("emulator_abc_posteriors_vaquita.RData")) {
  load("emulator_abc_posteriors_vaquita.RData")
  abc_summ <- summary(posteriors[["tol=0.010"]])$abc_rejection
  abc_est  <- abc_summ[, "Mean"]
}

cat(sprintf("\n%-12s %8s %12s %12s",
            "Param", "True", "MCMC Mean", "Grad. Opt."))
if (!is.null(abc_est)) cat(sprintf(" %12s", "ABC Mean"))
cat("\n")

for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %12.1f %12.1f",
              param_cols[i], true_vals[i],
              mcmc$point_estimate[i], opt$point_estimate[i]))
  if (!is.null(abc_est))
    cat(sprintf(" %12.1f", abc_est[i]))
  cat("\n")
}

cat(sprintf("\n--- Percent errors (MCMC) ---\n"))
for (i in seq_along(param_cols)) {
  err <- 100 * abs(mcmc$point_estimate[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("  %-12s: %.1f%%\n", param_cols[i], err))
}
mean_err <- mean(abs(mcmc$point_estimate - true_vals) / true_vals) * 100
cat(sprintf("  %-12s: %.1f%%\n", "Mean", mean_err))

# =====================================================================
# 6. Diagnostics
# =====================================================================
cat("\n=== 6. Diagnostics ===\n")
cat(sprintf("Chains: %d\n", mcmc$n_chains))
cat(sprintf("Acceptance rate: %.1f%% (mean across chains)\n", 100 * mcmc$acceptance_rate))
cat(sprintf("Per-chain acceptance: %s\n",
            paste(sprintf("%.0f%%", 100 * mcmc$chain_acceptance), collapse = " ")))
if (mcmc$acceptance_rate < 0.15 || mcmc$acceptance_rate > 0.35)
  cat("WARNING: acceptance rate outside [15%, 35%]\n")

cat(sprintf("ESS: %s\n",
            paste(sprintf("%s=%.0f", param_cols, mcmc$ess), collapse = " ")))
low_ess <- mcmc$ess[mcmc$ess < 200]
if (length(low_ess) > 0)
  cat(sprintf("WARNING: low ESS for: %s\n", paste(names(low_ess), collapse = ", ")))

if (!is.null(mcmc$rhat) && any(!is.na(mcmc$rhat))) {
  cat(sprintf("R-hat: %s\n",
              paste(sprintf("%s=%.3f", param_cols, mcmc$rhat), collapse = " ")))
  bad_rhat <- mcmc$rhat[!is.na(mcmc$rhat) & mcmc$rhat > 1.1]
  if (length(bad_rhat) > 0)
    cat(sprintf("WARNING: R-hat > 1.1 for: %s\n",
                paste(names(bad_rhat), collapse = ", ")))
}

# Check all samples within bounds
in_bounds <- TRUE
for (i in seq_along(param_cols)) {
  lo <- as.numeric(mcmc$prior_bounds$prior.1[i])
  hi <- as.numeric(mcmc$prior_bounds$prior.2[i])
  oob <- sum(mcmc$samples[, i] < lo | mcmc$samples[, i] > hi)
  if (oob > 0) {
    cat(sprintf("WARNING: %d samples out of bounds for %s\n",
                oob, param_cols[i]))
    in_bounds <- FALSE
  }
}
if (in_bounds)
  cat("All samples within prior bounds: OK\n")

# =====================================================================
# 7. Plot quick run to PDF
# =====================================================================
cat("\n=== 7. Saving quick-run plots ===\n")
pdf("emulator_mcmc_posterior_vaquita.pdf", width = 14, height = 12)
plot(mcmc, true_values = true_vals)
dev.off()
cat("Saved: emulator_mcmc_posterior_vaquita.pdf\n")

# =====================================================================
# 8. Production run (50K samples, 5K burn-in)
# =====================================================================
cat("\n=== 8. Production MCMC run (50K samples) ===\n")

t0 <- proc.time()
mcmc_long <- emulator.posterior(
  emulator,
  observed  = obs_raw,
  model     = model,
  n_samples = 50000L,
  burnin    = 5000L,
  thin      = 1L,
  n_chains  = 10L,
  sigma2    = mcmc$sigma2,   # reuse sigma2 from quick run
  adaptive  = TRUE,
  verbose   = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nProduction run completed in %.1f sec (%.1f min)\n",
            elapsed, elapsed / 60))

cat("\n=== 9. Production posterior summary ===\n")
print(summary(mcmc_long))

# Compare quick vs production
cat("\n=== 10. Quick vs Production comparison ===\n")
cat(sprintf("\n%-12s %8s %12s %12s %12s",
            "Param", "True", "Quick(5K)", "Prod(50K)", "Grad. Opt."))
if (!is.null(abc_est)) cat(sprintf(" %12s", "ABC Mean"))
cat("\n")
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %12.1f %12.1f %12.1f",
              param_cols[i], true_vals[i],
              mcmc$point_estimate[i],
              mcmc_long$point_estimate[i],
              opt$point_estimate[i]))
  if (!is.null(abc_est))
    cat(sprintf(" %12.1f", abc_est[i]))
  cat("\n")
}

cat(sprintf("\n%-12s %12s %12s\n", "Param", "Quick err%", "Prod err%"))
for (i in seq_along(param_cols)) {
  err_q <- 100 * abs(mcmc$point_estimate[i] - true_vals[i]) / true_vals[i]
  err_p <- 100 * abs(mcmc_long$point_estimate[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %11.1f%% %11.1f%%\n", param_cols[i], err_q, err_p))
}
mean_q <- mean(abs(mcmc$point_estimate - true_vals) / true_vals) * 100
mean_p <- mean(abs(mcmc_long$point_estimate - true_vals) / true_vals) * 100
cat(sprintf("%-12s %11.1f%% %11.1f%%\n", "Mean", mean_q, mean_p))

cat(sprintf("\n%-12s %12s %12s\n", "Diagnostic", "Quick(5K)", "Prod(50K)"))
cat(sprintf("%-12s %12.1f%% %11.1f%%\n", "Acceptance",
            100 * mcmc$acceptance_rate, 100 * mcmc_long$acceptance_rate))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %12.0f %12.0f    (ESS)\n",
              param_cols[i], mcmc$ess[i], mcmc_long$ess[i]))
}
if (!is.null(mcmc_long$rhat) && any(!is.na(mcmc_long$rhat))) {
  for (i in seq_along(param_cols)) {
    cat(sprintf("%-12s %12.3f %12.3f    (R-hat)\n",
                param_cols[i], mcmc$rhat[i], mcmc_long$rhat[i]))
  }
}

# =====================================================================
# 11. Plot production run to PDF
# =====================================================================
cat("\n=== 11. Saving production plots ===\n")
pdf("emulator_mcmc_posterior_vaquita_50K.pdf", width = 14, height = 12)
plot(mcmc_long, true_values = true_vals)
dev.off()
cat("Saved: emulator_mcmc_posterior_vaquita_50K.pdf\n")

# =====================================================================
# 12. Save results
# =====================================================================
save(mcmc, mcmc_long, true_vals, param_cols,
     file = "test_emulator_posterior_vaquita_results.RData")
cat("Saved: test_emulator_posterior_vaquita_results.RData\n")

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
