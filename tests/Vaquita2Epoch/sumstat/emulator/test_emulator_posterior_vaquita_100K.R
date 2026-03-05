#!/usr/bin/env Rscript
#
# Test emulator.MCMC() MCMC on Vaquita with 100K emulator.
# Compares 10K vs 100K emulator posteriors side by side.
#
# Usage (run from project root on segovia):
#   cd ~/Dropbox/github/PipeMaster
#   nohup Rscript tests/test_emulator_posterior_vaquita_100K.R > tests/test_emulator_posterior_vaquita_100K.log 2>&1 &
#

devtools::load_all()

cat("=====================================================\n")
cat("  Emulator MCMC Posterior — Vaquita 10K vs 100K\n")
cat("=====================================================\n\n")

param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

# =====================================================================
# 1. Load both emulators
# =====================================================================
cat("=== 1. Loading emulators ===\n")

emu_10K  <- load.emulator.result("tests/emulator_vaquita_results")
cat(sprintf("10K emulator: val_loss=%.6f\n", emu_10K$best_val_loss))

emu_100K <- load.emulator.result("tests/emulator_vaquita_100K_results")
cat(sprintf("100K emulator: val_loss=%.6f\n", emu_100K$best_val_loss))

# =====================================================================
# 2. Load observed stats + model
# =====================================================================
cat("\n=== 2. Loading observed stats ===\n")

data_dir <- file.path("tests", "data", "Vaquita2Epoch", "vaquita_100K")
load(file.path(data_dir, "model.RData"))
model <- Vaquita2Epoch_100K

pop_assign <- read.table(file.path(data_dir, "pop_assign_Vaquita2Epoch.txt"),
                         header = TRUE)
obs_ss <- obs.sumstat.ngs(
  model = model,
  path.to.phylip = file.path(data_dir, "phylip_Vaquita2Epoch.phy"),
  pop.assign = pop_assign
)

reftable <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                       header = TRUE, sep = "\t")
nuisance  <- c("mean.rate", "sd.rate")
stat_cols <- setdiff(colnames(reftable), c(param_cols, nuisance))
obs_raw   <- as.numeric(obs_ss[1, stat_cols])

cat(sprintf("Observed: %d summary statistics\n", length(obs_raw)))
cat(sprintf("True values: %s\n",
            paste(sprintf("%s=%.0f", names(true_vals), true_vals), collapse = " ")))

ptab <- get.prior.table(model)
cat("\nPrior bounds:\n")
print(ptab[ptab$Parameter %in% param_cols, ])

# =====================================================================
# 3. MCMC — 10K emulator (50K samples)
# =====================================================================
cat("\n=== 3. MCMC with 10K emulator (50K samples, 10 chains) ===\n")

t0 <- proc.time()
mcmc_10K <- emulator.MCMC(
  emu_10K,
  observed  = obs_raw,
  model     = model,
  n_samples = 50000L,
  burnin    = 5000L,
  thin      = 1L,
  n_chains  = 10L,
  adaptive  = TRUE,
  verbose   = TRUE
)
elapsed_10K <- (proc.time() - t0)[3]
cat(sprintf("\n10K MCMC completed in %.1f sec (%.1f min)\n",
            elapsed_10K, elapsed_10K / 60))

# =====================================================================
# 4. MCMC — 100K emulator (50K samples)
# =====================================================================
cat("\n=== 4. MCMC with 100K emulator (50K samples, 10 chains) ===\n")

t0 <- proc.time()
mcmc_100K <- emulator.MCMC(
  emu_100K,
  observed  = obs_raw,
  model     = model,
  n_samples = 50000L,
  burnin    = 5000L,
  thin      = 1L,
  n_chains  = 10L,
  adaptive  = TRUE,
  verbose   = TRUE
)
elapsed_100K <- (proc.time() - t0)[3]
cat(sprintf("\n100K MCMC completed in %.1f sec (%.1f min)\n",
            elapsed_100K, elapsed_100K / 60))

# =====================================================================
# 5. Summaries
# =====================================================================
cat("\n=== 5. Posterior summaries ===\n")

cat("\n--- 10K emulator ---\n")
print(summary(mcmc_10K))

cat("\n--- 100K emulator ---\n")
print(summary(mcmc_100K))

# =====================================================================
# 6. Comparison table
# =====================================================================
cat("\n=== 6. Comparison ===\n")

cat(sprintf("\n%-12s %8s %12s %12s\n",
            "Param", "True", "10K Mean", "100K Mean"))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %12.1f %12.1f\n",
              param_cols[i], true_vals[i],
              mcmc_10K$point_estimate[i],
              mcmc_100K$point_estimate[i]))
}

cat(sprintf("\n%-12s %12s %12s\n", "Param", "10K err%", "100K err%"))
for (i in seq_along(param_cols)) {
  err_10  <- 100 * abs(mcmc_10K$point_estimate[i]  - true_vals[i]) / true_vals[i]
  err_100 <- 100 * abs(mcmc_100K$point_estimate[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %11.1f%% %11.1f%%\n", param_cols[i], err_10, err_100))
}
mean_10  <- mean(abs(mcmc_10K$point_estimate  - true_vals) / true_vals) * 100
mean_100 <- mean(abs(mcmc_100K$point_estimate - true_vals) / true_vals) * 100
cat(sprintf("%-12s %11.1f%% %11.1f%%\n", "Mean", mean_10, mean_100))

# =====================================================================
# 7. Diagnostics
# =====================================================================
cat("\n=== 7. Diagnostics ===\n")

cat(sprintf("\n%-12s %12s %12s\n", "Diagnostic", "10K", "100K"))
cat(sprintf("%-12s %12.1f%% %11.1f%%\n", "Acceptance",
            100 * mcmc_10K$acceptance_rate, 100 * mcmc_100K$acceptance_rate))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %12.0f %12.0f    (ESS)\n",
              param_cols[i], mcmc_10K$ess[i], mcmc_100K$ess[i]))
}
if (!is.null(mcmc_100K$rhat) && any(!is.na(mcmc_100K$rhat))) {
  for (i in seq_along(param_cols)) {
    cat(sprintf("%-12s %12.3f %12.3f    (R-hat)\n",
                param_cols[i], mcmc_10K$rhat[i], mcmc_100K$rhat[i]))
  }
}

# =====================================================================
# 8. Plots
# =====================================================================
cat("\n=== 8. Saving plots ===\n")

# Individual plots
pdf("tests/emulator_mcmc_posterior_10K.pdf", width = 14, height = 12)
plot(mcmc_10K, true_values = true_vals)
dev.off()
cat("Saved: tests/emulator_mcmc_posterior_10K.pdf\n")

pdf("tests/emulator_mcmc_posterior_100K.pdf", width = 14, height = 12)
plot(mcmc_100K, true_values = true_vals)
dev.off()
cat("Saved: tests/emulator_mcmc_posterior_100K.pdf\n")

# Overlay comparison plot
pdf("tests/emulator_mcmc_comparison_10K_vs_100K.pdf", width = 14, height = 5)
par(mfrow = c(1, length(param_cols)), mar = c(4, 4, 3, 1))

prior_tab <- ptab[match(param_cols, ptab$Parameter), ]

for (i in seq_along(param_cols)) {
  p <- param_cols[i]
  lo <- as.numeric(prior_tab$prior.1[i])
  hi <- as.numeric(prior_tab$prior.2[i])

  s_10  <- mcmc_10K$samples[, i]
  s_100 <- mcmc_100K$samples[, i]

  d_10  <- density(s_10,  from = lo, to = hi)
  d_100 <- density(s_100, from = lo, to = hi)

  ylim <- c(0, max(d_10$y, d_100$y) * 1.1)

  plot(d_100, main = p, xlab = p, ylab = "Density",
       col = "red", lwd = 2, xlim = c(lo, hi), ylim = ylim)
  lines(d_10, col = "blue", lwd = 2, lty = 2)
  abline(v = true_vals[i], col = "black", lwd = 2, lty = 3)

  legend("topright",
         legend = c("100K emulator", "10K emulator", "True value"),
         col = c("red", "blue", "black"),
         lwd = c(2, 2, 2), lty = c(1, 2, 3), cex = 0.8)
}
dev.off()
cat("Saved: tests/emulator_mcmc_comparison_10K_vs_100K.pdf\n")

# =====================================================================
# 9. Save results
# =====================================================================
save(mcmc_10K, mcmc_100K, true_vals, param_cols,
     file = "tests/test_emulator_posterior_vaquita_100K_results.RData")
cat("Saved: tests/test_emulator_posterior_vaquita_100K_results.RData\n")

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
