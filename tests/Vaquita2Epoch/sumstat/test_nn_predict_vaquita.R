#!/usr/bin/env Rscript
#
# Run nn.predict() on the joint 3-param Vaquita model with conformal
# and mc_dropout methods. Generate posteriors and compare with true values.
#
# Usage: Rscript tests/test_nn_predict_vaquita.R
#
# Requires: keras, tensorflow

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
library(PipeMaster)

cat("=====================================================\n")
cat("  nn.predict() — Posteriors for Vaquita 3-param model\n")
cat("=====================================================\n\n")

# --- Paths ---
data_dir    <- "tests/data/Vaquita2Epoch/vaquita_100K"
results_dir <- "tests/tune_vaquita_results"

# --- Load data ---
reftable <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                       header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

cat(sprintf("Reftable: %d rows x %d cols\n", nrow(reftable), ncol(reftable)))

# --- Load observed data ---
load("tests/sumstats/vaquita_100K/segovia/shared_data.RData")
cat(sprintf("Loaded obs_raw (%d stats)\n\n", length(obs_raw)))

# --- Load joint model ---
joint_path <- file.path(results_dir, "joint_3param")
cat("Loading joint 3-param model...\n")
tune_result <- load.tune.result(joint_path)
cat(sprintf("val_loss: %.6f\n\n", tune_result$best_val_loss))

# --- Run nn.predict() ---
cat("Running nn.predict() with conformal + mc_dropout...\n\n")

t0 <- proc.time()
posterior <- nn.predict(
  tune.result = tune_result,
  observed    = obs_raw,
  reftable    = reftable,
  param.cols  = param_cols,
  type        = "sumstat",
  method      = c("conformal", "mc_dropout"),
  n_ensemble  = 5,
  n_mc        = 1000L,
  n_boot      = 0,
  max_epochs  = 500,
  cores       = 4L, gpus = 1L,
  seed        = 42,
  verbose     = TRUE
)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nnn.predict() completed in %.1f sec (%.1f min)\n\n", elapsed, elapsed / 60))

# --- Summary ---
cat("=====================================================\n")
cat("  RESULTS\n")
cat("=====================================================\n\n")

cat("--- Point estimate ---\n")
cat(sprintf("%-14s %8s %10s %10s\n", "Param", "True", "Estimate", "Error%"))
for (i in seq_along(param_cols)) {
  p <- param_cols[i]
  est <- posterior$point_estimate[i]
  err <- 100 * abs(est - true_vals[i]) / true_vals[i]
  cat(sprintf("%-14s %8.0f %10.0f %9.1f%%\n", p, true_vals[i], est, err))
}

# --- Per-method posteriors ---
methods <- c("conformal", "mc_dropout")
method_labels <- c(conformal = "Conformal", mc_dropout = "MC Dropout")

for (m in methods) {
  samples <- posterior[[m]]
  if (is.null(samples)) next

  cat(sprintf("\n--- %s posterior ---\n", method_labels[m]))
  cat(sprintf("%-14s %8s %10s %10s %10s %10s\n",
              "Param", "True", "Q5%", "Median", "Q95%", "Contains?"))
  for (i in seq_along(param_cols)) {
    p <- param_cols[i]
    s <- samples[, i]
    q05 <- quantile(s, 0.05)
    q50 <- quantile(s, 0.50)
    q95 <- quantile(s, 0.95)
    contains <- true_vals[i] >= q05 & true_vals[i] <= q95
    cat(sprintf("%-14s %8.0f %10.0f %10.0f %10.0f %10s\n",
                p, true_vals[i], q05, q50, q95, if (contains) "YES" else "NO"))
  }
}

# --- Save ---
save(posterior, true_vals, param_cols,
     file = file.path(results_dir, "posterior_results.RData"))
cat(sprintf("\nSaved: %s/posterior_results.RData\n", results_dir))

# --- Plot densities ---
pdf_file <- file.path(results_dir, "posterior_densities.pdf")
pdf(pdf_file, width = 12, height = 4 * length(param_cols))
par(mfrow = c(length(param_cols), 1), mar = c(4, 4, 3, 1))

for (i in seq_along(param_cols)) {
  p <- param_cols[i]
  all_vals <- c()
  colors <- c()
  labels <- c()

  # Collect samples from each method
  if (!is.null(posterior$conformal)) {
    all_vals <- c(all_vals, list(posterior$conformal[, i]))
    colors <- c(colors, "blue")
    labels <- c(labels, "Conformal")
  }
  if (!is.null(posterior$mc_dropout)) {
    all_vals <- c(all_vals, list(posterior$mc_dropout[, i]))
    colors <- c(colors, "red")
    labels <- c(labels, "MC Dropout")
  }

  # Get range
  all_flat <- unlist(all_vals)
  xlim <- range(all_flat[is.finite(all_flat)], na.rm = TRUE)

  # Plot
  plot(NULL, xlim = xlim, ylim = c(0, 1), xlab = p, ylab = "Density",
       main = sprintf("%s (true = %d)", p, true_vals[i]))

  max_d <- 0
  for (j in seq_along(all_vals)) {
    d <- density(all_vals[[j]], na.rm = TRUE)
    max_d <- max(max_d, max(d$y))
  }

  plot(NULL, xlim = xlim, ylim = c(0, max_d * 1.1), xlab = p, ylab = "Density",
       main = sprintf("%s (true = %d)", p, true_vals[i]))

  for (j in seq_along(all_vals)) {
    d <- density(all_vals[[j]], na.rm = TRUE)
    lines(d, col = colors[j], lwd = 2)
  }
  abline(v = true_vals[i], col = "black", lwd = 2, lty = 2)
  abline(v = posterior$point_estimate[i], col = "gray40", lwd = 1.5, lty = 3)
  legend("topright", legend = c(labels, "True value", "Point est."),
         col = c(colors, "black", "gray40"),
         lwd = c(rep(2, length(labels)), 2, 1.5),
         lty = c(rep(1, length(labels)), 2, 3), cex = 0.8)
}
dev.off()
cat(sprintf("Saved: %s\n", pdf_file))

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
