#!/usr/bin/env Rscript
#
# Step 2: Run nn.predict() on segovia (32 cores) using saved tune.nn result.
# Expects tune result saved at tests/hABC/data/Ubelmania/tune_result/
#
# Usage: Rscript tests/test_nn_predict.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

cat("========================================\n")
cat("  nn.predict() — Ubelmania pop 1\n")
cat("========================================\n\n")

# --- Load saved tune result ---
tune_dir <- "tests/hABC/data/Ubelmania/tune_result"
tune_result <- load.tune.result(tune_dir)
cat(sprintf("Best val_loss: %.6f\n", tune_result$best_val_loss))
cat(sprintf("Type: %s\n\n", tune_result$type))

# --- Load data ---
reftable <- read.table("tests/hABC/data/Ubelmania/abc_sims/SIMS_pop_1.txt",
                       header = TRUE, sep = "\t")
observed <- read.table("tests/hABC/data/Ubelmania/grid_search/obs_pop_1.txt",
                       header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")

cat(sprintf("Reftable: %d rows x %d cols\n", nrow(reftable), ncol(reftable)))
cat(sprintf("Observed: %d cols\n\n", ncol(observed)))

# --- nn.predict ---
t0 <- proc.time()
post <- nn.predict(tune_result, observed,
                   reftable = reftable,
                   param.cols = param_cols,
                   method = c("conformal", "bootstrap"),
                   n_boot = 20,
                   n_ensemble = 10,
                   max_epochs = 500,
                   cores = 32,
                   verbose = TRUE)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\nnn.predict completed in %.1f sec (%.1f min)\n\n", elapsed, elapsed / 60))

# --- Results ---
cat("=== Results ===\n")
cat("Point estimate:\n")
print(post$point_estimate)

cat(sprintf("\nConformal samples: %s\n",
            if (!is.null(post$conformal)) paste(dim(post$conformal), collapse = " x ") else "NULL"))
cat(sprintf("Bootstrap samples: %s\n",
            if (!is.null(post$bootstrap)) paste(dim(post$bootstrap), collapse = " x ") else "NULL"))

if (!is.null(post$conformal)) {
  cat("\nConformal summary:\n")
  print(apply(post$conformal, 2, summary))
}

if (!is.null(post$bootstrap)) {
  cat("\nBootstrap summary:\n")
  print(apply(post$bootstrap, 2, summary))
  n_na <- sum(is.na(post$bootstrap[, 1]))
  cat(sprintf("\nBootstrap NA rows: %d / %d\n", n_na, nrow(post$bootstrap)))
}

# --- Density plots ---
cat("\n=== Saving density plots ===\n")
pdf("tests/test_nn_predict_densities.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
for (j in 1:3) {
  pname <- post$param_names[j]
  if (!is.null(post$conformal)) {
    vals <- post$conformal[, j]
    vals <- vals[is.finite(vals)]
    plot(density(vals), main = pname, xlab = pname, ylab = "Density",
         col = "red", lwd = 2)
  }
  if (!is.null(post$bootstrap)) {
    vals <- post$bootstrap[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 2) {
      d2 <- density(vals)
      if (!is.null(post$conformal)) lines(d2, col = "blue", lwd = 2)
      else plot(d2, main = pname, xlab = pname, ylab = "Density", col = "blue", lwd = 2)
    }
  }
  abline(v = post$point_estimate[j], lty = 2, lwd = 2)
  legend("topright", legend = c("Conformal", "Bootstrap", "Point est."),
         col = c("red", "blue", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.8)
}
dev.off()
cat("Saved: tests/test_nn_predict_densities.pdf\n")

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
