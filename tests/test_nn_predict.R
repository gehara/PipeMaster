#!/usr/bin/env Rscript
#
# End-to-end test: tune.nn() + nn.predict() on Ubelmania pop 1 data
#

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

cat("========================================\n")
cat("  nn.predict() end-to-end test\n")
cat("========================================\n\n")

# --- Load data ---
reftable <- read.table("tests/hABC/data/Ubelmania/abc_sims/SIMS_pop_1.txt",
                       header = TRUE, sep = "\t")
observed <- read.table("tests/hABC/data/Ubelmania/grid_search/obs_pop_1.txt",
                       header = TRUE, sep = "\t")

param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")

cat(sprintf("Reftable: %d rows x %d cols\n", nrow(reftable), ncol(reftable)))
cat(sprintf("Observed: %d cols\n", ncol(observed)))
cat(sprintf("Params: %s\n\n", paste(param_cols, collapse = ", ")))

# --- Step 1: tune.nn (short run for testing) ---
cat("=== Step 1: tune.nn ===\n")
t0 <- proc.time()
tune_result <- tune.nn(reftable,
                       param.cols = param_cols,
                       type = "sumstat",
                       max_epochs = 50,
                       seed = 42)
t1 <- proc.time()
cat(sprintf("\ntune.nn completed in %.1f sec\n", (t1 - t0)[3]))
cat(sprintf("Best val_loss: %.6f\n", tune_result$best_val_loss))
cat(sprintf("Data stat_cols stored: %s\n", !is.null(tune_result$data$stat_cols)))
cat(sprintf("Type stored: %s\n\n", tune_result$type))

# --- Step 2: nn.predict ---
cat("=== Step 2: nn.predict (conformal + bootstrap) ===\n")
t0 <- proc.time()
post <- nn.predict(tune_result, observed,
                   reftable = reftable,
                   param.cols = param_cols,
                   method = c("conformal", "bootstrap"),
                   n_boot = 50,
                   n_ensemble = 10,
                   max_epochs = 500,
                   cores = 32,
                   verbose = TRUE)
t1 <- proc.time()
cat(sprintf("\nnn.predict completed in %.1f sec\n\n", (t1 - t0)[3]))

# --- Step 3: Verify results ---
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

# --- Step 4: Quick density plots ---
cat("\n=== Saving density plots ===\n")
pdf("tests/test_nn_predict_densities.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
for (j in 1:3) {
  pname <- post$param_names[j]

  # Conformal density
  if (!is.null(post$conformal)) {
    vals <- post$conformal[, j]
    vals <- vals[is.finite(vals)]
    d <- density(vals)
    plot(d, main = pname, xlab = pname, ylab = "Density",
         col = "red", lwd = 2)
  }

  # Bootstrap density
  if (!is.null(post$bootstrap)) {
    vals <- post$bootstrap[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 2) {
      d2 <- density(vals)
      if (!is.null(post$conformal)) {
        lines(d2, col = "blue", lwd = 2)
      } else {
        plot(d2, main = pname, xlab = pname, ylab = "Density",
             col = "blue", lwd = 2)
      }
    }
  }

  # Point estimate
  abline(v = post$point_estimate[j], lty = 2, lwd = 2)
  legend("topright", legend = c("Conformal", "Bootstrap", "Point est."),
         col = c("red", "blue", "black"), lty = c(1, 1, 2), lwd = 2, cex = 0.8)
}
dev.off()
cat("Saved: tests/test_nn_predict_densities.pdf\n")

cat("\n========================================\n")
cat("  TEST COMPLETE\n")
cat("========================================\n")
