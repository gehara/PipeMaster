#!/usr/bin/env Rscript
#
# Step 1: Run tune.nn() on GPU machine (lagarto), save result.
# Result is saved to tests/hABC/data/Ubelmania/tune_result/
# which syncs via Dropbox to segovia for nn.predict().
#
# Usage: Rscript tests/test_nn_tune.R

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

cat("========================================\n")
cat("  tune.nn() — Ubelmania pop 1\n")
cat("========================================\n\n")

# --- Load data ---
reftable <- read.table("tests/hABC/data/Ubelmania/abc_sims/SIMS_pop_1.txt",
                       header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")

cat(sprintf("Reftable: %d rows x %d cols\n", nrow(reftable), ncol(reftable)))
cat(sprintf("Params: %s\n\n", paste(param_cols, collapse = ", ")))

# --- tune.nn ---
t0 <- proc.time()
tune_result <- tune.nn(reftable,
                       param.cols = param_cols,
                       type = "sumstat",
                       max_epochs = 500,
                       seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\ntune.nn completed in %.1f sec (%.1f min)\n", elapsed, elapsed / 60))
cat(sprintf("Best val_loss: %.6f\n\n", tune_result$best_val_loss))

# --- Save result ---
out_dir <- "tests/hABC/data/Ubelmania/tune_result"
save.tune.result(tune_result, out_dir)

cat("\nDone. Now run nn.predict on segovia:\n")
cat("  Rscript tests/test_nn_predict.R\n")
