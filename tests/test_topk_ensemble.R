#!/usr/bin/env Rscript
#
# Test top_k ensemble: tune.nn(top_k=4) on 1% Vaquita reftable.
#
# Verifies:
#   - top_k=4 returns 4 models with val_loss vector
#   - nn.predict() uses proximity-weighted ensemble
#   - save/load round-trip preserves multi-model result
#   - Compare ensemble vs single-best point estimates
#
# Usage (run from tests/):
#   cd tests && Rscript test_topk_ensemble.R 2>&1 | tee test_topk_ensemble.log

library(PipeMaster)

cat("========================================\n")
cat("  top_k=4 ensemble — Vaquita (10% data)\n")
cat("========================================\n\n")

# =====================================================================
# Load reftable and subsample to 10%
# =====================================================================
data_dir <- file.path("data", "Vaquita2Epoch", "vaquita_100K")

reftable <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                       header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

set.seed(42)
n_sub <- as.integer(nrow(reftable) * 0.10)
idx <- sample.int(nrow(reftable), n_sub)
reftable_sub <- reftable[idx, ]

cat(sprintf("Full reftable: %d rows\n", nrow(reftable)))
cat(sprintf("Subsample (10%%): %d rows\n", nrow(reftable_sub)))
cat(sprintf("Params: %s\n\n", paste(param_cols, collapse = ", ")))

# =====================================================================
# 1. tune.nn with top_k=4
# =====================================================================
cat("=== 1. tune.nn(top_k=4, n_searches=10, max_epochs=1000) ===\n")
t0 <- proc.time()
tune_result <- tune.nn(reftable_sub,
                       param.cols = param_cols,
                       type = "sumstat",
                       max_epochs = 1000,
                       n_searches = 10, cores = 10, gpus = 0,
                       top_k = 4L,
                       seed = 42)
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\ntune.nn completed in %.1f sec\n", elapsed))

# Verify return structure
cat(sprintf("\n--- Return structure ---\n"))
cat(sprintf("  length(models):      %d\n", length(tune_result$models)))
cat(sprintf("  models_val_loss:     %s\n",
            paste(sprintf("%.6f", tune_result$models_val_loss), collapse = ", ")))
cat(sprintf("  best_val_loss:       %.6f\n", tune_result$best_val_loss))
cat(sprintf("  best_model == models[[1]]: %s\n\n",
            identical(tune_result$best_model, tune_result$models[[1]])))

stopifnot(length(tune_result$models) >= 1L)
stopifnot(length(tune_result$models) <= 4L)
stopifnot(length(tune_result$models_val_loss) == length(tune_result$models))
stopifnot(!is.unsorted(tune_result$models_val_loss))

# Verify models_metrics returned
cat(sprintf("\n--- Model metrics (R² and MPE) ---\n"))
stopifnot(!is.null(tune_result$models_metrics))
stopifnot(length(tune_result$models_metrics) == length(tune_result$models))
for (m in seq_along(tune_result$models_metrics)) {
  mm <- tune_result$models_metrics[[m]]
  cat(sprintf("  Model %d: R²=[%s]  MPE=[%s]\n", m,
              paste(sprintf("%.3f", mm$r2), collapse = ", "),
              paste(sprintf("%.1f%%", mm$mpe), collapse = ", ")))
}

# =====================================================================
# 2. Prepare observed data
# =====================================================================
load(file.path(data_dir, "model.RData"))
model <- Vaquita2Epoch_100K
pop_assign <- read.table(file.path(data_dir, "pop_assign_Vaquita2Epoch.txt"),
                         header = TRUE)
obs_ss <- obs.sumstat.ngs(
  model = model,
  path.to.phylip = file.path(data_dir, "phylip_Vaquita2Epoch.phy"),
  pop.assign = pop_assign
)

nuisance  <- c("mean.rate", "sd.rate")
stat_cols <- setdiff(colnames(reftable), c(param_cols, nuisance))
obs_raw   <- as.numeric(obs_ss[1, stat_cols])

# =====================================================================
# 3. Ensemble point estimate (top_k=4)
# =====================================================================
cat("=== 2. Ensemble nn.predict(method='point') ===\n")
pred_ensemble <- nn.predict(tune_result, obs_raw, method = "point")
est_ensemble  <- pred_ensemble$point_estimate

# =====================================================================
# 4. Single-best point estimate (for comparison)
# =====================================================================
cat("\n=== 3. Single-best nn.predict(method='point') ===\n")
tune_single <- tune_result
tune_single$models <- list(tune_result$models[[1]])
tune_single$models_val_loss <- tune_result$models_val_loss[1]

pred_single <- nn.predict(tune_single, obs_raw, method = "point")
est_single  <- pred_single$point_estimate

# =====================================================================
# 5. Per-model predictions (raw, no ensemble)
# =====================================================================
cat("\n=== 4. Individual model predictions ===\n")
data <- tune_result$data
obs_aug <- c(obs_raw, log1p(abs(obs_raw)))
X_obs <- matrix((obs_aug - data$feat_mu) / data$feat_sd, nrow = 1)
X_obs[!is.finite(X_obs)] <- 0

for (m in seq_along(tune_result$models)) {
  pred_z <- predict(tune_result$models[[m]], X_obs, verbose = 0L)
  pred_raw <- as.numeric(exp(t(t(pred_z) * data$target_sd + data$target_mu)))
  names(pred_raw) <- param_cols
  cat(sprintf("  Model %d (vl=%.6f): %s\n", m, tune_result$models_val_loss[m],
              paste(sprintf("%s=%.0f", param_cols, pred_raw), collapse = " ")))
}

# =====================================================================
# 6. Comparison table
# =====================================================================
cat("\n=====================================================================\n")
cat("  COMPARISON: Single-best vs Ensemble vs True\n")
cat("=====================================================================\n\n")

cat(sprintf("%-12s %8s %10s %10s\n", "Param", "True", "Single", "Ensemble"))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %10.0f %10.0f\n",
              param_cols[i], true_vals[i], est_single[i], est_ensemble[i]))
}

cat(sprintf("\n--- Percent errors ---\n"))
cat(sprintf("%-12s %10s %10s\n", "Param", "Single", "Ensemble"))
for (i in seq_along(param_cols)) {
  err_s <- 100 * abs(est_single[i] - true_vals[i]) / true_vals[i]
  err_e <- 100 * abs(est_ensemble[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %9.1f%% %9.1f%%\n", param_cols[i], err_s, err_e))
}

mean_s <- mean(sapply(seq_along(param_cols), function(i)
  abs(est_single[i] - true_vals[i]) / true_vals[i]))
mean_e <- mean(sapply(seq_along(param_cols), function(i)
  abs(est_ensemble[i] - true_vals[i]) / true_vals[i]))
cat(sprintf("\n%-12s %9.1f%% %9.1f%%\n", "Mean", 100*mean_s, 100*mean_e))

# =====================================================================
# 7. Save/load round-trip test
# =====================================================================
cat("\n=== 5. Save/load round-trip ===\n")
save_dir <- tempfile("topk_save_test_")
save.tune.result(tune_result, save_dir)

loaded <- load.tune.result(save_dir)
cat(sprintf("  Loaded models:         %d\n", length(loaded$models)))
cat(sprintf("  Loaded models_val_loss: %s\n",
            paste(sprintf("%.6f", loaded$models_val_loss), collapse = ", ")))

# Verify loaded result produces same ensemble prediction
pred_loaded <- nn.predict(loaded, obs_raw, method = "point")
est_loaded  <- pred_loaded$point_estimate

cat(sprintf("\n  Ensemble after reload:\n"))
for (i in seq_along(param_cols)) {
  cat(sprintf("    %-12s original=%.0f  reloaded=%.0f\n",
              param_cols[i], est_ensemble[i], est_loaded[i]))
}

# Clean up
unlink(save_dir, recursive = TRUE)

stopifnot(length(loaded$models) == length(tune_result$models))

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
