#!/usr/bin/env Rscript
#
# Test tune.nn() on the Vaquita 100K sumstat dataset (3 demographic params).
#
# Runs three models on the same train/val split and compares point estimates
# against known true values:
#   1. Native R Hyperband (tune.nn) — parallel multi-search
#   2. Baseline architecture — fixed conservative HPs (256-128-64, lr=1e-3)
#   3. Keras Tuner Hyperband — best HPs found by Python keras-tuner on segovia
#
# The Vaquita2Epoch model is a 1-population bottleneck with 3 parameters:
#   Ne0.pop1     = 2807  (current Ne)
#   Ne1.pop1     = 4485  (ancestral Ne)
#   t.Ne1.pop1   = 2162  (bottleneck time in generations)
#
# Usage:
#   Rscript tests/test_nn_tune_vaquita.R            # requires installed PipeMaster
#   Rscript tests/test_nn_tune_vaquita.R 2>&1 | tee tests/test_nn_tune_vaquita.log
#
# Requires: keras, tensorflow, PipeMaster (installed)

library(PipeMaster)

cat("========================================\n")
cat("  tune.nn() — Vaquita2Epoch (sumstat)\n")
cat("========================================\n\n")

# =====================================================================
# Load reftable: 100K coalescent simulations, 24 summary statistics
# Columns: 3 demographic params + 2 nuisance (mean.rate, sd.rate) + 24 sumstats
# =====================================================================
data_dir <- "tests/data/Vaquita2Epoch/vaquita_100K"

reftable <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                       header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

cat(sprintf("Reftable: %d rows x %d cols\n", nrow(reftable), ncol(reftable)))
cat(sprintf("Params: %s\n\n", paste(param_cols, collapse = ", ")))

# =====================================================================
# 1. Native R Hyperband (tune.nn)
#
# Launches n_searches independent Hyperband searches as parallel Rscript
# workers. Each search runs 4 brackets (s_max=3 with eta=3):
#   Bracket 3: 27 configs x 19 epochs  -> prune to 9 -> 3 -> 1 @ 500 ep
#   Bracket 2: 18 configs x 56 epochs  -> prune to 6 -> 2 @ 500 ep
#   Bracket 1: 12 configs x 167 epochs -> prune to 4 @ 500 ep
#   Bracket 0:  8 configs x 500 epochs (full training, no pruning)
#
# After all searches, the best config is retrained to max_epochs.
# Each worker gets floor(physical_cores / n_workers) TF intra-op threads
# (greedy=TRUE by default). Set greedy=FALSE for 1 thread per worker.
#
# GPU: all workers assigned to GPU 0 (gpu.threshold controls how many
# workers share a GPU before overflow to CPU-only).
# =====================================================================
cat("=== 1. Native R Hyperband (tune.nn) ===\n")
t0 <- proc.time()
tune_result <- tune.nn(reftable,
                       param.cols = param_cols,
                       type = "sumstat",
                       max_epochs = 500,
                       n_searches = 4L, cores = 4L, gpus = 1L,
                       seed = 42)
elapsed_tune <- (proc.time() - t0)[3]
cat(sprintf("\ntune.nn completed in %.1f sec (%.1f min)\n", elapsed_tune, elapsed_tune / 60))
cat(sprintf("Best val_loss: %.6f\n\n", tune_result$best_val_loss))

# =====================================================================
# Prepare observed data for prediction
#
# Load the Vaquita phylip alignment, compute observed summary statistics
# via obs.sumstat.ngs(), then extract the same stat columns used by the
# reftable (excluding params and nuisance columns).
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

nuisance <- c("mean.rate", "sd.rate")
stat_cols <- setdiff(colnames(reftable), c(param_cols, nuisance))
obs_raw <- as.numeric(obs_ss[1, stat_cols])

# =====================================================================
# Point prediction from tune_result using nn.predict(method = "point")
#
# Fast point estimate — no reftable or param.cols needed.
# Internally handles feature augmentation, Z-scoring, and inverse
# transform via the same pipeline used by the other methods.
# =====================================================================
pred_point <- nn.predict(tune_result, obs_raw, method = "point")
pred_tuned <- pred_point$point_estimate

# Prepare X_obs for baseline and KT models (same augmentation pipeline)
obs_aug <- c(obs_raw, log1p(abs(obs_raw)))
X_obs <- matrix((obs_aug - tune_result$data$feat_mu) / tune_result$data$feat_sd, nrow = 1)
X_obs[!is.finite(X_obs)] <- 0

# =====================================================================
# 2. Baseline architecture
#
# Conservative fixed HPs: 3-stage ResNet (256-128-64), no dropout,
# standard lr=1e-3, batch_size=512, huber_delta=1.0.
# Trained on the same data split as tune_result for fair comparison.
# Early stopping (patience=30) + ReduceLROnPlateau (patience=15).
# =====================================================================
cat("=== 2. Training baseline architecture ===\n")
baseline_hp <- list(
  units_1       = 256,
  n_resblocks_1 = 2,
  units_2       = 128,
  n_resblocks_2 = 1,
  units_3       = 64,
  learning_rate = 0.001,
  l2_reg        = 1e-4,
  dropout       = 0,
  use_dropout   = FALSE,
  batch_size    = 512,
  huber_delta   = 1.0
)

tensorflow::tf$random$set_seed(42L)
baseline_model <- PipeMaster:::.build.resnet(baseline_hp, tune_result$data)

t0 <- proc.time()
baseline_history <- baseline_model |> keras::fit(
  x = tune_result$data$X_train, y = tune_result$data$Y_train,
  validation_data = list(tune_result$data$X_val, tune_result$data$Y_val),
  epochs     = 500L,
  batch_size = 512L,
  callbacks  = list(
    keras::callback_early_stopping(monitor = "val_loss", patience = 30L,
                                   restore_best_weights = TRUE),
    keras::callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,
                                         factor = 0.5, min_lr = 1e-6, verbose = 0L)
  ),
  verbose = 0L
)
elapsed_base <- (proc.time() - t0)[3]

base_vl <- unlist(baseline_history$metrics$val_loss)
if (is.null(base_vl)) base_vl <- unlist(baseline_history$history$val_loss)
baseline_val_loss <- min(base_vl[is.finite(base_vl)])

pred_z_base <- predict(baseline_model, X_obs, verbose = 0L)
pred_baseline <- as.numeric(exp(t(t(pred_z_base) * tune_result$data$target_sd + tune_result$data$target_mu)))

cat(sprintf("Baseline val_loss: %.6f (%.1f sec)\n\n", baseline_val_loss, elapsed_base))

# =====================================================================
# 3. Keras Tuner Hyperband architecture
#
# Best HPs found by Python keras-tuner Hyperband on segovia (64 cores).
# Wider first layer (512), shallower (1 resblock), with dropout=0.05.
# Trained here from scratch on the same data split for fair comparison.
# =====================================================================
cat("=== 3. Training Keras Tuner best architecture ===\n")
kt_hp <- list(
  units_1       = 512,
  n_resblocks_1 = 1,
  units_2       = 64,
  n_resblocks_2 = 0,
  units_3       = 128,
  learning_rate = 0.000401,
  l2_reg        = 1.21e-5,
  dropout       = 0.05,
  use_dropout   = TRUE,
  batch_size    = 256,
  huber_delta   = 0.5
)

tensorflow::tf$random$set_seed(42L)
kt_model <- PipeMaster:::.build.resnet(kt_hp, tune_result$data)

t0 <- proc.time()
kt_history <- kt_model |> keras::fit(
  x = tune_result$data$X_train, y = tune_result$data$Y_train,
  validation_data = list(tune_result$data$X_val, tune_result$data$Y_val),
  epochs     = 500L,
  batch_size = 256L,
  callbacks  = list(
    keras::callback_early_stopping(monitor = "val_loss", patience = 30L,
                                   restore_best_weights = TRUE),
    keras::callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,
                                         factor = 0.5, min_lr = 1e-6, verbose = 0L)
  ),
  verbose = 0L
)
elapsed_kt <- (proc.time() - t0)[3]

kt_vl <- unlist(kt_history$metrics$val_loss)
if (is.null(kt_vl)) kt_vl <- unlist(kt_history$history$val_loss)
kt_val_loss <- min(kt_vl[is.finite(kt_vl)])

pred_z_kt <- predict(kt_model, X_obs, verbose = 0L)
pred_kt <- as.numeric(exp(t(t(pred_z_kt) * tune_result$data$target_sd + tune_result$data$target_mu)))

cat(sprintf("Keras Tuner val_loss: %.6f (%.1f sec)\n\n", kt_val_loss, elapsed_kt))

# =====================================================================
# Comparison: all three models side by side
#
# Compare hyperparameters, validation loss, point estimates, and errors.
# The "winner" for each parameter is the model closest to the true value.
# =====================================================================
cat("=====================================================================\n")
cat("  COMPARISON: Baseline vs Keras Tuner vs Native Hyperband\n")
cat("=====================================================================\n\n")

cat("--- Hyperparameters ---\n")
hp_names <- c("units_1", "n_resblocks_1", "units_2", "n_resblocks_2", "units_3",
              "learning_rate", "batch_size", "l2_reg", "huber_delta",
              "use_dropout", "dropout")
cat(sprintf("%-20s %12s %12s %12s\n", "", "Baseline", "KerasTuner", "NativeHB"))
for (nm in hp_names) {
  cat(sprintf("%-20s %12s %12s %12s\n", nm,
              as.character(baseline_hp[[nm]]),
              as.character(kt_hp[[nm]]),
              as.character(tune_result$best_hp[[nm]])))
}

cat(sprintf("\n--- Validation loss ---\n"))
cat(sprintf("  Baseline:      %.6f\n", baseline_val_loss))
cat(sprintf("  Keras Tuner:   %.6f\n", kt_val_loss))
cat(sprintf("  Native HB:     %.6f\n", tune_result$best_val_loss))

cat(sprintf("\n--- Point estimates vs true values ---\n"))
cat(sprintf("%-12s %8s %10s %10s %10s\n", "Param", "True", "Baseline", "KerasTuner", "NativeHB"))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %10.0f %10.0f %10.0f\n",
              param_cols[i], true_vals[i],
              pred_baseline[i], pred_kt[i], pred_tuned[i]))
}

cat(sprintf("\n--- Absolute errors ---\n"))
cat(sprintf("%-12s %10s %10s %10s\n", "Param", "Baseline", "KerasTuner", "NativeHB"))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %10.0f %10.0f %10.0f\n",
              param_cols[i],
              abs(pred_baseline[i] - true_vals[i]),
              abs(pred_kt[i] - true_vals[i]),
              abs(pred_tuned[i] - true_vals[i])))
}

cat(sprintf("\n--- Percent errors ---\n"))
cat(sprintf("%-12s %10s %10s %10s\n", "Param", "Baseline", "KerasTuner", "NativeHB"))
for (i in seq_along(param_cols)) {
  err_base <- 100 * abs(pred_baseline[i] - true_vals[i]) / true_vals[i]
  err_kt   <- 100 * abs(pred_kt[i] - true_vals[i]) / true_vals[i]
  err_hb   <- 100 * abs(pred_tuned[i] - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %9.1f%% %9.1f%% %9.1f%%\n",
              param_cols[i], err_base, err_kt, err_hb))
}

# Mean percent error across all 3 parameters
mean_base <- mean(sapply(seq_along(param_cols), function(i)
  abs(pred_baseline[i] - true_vals[i]) / true_vals[i]))
mean_kt <- mean(sapply(seq_along(param_cols), function(i)
  abs(pred_kt[i] - true_vals[i]) / true_vals[i]))
mean_hb <- mean(sapply(seq_along(param_cols), function(i)
  abs(pred_tuned[i] - true_vals[i]) / true_vals[i]))

cat(sprintf("\n%-12s %9.1f%% %9.1f%% %9.1f%%\n", "Mean", 100*mean_base, 100*mean_kt, 100*mean_hb))

# Best model per parameter (lowest absolute error)
cat(sprintf("\n--- Winner per parameter ---\n"))
for (i in seq_along(param_cols)) {
  errs <- c(Baseline = abs(pred_baseline[i] - true_vals[i]),
            KerasTuner = abs(pred_kt[i] - true_vals[i]),
            NativeHB = abs(pred_tuned[i] - true_vals[i]))
  cat(sprintf("%-12s -> %s\n", param_cols[i], names(which.min(errs))))
}

# =====================================================================
# Save results for later analysis
# =====================================================================
save(tune_result, baseline_hp, kt_hp,
     pred_baseline, pred_kt, pred_tuned,
     baseline_val_loss, kt_val_loss, true_vals,
     file = "tests/test_nn_tune_vaquita_results.RData")
cat("\nSaved: tests/test_nn_tune_vaquita_results.RData\n")

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
