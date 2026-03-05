#!/usr/bin/env Rscript
#
# Test train.emulator() + emulator.optimize() on the Vaquita dataset.
#
# Uses a 10K subset of the 100K reftable — the forward model (params → stats)
# is smoother and lower-dimensional than the inverse, so fewer training
# samples should suffice.
#
# The Vaquita2Epoch model is a 1-population bottleneck with 3 parameters:
#   Ne0.pop1     = 2807  (current Ne)
#   Ne1.pop1     = 4485  (ancestral Ne)
#   t.Ne1.pop1   = 2162  (bottleneck time in generations)
#
# Usage (run from tests/):
#   cd tests && Rscript test_emulator_vaquita.R
#   cd tests && Rscript test_emulator_vaquita.R 2>&1 | tee test_emulator_vaquita.log
#

devtools::load_all("..", quiet = TRUE)

cat("========================================\n")
cat("  Neural Emulator — Vaquita2Epoch\n")
cat("========================================\n\n")

# =====================================================================
# Load reftable and subsample to 10K
# =====================================================================
data_dir <- file.path("data", "Vaquita2Epoch", "vaquita_100K")

reftable_full <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                            header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

# Subsample 10K rows for the emulator (forward model needs fewer samples)
set.seed(42)
n_sub <- 10000L
idx <- sample(nrow(reftable_full), n_sub)
reftable <- reftable_full[idx, ]

cat(sprintf("Full reftable: %d rows x %d cols\n", nrow(reftable_full), ncol(reftable_full)))
cat(sprintf("Emulator subset: %d rows\n", nrow(reftable)))
cat(sprintf("Params: %s\n\n", paste(param_cols, collapse = ", ")))

# =====================================================================
# 1. Train emulator (params → stats)
# =====================================================================
cat("=== 1. Training neural emulator ===\n")
t0 <- proc.time()
emulator <- train.emulator(
  reftable,
  param.cols  = param_cols,
  max_epochs  = 500,
  n_searches  = 4,
  top_k       = 3,
  cores       = 4,
  gpus        = 0,
  seed        = 42,
  verbose     = TRUE
)
elapsed_train <- (proc.time() - t0)[3]
cat(sprintf("\ntrain.emulator completed in %.1f sec (%.1f min)\n",
            elapsed_train, elapsed_train / 60))
cat(sprintf("Best val_loss: %.6f\n\n", emulator$best_val_loss))

# Save emulator immediately so training isn't lost if optimize crashes
save.emulator.result(emulator, "emulator_vaquita_results")

# =====================================================================
# 2. Evaluate emulator quality on validation set
# =====================================================================
cat("=== 2. Emulator validation metrics ===\n")
if (!is.null(emulator$models_metrics) && length(emulator$models_metrics) > 0L) {
  m <- emulator$models_metrics[[1L]]
  cat(sprintf("%-20s %8s %8s\n", "Stat", "R²", "MPE(%)"))
  stat_names <- names(m$r2)
  for (i in seq_along(stat_names))
    cat(sprintf("%-20s %8.4f %8.1f\n", stat_names[i], m$r2[i], m$mpe[i]))
  cat(sprintf("\n%-20s %8.4f %8.1f\n", "Mean",
              mean(m$r2, na.rm = TRUE), mean(m$mpe, na.rm = TRUE)))
}

# =====================================================================
# 3. Prepare observed summary statistics
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

cat(sprintf("\nObserved: %d summary statistics\n", length(obs_raw)))

# =====================================================================
# 4. Gradient-based optimization
# =====================================================================
cat("\n=== 3. Gradient-based optimization ===\n")

# Extract prior bounds from model
prior_tab <- get.prior.table(model)
cat("Prior bounds:\n")
print(prior_tab[prior_tab$Parameter %in% param_cols, ])

t0 <- proc.time()
opt <- emulator.optimize(
  emulator,
  observed      = obs_raw,
  model         = model,
  n_starts      = 100,
  n_steps       = 2000,
  learning_rate = 0.01,
  verbose       = TRUE
)
elapsed_opt <- (proc.time() - t0)[3]
cat(sprintf("\nemulator.optimize completed in %.1f sec\n", elapsed_opt))

# =====================================================================
# 5. ABC rejection posterior via emulator
# =====================================================================
cat("\n=== 4. ABC rejection posterior ===\n")

t0 <- proc.time()
posterior <- emulator.ABC(
  emulator,
  observed  = obs_raw,
  model     = model,
  reftable  = reftable,
  n_samples = 1e6,
  tol       = 0.01,
  distance  = "euclidean",
  verbose   = TRUE
)
elapsed_abc <- (proc.time() - t0)[3]
cat(sprintf("\nemulator.ABC completed in %.1f sec\n", elapsed_abc))

cat("\n--- Posterior summary ---\n")
print(summary(posterior))

# Save posterior density plot
pdf("emulator_abc_posterior_vaquita.pdf", width = 12, height = 5)
plot(posterior, show_prior = TRUE, show_point_est = TRUE)
dev.off()
cat("Saved: emulator_abc_posterior_vaquita.pdf\n")

# Compare posterior mean/median to true values and gradient-based estimate
abc_summ <- summary(posterior)$abc_rejection
cat(sprintf("\n%-12s %8s %12s %12s %12s\n",
            "Param", "True", "ABC Mean", "ABC Median", "Grad. Opt."))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %12.1f %12.1f %12.1f\n",
              param_cols[i], true_vals[i],
              abc_summ[i, "Mean"], abc_summ[i, "Median"],
              opt$point_estimate[i]))
}

# =====================================================================
# 6. Compare to true values and nn.predict
# =====================================================================
cat("\n=== 5. Comparison ===\n")

# Load nn.predict point estimate for comparison
# Use load.tune.result() (not load()) — keras models need load_model_tf()
tune_dirs <- c("topk_ensemble_results", "tune_vaquita_results")
tune_loaded <- FALSE
for (td in tune_dirs) {
  if (dir.exists(td) && file.exists(file.path(td, "tune_result.rds"))) {
    cat(sprintf("Loading tune result from %s/ ...\n", td))
    tune_result <- load.tune.result(td)
    pred_nn <- nn.predict(tune_result, obs_raw, method = "point")
    nn_est  <- pred_nn$point_estimate
    tune_loaded <- TRUE
    break
  }
}
if (!tune_loaded) {
  cat("(No saved tune_result found, running quick tune.nn for comparison...)\n")
  tune_result <- tune.nn(reftable_full, param.cols = param_cols, type = "sumstat",
                         max_epochs = 500, n_searches = 4, cores = 4, seed = 42)
  pred_nn <- nn.predict(tune_result, obs_raw, method = "point")
  nn_est  <- pred_nn$point_estimate
}

emu_est <- opt$point_estimate

cat(sprintf("\n%-12s %8s %12s %12s\n", "Param", "True", "Emulator", "nn.predict"))
for (i in seq_along(param_cols)) {
  cat(sprintf("%-12s %8.0f %12.1f %12.1f\n",
              param_cols[i], true_vals[i], emu_est[i], nn_est[i]))
}

cat(sprintf("\n--- Percent errors ---\n"))
cat(sprintf("%-12s %12s %12s\n", "Param", "Emulator", "nn.predict"))
for (i in seq_along(param_cols)) {
  err_emu <- 100 * abs(emu_est[i] - true_vals[i]) / true_vals[i]
  err_nn  <- 100 * abs(nn_est[i]  - true_vals[i]) / true_vals[i]
  cat(sprintf("%-12s %11.1f%% %11.1f%%\n", param_cols[i], err_emu, err_nn))
}

mean_emu <- mean(abs(emu_est - true_vals) / true_vals) * 100
mean_nn  <- mean(abs(nn_est  - true_vals) / true_vals) * 100
cat(sprintf("%-12s %11.1f%% %11.1f%%\n", "Mean", mean_emu, mean_nn))

# =====================================================================
# 7. Loss distribution summary
# =====================================================================
cat(sprintf("\n--- Optimization loss distribution ---\n"))
cat(sprintf("  Min:    %.6f\n", min(opt$losses)))
cat(sprintf("  Median: %.6f\n", median(opt$losses)))
cat(sprintf("  Mean:   %.6f\n", mean(opt$losses)))
cat(sprintf("  Max:    %.6f\n", max(opt$losses)))

# Top-10 optima
n_top <- max(1L, floor(length(opt$losses) * 0.1))
top_idx <- order(opt$losses)[1:n_top]
cat(sprintf("\n--- Top-%d optima ---\n", n_top))
cat(sprintf("%-5s %12s %12s %12s %12s\n", "Rank",
            param_cols[1], param_cols[2], param_cols[3], "Loss"))
for (k in seq_len(min(n_top, 10))) {
  j <- top_idx[k]
  cat(sprintf("%-5d %12.1f %12.1f %12.1f %12.6f\n",
              k, opt$optima[j, 1], opt$optima[j, 2], opt$optima[j, 3],
              opt$losses[j]))
}

# =====================================================================
# 8. Save results
# =====================================================================
save(opt, posterior, obs_raw, true_vals,
     file = "test_emulator_vaquita_results.RData")
cat("Saved: test_emulator_vaquita_results.RData\n")

cat("\n========================================\n")
cat("  DONE\n")
cat("========================================\n")
