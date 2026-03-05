#!/usr/bin/env Rscript
#
# Test emulator.ABC() ABC rejection with three tolerances.
# Uses saved emulator from emulator_vaquita_results/.
#
# Usage: cd tests && Rscript test_emulator_abc_tol.R

devtools::load_all("..", quiet = TRUE)

cat("=== Loading saved emulator ===\n")
emulator <- load.emulator.result("emulator_vaquita_results")

# --- Load data ---
data_dir <- file.path("data", "Vaquita2Epoch", "vaquita_100K")
reftable_full <- read.table(file.path(data_dir, "SIMS_sumstat_vaq.txt"),
                            header = TRUE, sep = "\t")
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
true_vals  <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)

set.seed(42)
reftable <- reftable_full[sample(nrow(reftable_full), 10000L), ]

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
names(obs_raw) <- stat_cols

# Save observed stats for later usage
save(obs_raw, stat_cols, true_vals, param_cols,
     file = "obs_vaquita_100Kloci.RData")
cat("Saved: obs_vaquita_100Kloci.RData\n")
cat(sprintf("Observed: %d summary statistics\n\n", length(obs_raw)))

# === Single pass: sample, predict, compute distances ===
data       <- emulator$data
feat_mu    <- data$feat_mu;  feat_sd  <- data$feat_sd
target_mu  <- data$target_mu; target_sd <- data$target_sd
n_params   <- length(param_cols)
n_stats    <- length(stat_cols)

ptab <- get.prior.table(model)
idx  <- match(param_cols, ptab$Parameter)
lo   <- as.numeric(ptab$prior.1[idx])
hi   <- as.numeric(ptab$prior.2[idx])

n_samples <- 10000000L
cat(sprintf("Sampling %s theta from prior...\n",
            format(n_samples, big.mark = ",")))
theta_orig <- matrix(NA_real_, nrow = n_samples, ncol = n_params)
for (j in seq_len(n_params))
  theta_orig[, j] <- runif(n_samples, lo[j], hi[j])
colnames(theta_orig) <- param_cols

# Transform to emulator input space
theta_z <- t((t(log(theta_orig)) - feat_mu) / feat_sd)

# Ensemble prediction in batches
models     <- emulator$models
val_losses <- emulator$models_val_loss
weights    <- (1 / val_losses) / sum(1 / val_losses)
n_models   <- length(models)
batch_size <- 500000L
n_batches  <- ceiling(n_samples / batch_size)

pred_z <- matrix(0, nrow = n_samples, ncol = n_stats)
cat("Predicting through ensemble...")
t0 <- proc.time()
for (b in seq_len(n_batches)) {
  i_start <- (b - 1L) * batch_size + 1L
  i_end   <- min(b * batch_size, n_samples)
  X_batch <- theta_z[i_start:i_end, , drop = FALSE]
  batch_pred <- matrix(0, nrow = i_end - i_start + 1L, ncol = n_stats)
  for (m in seq_len(n_models)) {
    p <- as.matrix(predict(models[[m]], X_batch, verbose = 0L))
    batch_pred <- batch_pred + weights[m] * p
  }
  pred_z[i_start:i_end, ] <- batch_pred
  if (b %% 5L == 0L || b == n_batches)
    cat(sprintf(" [%d/%d]", b, n_batches))
}
elapsed <- (proc.time() - t0)[3]
cat(sprintf(" done (%.1f sec)\n", elapsed))

# Inverse transform predictions
pred_orig <- .inv.transform.emulator(pred_z, target_mu, target_sd)

# MAD-normalized distances
stat_mad <- apply(pred_orig, 2, function(x) {
  m <- median(x, na.rm = TRUE)
  v <- median(abs(x - m), na.rm = TRUE)
  if (v == 0) sd(x, na.rm = TRUE) else v
})
stat_mad[stat_mad == 0 | !is.finite(stat_mad)] <- 1
resid <- t((t(pred_orig) - obs_raw) / stat_mad)
dists <- sqrt(rowSums(resid^2))

# Sort once
ord <- order(dists)

# === Apply three tolerances ===
tols <- c(0.01, 0.005, 0.001)
posteriors <- list()

for (tol in tols) {
  n_accept <- max(1L, floor(n_samples * tol))
  acc_idx  <- ord[seq_len(n_accept)]
  accepted <- theta_orig[acc_idx, , drop = FALSE]
  acc_dist <- dists[acc_idx]

  # Weighted point estimate
  w <- 1 / (acc_dist + 1e-12)
  w <- w / sum(w)
  pt <- colSums(accepted * w)
  names(pt) <- param_cols

  # Prior samples from reftable
  prior_samples <- as.matrix(reftable[, param_cols, drop = FALSE])

  post <- list(
    point_estimate = pt,
    conformal = NULL, bootstrap = NULL, mc_dropout = NULL,
    quantile = NULL, q_probs = NULL,
    abc_rejection = accepted,
    prior = prior_samples,
    param_names = param_cols
  )
  class(post) <- "nn.posterior"

  label <- sprintf("tol=%.3f", tol)
  posteriors[[label]] <- post

  cat(sprintf("\n--- %s (n=%s) ---\n", label,
              format(n_accept, big.mark = ",")))
  cat(sprintf("  Distance range: [%.4f, %.4f]\n",
              min(acc_dist), max(acc_dist)))
  s <- summary(post)$abc_rejection
  cat(sprintf("  %-12s %8s %10s %10s %10s %10s\n",
              "Param", "True", "Mean", "Median", "2.5%", "97.5%"))
  for (i in seq_along(param_cols))
    cat(sprintf("  %-12s %8.0f %10.1f %10.1f %10.1f %10.1f\n",
                param_cols[i], true_vals[i], s[i, "Mean"], s[i, "Median"],
                s[i, "2.5%"], s[i, "97.5%"]))
  for (i in seq_along(param_cols)) {
    err <- 100 * abs(s[i, "Mean"] - true_vals[i]) / true_vals[i]
    cat(sprintf("  %-12s error: %.1f%%\n", param_cols[i], err))
  }
}

# === Comparison PDF ===
cat("\n=== Generating comparison PDF ===\n")
tol_labels <- names(posteriors)
tol_cols   <- c("steelblue", "firebrick", "darkgreen")

pdf("emulator_abc_tol_comparison_vaquita.pdf", width = 14, height = 10)
par(mfrow = c(2, n_params), mar = c(4, 4, 3, 1))

# Row 1: posteriors at three tolerances (+ prior)
for (j in seq_len(n_params)) {
  pname <- param_cols[j]
  prior_vals <- as.matrix(reftable[, param_cols, drop = FALSE])[, j]
  pr <- range(prior_vals, na.rm = TRUE)

  d_prior <- density(prior_vals, from = pr[1], to = pr[2])
  d_posts <- list()
  for (k in seq_along(tols)) {
    vals <- posteriors[[k]]$abc_rejection[, j]
    d_posts[[k]] <- density(vals, from = pr[1], to = pr[2])
  }
  ymax <- max(d_prior$y, sapply(d_posts, function(d) max(d$y)))

  plot(d_prior, main = pname, xlab = pname, ylab = "Density",
       col = "grey50", lwd = 2, lty = 3, xlim = pr,
       ylim = c(0, ymax * 1.1))
  for (k in seq_along(tols))
    lines(d_posts[[k]], col = tol_cols[k], lwd = 2)
  abline(v = true_vals[j], lty = 2, lwd = 2, col = "black")
  legend("topright",
         legend = c("prior", tol_labels, "true"),
         col = c("grey50", tol_cols, "black"),
         lty = c(3, rep(1, length(tols)), 2),
         lwd = 2, cex = 0.75)
}

# Row 2: zoomed posteriors only
for (j in seq_len(n_params)) {
  pname <- param_cols[j]
  vals_wide <- posteriors[[1]]$abc_rejection[, j]
  xr <- quantile(vals_wide, c(0.005, 0.995))

  d_posts <- list()
  for (k in seq_along(tols)) {
    vals <- posteriors[[k]]$abc_rejection[, j]
    d_posts[[k]] <- density(vals, from = xr[1], to = xr[2])
  }
  ymax <- max(sapply(d_posts, function(d) max(d$y)))

  plot(d_posts[[1]], main = paste(pname, "(zoomed)"), xlab = pname,
       ylab = "Density", col = tol_cols[1], lwd = 2,
       xlim = as.numeric(xr), ylim = c(0, ymax * 1.1))
  for (k in 2:length(tols))
    lines(d_posts[[k]], col = tol_cols[k], lwd = 2)
  abline(v = true_vals[j], lty = 2, lwd = 2, col = "black")
  legend("topright",
         legend = c(tol_labels, "true"),
         col = c(tol_cols, "black"),
         lty = c(rep(1, length(tols)), 2),
         lwd = 2, cex = 0.75)
}
dev.off()
cat("Saved: emulator_abc_tol_comparison_vaquita.pdf\n")

# Save results
save(posteriors, obs_raw, true_vals, param_cols, stat_cols,
     file = "emulator_abc_posteriors_vaquita.RData")
cat("Saved: emulator_abc_posteriors_vaquita.RData\n")
cat("\nDone.\n")
