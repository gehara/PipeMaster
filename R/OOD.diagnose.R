# ============================================================================
# Out-of-Distribution Diagnostic for Neural Network Inference
#
# OOD.diagnose() — OOD checks for train.emulator() and tune.nn() results,
#                  or standalone reftable + observed comparison
# ============================================================================

#' Out-of-Distribution Diagnostic for Neural Network Inference
#'
#' Evaluates whether observed summary statistics fall within the distribution
#' of the training simulations, and whether the model's predictions are
#' reliable. Combines four complementary checks. Works with both
#' \code{train.emulator()} (forward model θ→S) and \code{tune.nn()} (inverse
#' model S→θ) results of any architecture type (sumstat, sfs1d, sfs2d).
#'
#' Can also be called with \code{trained.nn = NULL} to run checks 1-3
#' (Mahalanobis, PCA, percentiles) using only the reftable and observed
#' statistics, without a trained model. In this mode, stat columns are
#' inferred from the observed data.frame column names or matched against
#' the reftable.
#'
#' @param trained.nn list or NULL — output from \code{train.emulator()} or
#'   \code{tune.nn()}. If NULL, only checks 1-3 are run (no model disagreement).
#' @param observed numeric vector or 1-row data.frame — observed summary statistics.
#'   If a data.frame with named columns, those names are used to identify stat
#'   columns in the reftable.
#' @param reftable data.frame — the reference table used for training (simulated stats).
#' @param stat.cols character vector or NULL — column names of summary statistics
#'   in the reftable. Required when \code{trained.nn = NULL} and \code{observed}
#'   is a numeric vector. Ignored when \code{trained.nn} is provided (stat columns
#'   are taken from the trained model).
#' @param theta numeric vector or NULL — optimized parameter estimate (from
#'   \code{emulator.optimize()}). If provided, model disagreement is evaluated
#'   at this point. Only used for emulator-type results; for \code{tune.nn()}
#'   results, model disagreement is evaluated at the observed stats directly.
#' @param pca.var numeric — cumulative variance threshold for PCA (default 0.95).
#' @param alpha numeric — significance level for flagging (default 0.01).
#' @param plot logical — produce diagnostic plots (default TRUE).
#'
#' @return A list of class \code{"OOD_diagnostic"} with:
#' \describe{
#'   \item{mahalanobis}{list with d2, p_value, pass (logical)}
#'   \item{pca}{list with d2, p_value, n_pcs, pass}
#'   \item{percentiles}{data.frame with stat, percentile, outlier (logical)}
#'   \item{model_disagreement}{list with mean_cv, per_stat_cv, reliable (logical).
#'     NULL if \code{trained.nn} is NULL or fewer than 2 models.}
#'   \item{overall}{character: \code{"pass"}, \code{"warn"}, or \code{"fail"}}
#' }
#'
#' @section Interpreting the diagnostic plots:
#'
#' The function produces a 2x2 panel of diagnostic plots:
#'
#' \strong{Panel 1 — PCA: observed vs simulated.}
#' Grey dots are simulations projected onto PC1 vs PC2; the red triangle is the
#' observed data. The observed point should sit inside the cloud. If it falls at
#' the edge or outside, the simulations do not cover the observed data well,
#' suggesting the prior range may need widening.
#'
#' \strong{Panel 2 — Per-stat percentiles.}
#' Each bar shows where the observed value ranks among simulations for that
#' statistic (0.5 = median). Red dashed lines mark the \code{alpha} and
#' \code{1 - alpha} thresholds. Bars near 0 or 1 (flagged in red) indicate
#' statistics where the observed value is extreme relative to simulations.
#' A few marginal outliers may be acceptable, but many suggest the model
#' cannot reproduce the observed data.
#'
#' \strong{Panel 3 — Mahalanobis distance.}
#' Grey histogram: Mahalanobis distances of all simulations from the simulation
#' centroid. Blue curve: theoretical chi-squared reference distribution.
#' Red vertical line: observed distance. The observed line should fall within
#' the bulk of the histogram, not in the right tail. Note that the histogram
#' may deviate from the chi-squared curve if the simulation distribution is
#' non-Gaussian.
#'
#' \strong{Panel 4 — Model disagreement (CV).}
#' Only shown when \code{trained.nn} is provided with multiple top-K models.
#' Each bar is the coefficient of variation across models when predicting at
#' the optimized theta (emulator) or observed stats (tune.nn). The red dashed
#' line at 0.10 (10\%) is the reliability threshold. Low CV indicates strong
#' model agreement; bars above the line flag statistics or parameters where the
#' neural network is uncertain.
#'
#' @section Overall assessment:
#'
#' The overall result is determined by combining all checks:
#' \describe{
#'   \item{PASS}{All checks pass.}
#'   \item{WARN}{Exactly one check fails.}
#'   \item{FAIL}{Two or more checks fail.}
#' }
#'
#' @export
OOD.diagnose <- function(trained.nn = NULL, observed, reftable,
                        stat.cols = NULL,
                        theta = NULL,
                        pca.var = 0.95,
                        alpha = 0.01,
                        plot = TRUE) {

  # --- Detect model type and determine stat columns ---
  has_model <- !is.null(trained.nn)

  if (has_model) {
    model_type <- trained.nn$type  # "emulator", "sumstat", "sfs1d", "sfs2d"
    if (is.null(model_type)) model_type <- "emulator"
    is_emulator <- model_type == "emulator"
    stat_cols <- trained.nn$data$stat_cols
  } else {
    model_type <- "none"
    is_emulator <- FALSE

    # Infer stat columns
    if (!is.null(stat.cols)) {
      stat_cols <- stat.cols
    } else if (is.data.frame(observed)) {
      stat_cols <- colnames(observed)
    } else {
      stop("When trained.nn is NULL, either 'stat.cols' must be provided or 'observed' must be a named data.frame.")
    }
  }

  # --- Prepare observed ---
  if (is.data.frame(observed) || is.matrix(observed)) {
    if (all(stat_cols %in% colnames(observed)))
      observed <- as.numeric(observed[1, stat_cols])
    else
      observed <- as.numeric(observed[1, ])
  }
  obs_raw <- observed  # keep a copy for check 4

  # Extract simulated stats from reftable
  sim_cols <- intersect(stat_cols, colnames(reftable))
  if (length(sim_cols) == 0)
    stop("No stat columns found in reftable. Check stat.cols or observed column names.")
  S_sim <- as.matrix(reftable[, sim_cols, drop = FALSE])

  # Remove non-finite rows
  good <- apply(S_sim, 1, function(x) all(is.finite(x)))
  S_sim <- S_sim[good, , drop = FALSE]

  n_stats <- ncol(S_sim)
  n_sims <- nrow(S_sim)

  if (length(observed) != n_stats)
    stop(sprintf("observed has %d elements but reftable has %d stat columns.",
                 length(observed), n_stats))

  results <- list()

  # =========================================================================
  # CHECK 1: Mahalanobis distance from simulation cloud
  # =========================================================================
  mu_sim <- colMeans(S_sim)
  Sigma_sim <- cov(S_sim)

  # Regularize if needed
  cn <- kappa(Sigma_sim, exact = FALSE)
  if (cn > 1e10) {
    Sigma_sim <- Sigma_sim + diag(1e-6 * max(diag(Sigma_sim)), n_stats)
  }

  d2_maha <- as.numeric(mahalanobis(matrix(observed, nrow = 1), mu_sim, Sigma_sim))
  p_maha <- 1 - pchisq(d2_maha, df = n_stats)

  results$mahalanobis <- list(
    d2 = d2_maha,
    p_value = p_maha,
    pass = p_maha >= alpha
  )

  # =========================================================================
  # CHECK 2: PCA-based distance (robust in high dimensions)
  # =========================================================================
  # Drop constant/zero-variance columns that would break scale. = TRUE
  col_sd <- apply(S_sim, 2, sd, na.rm = TRUE)
  pca_keep <- col_sd > 0 & is.finite(col_sd)
  S_sim_pca <- S_sim[, pca_keep, drop = FALSE]
  obs_pca   <- observed[pca_keep]

  if (sum(!pca_keep) > 0)
    cat(sprintf("PipeMaster:: OOD.diagnose: dropped %d zero-variance stats for PCA (%d remaining)\n",
                sum(!pca_keep), sum(pca_keep)))

  pca_fit <- prcomp(S_sim_pca, center = TRUE, scale. = TRUE)
  var_explained <- cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2))
  n_pcs <- which(var_explained >= pca.var)[1]
  if (is.na(n_pcs)) n_pcs <- n_stats  # fallback

  scores_sim <- pca_fit$x[, 1:n_pcs, drop = FALSE]
  obs_mat <- matrix(obs_pca, nrow = 1)
  colnames(obs_mat) <- colnames(S_sim_pca)
  score_obs <- predict(pca_fit, obs_mat)[, 1:n_pcs, drop = FALSE]

  mu_pca <- colMeans(scores_sim)
  Sigma_pca <- cov(scores_sim)

  cn_pca <- kappa(Sigma_pca, exact = FALSE)
  if (cn_pca > 1e10) {
    Sigma_pca <- Sigma_pca + diag(1e-6 * max(diag(Sigma_pca)), n_pcs)
  }

  d2_pca <- as.numeric(mahalanobis(score_obs, mu_pca, Sigma_pca))
  p_pca <- 1 - pchisq(d2_pca, df = n_pcs)

  results$pca <- list(
    d2 = d2_pca,
    p_value = p_pca,
    n_pcs = n_pcs,
    var_explained = var_explained[n_pcs],
    pass = p_pca >= alpha
  )

  # =========================================================================
  # CHECK 3: Per-stat percentile
  #
  # Classifies each stat as:
  #   "ok"            — observed within [alpha, 1-alpha] range
  #   "outlier"       — observed is extreme relative to meaningful sim variation
  #   "uninformative" — zero or near-zero variance in sims (and obs matches);
  #                     these carry no discriminating power and should not
  #                     count as outliers
  # =========================================================================
  sim_var <- apply(S_sim, 2, var)

  percentiles <- sapply(seq_len(n_stats), function(k) {
    mean(S_sim[, k] <= observed[k])
  })
  names(percentiles) <- sim_cols

  # Classify: uninformative when sim variance is negligible and obs is within
  # sim range (i.e. both are effectively constant at the same value)
  sim_range <- apply(S_sim, 2, function(x) diff(range(x)))
  obs_in_range <- vapply(seq_len(n_stats), function(k) {
    observed[k] >= min(S_sim[, k]) & observed[k] <= max(S_sim[, k])
  }, logical(1))
  uninformative <- (sim_var < .Machine$double.eps) |
                   (sim_range < .Machine$double.eps & obs_in_range)

  extreme <- percentiles < alpha | percentiles > (1 - alpha)
  outlier <- extreme & !uninformative

  reason <- rep("ok", n_stats)
  reason[uninformative] <- "uninformative"
  reason[outlier]       <- "outlier"

  n_outliers      <- sum(outlier)
  n_uninformative <- sum(uninformative)

  results$percentiles <- data.frame(
    stat = sim_cols,
    percentile = round(percentiles, 4),
    outlier = outlier,
    reason = reason,
    stringsAsFactors = FALSE
  )

  # =========================================================================
  # CHECK 4: Model disagreement (requires trained.nn with top-K models)
  #   Emulator (θ→S): predict stats at theta, compare across models
  #   tune.nn  (S→θ): predict params from observed, compare across models
  # =========================================================================
  results$model_disagreement <- NULL

  if (has_model) {
    models <- trained.nn$models
    run_check4 <- !is.null(models) && length(models) > 1

    # For emulator type, also require theta
    if (is_emulator && is.null(theta)) run_check4 <- FALSE

    if (run_check4) {
      feat_mu   <- trained.nn$data$feat_mu
      feat_sd   <- trained.nn$data$feat_sd
      target_mu <- trained.nn$data$target_mu
      target_sd <- trained.nn$data$target_sd

      if (is_emulator) {
        # Forward model: predict stats at theta
        param_cols <- trained.nn$param_cols
        theta_log <- log(theta[param_cols])
        theta_z <- (theta_log - feat_mu) / feat_sd
        X_input <- matrix(theta_z, nrow = 1)

        preds <- sapply(models, function(m) {
          p_z <- as.numeric(.predict.nn(m, X_input))
          p_t <- p_z * target_sd + target_mu
          sign(p_t) * expm1(abs(p_t))
        })
        pred_names <- sim_cols
      } else {
        # Inverse model: predict params from observed stats
        X_input <- .prep.observed(obs_raw, trained.nn$data, model_type, trained.nn$sfs.dims)

        preds <- sapply(models, function(m) {
          p_z <- as.numeric(.predict.nn(m, X_input))
          as.numeric(.inv.transform(matrix(p_z, nrow = 1), target_mu, target_sd))
        })
        pred_names <- names(target_mu)
        if (is.null(pred_names))
          pred_names <- paste0("param_", seq_len(nrow(preds)))
      }
      # preds: n_outputs x n_models

      pred_mean <- rowMeans(preds)
      pred_sd <- apply(preds, 1, sd)

      # Coefficient of variation (relative disagreement)
      cv <- ifelse(abs(pred_mean) > 1e-12, pred_sd / abs(pred_mean), pred_sd)
      names(cv) <- pred_names
      mean_cv <- mean(cv)

      results$model_disagreement <- list(
        mean_cv = mean_cv,
        per_stat_cv = cv,
        reliable = mean_cv < 0.1  # <10% CV is considered reliable
      )
    }
  }

  # =========================================================================
  # OVERALL ASSESSMENT
  # =========================================================================
  checks <- c(
    results$mahalanobis$pass,
    results$pca$pass,
    n_outliers <= ceiling(0.1 * n_stats),  # allow up to 10% outlier stats
    if (!is.null(results$model_disagreement)) results$model_disagreement$reliable else TRUE
  )

  if (all(checks)) {
    results$overall <- "pass"
  } else if (sum(!checks) <= 1) {
    results$overall <- "warn"
  } else {
    results$overall <- "fail"
  }

  # =========================================================================
  # PRINT SUMMARY
  # =========================================================================
  cat(sprintf("PipeMaster:: OOD.diagnose (%s)\n", model_type))
  cat(sprintf("  Mahalanobis: d2=%.2f, p=%.4f [%s]\n",
              d2_maha, p_maha, ifelse(results$mahalanobis$pass, "PASS", "FAIL")))
  cat(sprintf("  PCA (%d PCs, %.0f%% var): d2=%.2f, p=%.4f [%s]\n",
              n_pcs, var_explained[n_pcs] * 100, d2_pca, p_pca,
              ifelse(results$pca$pass, "PASS", "FAIL")))
  cat(sprintf("  Per-stat: %d outliers, %d uninformative, %d ok (of %d)\n",
              n_outliers, n_uninformative, n_stats - n_outliers - n_uninformative, n_stats))
  if (n_outliers > 0) {
    bad <- results$percentiles[results$percentiles$reason == "outlier", ]
    for (r in seq_len(nrow(bad)))
      cat(sprintf("    %s: percentile=%.3f\n", bad$stat[r], bad$percentile[r]))
  }
  if (!is.null(results$model_disagreement)) {
    cat(sprintf("  Model disagreement: mean CV=%.3f [%s]\n",
                mean_cv, ifelse(results$model_disagreement$reliable, "RELIABLE", "UNCERTAIN")))
  }
  cat(sprintf("  Overall: %s\n", toupper(results$overall)))

  # =========================================================================
  # PLOTS
  # =========================================================================
  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

    # Panel 1: PCA scatter (PC1 vs PC2) with observed point
    plot(scores_sim[, 1], scores_sim[, 2],
         col = rgb(0.5, 0.5, 0.5, 0.2), pch = 16, cex = 0.5,
         xlab = "PC1", ylab = "PC2", main = "PCA: observed vs simulated")
    points(score_obs[1], score_obs[2], col = "red", pch = 17, cex = 2)
    legend("topright", legend = c("simulated", "observed"),
           col = c("grey50", "red"), pch = c(16, 17), cex = 0.8)

    # Panel 2: Per-stat percentile barplot
    bar_cols <- ifelse(reason == "uninformative", "grey70",
                       ifelse(reason == "outlier", "firebrick", "steelblue"))
    barplot(percentiles, col = bar_cols, border = NA, las = 2, cex.names = 0.5,
            main = "Per-stat percentiles", ylab = "percentile", ylim = c(0, 1))
    abline(h = c(alpha, 1 - alpha), lty = 2, col = "red")
    legend("bottomright",
           legend = c("ok", "outlier", "uninformative"),
           fill = c("steelblue", "firebrick", "grey70"),
           cex = 0.7, bty = "n")

    # Panel 3: Mahalanobis distance vs chi-sq reference
    d2_sims <- mahalanobis(S_sim, mu_sim, Sigma_sim)
    hist(d2_sims, breaks = 50, col = "grey80", border = "white",
         main = "Mahalanobis distance", xlab = "d2", probability = TRUE)
    abline(v = d2_maha, col = "red", lwd = 2)
    x_seq <- seq(0, max(c(d2_sims, d2_maha)) * 1.1, length.out = 200)
    lines(x_seq, dchisq(x_seq, df = n_stats), col = "blue", lwd = 1.5)
    legend("topright",
           legend = c("simulations", "observed", paste0("chi2(", n_stats, ")")),
           col = c("grey50", "red", "blue"), lwd = c(5, 2, 1.5), cex = 0.7)

    # Panel 4: Model disagreement (if available)
    if (!is.null(results$model_disagreement)) {
      cv_vals <- results$model_disagreement$per_stat_cv
      cols_cv <- ifelse(cv_vals > 0.1, "firebrick", "steelblue")
      barplot(cv_vals, col = cols_cv, border = NA, las = 2, cex.names = 0.5,
              main = "Model disagreement (CV)", ylab = "CV")
      abline(h = 0.1, lty = 2, col = "red")
    }
  }

  class(results) <- "OOD_diagnostic"
  invisible(results)
}
