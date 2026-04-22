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
#' (per-stat percentile, NN density full-space, NN density PCA-space)
#' using only the reftable and observed statistics, without a trained model.
#' In this mode, stat columns are inferred from the observed data.frame
#' column names or matched against the reftable.
#'
#' @section Methods:
#'
#' The per-stat support check is the primary model-fit signal. Each stat is
#' classified as:
#' \itemize{
#'   \item \code{ok}: obs is inside the sim range (the model can produce it
#'     at least once across the prior predictive)
#'   \item \code{outlier}: obs is strictly outside the sim range — the model
#'     cannot produce this observation under the prior
#'   \item \code{uninformative}: sims are constant AND obs matches that
#'     constant (no discriminating information)
#' }
#' The fraction of informative stats flagged as outliers determines the
#' per-stat tier:
#' \itemize{
#'   \item pass: outlier fraction < 10\%
#'   \item warn: 10-25\%
#'   \item fail: > 25\%
#' }
#'
#' Secondary checks test density compatibility on the OOD-filtered stats only
#' (zero-variance and per-stat outliers excluded). Both use a nearest-neighbor
#' \emph{typicality} test against an empirical null: the distribution of
#' leave-one-out sim-to-sim NN distances. The observed NN distance is compared
#' to this null via one-sided empirical p-value. Distribution-free, no Gaussian
#' or chi-square assumption, and robust to correlated high-dimensional stats
#' where Mahalanobis distance is unreliable.
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
#' @param param.cols character vector or NULL — parameter columns in the reftable.
#'   When provided, an additional PLS-space NN density check is computed (PLS
#'   components maximize covariance between stats and params, aligning the
#'   projection with the directions ABC/NN use for inference). When
#'   \code{trained.nn} has \code{param_cols}, this is auto-inferred. When NULL
#'   and not inferable, the function runs on PCA alone.
#' @param theta numeric vector or NULL — optimized parameter estimate (from
#'   \code{emulator.optimize()}). If provided, model disagreement is evaluated
#'   at this point. Only used for emulator-type results; for \code{tune.nn()}
#'   results, model disagreement is evaluated at the observed stats directly.
#' @param pca.var numeric — cumulative variance threshold for PCA (default 0.95).
#' @param pls.n.comp integer — number of PLS components to fit (default 10,
#'   capped at \code{n_params + 5} and \code{n_stats - 1}).
#' @param alpha numeric — significance level for the NN density tests
#'   (default 0.01). Per-stat outliers use a strict support-based criterion
#'   (obs outside sim range) and do not depend on alpha.
#' @param plot logical — produce diagnostic plots (default TRUE).
#'
#' @return A list of class \code{"OOD_diagnostic"} with:
#' \describe{
#'   \item{per_stat}{list with n_outliers, n_uninformative, n_stats,
#'     outlier_frac, tier ("pass"/"warn"/"fail"), pass (logical)}
#'   \item{percentiles}{data.frame with stat, percentile, outlier, reason}
#'   \item{pca_all}{list with d_obs, p_value, pass, null_distribution, n_pcs,
#'     var_explained, n_stats, scores, score_obs (NN density in PCA of all
#'     informative stats — prior-predictive coverage)}
#'   \item{pca_filtered}{same structure as pca_all, but PCA fit on filtered
#'     stats only (zero-var + per-stat outliers removed).}
#'   \item{pls_all, pls_filtered}{same structure as pca_* but using PLS
#'     projection aligned with params. Only present when \code{param.cols}
#'     is available. Primary secondary check when present.}
#'   \item{model_disagreement}{list with mean_cv, per_stat_cv, reliable (logical).
#'     NULL if \code{trained.nn} is NULL or fewer than 2 models.}
#'   \item{overall}{character: \code{"pass"}, \code{"warn"}, or \code{"fail"}}
#' }
#'
#' @section Interpreting the diagnostic plots:
#'
#' A 2x3 panel layout with two rigorous NN density tests (cols 2) and
#' complementary views (cols 1 and 3):
#'
#' \strong{Panel 1 (row 1 col 1) — NN distance vs PC1.}
#' X-axis: PC1 score in the all-stats PCA basis (dominant direction of sim
#' variation). Y-axis: leave-one-out sim-to-sim nearest-neighbor distance
#' in the \emph{full} PC space (all retained PCs, not just PC1-PC2). Sims
#' form a horizontal band at their typical NN distance; obs (red triangle)
#' appears at its PC1 location with its actual NN distance to the cloud. A
#' regular 2D PC1-PC2 scatter can hide extremity along deep PCs because
#' those PCs aren't visualized; this plot solves that because the Y-axis
#' aggregates NN distance across all PCs. Log-scale Y when sim/obs NN
#' range exceeds 100x.
#'
#' \strong{Panel 2 (row 1 col 2) — NN distance (all stats).}
#' Leave-one-out sim-to-sim nearest-neighbor distances in the all-stats PCA
#' space (grey histogram) with the observed NN distance marked (red line).
#' Obs should sit in the bulk, not the right tail.
#'
#' \strong{Panel 3 (row 1 col 3) — Percentile meta-histogram.}
#' Distribution of mid-rank per-stat percentiles across informative stats.
#' Percentile = P(sim < obs) + 0.5 * P(sim == obs). Under perfect fit, the
#' distribution is uniform(0,1) and the histogram is flat. Heavy tails at 0
#' or 1 indicate systematic shift of obs relative to sims. Grey dashed line
#' is the expected uniform count per bin. Visualization only — does not
#' define outliers.
#'
#' \strong{Panel 4 (row 2 col 1) — Outlier stats.}
#' Horizontal bars for the stats where obs is strictly outside sim support.
#' Red bar (right) = obs > max(sim). Blue bar (left) = obs < min(sim).
#' Bar length encodes the mid-rank percentile (0 or 1 by construction since
#' obs is out of range).
#'
#' \strong{Panel 5 (row 2 col 2) — NN distance (filtered stats).}
#' Same method as Panel 2 but with PCA fit on filtered stats only (outliers
#' and zero-var removed). The secondary check that feeds the overall verdict.
#'
#' \strong{Panel 6 (row 2 col 3) — Model disagreement (CV).}
#' Only shown when \code{trained.nn} is provided with multiple top-K models.
#'
#' Use \code{OOD.outliers()} to inspect the full sim distribution for each
#' outlier stat individually.
#'
#' @section Overall assessment:
#'
#' The per-stat tier drives the overall verdict; the secondary checks (NN
#' density + model disagreement) can downgrade a "pass" to "warn", or a
#' "warn" to "fail", but cannot upgrade a "fail":
#' \describe{
#'   \item{per-stat = fail}{overall = FAIL (secondary checks ignored)}
#'   \item{per-stat = warn}{overall = WARN if secondary checks pass, else FAIL}
#'   \item{per-stat = pass}{overall = PASS if secondary checks pass, else WARN}
#' }
#'
#' @export
OOD.diagnose <- function(trained.nn = NULL, observed, reftable,
                        stat.cols = NULL, param.cols = NULL,
                        theta = NULL,
                        pca.var = 0.95,
                        pls.n.comp = 10L,
                        alpha = 0.01,
                        plot = TRUE) {

  # --- Detect model type and determine stat + param columns ---
  has_model <- !is.null(trained.nn)

  if (has_model) {
    model_type <- trained.nn$type  # "emulator", "sumstat", "sfs1d", "sfs2d"
    if (is.null(model_type)) model_type <- "emulator"
    is_emulator <- model_type == "emulator"
    stat_cols <- trained.nn$data$stat_cols
    if (is.null(param.cols) && !is.null(trained.nn$param_cols))
      param.cols <- trained.nn$param_cols
  } else {
    model_type <- "none"
    is_emulator <- FALSE

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
  obs_raw <- observed

  sim_cols <- intersect(stat_cols, colnames(reftable))
  if (length(sim_cols) == 0)
    stop("No stat columns found in reftable. Check stat.cols or observed column names.")
  S_sim <- as.matrix(reftable[, sim_cols, drop = FALSE])

  # Extract params (for PLS check); NULL if unavailable
  P_sim <- NULL
  if (!is.null(param.cols)) {
    missing_params <- setdiff(param.cols, colnames(reftable))
    if (length(missing_params) > 0) {
      warning(sprintf("param.cols not found in reftable: %s — PLS check skipped",
                      paste(missing_params, collapse = ", ")))
      param.cols <- NULL
    } else {
      P_sim <- as.matrix(reftable[, param.cols, drop = FALSE])
    }
  }

  good <- apply(S_sim, 1, function(x) all(is.finite(x)))
  if (!is.null(P_sim))
    good <- good & apply(P_sim, 1, function(x) all(is.finite(x)))
  S_sim <- S_sim[good, , drop = FALSE]
  if (!is.null(P_sim)) P_sim <- P_sim[good, , drop = FALSE]

  n_stats <- ncol(S_sim)
  n_sims <- nrow(S_sim)

  if (length(observed) != n_stats)
    stop(sprintf("observed has %d elements but reftable has %d stat columns.",
                 length(observed), n_stats))

  results <- list()

  # =========================================================================
  # CHECK 1: Per-stat support (primary model-fit signal)
  # Outlier = obs strictly outside sim range (model cannot produce obs).
  # =========================================================================
  sim_var   <- apply(S_sim, 2, var)
  sim_range <- apply(S_sim, 2, function(x) diff(range(x)))
  col_sd    <- apply(S_sim, 2, sd, na.rm = TRUE)

  # Mid-rank percentile (for visualization only — ties contribute 0.5).
  # Outlier classification below uses strict support-based criterion, not
  # a percentile threshold.
  percentiles <- sapply(seq_len(n_stats), function(k) {
    mean(S_sim[, k] < observed[k]) + 0.5 * mean(S_sim[, k] == observed[k])
  })
  names(percentiles) <- sim_cols

  sim_min <- apply(S_sim, 2, min)
  sim_max <- apply(S_sim, 2, max)
  obs_in_range <- observed >= sim_min & observed <= sim_max

  # Classification (support-based):
  #   uninformative : sim is constant AND obs matches that constant (no info)
  #   outlier       : obs is strictly outside the sim support (impossible
  #                   under the prior predictive — the model cannot produce
  #                   this observation)
  #   ok            : obs is inside the sim support, sims not constant
  sim_constant  <- sim_var < .Machine$double.eps | sim_range < .Machine$double.eps
  uninformative <- sim_constant & obs_in_range
  outlier       <- !obs_in_range

  reason <- rep("ok", n_stats)
  reason[uninformative] <- "uninformative"
  reason[outlier]       <- "outlier"

  n_outliers      <- sum(outlier)
  n_uninformative <- sum(uninformative)
  n_informative   <- n_stats - n_uninformative
  outlier_frac    <- if (n_informative > 0) n_outliers / n_informative else 0

  per_stat_tier <- if (outlier_frac < 0.10) "pass"
                   else if (outlier_frac < 0.25) "warn"
                   else "fail"

  results$percentiles <- data.frame(
    stat = sim_cols, percentile = round(percentiles, 4),
    outlier = outlier, reason = reason, stringsAsFactors = FALSE
  )
  results$per_stat <- list(
    n_outliers = n_outliers, n_uninformative = n_uninformative,
    n_stats = n_stats, outlier_frac = outlier_frac, tier = per_stat_tier,
    pass = per_stat_tier != "fail"
  )

  # Filter for the density checks: drop zero-var and per-stat outliers
  keep_stats <- col_sd > 0 & is.finite(col_sd) & !outlier & !uninformative
  n_kept <- sum(keep_stats)
  pca_keep_all <- col_sd > 0 & is.finite(col_sd)

  # =========================================================================
  # CHECK 2 & 3: NN density in PCA space (all-informative and filtered)
  #
  # Same method on both stat sets: PCA with scale=TRUE, project obs,
  # compute Euclidean NN distance in PC-score space, compare to the
  # leave-one-out sim-to-sim NN distribution (empirical null, no Gaussian
  # assumption). Two runs differ only in which stats go into the PCA.
  # =========================================================================
  .pca_nn <- function(keep_mask) {
    if (sum(keep_mask) < 2) {
      return(list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                  null_distribution = numeric(0), n_pcs = NA_integer_,
                  var_explained = NA_real_, n_stats = sum(keep_mask),
                  scores = NULL, score_obs = NULL))
    }
    S <- S_sim[, keep_mask, drop = FALSE]
    o <- observed[keep_mask]
    pca_fit <- prcomp(S, center = TRUE, scale. = TRUE)
    var_cum <- cumsum(pca_fit$sdev^2 / sum(pca_fit$sdev^2))
    n_pcs <- which(var_cum >= pca.var)[1]
    if (is.na(n_pcs)) n_pcs <- ncol(S)
    n_pcs <- max(n_pcs, 2)

    scores <- pca_fit$x[, 1:n_pcs, drop = FALSE]
    obs_mat <- matrix(o, nrow = 1); colnames(obs_mat) <- colnames(S)
    score_obs <- predict(pca_fit, obs_mat)[, 1:n_pcs, drop = FALSE]

    d_obs <- min(sqrt(rowSums(sweep(scores, 2, score_obs[1, ], "-")^2)))
    D <- as.matrix(dist(scores)); diag(D) <- Inf
    d_loo <- apply(D, 1, min)
    p <- mean(d_loo >= d_obs)

    list(d_obs = d_obs, p_value = p, pass = p >= alpha,
         null_distribution = d_loo, n_pcs = n_pcs,
         var_explained = var_cum[n_pcs], n_stats = ncol(S),
         scores = scores, score_obs = score_obs,
         pca_fit = pca_fit, keep_mask = keep_mask)
  }

  results$pca_all      <- .pca_nn(pca_keep_all)
  results$pca_filtered <- .pca_nn(keep_stats)
  if (n_kept < 2)
    warning("Fewer than 2 informative non-outlier stats; filtered NN density check skipped.")

  # Projection of obs into the all-stats PCA basis with outlier-stat values
  # replaced by their sim means. Visualization only: shows where obs would
  # sit if outlier stats were neutralized. Same basis as results$pca_all so
  # the two points are directly comparable on the scatter.
  results$pca_all$score_obs_projected <- NULL
  if (!is.null(results$pca_all$scores) && sum(outlier) > 0) {
    obs_proj <- observed
    obs_proj[outlier] <- colMeans(S_sim)[outlier]
    obs_proj_mat <- matrix(obs_proj[pca_keep_all], nrow = 1)
    colnames(obs_proj_mat) <-
      colnames(S_sim[, pca_keep_all, drop = FALSE])
    n_pcs_all <- results$pca_all$n_pcs
    results$pca_all$score_obs_projected <-
      predict(results$pca_all$pca_fit,
              obs_proj_mat)[, 1:n_pcs_all, drop = FALSE]
  }

  # =========================================================================
  # CHECK 2b & 3b: NN density in PLS space (all-informative and filtered)
  #
  # Only computed when param.cols is provided (or inferred from trained.nn).
  # PLS components maximize cov(stats, params) — directions aligned with
  # what ABC/NN use for inference. NN density in PLS space answers "is obs
  # typical along the param-informative directions of sim variation?"
  # =========================================================================
  # Normalize params: sign-safe log-transform (keeps sign, handles zeros),
  # then z-score. Prevents high-scale params (e.g., Ne) from dominating PLS
  # over low-scale ones (e.g., mig).
  P_sim_norm <- NULL
  if (!is.null(P_sim)) {
    Pl <- sign(P_sim) * log1p(abs(P_sim))
    pl_mu <- colMeans(Pl)
    pl_sd <- apply(Pl, 2, sd)
    pl_sd[pl_sd == 0 | !is.finite(pl_sd)] <- 1
    P_sim_norm <- sweep(sweep(Pl, 2, pl_mu, "-"), 2, pl_sd, "/")
  }

  .pls_nn <- function(keep_mask) {
    if (sum(keep_mask) < 2)
      return(list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                  null_distribution = numeric(0), n_comp = NA_integer_,
                  n_stats = sum(keep_mask),
                  scores = NULL, score_obs = NULL, pls_fit = NULL,
                  keep_mask = keep_mask))
    S <- S_sim[, keep_mask, drop = FALSE]
    o <- observed[keep_mask]
    n_comp <- min(pls.n.comp, sum(keep_mask) - 1, ncol(P_sim_norm) + 5)
    pls_fit <- pls.fit(stats = S, params = P_sim_norm,
                       n.comp = n_comp, scale = TRUE)
    n_comp_eff <- pls_fit$n.comp  # pls.fit may drop zero-var cols
    scores     <- pls_fit$scores
    score_obs  <- pls.project(pls_fit, matrix(o, nrow = 1,
                                              dimnames = list(NULL, colnames(S))))

    d_obs <- min(sqrt(rowSums(sweep(scores, 2, score_obs[1, ], "-")^2)))
    D <- as.matrix(dist(scores)); diag(D) <- Inf
    d_loo <- apply(D, 1, min)
    p <- mean(d_loo >= d_obs)

    list(d_obs = d_obs, p_value = p, pass = p >= alpha,
         null_distribution = d_loo, n_comp = n_comp_eff,
         n_stats = ncol(S),
         scores = scores, score_obs = score_obs,
         pls_fit = pls_fit, keep_mask = keep_mask)
  }

  if (!is.null(P_sim)) {
    results$pls_all      <- .pls_nn(pca_keep_all)
    results$pls_filtered <- .pls_nn(keep_stats)

    # Projection: obs with outlier stats replaced by sim means, in PLS basis
    results$pls_all$score_obs_projected <- NULL
    if (!is.null(results$pls_all$scores) && sum(outlier) > 0) {
      obs_proj <- observed
      obs_proj[outlier] <- colMeans(S_sim)[outlier]
      obs_proj_mat <- matrix(obs_proj[pca_keep_all], nrow = 1,
                             dimnames = list(NULL,
                               colnames(S_sim[, pca_keep_all, drop = FALSE])))
      results$pls_all$score_obs_projected <-
        pls.project(results$pls_all$pls_fit, obs_proj_mat)
    }
  } else {
    results$pls_all <- NULL
    results$pls_filtered <- NULL
  }

  # =========================================================================
  # CHECK 4: Model disagreement (requires trained.nn with top-K models)
  # =========================================================================
  results$model_disagreement <- NULL

  if (has_model) {
    models <- trained.nn$models
    run_check4 <- !is.null(models) && length(models) > 1
    if (is_emulator && is.null(theta)) run_check4 <- FALSE

    if (run_check4) {
      feat_mu   <- trained.nn$data$feat_mu
      feat_sd   <- trained.nn$data$feat_sd
      target_mu <- trained.nn$data$target_mu
      target_sd <- trained.nn$data$target_sd

      if (is_emulator) {
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
        X_input <- .prep.observed(obs_raw, trained.nn$data, model_type, trained.nn$sfs.dims)

        preds <- sapply(models, function(m) {
          p_z <- as.numeric(.predict.nn(m, X_input))
          as.numeric(.inv.transform(matrix(p_z, nrow = 1), target_mu, target_sd))
        })
        pred_names <- names(target_mu)
        if (is.null(pred_names))
          pred_names <- paste0("param_", seq_len(nrow(preds)))
      }

      pred_mean <- rowMeans(preds)
      pred_sd <- apply(preds, 1, sd)
      cv <- ifelse(abs(pred_mean) > 1e-12, pred_sd / abs(pred_mean), pred_sd)
      names(cv) <- pred_names
      mean_cv <- mean(cv)

      results$model_disagreement <- list(
        mean_cv = mean_cv, per_stat_cv = cv, reliable = mean_cv < 0.1
      )
    }
  }

  # =========================================================================
  # OVERALL ASSESSMENT
  #   per-stat tier is primary; secondary checks can downgrade but not upgrade
  # =========================================================================
  pca_filt_pass <- isTRUE(results$pca_filtered$pass)
  pls_filt_pass <- if (!is.null(results$pls_filtered))
                     isTRUE(results$pls_filtered$pass) else TRUE
  md_pass       <- if (!is.null(results$model_disagreement))
                     results$model_disagreement$reliable else TRUE
  secondary_ok  <- pca_filt_pass && pls_filt_pass && md_pass

  results$overall <- if (per_stat_tier == "fail") {
    "fail"
  } else if (per_stat_tier == "warn") {
    if (secondary_ok) "warn" else "fail"
  } else {
    if (secondary_ok) "pass" else "warn"
  }

  # =========================================================================
  # PRINT SUMMARY
  # =========================================================================
  cat(sprintf("PipeMaster:: OOD.diagnose (%s)\n", model_type))
  cat(sprintf("  Stats: %d total  =  %d informative (%d ok + %d outlier)  +  %d uninformative (sim=obs constant)\n",
              n_stats, n_informative, n_informative - n_outliers, n_outliers,
              n_uninformative))
  cat(sprintf("  Per-stat tier: outlier_frac = %d/%d = %.1f%% [%s]\n",
              n_outliers, n_informative, 100 * outlier_frac,
              toupper(per_stat_tier)))
  if (n_outliers > 0 && n_outliers <= 20) {
    bad <- results$percentiles[results$percentiles$reason == "outlier", ]
    for (r in seq_len(nrow(bad)))
      cat(sprintf("    %s: percentile=%.3f\n", bad$stat[r], bad$percentile[r]))
  } else if (n_outliers > 20) {
    cat(sprintf("    (%d outliers — inspect sim distributions with OOD.outliers())\n",
                n_outliers))
  }
  if (!is.na(results$pca_all$pass))
    cat(sprintf("  NN density (PCA all,      %3d stats, %d PCs, %.0f%% var): d_obs=%.3f, p=%.4f [%s]\n",
                results$pca_all$n_stats, results$pca_all$n_pcs,
                results$pca_all$var_explained * 100,
                results$pca_all$d_obs, results$pca_all$p_value,
                ifelse(results$pca_all$pass, "PASS", "FAIL")))
  if (!is.na(results$pca_filtered$pass))
    cat(sprintf("  NN density (PCA filtered, %3d stats, %d PCs, %.0f%% var): d_obs=%.3f, p=%.4f [%s]\n",
                results$pca_filtered$n_stats, results$pca_filtered$n_pcs,
                results$pca_filtered$var_explained * 100,
                results$pca_filtered$d_obs, results$pca_filtered$p_value,
                ifelse(results$pca_filtered$pass, "PASS", "FAIL")))
  if (!is.null(results$pls_all) && !is.na(results$pls_all$pass))
    cat(sprintf("  NN density (PLS all,      %3d stats, %d comps): d_obs=%.3f, p=%.4f [%s]\n",
                results$pls_all$n_stats, results$pls_all$n_comp,
                results$pls_all$d_obs, results$pls_all$p_value,
                ifelse(results$pls_all$pass, "PASS", "FAIL")))
  if (!is.null(results$pls_filtered) && !is.na(results$pls_filtered$pass))
    cat(sprintf("  NN density (PLS filtered, %3d stats, %d comps): d_obs=%.3f, p=%.4f [%s]\n",
                results$pls_filtered$n_stats, results$pls_filtered$n_comp,
                results$pls_filtered$d_obs, results$pls_filtered$p_value,
                ifelse(results$pls_filtered$pass, "PASS", "FAIL")))
  if (!is.null(results$model_disagreement)) {
    cat(sprintf("  Model disagreement: mean CV=%.3f [%s]\n",
                results$model_disagreement$mean_cv,
                ifelse(results$model_disagreement$reliable, "RELIABLE", "UNCERTAIN")))
  }
  cat(sprintf("  Overall: %s\n", toupper(results$overall)))

  # =========================================================================
  # PLOTS (2x3 layout)
  #   Row 1: PCA scatter (obs before/after)   NN hist all         per-stat meta-hist
  #   Row 2: Top-N outlier bars               NN hist filtered    model disagreement
  # =========================================================================
  # Pick PLS when available (more relevant to ABC/SML), fall back to PCA
  if (!is.null(results$pls_all)) {
    density_all      <- results$pls_all
    density_filtered <- results$pls_filtered
    density_method   <- "PLS"
    density_axis_n   <- results$pls_all$n_comp
    density_x_label  <- "PLS1 (sim score)"
  } else {
    density_all      <- results$pca_all
    density_filtered <- results$pca_filtered
    density_method   <- "PCA"
    density_axis_n   <- results$pca_all$n_pcs
    density_x_label  <- "PC1 (sim score)"
  }

  if (plot) {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

    .plot_nn_hist <- function(nn_res, label, method, n_axis) {
      if (is.na(nn_res$pass)) {
        plot.new(); title(main = sprintf("NN distance (%s) — unavailable", label))
        return(invisible(NULL))
      }
      hist(nn_res$null_distribution, breaks = 50,
           col = "grey80", border = "white", probability = TRUE,
           main = sprintf("NN distance %s (%s, %d dims, p=%.3f) [%s]",
                          label, method, n_axis, nn_res$p_value,
                          ifelse(nn_res$pass, "PASS", "FAIL")),
           xlab = sprintf("nearest-neighbor distance (%s-space)", method),
           xlim = range(c(nn_res$null_distribution, nn_res$d_obs)))
      abline(v = nn_res$d_obs, col = "red", lwd = 2)
      legend("topright", legend = c("sim LOO NN", "observed NN"),
             col = c("grey50", "red"), lwd = c(5, 2), cex = 0.7)
    }

    # Row 1 col 1: NN distance vs first component of projection basis.
    #   X-axis: first component score (PLS1 if available, else PC1)
    #   Y-axis: LOO NN distance to nearest sim in full component space
    #
    # Sims form a horizontal band at their typical NN distance. Obs's Y
    # value captures anomaly across all components — no hiding on deep
    # components not plotted. Arrow to projected obs shows where obs
    # would sit if outlier stats were neutralized. Log-scale Y when
    # obs / sim NN range exceeds 100x.
    if (!is.null(density_all$scores) && !is.na(density_all$pass)) {
      y_sim <- density_all$null_distribution
      y_obs <- density_all$d_obs
      # For the projected obs, compute its full-space NN distance so it
      # appears at the honest height on the Y-axis
      y_proj <- NULL; x_proj <- NULL
      if (!is.null(density_all$score_obs_projected)) {
        y_proj <- min(sqrt(rowSums(sweep(density_all$scores, 2,
          density_all$score_obs_projected[1, ], "-")^2)))
        x_proj <- density_all$score_obs_projected[1]
      }
      y_all <- c(y_sim, y_obs, y_proj)
      y_range_ratio <- if (min(y_sim) > 0) max(y_all) / min(y_sim) else Inf
      use_log_y <- is.finite(y_range_ratio) && y_range_ratio > 100
      plot(density_all$scores[, 1], y_sim,
           col = rgb(0.5, 0.5, 0.5, 0.3), pch = 16, cex = 0.5,
           log = if (use_log_y) "y" else "",
           ylim = range(y_all),
           xlab = density_x_label,
           ylab = sprintf("NN distance to nearest sim (full %d-%s space)",
                          density_axis_n, density_method),
           main = sprintf("NN vs %s1 (%d informative of %d stats, %s basis)",
                          density_method, density_all$n_stats, n_stats,
                          density_method))
      points(density_all$score_obs[1], y_obs, col = "red", pch = 17, cex = 2)
      abline(h = y_obs, col = "red", lty = 3, lwd = 0.8)
      if (!is.null(y_proj)) {
        arrows(density_all$score_obs[1], y_obs, x_proj, y_proj,
               col = "grey40", lwd = 1.5, length = 0.1)
        points(x_proj, y_proj, col = "orange", pch = 2, cex = 2, lwd = 2)
        legend("topright",
               legend = c("simulated (LOO NN)", "observed NN",
                          "obs projected\n(outliers neutralized)"),
               col = c("grey50", "red", "orange"),
               pch = c(16, 17, 2), cex = 0.7, bty = "n")
      } else {
        legend("topright",
               legend = c("simulated (LOO NN)", "observed NN"),
               col = c("grey50", "red"), pch = c(16, 17), cex = 0.8)
      }
    } else {
      plot.new(); title(main = sprintf("NN vs %s1 — unavailable", density_method))
    }

    # Row 1 col 2: NN distance (all)
    .plot_nn_hist(density_all, "all", density_method, density_axis_n)

    # Row 1 col 3: per-stat meta-histogram (visualization only; outliers
    # are defined by support, not by percentile bins)
    informative_pct <- percentiles[!uninformative]
    bin_width <- 0.05
    hist_breaks <- seq(0, 1, by = bin_width)
    h <- hist(informative_pct, breaks = hist_breaks, plot = FALSE)
    unif_count <- length(informative_pct) * bin_width
    plot(h, col = "steelblue", border = "white",
         main = sprintf("Percentile distribution (%d informative: %d ok + %d outlier; %d uninformative) [%s]",
                        length(informative_pct), length(informative_pct) - n_outliers,
                        n_outliers, n_uninformative, toupper(per_stat_tier)),
         xlab = "obs mid-rank percentile in sim distribution", ylab = "n stats")
    abline(h = unif_count, col = "grey40", lty = 2)
    legend("top", legend = "uniform (perfect fit)",
           col = "grey40", lty = 2, cex = 0.7, bty = "n")

    # Row 2 col 1: Outlier bars (obs strictly outside sim support).
    # Percentile is 0 or 1 by construction, so we encode direction (above
    # vs below) and rank stats in the order they appear in the result.
    if (n_outliers > 0) {
      bad <- results$percentiles[results$percentiles$reason == "outlier", ]
      top_n <- min(25, nrow(bad))
      bad_top <- bad[seq_len(top_n), , drop = FALSE]
      bad_top <- bad_top[rev(seq_len(top_n)), , drop = FALSE]
      bar_dir_col <- ifelse(bad_top$percentile >= 0.5, "firebrick", "royalblue")
      bar_heights <- ifelse(bad_top$percentile >= 0.5, 1, -1)
      old_mar <- par("mar")
      par(mar = c(4, 9, 3, 1))
      barplot(bar_heights, names.arg = bad_top$stat,
              horiz = TRUE, las = 1, cex.names = 0.6,
              col = bar_dir_col, border = NA,
              xlim = c(-1.1, 1.1),
              xlab = "direction (obs < min sim  |  obs > max sim)",
              main = sprintf("Outliers: %d of %d (red: obs>max sim, blue: obs<min sim)",
                             top_n, nrow(bad)))
      abline(v = 0, col = "grey60", lty = 1)
      par(mar = old_mar)
    } else {
      plot.new(); title(main = "No outliers")
    }

    # Row 2 col 2: NN distance (filtered)
    n_axis_filt <- if (density_method == "PLS")
                     density_filtered$n_comp else density_filtered$n_pcs
    .plot_nn_hist(density_filtered, "filtered", density_method, n_axis_filt)

    # Row 2 col 3: Model disagreement (if available)
    if (!is.null(results$model_disagreement)) {
      cv_vals <- results$model_disagreement$per_stat_cv
      cols_cv <- ifelse(cv_vals > 0.1, "firebrick", "steelblue")
      barplot(cv_vals, col = cols_cv, border = NA, las = 2, cex.names = 0.5,
              main = "Model disagreement (CV)", ylab = "CV")
      abline(h = 0.1, lty = 2, col = "red")
    } else {
      plot.new()
    }
  }

  class(results) <- "OOD_diagnostic"
  invisible(results)
}


#' Diagnose why the "obs projected (outliers neutralized)" arrow lands where
#' it does in the PCA scatter of \code{OOD.diagnose()}.
#'
#' Runs three complementary checks:
#' \enumerate{
#'   \item Full-PC NN distance: does the projection actually bring obs
#'     closer to the sim cloud when all retained PCs are considered, or
#'     only in PC1-PC2? Also shows cumulative distance vs number of PCs.
#'   \item Per-PC obs extremity: absolute z-score of obs on each PC
#'     (score / sim PC sd). If obs is extreme on high PCs that are not
#'     visualized in the scatter, the 2D picture under-represents misfit.
#'   \item Top outlier-stat loadings on PC1 and PC2: which specific outlier
#'     stats are pulling obs along PC1 and PC2, driving the visible arrow
#'     direction when they are neutralized.
#' }
#'
#' @param ood_result output of \code{OOD.diagnose()}
#' @param pdf_file character or NULL — write the 2x2 plot to this PDF.
#' @param top_n integer — how many top stats to show in the loading plots.
#' @param pdf_width,pdf_height numeric — PDF dimensions if pdf_file given.
#'
#' @return Invisibly a list with:
#' \describe{
#'   \item{nn_obs_full}{Full-PC NN distance from obs to nearest sim}
#'   \item{nn_proj_full}{Full-PC NN distance from projected obs to nearest sim}
#'   \item{cumulative_distance}{matrix with columns obs, projected — cumulative
#'     Euclidean distance to nearest sim as PCs are added one by one}
#'   \item{obs_pc_z}{numeric — |score_obs[k]| / sd(scores[,k]) per PC}
#'   \item{loadings_pc1, loadings_pc2}{data.frames of outlier-stat loadings,
#'     sorted by |loading| descending}
#' }
#' @export
OOD.projection.diagnose <- function(ood_result, pdf_file = NULL,
                                    top_n = 15,
                                    pdf_width = 12, pdf_height = 10) {
  if (!inherits(ood_result, "OOD_diagnostic"))
    stop("ood_result must be the output of OOD.diagnose().")

  pa <- ood_result$pca_all
  if (is.null(pa$scores) || is.null(pa$pca_fit))
    stop("ood_result$pca_all$scores / pca_fit missing — rerun OOD.diagnose().")
  if (is.null(pa$score_obs_projected))
    stop("No projected obs in result — no outlier stats, nothing to diagnose.")

  scores     <- pa$scores
  score_obs  <- as.numeric(pa$score_obs[1, ])
  score_proj <- as.numeric(pa$score_obs_projected[1, ])
  n_pcs      <- ncol(scores)
  pca_fit    <- pa$pca_fit

  # --- (A) Full-PC NN distance: obs vs projected ---
  nn_obs_full  <- min(sqrt(rowSums(sweep(scores, 2, score_obs,  "-")^2)))
  nn_proj_full <- min(sqrt(rowSums(sweep(scores, 2, score_proj, "-")^2)))

  # Cumulative distance as PCs are added: for each k in 1..n_pcs, compute
  # NN distance using only the first k PCs.
  cumdist <- matrix(NA_real_, nrow = n_pcs, ncol = 2,
                    dimnames = list(NULL, c("obs", "projected")))
  for (k in seq_len(n_pcs)) {
    cumdist[k, 1] <- min(sqrt(rowSums(
      sweep(scores[, 1:k, drop = FALSE], 2, score_obs[1:k],  "-")^2)))
    cumdist[k, 2] <- min(sqrt(rowSums(
      sweep(scores[, 1:k, drop = FALSE], 2, score_proj[1:k], "-")^2)))
  }

  # --- (B) per-PC obs extremity ---
  pc_sd <- apply(scores, 2, sd)
  obs_pc_z  <- abs(score_obs)  / pc_sd
  proj_pc_z <- abs(score_proj) / pc_sd

  # --- (C) outlier-stat loadings on PC1 and PC2 ---
  rotation <- pca_fit$rotation
  stat_names_pca <- rownames(rotation)
  outlier_stats <- ood_result$percentiles$stat[
    ood_result$percentiles$reason == "outlier"]
  outlier_in_pca <- intersect(outlier_stats, stat_names_pca)

  get_load_df <- function(pc_idx) {
    loads <- rotation[outlier_in_pca, pc_idx]
    df <- data.frame(stat = names(loads), loading = as.numeric(loads),
                     stringsAsFactors = FALSE)
    df <- df[order(-abs(df$loading)), ]
    head(df, top_n)
  }
  load_pc1 <- get_load_df(1)
  load_pc2 <- get_load_df(2)

  # --- plotting ---
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

  # Panel 1: cumulative NN distance vs #PCs
  plot(seq_len(n_pcs), cumdist[, 1], type = "b", col = "red", pch = 16,
       ylim = range(cumdist, na.rm = TRUE),
       xlab = "number of PCs included", ylab = "NN distance to nearest sim",
       main = sprintf("Full-PC NN: obs=%.3f  projected=%.3f",
                      nn_obs_full, nn_proj_full))
  lines(seq_len(n_pcs), cumdist[, 2], type = "b", col = "orange", pch = 2)
  legend("bottomright",
         legend = c("obs", "obs projected (outliers neutralized)"),
         col = c("red", "orange"), pch = c(16, 2), lty = 1, cex = 0.8)

  # Panel 2: per-PC obs z-score
  bar_mat <- rbind(obs = obs_pc_z, projected = proj_pc_z)
  barplot(bar_mat, beside = TRUE,
          names.arg = paste0("PC", seq_len(n_pcs)),
          col = c("red", "orange"), border = NA,
          main = sprintf("Per-PC obs extremity (|score| / PC sd)  [n_pcs=%d]",
                         n_pcs),
          ylab = "z-score", las = 1, cex.names = 0.8)
  abline(h = c(2, 3), lty = 2, col = "grey60")
  legend("topright", legend = c("obs", "projected"),
         fill = c("red", "orange"), cex = 0.8, bty = "n")

  # Panel 3: top outlier-stat loadings on PC1
  if (nrow(load_pc1) > 0) {
    par(mar = c(4, 9, 3, 1))
    load_pc1_plot <- load_pc1[rev(seq_len(nrow(load_pc1))), ]
    barplot(load_pc1_plot$loading,
            names.arg = load_pc1_plot$stat, horiz = TRUE, las = 1,
            col = ifelse(load_pc1_plot$loading > 0, "firebrick", "royalblue"),
            border = NA, cex.names = 0.7,
            main = sprintf("Top %d outlier-stat loadings on PC1",
                           nrow(load_pc1)),
            xlab = "PC1 loading")
    abline(v = 0, col = "grey60")
    par(mar = c(4, 4, 3, 1))
  } else {
    plot.new(); title(main = "No outlier stats in PCA basis")
  }

  # Panel 4: top outlier-stat loadings on PC2
  if (nrow(load_pc2) > 0) {
    par(mar = c(4, 9, 3, 1))
    load_pc2_plot <- load_pc2[rev(seq_len(nrow(load_pc2))), ]
    barplot(load_pc2_plot$loading,
            names.arg = load_pc2_plot$stat, horiz = TRUE, las = 1,
            col = ifelse(load_pc2_plot$loading > 0, "firebrick", "royalblue"),
            border = NA, cex.names = 0.7,
            main = sprintf("Top %d outlier-stat loadings on PC2",
                           nrow(load_pc2)),
            xlab = "PC2 loading")
    abline(v = 0, col = "grey60")
    par(mar = c(4, 4, 3, 1))
  } else {
    plot.new()
  }

  invisible(list(
    nn_obs_full = nn_obs_full, nn_proj_full = nn_proj_full,
    cumulative_distance = cumdist,
    obs_pc_z = obs_pc_z, proj_pc_z = proj_pc_z,
    loadings_pc1 = load_pc1, loadings_pc2 = load_pc2
  ))
}


#' Plot simulated distributions for OOD outlier statistics
#'
#' Small-multiples diagnostic: one histogram per outlier statistic showing the
#' simulation distribution with the observed value overlaid as a vertical line.
#' Complements the meta-histogram in \code{OOD.diagnose()} by letting the user
#' inspect each problematic stat in detail.
#'
#' Outliers are sorted by extremity (distance from median percentile, i.e.
#' \code{max(percentile, 1 - percentile)}) and truncated to \code{max_stats}.
#'
#' @param ood_result output of \code{OOD.diagnose()} (class \code{OOD_diagnostic})
#' @param observed numeric vector or 1-row data.frame of observed stats, same
#'   as passed to \code{OOD.diagnose()}
#' @param reftable data.frame of simulated stats (the reference table)
#' @param max_stats integer — maximum number of outlier stats to plot
#'   (default 24). Most extreme first.
#' @param ncol integer — number of histograms per row in the grid (default 4).
#' @param pdf_file character or NULL — if provided, writes the plot to this
#'   PDF file. Otherwise plots to the active device.
#' @param pdf_width,pdf_height numeric — PDF dimensions if \code{pdf_file} is given.
#'
#' @return Invisibly returns the outlier data.frame (sorted by extremity, truncated).
#' @export
OOD.outliers <- function(ood_result, observed, reftable,
                         max_stats = 24, ncol = 4,
                         pdf_file = NULL,
                         pdf_width = 12, pdf_height = 9) {
  if (!inherits(ood_result, "OOD_diagnostic"))
    stop("ood_result must be the output of OOD.diagnose().")

  outliers <- ood_result$percentiles[ood_result$percentiles$reason == "outlier", ]
  if (nrow(outliers) == 0) {
    message("No outlier stats to plot.")
    return(invisible(NULL))
  }

  outliers$extremity <- pmax(outliers$percentile, 1 - outliers$percentile)
  outliers <- outliers[order(-outliers$extremity), , drop = FALSE]
  n_total <- nrow(outliers)
  truncated <- n_total > max_stats
  if (truncated) outliers <- outliers[seq_len(max_stats), , drop = FALSE]
  n_plots <- nrow(outliers)

  stat_cols <- ood_result$percentiles$stat
  if (is.data.frame(observed) || is.matrix(observed)) {
    if (all(stat_cols %in% colnames(observed)))
      obs_vec <- setNames(as.numeric(observed[1, stat_cols]), stat_cols)
    else
      obs_vec <- setNames(as.numeric(observed[1, ]), stat_cols)
  } else {
    obs_vec <- setNames(as.numeric(observed), stat_cols)
  }

  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    on.exit(dev.off(), add = TRUE)
  }

  nrow_plot <- ceiling(n_plots / ncol)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow_plot, ncol), mar = c(3, 3, 2.5, 1),
      oma = c(0, 0, if (truncated) 2 else 0, 0),
      mgp = c(1.8, 0.6, 0))

  for (i in seq_len(n_plots)) {
    stat <- outliers$stat[i]
    pct  <- outliers$percentile[i]
    sim_vals <- reftable[[stat]]
    if (is.null(sim_vals)) {
      plot.new(); title(main = sprintf("%s (not in reftable)", stat), cex.main = 0.8)
      next
    }
    sim_vals <- sim_vals[is.finite(sim_vals)]
    obs_val  <- obs_vec[stat]

    obs_col <- if (pct > 0.5) "firebrick" else "royalblue"

    hist(sim_vals, breaks = 40, col = "grey85", border = "white",
         main = sprintf("%s\npct=%.3f", stat, pct),
         xlab = "", ylab = "", cex.main = 0.8,
         xlim = range(c(sim_vals, obs_val), na.rm = TRUE))
    abline(v = obs_val, col = obs_col, lwd = 2)
  }

  if (truncated) {
    mtext(sprintf("Showing top %d of %d outliers (sorted by extremity)",
                  max_stats, n_total),
          outer = TRUE, cex = 1)
  }

  invisible(outliers)
}
