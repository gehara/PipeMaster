# ============================================================================
# Out-of-Distribution Diagnostics for PipeMaster
#
# Two-tier API:
#   OOD.pretrain()  — prior-predictive coverage (no trained NN required)
#                     per-stat support + NN density in PCA + NN density in PLS
#   OOD.posttrain() — model-fit + posterior fidelity (requires trained NN)
#                     NN-latent NN density + ensemble disagreement
#                     Optionally consumes an OOD.pretrain() result to reuse
#                     the per-stat / PCA / PLS context for verdict and plots.
#
# Forensics for either tier: OOD.projection.diagnose() in R/OOD.projection.R
# (class-dispatched on OOD_pretrain vs OOD_posttrain).
#
# Both public functions share the internal NN-density / projection helpers
# defined below so the math is identical across tiers.
# ============================================================================


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Leave-one-out NN distance via kd-tree. k=2 trick: col 1 is self (d=0),
# col 2 is the nearest *other* sim. O(N log N) memory; replaces an N x N
# as.matrix(dist()) that OOMs at N >= ~30K.
.ood.knn.loo <- function(scores) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("OOD diagnostics require the 'RANN' package. ",
         "Install with install.packages('RANN').")
  RANN::nn2(scores, scores, k = 2L)$nn.dists[, 2]
}

.ood.knn.obs <- function(scores, score_obs) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("OOD diagnostics require the 'RANN' package. ",
         "Install with install.packages('RANN').")
  RANN::nn2(scores, score_obs, k = 1L)$nn.dists[1, 1]
}

# k-NN average distance (LOO, per sim). Returns vector of length nrow(scores)
# containing the mean distance from each sim to its k nearest OTHER sims.
# Sun et al. 2022 (ICLR) "Out-of-Distribution Detection with Deep Nearest
# Neighbors" -- k-NN average is the current standard non-Gaussian OOD metric.
.ood.knn.k.loo <- function(scores, k = 10L) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("OOD diagnostics require the 'RANN' package.")
  k <- as.integer(k)
  k_eff <- min(k, nrow(scores) - 1L)
  if (k_eff < 1L) return(rep(NA_real_, nrow(scores)))
  nn <- RANN::nn2(scores, scores, k = k_eff + 1L)$nn.dists
  rowMeans(nn[, 2:(k_eff + 1L), drop = FALSE])
}

# k-NN average distance from obs to k nearest sims.
.ood.knn.k.obs <- function(scores, score_obs, k = 10L) {
  if (!requireNamespace("RANN", quietly = TRUE))
    stop("OOD diagnostics require the 'RANN' package.")
  k <- as.integer(k)
  k_eff <- min(k, nrow(scores))
  if (k_eff < 1L) return(NA_real_)
  nn <- RANN::nn2(scores, score_obs, k = k_eff)$nn.dists
  mean(nn[1, ])
}

# Pooled within-class covariance matrix (weighted average of per-class S_k).
# Matches the construction used in Lee et al. 2018 (NeurIPS) for Mahalanobis
# OOD detection in penultimate-layer feature space.
.ood.pooled.cov <- function(scores, class_id) {
  uniq <- sort(unique(class_id))
  d <- ncol(scores)
  total <- matrix(0, d, d)
  total_w <- 0
  for (k in uniq) {
    rows <- which(class_id == k)
    if (length(rows) < 2L) next
    S_k <- cov(scores[rows, , drop = FALSE])
    w   <- length(rows) - 1L
    total   <- total + w * S_k
    total_w <- total_w + w
  }
  if (total_w == 0) diag(1, d, d) else (total / total_w)
}

# Mahalanobis distance (LOO, per sim) from each sim to its OWN class mean,
# using the pooled within-class covariance. Equivalent in spirit to Lee
# et al. 2018 (NeurIPS) "A Simple Unified Framework for Detecting OOD
# Samples and Adversarial Attacks", but expressed as the LOO empirical
# null distribution for the per-class geometric verdict.
.ood.mahalanobis.loo <- function(scores, class_id) {
  uniq <- sort(unique(class_id))
  Sigma <- .ood.pooled.cov(scores, class_id)
  Sigma <- Sigma + diag(1e-6, ncol(Sigma))
  Sigma_inv <- tryCatch(solve(Sigma),
                        error = function(e) MASS::ginv(Sigma))
  d_loo <- rep(NA_real_, nrow(scores))
  for (k in uniq) {
    rows <- which(class_id == k)
    if (length(rows) < 2L) next
    mu_k <- colMeans(scores[rows, , drop = FALSE])
    centered <- sweep(scores[rows, , drop = FALSE], 2, mu_k, "-")
    d2 <- rowSums((centered %*% Sigma_inv) * centered)
    d_loo[rows] <- sqrt(pmax(d2, 0))
  }
  d_loo
}

# Mahalanobis distance from obs to each class mean, using pooled Sigma.
.ood.mahalanobis.obs.per.class <- function(scores, score_obs, class_id) {
  uniq <- sort(unique(class_id))
  Sigma <- .ood.pooled.cov(scores, class_id)
  Sigma <- Sigma + diag(1e-6, ncol(Sigma))
  Sigma_inv <- tryCatch(solve(Sigma),
                        error = function(e) MASS::ginv(Sigma))
  d_out <- setNames(rep(NA_real_, length(uniq)), as.character(uniq))
  for (i in seq_along(uniq)) {
    k <- uniq[i]
    rows <- which(class_id == k)
    if (length(rows) < 2L) next
    mu_k <- colMeans(scores[rows, , drop = FALSE])
    delta <- as.numeric(score_obs[1, ]) - mu_k
    d2 <- sum((delta %*% Sigma_inv) * delta)
    d_out[i] <- sqrt(max(d2, 0))
  }
  d_out
}

# Per-method geometric PASS / WARN / FAIL verdict.
#
# Pure descriptive geometry vs. the empirical sim-to-sim LOO-NN cloud:
#   PASS = obs is INSIDE  the cloud  (d_obs <= Q95(d_loo))
#   WARN = obs is AT  THE EDGE       (Q95 < d_obs <= max(d_loo))
#   FAIL = obs is OUTSIDE the cloud  (d_obs >  max(d_loo))
#
# No combination with alpha, no interpretation. The user reads it and
# decides what it means for their problem.
.ood.geometric.verdict <- function(d_obs, d_loo, edge.quantile = 0.95) {
  if (!is.finite(d_obs) || length(d_loo) < 2 || !any(is.finite(d_loo)))
    return(NA_character_)
  d_loo <- d_loo[is.finite(d_loo)]
  q_edge <- as.numeric(quantile(d_loo, probs = edge.quantile, na.rm = TRUE))
  d_max  <- max(d_loo)
  if (d_obs <= q_edge) "PASS"
  else if (d_obs <= d_max) "WARN"
  else "FAIL"
}

# Per-stat support classification + mid-rank percentiles.
# Outlier = obs strictly outside sim range (model cannot produce obs).
# Uninformative = sim is constant AND obs matches that constant (no info).
.ood.per.stat <- function(S_sim, observed, sim_cols) {
  n_stats <- ncol(S_sim)
  sim_var   <- apply(S_sim, 2, var)
  sim_range <- apply(S_sim, 2, function(x) diff(range(x)))
  col_sd    <- apply(S_sim, 2, sd, na.rm = TRUE)

  percentiles <- vapply(seq_len(n_stats), function(k) {
    mean(S_sim[, k] < observed[k]) + 0.5 * mean(S_sim[, k] == observed[k])
  }, numeric(1))
  names(percentiles) <- sim_cols

  sim_min <- apply(S_sim, 2, min)
  sim_max <- apply(S_sim, 2, max)
  obs_in_range  <- observed >= sim_min & observed <= sim_max
  sim_constant  <- sim_var < .Machine$double.eps |
                   sim_range < .Machine$double.eps
  uninformative <- sim_constant & obs_in_range
  outlier       <- !obs_in_range

  reason <- rep("ok", n_stats)
  reason[uninformative] <- "uninformative"
  reason[outlier]       <- "outlier"

  n_outliers      <- sum(outlier)
  n_uninformative <- sum(uninformative)
  n_informative   <- n_stats - n_uninformative
  outlier_frac    <- if (n_informative > 0) n_outliers / n_informative else 0

  tier <- if (outlier_frac < 0.10) "pass"
          else if (outlier_frac < 0.25) "warn"
          else "fail"

  list(
    percentiles = data.frame(
      stat = sim_cols, percentile = round(percentiles, 4),
      outlier = outlier, reason = reason, stringsAsFactors = FALSE
    ),
    summary = list(
      n_outliers = n_outliers, n_uninformative = n_uninformative,
      n_stats = n_stats, outlier_frac = outlier_frac, tier = tier,
      pass = tier != "fail"
    ),
    outlier_mask = outlier,
    uninformative_mask = uninformative,
    keep_mask_filtered = col_sd > 0 & is.finite(col_sd) &
                          !outlier & !uninformative,
    keep_mask_all = col_sd > 0 & is.finite(col_sd),
    col_sd = col_sd
  )
}

# Sign-safe log + z-score for params (used by PLS).
.ood.normalize.params <- function(P_sim) {
  Pl <- sign(P_sim) * log1p(abs(P_sim))
  pl_mu <- colMeans(Pl)
  pl_sd <- apply(Pl, 2, sd)
  pl_sd[pl_sd == 0 | !is.finite(pl_sd)] <- 1
  sweep(sweep(Pl, 2, pl_mu, "-"), 2, pl_sd, "/")
}

# NN density in PCA space. Used for both "all" and "filtered" runs.
.ood.pca.nn <- function(S_sim, observed, keep_mask, pca.var, alpha) {
  if (sum(keep_mask) < 2) {
    return(list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                null_distribution = numeric(0), n_pcs = NA_integer_,
                var_explained = NA_real_, n_stats = sum(keep_mask),
                scores = NULL, score_obs = NULL,
                pca_fit = NULL, keep_mask = keep_mask))
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

  d_obs <- .ood.knn.obs(scores, score_obs)
  d_loo <- .ood.knn.loo(scores)
  p <- mean(d_loo >= d_obs)
  verdict <- .ood.geometric.verdict(d_obs, d_loo)

  list(d_obs = d_obs, p_value = p, pass = p >= alpha,
       verdict = verdict,
       d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
       d_loo_max = max(d_loo),
       null_distribution = d_loo, n_pcs = n_pcs,
       var_explained = var_cum[n_pcs], n_stats = ncol(S),
       scores = scores, score_obs = score_obs,
       pca_fit = pca_fit, keep_mask = keep_mask)
}

# NN density in PLS space.
.ood.pls.nn <- function(S_sim, P_sim_norm, observed, keep_mask,
                        pls.n.comp, alpha) {
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
  n_comp_eff <- pls_fit$n.comp
  scores     <- pls_fit$scores
  score_obs  <- pls.project(pls_fit,
                            matrix(o, nrow = 1,
                                   dimnames = list(NULL, colnames(S))))

  d_obs <- .ood.knn.obs(scores, score_obs)
  d_loo <- .ood.knn.loo(scores)
  p <- mean(d_loo >= d_obs)
  verdict <- .ood.geometric.verdict(d_obs, d_loo)

  list(d_obs = d_obs, p_value = p, pass = p >= alpha,
       verdict = verdict,
       d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
       d_loo_max = max(d_loo),
       null_distribution = d_loo, n_comp = n_comp_eff,
       n_stats = ncol(S),
       scores = scores, score_obs = score_obs,
       pls_fit = pls_fit, keep_mask = keep_mask)
}

# Project obs (with outlier stats replaced by sim means) into the all-stats
# basis so the panel-1 arrow shows where obs would sit if outliers were
# neutralized.
.ood.project.outlier.neutralized <- function(observed, S_sim, outlier_mask,
                                             keep_mask_all, basis_fit,
                                             type = c("pca", "pls"), n_pcs = NULL) {
  type <- match.arg(type)
  obs_proj <- observed
  obs_proj[outlier_mask] <- colMeans(S_sim)[outlier_mask]
  obs_proj_mat <- matrix(obs_proj[keep_mask_all], nrow = 1,
                         dimnames = list(NULL,
                           colnames(S_sim[, keep_mask_all, drop = FALSE])))
  if (type == "pca") {
    predict(basis_fit, obs_proj_mat)[, seq_len(n_pcs), drop = FALSE]
  } else {
    pls.project(basis_fit, obs_proj_mat)
  }
}

# Resolve observed to a numeric vector matching stat_cols.
.ood.prepare.observed <- function(observed, stat_cols) {
  if (is.data.frame(observed) || is.matrix(observed)) {
    if (all(stat_cols %in% colnames(observed)))
      as.numeric(observed[1, stat_cols])
    else
      as.numeric(observed[1, ])
  } else {
    as.numeric(observed)
  }
}

# Build the sim-stat matrix from reftable, drop rows with non-finite stats
# (and params, when provided). Returns S_sim, P_sim (or NULL), sim_cols.
.ood.prepare.sims <- function(reftable, stat_cols, param.cols = NULL) {
  sim_cols <- intersect(stat_cols, colnames(reftable))
  if (length(sim_cols) == 0)
    stop("No stat columns found in reftable. Check stat.cols or observed column names.")
  S_sim <- as.matrix(reftable[, sim_cols, drop = FALSE])

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

  list(S_sim = S_sim, P_sim = P_sim, sim_cols = sim_cols, n_dropped = sum(!good))
}

# Verdict: tier-based on per-stat outlier fraction; secondary checks (joint
# NN density + ensemble disagreement) can downgrade pass→warn or warn→fail
# but never upgrade.
.ood.combine.verdict <- function(per_stat_tier, secondary_ok) {
  if (per_stat_tier == "fail") return("fail")
  if (per_stat_tier == "warn") return(if (secondary_ok) "warn" else "fail")
  if (secondary_ok) "pass" else "warn"
}


# ---------------------------------------------------------------------------
# Public: OOD.pretrain() — prior-predictive coverage check
# ---------------------------------------------------------------------------

#' Pre-training Out-of-Distribution Diagnostic (Prior-Predictive Coverage)
#'
#' Evaluates whether the observed summary statistics are reachable under the
#' prior predictive distribution, before any model is trained. Operates on
#' the reftable (prior-predictive sims) and observed stats alone — no neural
#' network required.
#'
#' Three complementary checks:
#' \itemize{
#'   \item \strong{Per-stat support} (marginal): each stat is classified as
#'     \code{ok} (obs inside sim range), \code{outlier} (obs strictly outside
#'     sim range — model cannot produce this observation under the prior),
#'     or \code{uninformative} (sims constant + obs matches constant).
#'     Outlier fraction drives the per-stat tier:
#'     \code{pass} (< 10\%), \code{warn} (10-25\%), \code{fail} (> 25\%).
#'   \item \strong{NN density in PCA space} (joint linear, variance-aligned),
#'     all-stats and filtered (zero-variance + per-stat outliers removed).
#'     Observed nearest-neighbor distance to the sim cloud is compared to the
#'     empirical leave-one-out sim-to-sim NN distribution (one-sided p-value);
#'     distribution-free, no Gaussian assumption.
#'   \item \strong{NN density in PLS space} (joint linear, param-aligned),
#'     all-stats and filtered. PLS components maximize covariance between
#'     stats and params, aligning with the directions ABC/NN use for
#'     inference. Computed only when \code{param.cols} is provided.
#' }
#'
#' Use \code{OOD.projection.diagnose()} to drill into why obs lands where it
#' does in PCA/PLS space, and \code{OOD.outliers()} for per-outlier sim
#' distributions.
#'
#' @param observed numeric vector or 1-row data.frame — observed summary
#'   statistics. If a data.frame with named columns, those names are matched
#'   against the reftable.
#' @param reftable data.frame — the prior-predictive reference table.
#' @param stat.cols character vector — stat columns in the reftable.
#'   Required unless \code{observed} is a data.frame whose column names
#'   identify stats.
#' @param param.cols character vector or NULL — parameter columns in the
#'   reftable. When provided, the PLS NN density check is computed.
#' @param pca.var numeric — cumulative variance threshold for PCA (default 0.95).
#' @param pls.n.comp integer — number of PLS components to fit (default 10,
#'   capped at \code{n_params + 5} and \code{n_stats - 1}).
#' @param alpha numeric — significance level for the NN density tests
#'   (default 0.01). Per-stat outliers use a strict support-based criterion
#'   (obs outside sim range) and do not depend on alpha.
#' @param plot logical — produce diagnostic plots (default TRUE).
#'
#' @return A list of class \code{"OOD_pretrain"} with:
#' \describe{
#'   \item{per_stat}{list with n_outliers, n_uninformative, n_stats,
#'     outlier_frac, tier, pass}
#'   \item{percentiles}{data.frame with stat, percentile, outlier, reason}
#'   \item{pca_all, pca_filtered}{NN density in PCA space (all stats, and
#'     stats minus zero-var + per-stat outliers)}
#'   \item{pls_all, pls_filtered}{NN density in PLS space — present when
#'     \code{param.cols} is provided}
#'   \item{overall}{character: \code{"pass"}, \code{"warn"}, or \code{"fail"}}
#' }
#'
#' @section Plot layout:
#' 2x3 panel layout. Density panels prefer PLS over PCA when params are
#' available. Panels: NN-vs-component scatter, NN density (primary),
#' percentile meta-histogram, outlier bars, NN density (filtered),
#' empty (reserved for ensemble disagreement in \code{OOD.posttrain()}).
#'
#' @seealso \code{\link{OOD.posttrain}} for post-training checks (NN-latent
#'   density, ensemble disagreement). \code{\link{OOD.projection.diagnose}}
#'   for projection forensics. \code{\link{OOD.outliers}} for per-outlier
#'   sim distributions.
#'
#' @export
OOD.pretrain <- function(observed, reftable, stat.cols = NULL,
                         param.cols = NULL,
                         pca.var = 0.95,
                         pls.n.comp = 10L,
                         alpha = 0.01,
                         plot = TRUE,
                         verbose = TRUE) {

  # --- Resolve stat columns ---
  if (!is.null(stat.cols)) {
    stat_cols <- stat.cols
  } else if (is.data.frame(observed)) {
    stat_cols <- colnames(observed)
  } else {
    stop("stat.cols must be provided when observed is a numeric vector.")
  }

  # --- Prepare observed and sim matrices ---
  observed <- .ood.prepare.observed(observed, stat_cols)
  prep <- .ood.prepare.sims(reftable, stat_cols, param.cols)
  S_sim <- prep$S_sim; P_sim <- prep$P_sim; sim_cols <- prep$sim_cols
  n_stats <- ncol(S_sim)

  if (length(observed) != n_stats)
    stop(sprintf("observed has %d elements but reftable has %d stat columns.",
                 length(observed), n_stats))

  results <- list()

  # --- Per-stat support ---
  per <- .ood.per.stat(S_sim, observed, sim_cols)
  results$percentiles <- per$percentiles
  results$per_stat    <- per$summary

  # --- PCA NN density (all + filtered) ---
  results$pca_all      <- .ood.pca.nn(S_sim, observed, per$keep_mask_all,
                                       pca.var, alpha)
  results$pca_filtered <- .ood.pca.nn(S_sim, observed, per$keep_mask_filtered,
                                       pca.var, alpha)
  if (sum(per$keep_mask_filtered) < 2)
    warning("Fewer than 2 informative non-outlier stats; filtered NN density check skipped.")

  # Project obs (outliers neutralized) into all-stats PCA basis
  if (!is.null(results$pca_all$scores) && sum(per$outlier_mask) > 0) {
    results$pca_all$score_obs_projected <-
      .ood.project.outlier.neutralized(observed, S_sim, per$outlier_mask,
                                        per$keep_mask_all,
                                        results$pca_all$pca_fit, "pca",
                                        n_pcs = results$pca_all$n_pcs)
  }

  # --- PLS NN density (all + filtered), if params available ---
  if (!is.null(P_sim)) {
    P_sim_norm <- .ood.normalize.params(P_sim)
    results$pls_all      <- .ood.pls.nn(S_sim, P_sim_norm, observed,
                                         per$keep_mask_all, pls.n.comp, alpha)
    results$pls_filtered <- .ood.pls.nn(S_sim, P_sim_norm, observed,
                                         per$keep_mask_filtered, pls.n.comp, alpha)
    if (!is.null(results$pls_all$scores) && sum(per$outlier_mask) > 0) {
      results$pls_all$score_obs_projected <-
        .ood.project.outlier.neutralized(observed, S_sim, per$outlier_mask,
                                          per$keep_mask_all,
                                          results$pls_all$pls_fit, "pls")
    }
  } else {
    results$pls_all <- NULL
    results$pls_filtered <- NULL
  }

  # --- Overall verdict ---
  density_pass <- if (!is.null(results$pls_filtered) &&
                      !is.na(results$pls_filtered$pass)) {
    isTRUE(results$pls_filtered$pass)
  } else {
    isTRUE(results$pca_filtered$pass)
  }
  results$overall <- .ood.combine.verdict(per$summary$tier, density_pass)

  # Stash data for downstream (printing / plotting / posttrain reuse /
  # prior-suggestion forensics). P_sim + param.cols carried only when
  # params were supplied (PLS fired).
  results$.context <- list(
    S_sim = S_sim, observed = observed, sim_cols = sim_cols,
    n_stats = n_stats, alpha = alpha,
    keep_mask_all = per$keep_mask_all,
    keep_mask_filtered = per$keep_mask_filtered,
    outlier_mask = per$outlier_mask,
    uninformative_mask = per$uninformative_mask,
    P_sim = P_sim,
    param.cols = if (!is.null(P_sim)) param.cols else NULL
  )

  class(results) <- c("OOD_pretrain", "OOD_diagnostic")

  if (verbose) .ood.print.summary(results, model_type = "pretrain")
  if (plot) .ood.plot.pretrain(results)

  invisible(results)
}


# ---------------------------------------------------------------------------
# Public: OOD.posttrain() — model-fit + posterior fidelity
# ---------------------------------------------------------------------------

#' Post-training Out-of-Distribution Diagnostic (Model Fit & Posterior Fidelity)
#'
#' Evaluates whether a trained neural network's representation and ensemble
#' disagreement at the observed point indicate reliable inference. Operates
#' on the trained NN plus the reftable used for training. Optionally consumes
#' an \code{OOD.pretrain()} result so the prior-predictive checks (per-stat,
#' PCA, PLS) are not recomputed and the overall verdict integrates both tiers.
#'
#' Two complementary checks (in addition to anything inherited from
#' \code{pretrain}):
#' \itemize{
#'   \item \strong{NN-latent NN density}: penultimate-layer activations of
#'     the top-ranked ensemble model, z-scored per dim. The nonlinear joint
#'     manifold the NN actually uses for regression. Catches "curved" holes
#'     that linear PCA/PLS miss. Torch inverse models only.
#'   \item \strong{Ensemble disagreement}: per-param coefficient of variation
#'     across the ensemble at obs, compared to the empirical distribution of
#'     per-sim mean CVs (torch inverse models). High disagreement at obs
#'     relative to sims is a model-relative epistemic-uncertainty signal
#'     that the trained model is extrapolating at obs.
#' }
#'
#' @param trained.nn list — output from \code{train.emulator()} or
#'   \code{tune.nn()}.
#' @param observed numeric vector or 1-row data.frame — observed summary
#'   statistics.
#' @param reftable data.frame — the reference table used for training.
#' @param pretrain object of class \code{"OOD_pretrain"} or NULL — when
#'   provided, the per-stat/PCA/PLS checks are reused (faster + verdict
#'   integrates both tiers). When NULL, posttrain runs the pre-training
#'   checks itself for verdict context.
#' @param theta numeric vector or NULL — optimized parameter estimate (from
#'   \code{emulator.optimize()}). Required for ensemble disagreement on
#'   forward emulator models; ignored for inverse \code{tune.nn()} models.
#' @param pca.var numeric — PCA cumulative variance threshold (default 0.95).
#'   Used when \code{pretrain = NULL}.
#' @param pls.n.comp integer — number of PLS components (default 10L). Used
#'   when \code{pretrain = NULL}.
#' @param alpha numeric — significance level for NN density tests (default 0.01).
#' @param plot logical — produce diagnostic plots (default TRUE).
#'
#' @return A list of class \code{"OOD_posttrain"} with all fields from
#'   \code{OOD.pretrain()} plus:
#' \describe{
#'   \item{nn_latent}{NN density in the penultimate-layer representation
#'     of the top-ranked ensemble model. Torch inverse models only.}
#'   \item{model_disagreement}{per-param CV at obs and (torch ensembles)
#'     the empirical null of per-sim mean CV with one-sided p-value}
#'   \item{pretrain}{the \code{OOD_pretrain} result used for context
#'     (either passed in or computed internally)}
#'   \item{overall}{character: \code{"pass"}, \code{"warn"}, or \code{"fail"}}
#' }
#'
#' @section Verdict logic:
#' Per-stat tier (from pretrain) is primary. Secondary checks — joint NN
#' density (PLS-filtered or NN-latent, whichever is richer) and ensemble
#' disagreement — can downgrade pass→warn or warn→fail but never upgrade.
#'
#' @seealso \code{\link{OOD.pretrain}}, \code{\link{OOD.projection.diagnose}}.
#'
#' @export
OOD.posttrain <- function(trained.nn, observed, reftable,
                          pretrain = NULL,
                          theta = NULL,
                          pca.var = 0.95,
                          pls.n.comp = 10L,
                          alpha = 0.01,
                          plot = TRUE) {

  if (is.null(trained.nn))
    stop("trained.nn is required. Use OOD.pretrain() for pre-training checks.")

  model_type <- trained.nn$type
  if (is.null(model_type)) model_type <- "emulator"
  is_emulator <- model_type == "emulator"

  stat_cols  <- trained.nn$data$stat_cols
  param.cols <- trained.nn$param_cols

  observed <- .ood.prepare.observed(observed, stat_cols)
  obs_raw  <- observed

  # --- Pretrain context: reuse if supplied, else recompute ---
  if (!is.null(pretrain)) {
    if (!inherits(pretrain, "OOD_pretrain"))
      stop("pretrain must be the output of OOD.pretrain().")
    pre <- pretrain
    # Validate stat compatibility
    if (!setequal(pre$.context$sim_cols, intersect(stat_cols, colnames(reftable))))
      warning("pretrain stat columns differ from trained.nn stat columns — verdict may be inconsistent.")
  } else {
    pre <- OOD.pretrain(observed = observed, reftable = reftable,
                        stat.cols = stat_cols, param.cols = param.cols,
                        pca.var = pca.var, pls.n.comp = pls.n.comp,
                        alpha = alpha, plot = FALSE, verbose = FALSE)
  }

  # Inherit pretrain results into the posttrain return
  results <- pre[setdiff(names(pre), c("overall"))]
  S_sim       <- pre$.context$S_sim
  sim_cols    <- pre$.context$sim_cols
  outlier     <- pre$.context$outlier_mask
  keep_filt   <- pre$.context$keep_mask_filtered
  n_stats     <- pre$.context$n_stats
  n_kept      <- sum(keep_filt)

  # --- NN-latent density (torch inverse models only) ---
  results$nn_latent <- NULL
  if (!is_emulator && n_kept >= 2) {
    models <- trained.nn$models
    m1 <- if (!is.null(models) && length(models) >= 1) models[[1]] else NULL
    if (!is.null(m1) && inherits(m1, "nn_module")) {
      results$nn_latent <- tryCatch(
        .ood.nn.latent(m1, trained.nn, S_sim, observed, obs_raw, outlier,
                       model_type, alpha),
        error = function(e) {
          warning(sprintf("NN-latent check skipped: %s", conditionMessage(e)))
          NULL
        })
    }
  }

  # --- Ensemble disagreement ---
  results$model_disagreement <- .ood.ensemble.disagreement(
    trained.nn, S_sim, observed, obs_raw, sim_cols,
    model_type, is_emulator, theta, alpha)

  # --- Combined verdict ---
  if (!is.null(results$nn_latent) && !is.na(results$nn_latent$pass)) {
    density_pass <- isTRUE(results$nn_latent$pass)
  } else if (!is.null(results$pls_filtered) &&
             !is.na(results$pls_filtered$pass)) {
    density_pass <- isTRUE(results$pls_filtered$pass)
  } else {
    density_pass <- isTRUE(results$pca_filtered$pass)
  }
  md_pass <- if (!is.null(results$model_disagreement))
               results$model_disagreement$reliable else TRUE
  secondary_ok <- density_pass && md_pass

  results$overall <- .ood.combine.verdict(pre$per_stat$tier, secondary_ok)
  results$pretrain <- pre

  class(results) <- c("OOD_posttrain", "OOD_diagnostic")

  .ood.print.summary(results, model_type = model_type)
  if (plot) .ood.plot.posttrain(results)

  invisible(results)
}


# ---------------------------------------------------------------------------
# Internal: NN-latent density (post-training only)
# ---------------------------------------------------------------------------

.ood.nn.latent <- function(m1, trained.nn, S_sim, observed, obs_raw,
                           outlier_mask, model_type, alpha) {
  X_sim <- .prep.sims(S_sim, trained.nn$data, model_type, trained.nn$sfs.dims)
  X_obs <- .prep.observed(obs_raw, trained.nn$data, model_type,
                          trained.nn$sfs.dims)

  Z_sim <- .torch.penultimate(m1, X_sim)
  Z_obs <- .torch.penultimate(m1, X_obs)

  z_mu <- colMeans(Z_sim)
  z_sd <- apply(Z_sim, 2, sd)
  z_sd[z_sd == 0 | !is.finite(z_sd)] <- 1
  Z_sim_n <- sweep(sweep(Z_sim, 2, z_mu, "-"), 2, z_sd, "/")
  Z_obs_n <- sweep(sweep(Z_obs, 2, z_mu, "-"), 2, z_sd, "/")

  d_obs <- .ood.knn.obs(Z_sim_n, Z_obs_n)
  d_loo <- .ood.knn.loo(Z_sim_n)
  p <- mean(d_loo >= d_obs)
  verdict <- .ood.geometric.verdict(d_obs, d_loo)

  # Project obs (outliers neutralized) through the NN for honest panel-1 Y
  Z_obs_proj_n <- NULL
  if (sum(outlier_mask) > 0) {
    obs_proj <- observed
    obs_proj[outlier_mask] <- colMeans(S_sim)[outlier_mask]
    X_obs_proj <- .prep.observed(obs_proj, trained.nn$data, model_type,
                                  trained.nn$sfs.dims)
    Z_obs_proj <- .torch.penultimate(m1, X_obs_proj)
    Z_obs_proj_n <- sweep(sweep(Z_obs_proj, 2, z_mu, "-"), 2, z_sd, "/")
  }

  list(d_obs = d_obs, p_value = p, pass = p >= alpha,
       verdict = verdict,
       d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
       d_loo_max = max(d_loo),
       null_distribution = d_loo, n_dims = ncol(Z_sim),
       n_stats = ncol(S_sim),
       scores = Z_sim_n, score_obs = Z_obs_n,
       score_obs_projected = Z_obs_proj_n)
}


# ---------------------------------------------------------------------------
# Internal: ensemble disagreement (post-training only)
# ---------------------------------------------------------------------------

.ood.ensemble.disagreement <- function(trained.nn, S_sim, observed, obs_raw,
                                       sim_cols, model_type, is_emulator,
                                       theta, alpha) {
  models <- trained.nn$models
  if (is.null(models) || length(models) <= 1) return(NULL)
  if (is_emulator && is.null(theta)) return(NULL)

  feat_mu   <- trained.nn$data$feat_mu
  feat_sd   <- trained.nn$data$feat_sd
  target_mu <- trained.nn$data$target_mu
  target_sd <- trained.nn$data$target_sd

  is_torch_ensemble <- all(vapply(models,
                                  function(m) inherits(m, "nn_module"),
                                  logical(1)))

  # Predictions at obs
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
    X_input <- .prep.observed(obs_raw, trained.nn$data, model_type,
                              trained.nn$sfs.dims)
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

  # Empirical null (torch inverse only)
  null_mean_cv <- NULL
  p_cv <- NA_real_
  pass_cv <- NA
  if (is_torch_ensemble && !is_emulator) {
    null_mean_cv <- tryCatch({
      X_sim <- .prep.sims(S_sim, trained.nn$data, model_type,
                          trained.nn$sfs.dims)
      preds_sim <- vapply(models, function(m) {
        p_z <- .torch.predict(m, X_sim)
        as.matrix(.inv.transform(p_z, target_mu, target_sd))
      }, FUN.VALUE = matrix(0, nrow(S_sim), length(target_mu)))
      pm <- apply(preds_sim, c(1, 2), mean)
      ps <- apply(preds_sim, c(1, 2), sd)
      cv_sim <- ifelse(abs(pm) > 1e-12, ps / abs(pm), ps)
      rowMeans(cv_sim)
    }, error = function(e) {
      warning(sprintf("Ensemble null skipped: %s", conditionMessage(e)))
      NULL
    })
    if (!is.null(null_mean_cv) && length(null_mean_cv) > 0) {
      p_cv <- mean(null_mean_cv >= mean_cv)
      pass_cv <- p_cv >= alpha
    }
  }

  reliable <- if (!is.na(pass_cv)) isTRUE(pass_cv) else (mean_cv < 0.1)

  list(mean_cv = mean_cv, per_stat_cv = cv,
       null_mean_cv = null_mean_cv, p_value = p_cv,
       reliable = reliable)
}


# ---------------------------------------------------------------------------
# Internal: stdout summary printer
# ---------------------------------------------------------------------------

.ood.print.summary <- function(results, model_type) {
  ps  <- results$per_stat
  pct <- results$percentiles

  cat(sprintf("PipeMaster:: OOD.%s (%s)\n",
              if (inherits(results, "OOD_posttrain")) "posttrain" else "pretrain",
              model_type))
  cat(sprintf("  Stats: %d total  =  %d informative (%d ok + %d outlier)  +  %d uninformative\n",
              ps$n_stats,
              ps$n_stats - ps$n_uninformative,
              ps$n_stats - ps$n_uninformative - ps$n_outliers,
              ps$n_outliers, ps$n_uninformative))
  cat(sprintf("  Per-stat tier: outlier_frac = %d/%d = %.1f%% [%s]\n",
              ps$n_outliers, ps$n_stats - ps$n_uninformative,
              100 * ps$outlier_frac, toupper(ps$tier)))
  if (ps$n_outliers > 0 && ps$n_outliers <= 20) {
    bad <- pct[pct$reason == "outlier", ]
    for (r in seq_len(nrow(bad)))
      cat(sprintf("    %s: percentile=%.3f\n", bad$stat[r], bad$percentile[r]))
  } else if (ps$n_outliers > 20) {
    cat(sprintf("    (%d outliers — inspect with OOD.outliers())\n", ps$n_outliers))
  }

  .one_density <- function(label, x, units = "PCs") {
    if (is.null(x) || is.na(x$pass)) return()
    n_axis <- if (!is.null(x$n_pcs)) x$n_pcs
              else if (!is.null(x$n_comp)) x$n_comp
              else x$n_dims
    if (!is.null(x$var_explained)) {
      cat(sprintf("  NN density (%s, %3d stats, %d %s, %.0f%% var): d_obs=%.3f, p=%.4f [%s]\n",
                  label, x$n_stats, n_axis, units,
                  x$var_explained * 100, x$d_obs, x$p_value,
                  ifelse(x$pass, "PASS", "FAIL")))
    } else {
      cat(sprintf("  NN density (%s, %3d stats, %d %s): d_obs=%.3f, p=%.4f [%s]\n",
                  label, x$n_stats, n_axis, units,
                  x$d_obs, x$p_value, ifelse(x$pass, "PASS", "FAIL")))
    }
  }
  .one_density("PCA all     ", results$pca_all,      "PCs")
  .one_density("PCA filtered", results$pca_filtered, "PCs")
  .one_density("PLS all     ", results$pls_all,      "comps")
  .one_density("PLS filtered", results$pls_filtered, "comps")
  .one_density("NN-latent   ", results$nn_latent,    "dims")

  if (!is.null(results$model_disagreement)) {
    md <- results$model_disagreement
    if (!is.null(md$p_value) && !is.na(md$p_value)) {
      cat(sprintf("  Ensemble disagreement: mean CV=%.3f, p=%.4f vs sim null (%d sims) [%s]\n",
                  md$mean_cv, md$p_value, length(md$null_mean_cv),
                  ifelse(md$reliable, "RELIABLE", "UNCERTAIN")))
    } else {
      cat(sprintf("  Ensemble disagreement: mean CV=%.3f (no sim null) [%s]\n",
                  md$mean_cv, ifelse(md$reliable, "RELIABLE", "UNCERTAIN")))
    }
  }
  cat(sprintf("  Overall: %s\n", toupper(results$overall)))
}


# ---------------------------------------------------------------------------
# Internal: plot helpers (shared by pretrain and posttrain)
# ---------------------------------------------------------------------------

# Histogram of LOO sim NN distances with obs marked.
.ood.panel.nn.hist <- function(nn_res, label, method, n_axis) {
  if (is.null(nn_res) || is.na(nn_res$pass)) {
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

# NN-vs-component-1 scatter with obs and (optionally) projected obs.
.ood.panel.scatter <- function(density_all, density_method, density_axis_n,
                               density_x_label, n_stats) {
  if (is.null(density_all$scores) || is.na(density_all$pass)) {
    plot.new(); title(main = sprintf("NN vs %s1 — unavailable", density_method))
    return(invisible(NULL))
  }
  y_sim <- density_all$null_distribution
  y_obs <- density_all$d_obs
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
}

.ood.panel.percentile.hist <- function(percentiles, uninformative_mask,
                                       n_outliers, n_uninformative, tier) {
  informative_pct <- percentiles[!uninformative_mask]
  bin_width <- 0.05
  hist_breaks <- seq(0, 1, by = bin_width)
  h <- hist(informative_pct, breaks = hist_breaks, plot = FALSE)
  unif_count <- length(informative_pct) * bin_width
  plot(h, col = "steelblue", border = "white",
       main = sprintf("Percentile distribution (%d informative: %d ok + %d outlier; %d uninformative) [%s]",
                      length(informative_pct), length(informative_pct) - n_outliers,
                      n_outliers, n_uninformative, toupper(tier)),
       xlab = "obs mid-rank percentile in sim distribution", ylab = "n stats")
  abline(h = unif_count, col = "grey40", lty = 2)
  legend("top", legend = "uniform (perfect fit)",
         col = "grey40", lty = 2, cex = 0.7, bty = "n")
}

.ood.panel.outlier.bars <- function(percentiles_df, n_outliers) {
  if (n_outliers == 0) {
    plot.new(); title(main = "No outliers")
    return(invisible(NULL))
  }
  bad <- percentiles_df[percentiles_df$reason == "outlier", ]
  top_n <- min(25, nrow(bad))
  bad_top <- bad[seq_len(top_n), , drop = FALSE]
  bad_top <- bad_top[rev(seq_len(top_n)), , drop = FALSE]
  bar_dir_col <- ifelse(bad_top$percentile >= 0.5, "firebrick", "royalblue")
  bar_heights <- ifelse(bad_top$percentile >= 0.5, 1, -1)
  old_mar <- par("mar"); par(mar = c(4, 9, 3, 1))
  barplot(bar_heights, names.arg = bad_top$stat,
          horiz = TRUE, las = 1, cex.names = 0.6,
          col = bar_dir_col, border = NA,
          xlim = c(-1.1, 1.1),
          xlab = "direction (obs < min sim  |  obs > max sim)",
          main = sprintf("Outliers: %d of %d (red: obs>max sim, blue: obs<min sim)",
                         top_n, nrow(bad)))
  abline(v = 0, col = "grey60", lty = 1)
  par(mar = old_mar)
}

.ood.panel.ensemble <- function(md) {
  if (is.null(md)) {
    plot.new()
    return(invisible(NULL))
  }
  if (!is.null(md$null_mean_cv) && length(md$null_mean_cv) > 0) {
    verdict <- ifelse(md$reliable, "RELIABLE", "UNCERTAIN")
    hist(md$null_mean_cv, breaks = 50, col = "grey80", border = "white",
         probability = TRUE,
         main = sprintf("Ensemble CV (p=%.3f vs sim null) [%s]",
                        md$p_value, verdict),
         xlab = "mean CV across ensemble (per-sim)",
         xlim = range(c(md$null_mean_cv, md$mean_cv)))
    abline(v = md$mean_cv, col = "red", lwd = 2)
    legend("topright",
           legend = c("sim per-sim mean CV", "observed mean CV"),
           col = c("grey50", "red"), lwd = c(5, 2), cex = 0.7)
  } else {
    cv_vals <- md$per_stat_cv
    cols_cv <- ifelse(cv_vals > 0.1, "firebrick", "steelblue")
    barplot(cv_vals, col = cols_cv, border = NA, las = 2, cex.names = 0.5,
            main = "Ensemble disagreement (CV)", ylab = "CV")
    abline(h = 0.1, lty = 2, col = "red")
  }
}

# Pick the primary density basis: NN-latent > PLS > PCA. Returns a list with
# density_all, density_method, density_axis_n, density_x_label, plus the
# filtered counterpart (NN-latent has no filtered variant; falls back to PLS > PCA).
.ood.pick.basis <- function(results) {
  if (!is.null(results$nn_latent) && !is.na(results$nn_latent$pass)) {
    primary <- list(all = results$nn_latent, method = "NN-latent",
                    axis_n = results$nn_latent$n_dims,
                    x_label = "latent dim 1 (sim score)")
  } else if (!is.null(results$pls_all)) {
    primary <- list(all = results$pls_all, method = "PLS",
                    axis_n = results$pls_all$n_comp,
                    x_label = "PLS1 (sim score)")
  } else {
    primary <- list(all = results$pca_all, method = "PCA",
                    axis_n = results$pca_all$n_pcs,
                    x_label = "PC1 (sim score)")
  }
  if (!is.null(results$pls_filtered)) {
    filtered <- list(filtered = results$pls_filtered, method = "PLS",
                     axis_n = results$pls_filtered$n_comp)
  } else {
    filtered <- list(filtered = results$pca_filtered, method = "PCA",
                     axis_n = results$pca_filtered$n_pcs)
  }
  list(primary = primary, filtered = filtered)
}

# Per-tier verdict summary text panel for pretrain row 1 col 3
.ood.panel.verdict.summary <- function(results) {
  ps <- results$per_stat
  lines <- c()
  lines <- c(lines, sprintf("OVERALL:    %s", toupper(results$overall)))
  lines <- c(lines, "")
  lines <- c(lines, sprintf("Per-stat:   %s", toupper(ps$tier)))
  lines <- c(lines, sprintf("  %d outliers / %d informative (%.1f%%)",
                             ps$n_outliers,
                             ps$n_stats - ps$n_uninformative,
                             100 * ps$outlier_frac))
  lines <- c(lines, sprintf("  %d uninformative", ps$n_uninformative))
  lines <- c(lines, sprintf("  %d stats total", ps$n_stats))
  lines <- c(lines, "")
  fmt <- function(label, x) {
    if (is.null(x) || is.na(x$pass)) return(NULL)
    sprintf("%-12s d=%6.2f  p=%.4f  %s",
            label, x$d_obs, x$p_value,
            ifelse(x$pass, "[PASS]", "[FAIL]"))
  }
  lines <- c(lines, fmt("PCA all:",      results$pca_all))
  lines <- c(lines, fmt("PCA filtered:", results$pca_filtered))
  if (!is.null(results$pls_all)) {
    lines <- c(lines, fmt("PLS all:",      results$pls_all))
    lines <- c(lines, fmt("PLS filtered:", results$pls_filtered))
  } else {
    lines <- c(lines, "PLS:        not run")
    lines <- c(lines, "  (pass param.cols)")
  }
  lines <- c(lines, "")
  lines <- c(lines, "Tier rules:")
  lines <- c(lines, "  pass < 10%")
  lines <- c(lines, "  warn 10-25%")
  lines <- c(lines, "  fail > 25%")

  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title(main = "Verdict summary")
  text(0.02, 0.97, paste(lines, collapse = "\n"),
       adj = c(0, 1), family = "mono", cex = 0.85)
}

# Wrapper so each density-row can fail gracefully when its basis isn't available
.ood.plot.density.row <- function(label, scatter_basis, hist_all, hist_filtered) {
  if (is.null(scatter_basis$all) || is.null(scatter_basis$all$scores)) {
    plot.new(); title(main = sprintf("%s scatter — unavailable", label))
    plot.new(); title(main = sprintf("%s NN density (all) — unavailable", label))
    plot.new(); title(main = sprintf("%s NN density (filtered) — unavailable", label))
    return(invisible(NULL))
  }
  .ood.panel.scatter(scatter_basis$all, scatter_basis$method,
                     scatter_basis$axis_n, scatter_basis$x_label,
                     scatter_basis$all$n_stats)
  .ood.panel.nn.hist(hist_all, "all",
                     scatter_basis$method, scatter_basis$axis_n)
  .ood.panel.nn.hist(hist_filtered, "filtered",
                     scatter_basis$method, scatter_basis$axis_n)
}

# 3x3 plot for pretrain — rows: per-stat, PCA, PLS. Each row owns 3 panels.
.ood.plot.pretrain <- function(results) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  ctx <- results$.context
  ps  <- results$per_stat

  # Row 1: per-stat
  .ood.panel.percentile.hist(results$percentiles$percentile,
                              ctx$uninformative_mask,
                              ps$n_outliers, ps$n_uninformative, ps$tier)
  .ood.panel.outlier.bars(results$percentiles, ps$n_outliers)
  .ood.panel.verdict.summary(results)

  # Row 2: PCA
  pca_basis <- list(all = results$pca_all, method = "PCA",
                    axis_n = results$pca_all$n_pcs,
                    x_label = "PC1 (sim score)")
  .ood.plot.density.row("PCA", pca_basis,
                         results$pca_all, results$pca_filtered)

  # Row 3: PLS
  if (!is.null(results$pls_all)) {
    pls_basis <- list(all = results$pls_all, method = "PLS",
                      axis_n = results$pls_all$n_comp,
                      x_label = "PLS1 (sim score)")
    .ood.plot.density.row("PLS", pls_basis,
                           results$pls_all, results$pls_filtered)
  } else {
    plot.new(); title(main = "PLS scatter — pass param.cols to enable")
    plot.new(); title(main = "PLS NN density (all) — pass param.cols")
    plot.new(); title(main = "PLS NN density (filtered) — pass param.cols")
  }
}

# 2x3 plot for posttrain. Same layout; panel 6 shows ensemble disagreement.
.ood.plot.posttrain <- function(results) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

  basis <- .ood.pick.basis(results)
  ctx <- results$.context
  ps  <- results$per_stat

  .ood.panel.scatter(basis$primary$all, basis$primary$method,
                     basis$primary$axis_n, basis$primary$x_label, ctx$n_stats)
  .ood.panel.nn.hist(basis$primary$all, "all",
                     basis$primary$method, basis$primary$axis_n)
  .ood.panel.percentile.hist(results$percentiles$percentile,
                              ctx$uninformative_mask,
                              ps$n_outliers, ps$n_uninformative, ps$tier)
  .ood.panel.outlier.bars(results$percentiles, ps$n_outliers)
  .ood.panel.nn.hist(basis$filtered$filtered, "filtered",
                     basis$filtered$method, basis$filtered$axis_n)
  .ood.panel.ensemble(results$model_disagreement)
}


# ============================================================================
# OOD.pretrain.classify — Pre-training diagnostic for model selection
#
# Multi-model variant of OOD.pretrain(): given a stacked reftable labelled by
# `model_col` and an observed dataset, evaluates whether each candidate model's
# prior predictive distribution covers the observation, identifies stats that
# are universal outliers (rejected by every model — drop them from the
# classifier feature set), and visualizes the K classes in PCA + PLS-DA space.
#
# Output drives:
#   1. Feature filtering for tune.nn.classify() — drop universal-outlier stats
#   2. Visual diagnostic on whether obs lands inside any model's cloud
#   3. Per-model NN-density p-value: probability obs is plausible under model m
# ============================================================================

# --- Internal: per-model rejection table (K x n_stats logical) ---
.ood.classify.rejection <- function(S_sim, model_id_int, observed, stat_cols,
                                    class_names) {
  K <- length(class_names)
  n_stats <- length(stat_cols)
  rejection <- matrix(FALSE, nrow = K, ncol = n_stats,
                      dimnames = list(class_names, stat_cols))
  per_model_pct <- matrix(NA_real_, nrow = K, ncol = n_stats,
                          dimnames = list(class_names, stat_cols))
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) next
    S_k <- S_sim[rows, , drop = FALSE]
    sim_min <- apply(S_k, 2, min)
    sim_max <- apply(S_k, 2, max)
    rejection[k, ] <- observed < sim_min | observed > sim_max
    per_model_pct[k, ] <- vapply(seq_len(n_stats), function(j)
      mean(S_k[, j] < observed[j]) + 0.5 * mean(S_k[, j] == observed[j]),
      numeric(1))
  }
  list(rejection = rejection, percentiles = per_model_pct)
}

# --- Internal: combined PCA on stacked sims (z-scored) ---
.ood.classify.pca <- function(S_sim, observed, keep_mask, pca.var) {
  if (sum(keep_mask) < 2L)
    return(list(scores = NULL, score_obs = NULL, fit = NULL,
                n_pcs = NA_integer_, var_explained = NA_real_,
                n_stats = sum(keep_mask)))
  S <- S_sim[, keep_mask, drop = FALSE]
  obs_keep <- observed[keep_mask]
  fit <- prcomp(S, center = TRUE, scale. = TRUE)
  var_cum <- cumsum(fit$sdev^2 / sum(fit$sdev^2))
  n_pcs <- which(var_cum >= pca.var)[1]
  if (is.na(n_pcs)) n_pcs <- ncol(S)
  n_pcs <- max(n_pcs, 2L)
  scores <- fit$x[, seq_len(n_pcs), drop = FALSE]
  score_obs <- predict(fit,
    matrix(obs_keep, nrow = 1, dimnames = list(NULL, colnames(S))))[, seq_len(n_pcs),
                                                                    drop = FALSE]
  list(scores = scores, score_obs = score_obs, fit = fit,
       n_pcs = n_pcs, var_explained = var_cum[n_pcs], n_stats = ncol(S))
}

# --- Internal: PLS-DA via pls.fit() with one-hot dummies ---
.ood.classify.plsda <- function(S_sim, model_id_int, observed, keep_mask,
                                class_names, n.comp = NULL) {
  K <- length(class_names)
  if (sum(keep_mask) < 2L)
    return(list(scores = NULL, score_obs = NULL, fit = NULL,
                n_comp = NA_integer_, n_stats = sum(keep_mask)))
  S <- S_sim[, keep_mask, drop = FALSE]
  obs_keep <- observed[keep_mask]

  Y_dummy <- matrix(0, nrow = nrow(S), ncol = K,
                    dimnames = list(NULL, class_names))
  for (k in seq_len(K)) Y_dummy[model_id_int == k, k] <- 1
  # Center dummies so PLS-DA fits class indicators (not raw 0/1)
  Y_dummy <- sweep(Y_dummy, 2, colMeans(Y_dummy), "-")

  if (is.null(n.comp)) n.comp <- max(2L, K - 1L)
  n.comp <- min(n.comp, sum(keep_mask) - 1L)

  fit <- pls.fit(stats = S, params = Y_dummy, n.comp = n.comp,
                 scale = TRUE, max.rows = 50000L)
  # pls.fit may subsample rows for fitting; project ALL rows back so scores
  # aligns with model_id_int (length nrow(S_sim)).
  scores <- pls.project(fit, S)
  score_obs <- pls.project(fit, matrix(obs_keep, nrow = 1,
                                       dimnames = list(NULL, colnames(S))))
  list(scores = scores, score_obs = score_obs, fit = fit,
       n_comp = fit$n.comp, n_stats = ncol(S))
}

# --- Internal: per-model NN density. For each model k, distance from obs to
#     nearest sim of model k (in keep-mask stat space, z-scored), compared to
#     the LOO sim-to-sim NN distribution within model k. Returns d_obs and
#     one-sided p-value per model. ---
.ood.classify.per.model.nn <- function(S_sim, model_id_int, observed,
                                       keep_mask, class_names, alpha) {
  K <- length(class_names)
  if (sum(keep_mask) < 2L)
    return(setNames(replicate(K, list(d_obs = NA_real_, p_value = NA_real_,
                                       pass = NA, null_distribution = numeric(0)),
                              simplify = FALSE),
                    class_names))

  S <- S_sim[, keep_mask, drop = FALSE]
  o <- observed[keep_mask]
  # Z-score using global sim mu/sd to put models in the same metric
  mu <- colMeans(S); sd_ <- apply(S, 2, sd); sd_[sd_ == 0] <- 1
  Z   <- sweep(sweep(S, 2, mu, "-"), 2, sd_, "/")
  z_o <- (o - mu) / sd_

  out <- vector("list", K); names(out) <- class_names
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) {
      out[[k]] <- list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                       null_distribution = numeric(0))
      next
    }
    Z_k <- Z[rows, , drop = FALSE]
    z_o_mat <- matrix(z_o, nrow = 1L)
    d_obs <- .ood.knn.obs(Z_k, z_o_mat)
    d_loo <- .ood.knn.loo(Z_k)
    p <- mean(d_loo >= d_obs)
    verdict <- .ood.geometric.verdict(d_obs, d_loo)
    out[[k]] <- list(d_obs = d_obs, p_value = p, pass = p >= alpha,
                     verdict = verdict,
                     d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
                     d_loo_max = max(d_loo),
                     null_distribution = d_loo)
  }
  out
}

# Per-class NN density in PLS-DA space. Takes precomputed PLS-DA scores
# (from .ood.classify.plsda) + obs's projection, partitions by class,
# and computes the geometric verdict per class. Z-scoring is done in
# the global PLS-DA space so cross-class distances share a metric.
## Per-class Mahalanobis-to-centroid verdict in score space (PCA or PLS-DA).
## Complements the K=1 NN-density verdict — gives a centroid distance,
## which is more interpretable for "is obs at the cloud center". Reuses
## the standalone helpers .ood.mahalanobis.loo() and
## .ood.mahalanobis.obs.per.class() (Lee et al. 2018).
.ood.classify.per.model.mahalanobis <- function(scores, score_obs,
                                                model_id_int, class_names,
                                                alpha) {
  K <- length(class_names)
  if (is.null(scores) || nrow(scores) < 2L)
    return(setNames(replicate(K, list(d_obs = NA_real_, p_value = NA_real_,
                                       pass = NA, verdict = NA_character_,
                                       null_distribution = numeric(0)),
                              simplify = FALSE),
                    class_names))

  d_loo_all       <- .ood.mahalanobis.loo(scores, model_id_int)
  d_obs_per_class <- .ood.mahalanobis.obs.per.class(scores, score_obs,
                                                    model_id_int)

  out <- vector("list", K); names(out) <- class_names
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) {
      out[[k]] <- list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                       verdict = NA_character_,
                       null_distribution = numeric(0))
      next
    }
    d_loo <- d_loo_all[rows]
    d_loo <- d_loo[is.finite(d_loo)]
    d_obs <- d_obs_per_class[as.character(k)]
    if (is.na(d_obs)) d_obs <- d_obs_per_class[k]
    p <- mean(d_loo >= d_obs)
    verdict <- .ood.geometric.verdict(d_obs, d_loo)
    out[[k]] <- list(d_obs = as.numeric(d_obs), p_value = p, pass = p >= alpha,
                     verdict = verdict,
                     d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
                     d_loo_max = max(d_loo),
                     null_distribution = d_loo)
  }
  out
}

.ood.classify.per.model.plsda <- function(scores, score_obs,
                                          model_id_int, class_names, alpha) {
  K <- length(class_names)
  if (is.null(scores) || nrow(scores) < 2L)
    return(setNames(replicate(K, list(d_obs = NA_real_, p_value = NA_real_,
                                       pass = NA, verdict = NA_character_,
                                       null_distribution = numeric(0)),
                              simplify = FALSE),
                    class_names))
  mu <- colMeans(scores); sd_ <- apply(scores, 2, sd); sd_[sd_ == 0] <- 1
  Z   <- sweep(sweep(scores, 2, mu, "-"), 2, sd_, "/")
  z_o <- (as.numeric(score_obs[1, ]) - mu) / sd_

  out <- vector("list", K); names(out) <- class_names
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) {
      out[[k]] <- list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                       verdict = NA_character_,
                       null_distribution = numeric(0))
      next
    }
    Z_k <- Z[rows, , drop = FALSE]
    z_o_mat <- matrix(z_o, nrow = 1L)
    d_obs <- .ood.knn.obs(Z_k, z_o_mat)
    d_loo <- .ood.knn.loo(Z_k)
    p <- mean(d_loo >= d_obs)
    verdict <- .ood.geometric.verdict(d_obs, d_loo)
    out[[k]] <- list(d_obs = d_obs, p_value = p, pass = p >= alpha,
                     verdict = verdict,
                     d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
                     d_loo_max = max(d_loo),
                     null_distribution = d_loo)
  }
  out
}

#' Pre-training Diagnostic for Model Selection (Multi-Model OOD)
#'
#' Runs prior-predictive coverage checks across multiple candidate demographic
#' models simultaneously. For each model independently, identifies stats whose
#' observed value falls outside that model's sim support (per-model outliers).
#' Stats rejected by *every* candidate model (universal outliers) cannot help
#' the classifier and are flagged for removal from the feature set.
#'
#' Visualizes the four (or K) candidate models in a single PCA and a single
#' PLS-DA projection, both colored by model_id with the observation marked.
#'
#' @param reftable data.frame — combined reference table with sims from all
#'   candidate models stacked, and a label column (`model_col`) identifying
#'   which model each sim came from.
#' @param model_col character — name of the label column.
#' @param observed numeric vector or 1-row data.frame — observed summary stats.
#' @param stat.cols character vector or NULL — feature columns. If NULL, all
#'   non-param non-label columns are used.
#' @param subsample integer or NULL — per-model subsample for PCA/PLS-DA
#'   (default NULL = use all sims). Useful for very large reftables.
#' @param alpha numeric — significance level for per-model NN density p-values
#'   (default 0.01).
#' @param universal_fail_threshold numeric — fraction of universal outliers
#'   above which the diagnostic fails (default 0.25).
#' @param pca.var numeric — PCA cumulative variance threshold (default 0.95).
#' @param pls.n.comp integer or NULL — PLS-DA components (default K-1 or 2,
#'   whichever is larger).
#' @param seed integer — random seed for subsampling.
#' @param plot logical — produce diagnostic plots (default TRUE).
#' @param verbose logical — print summary (default TRUE).
#'
#' @return A list of class `c("OOD_pretrain_classify", "OOD_diagnostic")`
#'   with:
#' \describe{
#'   \item{class_names}{character vector of K model labels}
#'   \item{n_classes}{integer K}
#'   \item{n_stats}{integer total stats considered}
#'   \item{rejection}{K x n_stats logical matrix — TRUE if obs outside model k's support for stat j}
#'   \item{rejection_count}{integer vector length n_stats — number of models rejecting each stat}
#'   \item{universal_mask}{logical mask — stats rejected by all K models}
#'   \item{specific_mask}{logical mask — stats rejected by 1..K-1 models}
#'   \item{keep_mask}{logical mask — feature-set filter for tune.nn.classify (= !universal)}
#'   \item{stat_cols}{character — original stat column names}
#'   \item{percentiles_per_model}{K x n_stats matrix of mid-rank percentiles}
#'   \item{pca}{combined PCA: scores, score_obs, fit, n_pcs, var_explained}
#'   \item{pls_da}{PLS-DA: scores, score_obs, fit, n_comp}
#'   \item{per_model_nn}{list of K, each with d_obs / p_value / pass / null_distribution}
#'   \item{universal_frac}{fraction of stats rejected by all models}
#'   \item{verdict}{`"pass"` / `"warn"` / `"fail"`}
#' }
#'
#' @section Verdict tiers:
#' Based on `universal_frac`:
#' `pass` if < 10\%, `warn` if 10-`universal_fail_threshold`, `fail` otherwise.
#' Halt classifier training when `verdict == "fail"`.
#'
#' @section Plot layout:
#' 3x3 grid:
#' \itemize{
#'   \item Row 1: rejection-count histogram, per-model rejection bars,
#'                summary text + verdict
#'   \item Row 2: PCA scatter (PC1 vs PC2) colored by model + obs star,
#'                PCA NN density per model, per-model PCA NN p-values
#'   \item Row 3: PLS-DA scatter (Comp1 vs Comp2) colored by model + obs star,
#'                PLS-DA NN density per model, per-model PLS-DA NN p-values
#' }
#'
#' @seealso \code{\link{tune.nn.classify}} for downstream NN classifier;
#'   \code{\link{OOD.pretrain}} for the single-model variant.
#'
#' @export
OOD.pretrain.classify <- function(reftable, model_col, observed,
                                   stat.cols = NULL,
                                   subsample = NULL,
                                   alpha = 0.01,
                                   universal_fail_threshold = 0.25,
                                   pca.var = 0.95,
                                   pls.n.comp = NULL,
                                   seed = 42,
                                   plot = TRUE,
                                   verbose = TRUE) {

  if (!model_col %in% colnames(reftable))
    stop(sprintf("OOD.pretrain.classify: model_col '%s' not in reftable",
                 model_col))

  raw_labels  <- as.character(reftable[[model_col]])
  class_names <- sort(unique(raw_labels))
  K <- length(class_names)
  if (K < 2L) stop("OOD.pretrain.classify needs >= 2 classes; got ", K)
  model_id_int <- match(raw_labels, class_names)

  # Resolve stat columns
  if (is.null(stat.cols)) {
    param_pattern <- "^(Ne[0-9]*\\.|Ne\\.anc|t\\.Ne|alpha|join[0-9]+|mig[0-9]+|mean\\.rate|sd\\.rate|mean\\.rec\\.rate|sd\\.rec\\.rate|rec\\.rate)"
    param_cols <- grep(param_pattern, colnames(reftable), value = TRUE)
    stat_cols <- setdiff(colnames(reftable), c(model_col, param_cols))
  } else {
    stat_cols <- stat.cols
  }

  observed <- .ood.prepare.observed(observed, stat_cols)
  if (length(observed) != length(stat_cols))
    stop(sprintf("observed has %d elements but %d stat columns specified",
                 length(observed), length(stat_cols)))

  # Optional per-class subsample for PCA/PLS-DA scalability
  if (!is.null(subsample) && is.finite(subsample) && subsample > 0L) {
    set.seed(seed)
    keep_rows <- integer(0)
    for (k in seq_len(K)) {
      rows_k <- which(model_id_int == k)
      if (length(rows_k) > subsample) rows_k <- sample(rows_k, subsample)
      keep_rows <- c(keep_rows, rows_k)
    }
    reftable     <- reftable[keep_rows, , drop = FALSE]
    model_id_int <- model_id_int[keep_rows]
    raw_labels   <- raw_labels[keep_rows]
  }

  S_sim <- as.matrix(reftable[, stat_cols, drop = FALSE])
  good_rows <- apply(S_sim, 1, function(x) all(is.finite(x)))
  if (any(!good_rows)) {
    if (verbose)
      cat(sprintf("PipeMaster:: dropping %d rows with non-finite stats\n",
                  sum(!good_rows)))
    S_sim        <- S_sim[good_rows, , drop = FALSE]
    model_id_int <- model_id_int[good_rows]
    raw_labels   <- raw_labels[good_rows]
  }

  # --- 1. Per-model rejection table ---
  rej <- .ood.classify.rejection(S_sim, model_id_int, observed,
                                 stat_cols, class_names)
  rejection      <- rej$rejection
  per_model_pct  <- rej$percentiles
  rejection_count <- colSums(rejection)
  universal_mask  <- rejection_count == K
  specific_mask   <- rejection_count > 0L & rejection_count < K
  matched_mask    <- rejection_count == 0L
  keep_mask       <- !universal_mask
  universal_frac <- mean(universal_mask)

  verdict <- if (universal_frac < 0.10) "pass"
             else if (universal_frac < universal_fail_threshold) "warn"
             else "fail"

  # --- 2. Combined PCA (using kept stats only) ---
  col_sd <- apply(S_sim, 2, sd, na.rm = TRUE)
  pca_keep <- keep_mask & col_sd > 0 & is.finite(col_sd)
  pca_res <- .ood.classify.pca(S_sim, observed, pca_keep, pca.var)

  # --- 3. PLS-DA ---
  plsda_res <- .ood.classify.plsda(S_sim, model_id_int, observed,
                                    pca_keep, class_names, n.comp = pls.n.comp)

  # --- 4. Per-model NN density (PCA-space and PLS-DA-space optional) ---
  per_model_nn <- .ood.classify.per.model.nn(S_sim, model_id_int, observed,
                                              pca_keep, class_names, alpha)
  per_model_plsda <- .ood.classify.per.model.plsda(
    plsda_res$scores, plsda_res$score_obs,
    model_id_int, class_names, alpha
  )

  ## Page-2 complement: Mahalanobis-to-centroid verdicts in PCA and PLS-DA
  ## score spaces. NN-density (above) catches local-density coverage;
  ## Mahalanobis catches centroid distance. Both reported.
  per_model_pca_mahal   <- .ood.classify.per.model.mahalanobis(
    pca_res$scores, pca_res$score_obs,
    model_id_int, class_names, alpha
  )
  per_model_plsda_mahal <- .ood.classify.per.model.mahalanobis(
    plsda_res$scores, plsda_res$score_obs,
    model_id_int, class_names, alpha
  )

  results <- list(
    class_names           = class_names,
    n_classes             = K,
    n_stats               = length(stat_cols),
    stat_cols             = stat_cols,
    rejection             = rejection,
    rejection_count       = rejection_count,
    universal_mask        = universal_mask,
    specific_mask         = specific_mask,
    matched_mask          = matched_mask,
    keep_mask             = keep_mask,
    keep_stat_cols        = stat_cols[keep_mask],
    universal_frac        = universal_frac,
    verdict               = verdict,
    percentiles_per_model = per_model_pct,
    pca                   = pca_res,
    pls_da                = plsda_res,
    per_model_nn          = per_model_nn,
    per_model_plsda       = per_model_plsda,
    per_model_pca_mahal   = per_model_pca_mahal,
    per_model_plsda_mahal = per_model_plsda_mahal,
    .context              = list(
      S_sim        = S_sim,
      model_id_int = model_id_int,
      observed     = observed,
      pca_keep     = pca_keep,
      alpha        = alpha,
      universal_fail_threshold = universal_fail_threshold
    )
  )
  class(results) <- c("OOD_pretrain_classify", "OOD_diagnostic")

  if (verbose) .ood.classify.print.summary(results)
  if (plot)    .ood.classify.plot(results)

  invisible(results)
}


# ---------------------------------------------------------------------------
# Print summary for OOD_pretrain_classify
# ---------------------------------------------------------------------------

.ood.classify.print.summary <- function(results) {
  cat("Pre-training OOD diagnostic — model selection\n")
  cat(sprintf("  Classes (K=%d): %s\n",
              results$n_classes,
              paste(results$class_names, collapse = ", ")))
  cat(sprintf("  Stats: %d total | %d matched (no model rejects) | %d model-specific outliers | %d universal outliers\n",
              results$n_stats, sum(results$matched_mask),
              sum(results$specific_mask), sum(results$universal_mask)))
  cat(sprintf("  Universal outlier fraction: %.1f%%\n",
              100 * results$universal_frac))

  cat("  Per-model rejection counts (out of all stats):\n")
  for (k in seq_along(results$class_names))
    cat(sprintf("    %-20s: %d / %d stats outside support\n",
                results$class_names[k],
                sum(results$rejection[k, ]), results$n_stats))

  if (!is.null(results$pca$scores))
    cat(sprintf("  PCA: %d components, %.1f%% variance, %d kept stats\n",
                results$pca$n_pcs, 100 * results$pca$var_explained,
                results$pca$n_stats))
  if (!is.null(results$pls_da$scores))
    cat(sprintf("  PLS-DA: %d components, %d kept stats\n",
                results$pls_da$n_comp, results$pls_da$n_stats))

  cat("\n  Per-method geometric verdict (PASS=inside, WARN=edge, FAIL=outside the sim cloud):\n")
  cat(sprintf("    %-20s  %-7s  %-7s  %-7s  %-7s  %-7s\n",
              "class", "d_obs", "Q95", "max", "p_val", "verdict"))

  ## PCA (per-class)
  if (length(results$per_model_nn) > 0L &&
      !is.null(results$per_model_nn[[1]]$verdict)) {
    cat(sprintf("  -- PCA (%d PCs) --\n",
                if (!is.null(results$pca$n_pcs)) results$pca$n_pcs else NA))
    for (k in seq_along(results$class_names)) {
      nn <- results$per_model_nn[[k]]
      if (is.null(nn$d_obs) || is.na(nn$d_obs)) next
      cat(sprintf("    %-20s  %-7.3f  %-7.3f  %-7.3f  %-7.4f  %s\n",
                  results$class_names[k],
                  nn$d_obs, nn$d_loo_q95, nn$d_loo_max,
                  nn$p_value, nn$verdict))
    }
  }

  ## PLS-DA (per-class)
  if (!is.null(results$per_model_plsda) &&
      !is.null(results$per_model_plsda[[1]]$verdict) &&
      !is.na(results$per_model_plsda[[1]]$verdict)) {
    cat(sprintf("  -- PLS-DA (%d comps) --\n",
                if (!is.null(results$pls_da$n_comp)) results$pls_da$n_comp else NA))
    for (k in seq_along(results$class_names)) {
      nn <- results$per_model_plsda[[k]]
      if (is.null(nn$d_obs) || is.na(nn$d_obs)) next
      cat(sprintf("    %-20s  %-7.3f  %-7.3f  %-7.3f  %-7.4f  %s\n",
                  results$class_names[k],
                  nn$d_obs, nn$d_loo_q95, nn$d_loo_max,
                  nn$p_value, nn$verdict))
    }
  }

  ## Overall PCA (single, all classes pooled)
  if (!is.null(results$pca) && !is.null(results$pca$verdict)) {
    cat(sprintf("  -- PCA  (pooled, %d PCs, %.0f%% var) -- d_obs=%.3f  Q95=%.3f  max=%.3f  verdict=%s\n",
                results$pca$n_pcs, 100 * results$pca$var_explained,
                results$pca$d_obs, results$pca$d_loo_q95,
                results$pca$d_loo_max, results$pca$verdict))
  }
  if (!is.null(results$pls_da) && !is.null(results$pls_da$verdict)) {
    cat(sprintf("  -- PLS  (pooled, %d comps) -- d_obs=%.3f  Q95=%.3f  max=%.3f  verdict=%s\n",
                results$pls_da$n_comp,
                results$pls_da$d_obs, results$pls_da$d_loo_q95,
                results$pls_da$d_loo_max, results$pls_da$verdict))
  }

  if (sum(results$universal_mask) > 0L) {
    cat("\n  Universal outliers (drop from classifier feature set):\n")
    bad_stats <- results$stat_cols[results$universal_mask]
    if (length(bad_stats) <= 20L) {
      for (s in bad_stats) cat(sprintf("    %s\n", s))
    } else {
      cat(sprintf("    (%d stats — see results$stat_cols[results$universal_mask])\n",
                  length(bad_stats)))
    }
  }
}


# ---------------------------------------------------------------------------
# Plot for OOD_pretrain_classify (3x3 layout)
# ---------------------------------------------------------------------------

.ood.classify.plot <- function(results) {
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))

  K <- results$n_classes
  cols <- if (K <= 8L)
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")[seq_len(K)]
  else rainbow(K)
  names(cols) <- results$class_names

  # ---- Row 1 ----
  # 1: rejection-count histogram (0..K)
  rc <- factor(results$rejection_count, levels = 0:K)
  bp <- barplot(table(rc), col = c("steelblue", rep("orange", K - 1L), "red"),
                border = "white",
                main = sprintf("Stats by # models rejecting (K=%d)", K),
                xlab = "# models rejecting", ylab = "n stats")
  abline(h = 0, col = "grey50")

  # 2: per-model rejection bar
  rej_per_model <- rowSums(results$rejection)
  barplot(rej_per_model, col = cols, border = "white",
          main = "Outliers per model (out of n_stats)",
          ylab = "n outlier stats", las = 2,
          ylim = c(0, max(rej_per_model, 1) * 1.1))

  # 3: per-class verdict panel (geometric)
  plot.new(); title(main = "Per-class PCA-NN verdict")
  v_col_for <- function(v) switch(as.character(v),
                                   PASS = "darkgreen",
                                   WARN = "darkorange",
                                   FAIL = "red",
                                   "grey40")
  text(0.5, 0.95, "PASS=inside, WARN=edge, FAIL=outside the sim cloud",
       cex = 0.8, col = "grey40")
  text(0.5, 0.85, sprintf("Universal outliers: %d / %d (%.1f%%)",
                          sum(results$universal_mask), results$n_stats,
                          100 * results$universal_frac), cex = 0.95)
  K_local <- length(results$class_names)
  y_top <- 0.70
  y_step <- 0.5 / max(1, K_local)
  for (i in seq_len(K_local)) {
    cn <- results$class_names[i]
    e  <- results$per_model_nn[[i]]
    v  <- if (!is.null(e$verdict) && !is.na(e$verdict)) e$verdict else "--"
    yi <- y_top - (i - 1) * y_step
    text(0.5, yi,
         sprintf("%s : %s", cn, v),
         col = v_col_for(v), cex = 1.0, font = 2)
  }
  text(0.5, 0.10, sprintf("Drop %d universal-outlier stats from features",
                          sum(results$universal_mask)),
       cex = 0.8, col = "grey40")

  # ---- Row 2: PCA — Comp1 vs LOO-NN in full PC space (per class) ----
  if (!is.null(results$pca$scores)) {
    sc <- results$pca$scores
    so <- results$pca$score_obs
    nn_pca <- .ood.classify.scatter.data(sc, so,
                                         results$.context$model_id_int, K)
    cls_v <- cols[results$class_names[results$.context$model_id_int]]
    y_all <- c(nn_pca$sim_loo, nn_pca$obs_d)
    y_all <- y_all[is.finite(y_all) & y_all > 0]
    use_log_y <- length(y_all) > 0 &&
                 (max(y_all) / min(y_all)) > 100
    pct_pc1 <- 100 * results$pca$fit$sdev[1]^2 /
                     sum(results$pca$fit$sdev^2)
    plot(sc[, 1], nn_pca$sim_loo,
         col = adjustcolor(cls_v, 0.5), pch = 16, cex = 0.5,
         log = if (use_log_y) "y" else "",
         xlab = sprintf("PC1 (%.1f%%)", pct_pc1),
         ylab = sprintf("LOO-NN distance in full %d-PC space (z-scored)",
                        results$pca$n_pcs),
         main = sprintf("PCA: PC1 vs LOO-NN (%d sims, %d stats, %d PCs)",
                        nrow(sc), results$pca$n_stats, results$pca$n_pcs))
    # K obs markers, one per class, at PC1_obs vs obs-to-nearest in class k
    for (k in seq_len(K)) {
      cn <- results$class_names[k]
      if (!is.na(nn_pca$obs_d[k]))
        points(so[1, 1], nn_pca$obs_d[k],
               col = cols[cn], pch = 8, cex = 2.0, lwd = 2.0)
    }
    legend("topright", legend = c(results$class_names, "obs (per class)"),
           col = c(cols, "black"), pch = c(rep(16, K), 8),
           pt.cex = c(rep(1, K), 1.6), cex = 0.7, bty = "n")
  } else {
    plot.new(); title(main = "PCA — unavailable")
  }

  # 5: per-model PCA NN distance bars (d_obs)
  d_obs_pca <- vapply(results$per_model_nn, function(x) x$d_obs, numeric(1))
  d_obs_pca[is.na(d_obs_pca)] <- 0
  barplot(d_obs_pca, col = cols, border = "white", names.arg = results$class_names,
          cex.names = 0.7,
          main = "PCA-space d_obs to nearest sim",
          ylab = "z-scored Euclidean distance", las = 2)

  # 6: per-class PCA verdict + percentile of obs in sim NN distribution
  .ood.verdict.panel <- function(nn_list, class_names, title_main) {
    plot.new(); title(main = title_main)
    v_col <- function(v) switch(as.character(v),
                                 PASS = "darkgreen", WARN = "darkorange",
                                 FAIL = "red", "grey60")
    text(0.5, 0.92,
         "PASS=inside  WARN=edge  FAIL=outside",
         cex = 0.75, col = "grey40")
    text(0.18, 0.78, "class",   adj = 0, cex = 0.9, font = 2)
    text(0.62, 0.78, "verdict", adj = 0, cex = 0.9, font = 2)
    text(0.85, 0.78, "pct",     adj = 1, cex = 0.9, font = 2)
    K_local <- length(class_names)
    y_top  <- 0.62
    y_step <- 0.5 / max(1, K_local)
    for (i in seq_len(K_local)) {
      nn <- nn_list[[i]]
      yi <- y_top - (i - 1) * y_step
      v  <- if (!is.null(nn$verdict) && !is.na(nn$verdict)) nn$verdict else "--"
      p  <- if (!is.null(nn$p_value) && !is.na(nn$p_value)) nn$p_value else NA_real_
      pct_label <- if (is.na(p)) "--"
                   else if (v == "FAIL") ">100"
                   else sprintf("%.1f", 100 * (1 - p))
      text(0.18, yi, class_names[i], adj = 0, cex = 0.95)
      text(0.62, yi, v,
           adj = 0, cex = 1.0, font = 2, col = v_col(v))
      text(0.85, yi, pct_label, adj = 1, cex = 0.95, col = "grey20")
    }
  }
  .ood.verdict.panel(results$per_model_nn, results$class_names,
                     "PCA per-class verdict")

  # ---- Row 3: PLS-DA — Comp1 vs LOO-NN in full PLS space (per class) ----
  if (!is.null(results$pls_da$scores) && ncol(results$pls_da$scores) >= 2L) {
    sc <- results$pls_da$scores
    so <- results$pls_da$score_obs
    nn_pls <- .ood.classify.scatter.data(sc, so,
                                         results$.context$model_id_int, K)
    cls_v <- cols[results$class_names[results$.context$model_id_int]]
    y_all <- c(nn_pls$sim_loo, nn_pls$obs_d)
    y_all <- y_all[is.finite(y_all) & y_all > 0]
    use_log_y <- length(y_all) > 0 &&
                 (max(y_all) / min(y_all)) > 100
    plot(sc[, 1], nn_pls$sim_loo,
         col = adjustcolor(cls_v, 0.5), pch = 16, cex = 0.5,
         log = if (use_log_y) "y" else "",
         xlab = "PLS Comp 1",
         ylab = sprintf("LOO-NN distance in full %d-comp space (z-scored)",
                        results$pls_da$n_comp),
         main = sprintf("PLS-DA: Comp1 vs LOO-NN (%d comps, %d stats)",
                        results$pls_da$n_comp, results$pls_da$n_stats))
    for (k in seq_len(K)) {
      cn <- results$class_names[k]
      if (!is.na(nn_pls$obs_d[k]))
        points(so[1, 1], nn_pls$obs_d[k],
               col = cols[cn], pch = 8, cex = 2.0, lwd = 2.0)
    }
    legend("topright", legend = c(results$class_names, "obs (per class)"),
           col = c(cols, "black"), pch = c(rep(16, K), 8),
           pt.cex = c(rep(1, K), 1.6), cex = 0.7, bty = "n")

    # 8: per-model PLS-DA d_obs and 9: p-values (reuse nn_pls computation)
    d_obs_pls <- nn_pls$obs_d; d_obs_pls[is.na(d_obs_pls)] <- 0
    barplot(d_obs_pls, col = cols, border = "white",
            names.arg = results$class_names,
            cex.names = 0.7,
            main = "PLS-DA-space d_obs to nearest sim",
            ylab = "z-scored Euclidean distance", las = 2)

    # 9: per-class PLS-DA verdict + percentile
    .ood.verdict.panel(results$per_model_plsda, results$class_names,
                       "PLS-DA per-class verdict")
  } else {
    plot.new(); title(main = "PLS-DA — unavailable")
    plot.new(); plot.new()
  }

  # ============================================================================
  # Page 2 — Mahalanobis-to-centroid verdicts + 2D (Comp1 vs Comp2) scatters.
  # Complements page 1's K=1 NN-density view (local density) with a centroid-
  # distance view (Mahalanobis), and exposes the second PC/PLS axis the
  # page-1 scatter hides.
  # ============================================================================
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

  ## --- Row 1: PCA Comp1 vs Comp2 + per-class Mahalanobis ---
  if (!is.null(results$pca$scores) && ncol(results$pca$scores) >= 2L) {
    sc  <- results$pca$scores
    so  <- results$pca$score_obs
    mid <- results$.context$model_id_int
    cls_v <- cols[results$class_names[mid]]
    pct_pc1 <- 100 * results$pca$fit$sdev[1]^2 / sum(results$pca$fit$sdev^2)
    pct_pc2 <- 100 * results$pca$fit$sdev[2]^2 / sum(results$pca$fit$sdev^2)
    plot(sc[, 1], sc[, 2],
         col = adjustcolor(cls_v, 0.4), pch = 16, cex = 0.5,
         xlab = sprintf("PC1 (%.1f%%)", pct_pc1),
         ylab = sprintf("PC2 (%.1f%%)", pct_pc2),
         main = sprintf("PCA: PC1 vs PC2 (%.1f%% var)", pct_pc1 + pct_pc2))
    ## per-class centroids
    for (k in seq_len(K)) {
      rows <- which(mid == k)
      if (length(rows) < 1L) next
      cn <- results$class_names[k]
      points(mean(sc[rows, 1]), mean(sc[rows, 2]),
             col = "black", pch = 21, bg = cols[cn],
             cex = 2.4, lwd = 1.4)
    }
    ## obs (single point — same in all classes)
    points(so[1, 1], so[1, 2], pch = 8, col = "red",
           cex = 2.6, lwd = 2.5)
    legend("topright", legend = c(results$class_names, "centroid", "obs"),
           col = c(cols, "black", "red"),
           pch = c(rep(16, K), 21, 8),
           pt.bg = c(rep(NA, K), "grey80", NA),
           pt.cex = c(rep(1, K), 1.6, 1.6),
           cex = 0.7, bty = "n")
  } else {
    plot.new(); title(main = "PCA Comp1 vs Comp2 — unavailable")
  }

  ## per-class Mahalanobis-to-centroid bar (PCA space)
  d_mahal_pca <- vapply(results$per_model_pca_mahal,
                        function(x) x$d_obs, numeric(1))
  d_mahal_pca[is.na(d_mahal_pca)] <- 0
  barplot(d_mahal_pca, col = cols, border = "white",
          names.arg = results$class_names, cex.names = 0.7,
          main = "PCA Mahalanobis: obs to centroid",
          ylab = "Mahalanobis distance", las = 2)

  ## per-class Mahalanobis verdict
  .ood.verdict.panel(results$per_model_pca_mahal, results$class_names,
                     "PCA per-class Mahalanobis verdict")

  ## --- Row 2: PLS-DA Comp1 vs Comp2 + per-class Mahalanobis ---
  if (!is.null(results$pls_da$scores) && ncol(results$pls_da$scores) >= 2L) {
    sc  <- results$pls_da$scores
    so  <- results$pls_da$score_obs
    mid <- results$.context$model_id_int
    cls_v <- cols[results$class_names[mid]]
    plot(sc[, 1], sc[, 2],
         col = adjustcolor(cls_v, 0.4), pch = 16, cex = 0.5,
         xlab = "PLS Comp 1", ylab = "PLS Comp 2",
         main = sprintf("PLS-DA: Comp1 vs Comp2 (%d comps)",
                        results$pls_da$n_comp))
    for (k in seq_len(K)) {
      rows <- which(mid == k)
      if (length(rows) < 1L) next
      cn <- results$class_names[k]
      points(mean(sc[rows, 1]), mean(sc[rows, 2]),
             col = "black", pch = 21, bg = cols[cn],
             cex = 2.4, lwd = 1.4)
    }
    points(so[1, 1], so[1, 2], pch = 8, col = "red",
           cex = 2.6, lwd = 2.5)
    legend("topright", legend = c(results$class_names, "centroid", "obs"),
           col = c(cols, "black", "red"),
           pch = c(rep(16, K), 21, 8),
           pt.bg = c(rep(NA, K), "grey80", NA),
           pt.cex = c(rep(1, K), 1.6, 1.6),
           cex = 0.7, bty = "n")
  } else {
    plot.new(); title(main = "PLS-DA Comp1 vs Comp2 — unavailable")
  }

  d_mahal_pls <- vapply(results$per_model_plsda_mahal,
                        function(x) x$d_obs, numeric(1))
  d_mahal_pls[is.na(d_mahal_pls)] <- 0
  barplot(d_mahal_pls, col = cols, border = "white",
          names.arg = results$class_names, cex.names = 0.7,
          main = "PLS-DA Mahalanobis: obs to centroid",
          ylab = "Mahalanobis distance", las = 2)

  .ood.verdict.panel(results$per_model_plsda_mahal, results$class_names,
                     "PLS-DA per-class Mahalanobis verdict")
}


# --- Internal: per-class LOO-NN distance in score space (PCA / PLS-DA / etc.).
#     Returns sim_loo (per-row LOO-NN distance to nearest same-class sim,
#     numeric of length nrow(scores), NA outside any class) and obs_d
#     (length-K vector of obs-to-nearest distance per class), plus p_value
#     (length K, one-sided rank test of obs vs LOO null). All distances in
#     z-scored full-score space, so cross-class comparisons share a metric. ---
.ood.classify.scatter.data <- function(scores, score_obs, model_id_int, K) {
  n <- nrow(scores)
  mu  <- colMeans(scores)
  sd_ <- apply(scores, 2, sd); sd_[sd_ == 0] <- 1
  Z   <- sweep(sweep(scores, 2, mu, "-"), 2, sd_, "/")
  z_o <- (as.numeric(score_obs[1, ]) - mu) / sd_
  z_o_mat <- matrix(z_o, nrow = 1L)

  sim_loo <- rep(NA_real_, n)
  obs_d   <- rep(NA_real_, K)
  p_value <- rep(NA_real_, K)

  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) next
    Z_k <- Z[rows, , drop = FALSE]
    d_loo <- .ood.knn.loo(Z_k)
    d_obs <- .ood.knn.obs(Z_k, z_o_mat)
    sim_loo[rows] <- d_loo
    obs_d[k]   <- d_obs
    p_value[k] <- mean(d_loo >= d_obs)
  }
  list(sim_loo = sim_loo, obs_d = obs_d, p_value = p_value)
}


#' @export
print.OOD_pretrain_classify <- function(x, ...) {
  .ood.classify.print.summary(x)
  invisible(x)
}

#' @export
plot.OOD_pretrain_classify <- function(x, ...) {
  .ood.classify.plot(x)
  invisible(NULL)
}

#' @export
summary.OOD_pretrain_classify <- function(object, ...) {
  out <- data.frame(
    model = object$class_names,
    n_outliers = rowSums(object$rejection),
    pca_d_obs = vapply(object$per_model_nn, function(x)
      if (is.null(x$d_obs)) NA_real_ else x$d_obs, numeric(1)),
    pca_p = vapply(object$per_model_nn, function(x)
      if (is.null(x$p_value)) NA_real_ else x$p_value, numeric(1)),
    stringsAsFactors = FALSE
  )
  cat(sprintf("Verdict: %s | universal frac = %.1f%% | %d stats dropped from classifier\n\n",
              toupper(object$verdict), 100 * object$universal_frac,
              sum(object$universal_mask)))
  print(out)
  invisible(out)
}


# ============================================================================
# OOD.posttrain.classify — Post-training diagnostic for model selection
#
# Companion to OOD.pretrain.classify(): uses the trained classifier's
# penultimate-layer activations (NN-latent space) to ask whether the
# observation lies inside any class's learned manifold, and whether it
# lies inside the predicted class specifically.
# ============================================================================

# --- Internal: per-class NN-latent density. Globally z-score the latent so
#     cross-class distances live in the same metric, then run LOO-NN within
#     each class subset. Mirrors .ood.classify.per.model.nn but on latent
#     coordinates, not raw stats. ---
.ood.classify.latent.nn <- function(Z_sim, model_id_int, Z_obs,
                                    class_names, alpha, knn.k = 10L) {
  K <- length(class_names)
  mu <- colMeans(Z_sim)
  sd_ <- apply(Z_sim, 2, sd)
  sd_[sd_ == 0 | !is.finite(sd_)] <- 1
  Z_n   <- sweep(sweep(Z_sim, 2, mu, "-"), 2, sd_, "/")
  Z_o_n <- sweep(sweep(Z_obs, 2, mu, "-"), 2, sd_, "/")

  out <- vector("list", K); names(out) <- class_names
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) {
      out[[k]] <- list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                       verdict = NA_character_,
                       null_distribution = numeric(0),
                       knn.k = as.integer(knn.k))
      next
    }
    Z_k <- Z_n[rows, , drop = FALSE]
    # k-NN average distance (Sun et al. 2022, ICLR).
    d_obs <- .ood.knn.k.obs(Z_k, Z_o_n, k = knn.k)
    d_loo <- .ood.knn.k.loo(Z_k,        k = knn.k)
    p <- mean(d_loo >= d_obs)
    verdict <- .ood.geometric.verdict(d_obs, d_loo)
    out[[k]] <- list(d_obs = d_obs, p_value = p, pass = p >= alpha,
                     verdict = verdict,
                     d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
                     d_loo_max = max(d_loo),
                     null_distribution = d_loo,
                     knn.k = as.integer(knn.k))
  }
  list(per_class = out, scores = Z_n, score_obs = Z_o_n,
       n_dims = ncol(Z_sim), knn.k = as.integer(knn.k))
}

# Per-class Mahalanobis distance in penultimate latent space.
# Distance from sim/obs to each class mean using pooled within-class Sigma
# (Lee et al. 2018, NeurIPS). Per-sim LOO distance to its OWN class mean
# forms the empirical null for the geometric verdict.
.ood.classify.latent.mahalanobis <- function(Z_sim, model_id_int, Z_obs,
                                              class_names, alpha) {
  K <- length(class_names)
  mu <- colMeans(Z_sim)
  sd_ <- apply(Z_sim, 2, sd)
  sd_[sd_ == 0 | !is.finite(sd_)] <- 1
  Z_n   <- sweep(sweep(Z_sim, 2, mu, "-"), 2, sd_, "/")
  Z_o_n <- sweep(sweep(Z_obs, 2, mu, "-"), 2, sd_, "/")

  d_loo_all <- .ood.mahalanobis.loo(Z_n, model_id_int)
  d_obs_per_class <- .ood.mahalanobis.obs.per.class(
    Z_n, Z_o_n, model_id_int)

  out <- vector("list", K); names(out) <- class_names
  for (k in seq_len(K)) {
    rows <- which(model_id_int == k)
    if (length(rows) < 2L) {
      out[[k]] <- list(d_obs = NA_real_, p_value = NA_real_, pass = NA,
                       verdict = NA_character_,
                       null_distribution = numeric(0))
      next
    }
    d_loo <- d_loo_all[rows]
    d_loo <- d_loo[is.finite(d_loo)]
    d_obs <- d_obs_per_class[as.character(k)]
    if (is.na(d_obs)) d_obs <- d_obs_per_class[k]
    p <- mean(d_loo >= d_obs)
    verdict <- .ood.geometric.verdict(d_obs, d_loo)
    out[[k]] <- list(d_obs = as.numeric(d_obs), p_value = p, pass = p >= alpha,
                     verdict = verdict,
                     d_loo_q95 = as.numeric(quantile(d_loo, 0.95, na.rm = TRUE)),
                     d_loo_max = max(d_loo),
                     null_distribution = d_loo)
  }
  list(per_class = out)
}


#' Post-training Out-of-Distribution Diagnostic for Model Selection
#'
#' Companion to \code{\link{OOD.pretrain.classify}}. Uses the trained
#' classifier's penultimate-layer (NN-latent) representation — the space
#' in which the classifier discriminates between models — to check
#' whether the observation lies inside any class's learned manifold,
#' and inside the predicted class specifically.
#'
#' Three outcomes:
#' \itemize{
#'   \item \strong{pass}: obs lies inside the predicted class's latent cloud.
#'   \item \strong{warn}: obs lies inside some class but not the predicted one
#'         (classifier picked the wrong manifold; prediction unreliable).
#'   \item \strong{fail}: obs lies outside every class's latent cloud
#'         (extrapolation; classifier output is unreliable).
#' }
#'
#' @param classifier Object returned by \code{\link{tune.nn.classify}}.
#' @param observed Numeric vector or 1-row data.frame/matrix of observed
#'   summary statistics. Column names must match the classifier's
#'   \code{feature.cols} (pre-augmentation).
#' @param reftable A reference table (training or held-out) containing
#'   the classifier's feature columns plus a label column.
#' @param model_col Name of the label column in \code{reftable}.
#' @param pretrain Optional \code{OOD_pretrain_classify} result to inherit
#'   for the combined verdict.
#' @param alpha Significance level for per-class NN density p-values
#'   (default 0.01).
#' @param subsample Per-class subsample size for the latent NN check
#'   (default NULL = use all sims).
#' @param seed Random seed for subsampling.
#' @param plot Logical; produce diagnostic plot.
#' @param verbose Logical; print summary.
#'
#' @return A list of class \code{c("OOD_posttrain_classify", "OOD_diagnostic")}
#'   with fields:
#' \describe{
#'   \item{class_names, n_classes, n_features}{classifier metadata}
#'   \item{prediction}{full output of \code{nn.predict.classify}}
#'   \item{predicted_class, predicted_prob, margin}{top-1 class, its prob,
#'     and \code{P(top1) - P(top2)}}
#'   \item{latent_nn}{list of K, each with \code{d_obs}, \code{p_value},
#'     \code{pass}, \code{null_distribution}}
#'   \item{latent_scores, latent_score_obs, latent_class_id}{z-scored
#'     latent coordinates for sims and obs}
#'   \item{n_latent_dims}{penultimate-layer width}
#'   \item{obs_in_predicted_class}{\code{latent_nn[[predicted]]$pass}}
#'   \item{obs_in_any_class}{any class p >= alpha}
#'   \item{pretrain}{inherited \code{OOD_pretrain_classify} (or NULL)}
#'   \item{verdict}{combined pass/warn/fail}
#' }
#'
#' @section Plot layout:
#' 2x3 grid:
#' \itemize{
#'   \item Row 1: predicted P(class) bars; per-class latent NN p-value bars;
#'         verdict text panel
#'   \item Row 2: latent scatter (Z1 vs Z2) colored by class with obs star;
#'         per-class LOO-NN null distributions with d_obs marked;
#'         d_obs / median d_loo ratio per class
#' }
#'
#' @seealso \code{\link{OOD.pretrain.classify}} for the pre-training variant;
#'   \code{\link{tune.nn.classify}} for the classifier;
#'   \code{\link{OOD.posttrain}} for the regression analogue.
#'
#' @export
OOD.posttrain.classify <- function(classifier, observed, reftable,
                                    model_col = "model",
                                    pretrain = NULL,
                                    alpha = 0.01,
                                    subsample = NULL,
                                    seed = 42,
                                    plot = TRUE,
                                    verbose = TRUE) {
  if (!inherits(classifier, "model_selection_classifier"))
    stop("classifier must come from tune.nn.classify()")
  if (!model_col %in% colnames(reftable))
    stop(sprintf("OOD.posttrain.classify: model_col '%s' not in reftable",
                 model_col))
  if (!requireNamespace("torch", quietly = TRUE))
    stop("OOD.posttrain.classify requires the 'torch' package.")

  feat_cols   <- classifier$feature.cols
  feat_mu     <- classifier$feat_mu
  feat_sd     <- classifier$feat_sd
  K           <- classifier$n_classes
  class_names <- classifier$class_names
  device <- if (!is.null(classifier$device)) classifier$device else "cpu"

  # --- Resolve labels and subsample ---
  raw_labels <- as.character(reftable[[model_col]])
  miss_lbl <- setdiff(unique(raw_labels), class_names)
  if (length(miss_lbl) > 0L)
    stop("reftable contains labels not seen by classifier: ",
         paste(miss_lbl, collapse = ", "))
  model_id_int <- match(raw_labels, class_names)

  set.seed(seed)
  if (!is.null(subsample)) {
    counts <- table(model_id_int)
    if (any(counts > subsample)) {
      keep_idx <- unlist(lapply(seq_len(K), function(k) {
        r <- which(model_id_int == k)
        if (length(r) <= subsample) r else sample(r, subsample)
      }))
      reftable <- reftable[keep_idx, , drop = FALSE]
      model_id_int <- model_id_int[keep_idx]
    }
  }

  # --- Resolve observed (named numeric vector or 1-row df/matrix) ---
  if (is.data.frame(observed) || is.matrix(observed)) {
    if (nrow(observed) > 1L)
      stop("observed must be a single row")
    obs_named <- setNames(as.numeric(observed[1, ]), colnames(observed))
  } else {
    if (is.null(names(observed)))
      stop("observed must be a named numeric vector")
    obs_named <- setNames(as.numeric(observed), names(observed))
  }
  miss_feat <- setdiff(feat_cols, names(obs_named))
  if (length(miss_feat) > 0L)
    stop("observed missing feature columns: ",
         paste(miss_feat, collapse = ", "))
  obs_raw <- obs_named[feat_cols]

  # --- Build feature matrix matching training-time augmentation ---
  S_raw <- as.matrix(reftable[, feat_cols, drop = FALSE])
  bad <- apply(S_raw, 1, function(x) any(!is.finite(x)))
  if (any(bad)) {
    if (verbose)
      cat(sprintf("PipeMaster:: dropping %d rows with non-finite features\n",
                  sum(bad)))
    S_raw <- S_raw[!bad, , drop = FALSE]
    model_id_int <- model_id_int[!bad]
  }
  S_aug <- cbind(S_raw, log1p(abs(S_raw)))
  S_z   <- sweep(sweep(S_aug, 2, feat_mu, "-"), 2, feat_sd, "/")

  obs_aug <- c(obs_raw, log1p(abs(obs_raw)))
  obs_z   <- matrix((as.numeric(obs_aug) - feat_mu) / feat_sd, nrow = 1L)

  # --- Extract penultimate-layer activations ---
  if (verbose)
    cat(sprintf("PipeMaster:: extracting latent features (%d sims, K=%d)\n",
                nrow(S_z), K))
  dev <- torch::torch_device(device)
  Z_sim <- .torch.penultimate(classifier$model, S_z, device = dev)
  Z_obs <- .torch.penultimate(classifier$model, obs_z, device = dev)

  # --- Per-class latent k-NN density (Sun 2022) + Mahalanobis (Lee 2018) ---
  latent <- .ood.classify.latent.nn(Z_sim, model_id_int, Z_obs,
                                    class_names, alpha, knn.k = 10L)
  mahal  <- .ood.classify.latent.mahalanobis(Z_sim, model_id_int, Z_obs,
                                              class_names, alpha)

  # --- Predicted class via the classifier ---
  obs_df <- as.data.frame(matrix(obs_raw, nrow = 1L,
                                 dimnames = list(NULL, feat_cols)))
  prediction <- nn.predict.classify(classifier, obs_df,
                                    mc_samples = 0L, verbose = FALSE)
  pred_class  <- prediction$best_model
  pred_prob   <- prediction$best_prob
  prob_sorted <- sort(prediction$prob_mean, decreasing = TRUE)
  margin <- if (length(prob_sorted) >= 2L)
              unname(prob_sorted[1] - prob_sorted[2]) else NA_real_

  # --- Bayes factors of predicted class vs each other class ---
  # Under uniform model prior, softmax_k ~ P(M_k | obs); BF_ij = p_i / p_j.
  # Returns named vector (per class) of BF(predicted vs class). Self BF is 1.
  pred_p <- prediction$prob_mean[pred_class]
  bf_vec <- sapply(class_names, function(cn) {
    p_k <- prediction$prob_mean[cn]
    if (is.null(p_k) || !is.finite(p_k) || p_k == 0) Inf
    else as.numeric(pred_p / p_k)
  })
  names(bf_vec) <- class_names
  log10_bf_vec <- log10(bf_vec)
  log10_bf_runner_up <- if (length(prob_sorted) >= 2L)
    log10(prob_sorted[1] / prob_sorted[2]) else NA_real_
  jeffreys_label <- function(lbf) {
    if (!is.finite(lbf)) return("--")
    if (abs(lbf) > 2)   return("decisive")
    if (abs(lbf) > 1.5) return("very strong")
    if (abs(lbf) > 1)   return("strong")
    if (abs(lbf) > 0.5) return("substantial")
    "weak"
  }
  bf_runner_up_label <- jeffreys_label(log10_bf_runner_up)

  # --- Verdict ---
  # "Obs in class" = obs is at most at the edge of the cloud, i.e. NOT FAIL.
  # PASS = inside  -> YES,  WARN = at the edge -> YES,  FAIL = outside -> NO.
  .in_cloud <- function(e) {
    v <- e$verdict
    !is.null(v) && !is.na(v) && v %in% c("PASS", "WARN")
  }
  pred_pass <- .in_cloud(latent$per_class[[pred_class]])
  any_pass  <- any(vapply(latent$per_class, .in_cloud, logical(1)))
  verdict_secondary <- if (pred_pass) "pass"
                       else if (any_pass) "warn"
                       else "fail"

  pre_tier <- if (!is.null(pretrain) && !is.null(pretrain$verdict))
                pretrain$verdict else "pass"
  combined <- if (pre_tier == "fail" || verdict_secondary == "fail") "fail"
              else if (pre_tier == "warn" || verdict_secondary == "warn") "warn"
              else "pass"

  results <- list(
    class_names      = class_names,
    n_classes        = K,
    n_features       = length(feat_cols),
    prediction       = prediction,
    predicted_class  = pred_class,
    predicted_prob   = pred_prob,
    margin           = margin,
    bayes_factor     = bf_vec,
    log10_bayes_factor = log10_bf_vec,
    log10_bf_runner_up = log10_bf_runner_up,
    bf_runner_up_label = bf_runner_up_label,
    latent_nn        = latent$per_class,
    latent_mahal     = mahal$per_class,
    knn.k            = latent$knn.k,
    latent_scores    = latent$scores,
    latent_score_obs = latent$score_obs,
    latent_class_id  = model_id_int,
    n_latent_dims    = latent$n_dims,
    obs_in_predicted_class = pred_pass,
    obs_in_any_class       = any_pass,
    pretrain         = pretrain,
    verdict          = combined,
    .context = list(alpha = alpha)
  )
  class(results) <- c("OOD_posttrain_classify", "OOD_diagnostic")

  if (verbose) .ood.posttrain.classify.print.summary(results)
  if (plot)    .ood.posttrain.classify.plot(results)

  invisible(results)
}


# --- Internal: print summary for OOD_posttrain_classify ---
.ood.posttrain.classify.print.summary <- function(x) {
  cat("\nOOD.posttrain.classify summary\n")
  cat(sprintf("  Predicted class : %s (P=%.3f, margin=%.3f)\n",
              x$predicted_class, x$predicted_prob, x$margin))
  if (!is.null(x$log10_bf_runner_up) && is.finite(x$log10_bf_runner_up))
    cat(sprintf("  log10 BF (pred vs runner-up) : %.2f  [%s]\n",
                x$log10_bf_runner_up,
                if (!is.null(x$bf_runner_up_label)) x$bf_runner_up_label else "--"))
  cat(sprintf("  Latent dims     : %d (penultimate layer)\n",
              x$n_latent_dims))

  # Bayes factor table: predicted vs each other class
  if (!is.null(x$bayes_factor)) {
    cat(sprintf("\n  Bayes factors (predicted=%s vs each class):\n",
                x$predicted_class))
    cat(sprintf("    %-20s  %-10s  %-10s  %-12s\n",
                "class", "BF", "log10 BF", "Jeffreys"))
    jeffreys_label <- function(lbf) {
      if (!is.finite(lbf)) return("--")
      if (abs(lbf) > 2)   return("decisive")
      if (abs(lbf) > 1.5) return("very strong")
      if (abs(lbf) > 1)   return("strong")
      if (abs(lbf) > 0.5) return("substantial")
      "weak"
    }
    for (cn in x$class_names) {
      if (cn == x$predicted_class) next
      bf <- x$bayes_factor[[cn]]
      lbf <- x$log10_bayes_factor[[cn]]
      bf_str <- if (!is.finite(bf)) "Inf" else
                if (bf >= 100) sprintf("%.1f", bf) else
                sprintf("%.2f", bf)
      lbf_str <- if (!is.finite(lbf)) "Inf" else sprintf("%.2f", lbf)
      cat(sprintf("    %-20s  %-10s  %-10s  %-12s\n",
                  cn, bf_str, lbf_str, jeffreys_label(lbf)))
    }
  }
  cat(sprintf("\n  Per-class latent verdicts (PASS=inside, WARN=edge, FAIL=outside):\n"))
  cat(sprintf("    k-NN k=%d (Sun 2022) + Mahalanobis (Lee 2018)\n\n",
              if (!is.null(x$knn.k)) x$knn.k else 10L))
  cat(sprintf("    %-18s | %-7s  %-7s  %-7s | %-7s  %-7s  %-7s\n",
              "class",
              "k-NN d", "k-NN pct", "k-NN",
              "M d",    "M pct",    "M"))
  for (cn in x$class_names) {
    e_nn <- x$latent_nn[[cn]]
    e_m  <- x$latent_mahal[[cn]]
    if (is.null(e_nn$d_obs) || is.na(e_nn$d_obs)) next
    star <- if (cn == x$predicted_class) " *" else "  "
    v_nn <- if (!is.null(e_nn$verdict)) e_nn$verdict else "--"
    v_m  <- if (!is.null(e_m$verdict))  e_m$verdict  else "--"
    pct_nn <- if (is.na(e_nn$p_value)) "--"
              else if (v_nn == "FAIL") ">100"
              else sprintf("%.1f", 100 * (1 - e_nn$p_value))
    pct_m  <- if (is.null(e_m$p_value) || is.na(e_m$p_value)) "--"
              else if (v_m == "FAIL") ">100"
              else sprintf("%.1f", 100 * (1 - e_m$p_value))
    cat(sprintf("    %-16s%s | %-7.3f  %-7s  %-7s | %-7.3f  %-7s  %-7s\n",
                cn, star,
                e_nn$d_obs, pct_nn, v_nn,
                if (is.null(e_m$d_obs) || is.na(e_m$d_obs)) NA_real_ else e_m$d_obs,
                pct_m, v_m))
  }
  cat("\n    (k-NN d = mean distance to k nearest sims; M d = Mahalanobis)\n")
  cat("    * = predicted class\n")
  cat(sprintf("\n  Obs in predicted class : %s\n",
              ifelse(x$obs_in_predicted_class, "YES", "NO")))
  cat(sprintf("  Obs in any class       : %s\n",
              ifelse(x$obs_in_any_class, "YES", "NO")))
  cat("\n")
}


# --- Internal: 2x3 plot for OOD_posttrain_classify ---
.ood.posttrain.classify.plot <- function(x) {
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))

  K <- x$n_classes
  cols <- if (K <= 8L)
    c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
      "#9467bd", "#8c564b", "#e377c2", "#7f7f7f")[seq_len(K)]
  else rainbow(K)
  names(cols) <- x$class_names
  alpha <- x$.context$alpha

  # (1,1) Predicted P(class) bars
  pm <- x$prediction$prob_mean[x$class_names]
  bar_col <- ifelse(names(pm) == x$predicted_class,
                    cols[x$predicted_class],
                    adjustcolor(cols[names(pm)], 0.5))
  border_col <- ifelse(names(pm) == x$predicted_class, "black", "white")
  bp <- barplot(pm, col = bar_col, border = border_col,
                ylim = c(0, max(pm) * 1.18),
                main = sprintf("Predicted P(class) | margin=%.3f", x$margin),
                ylab = "P(class | obs)", las = 2)
  text(bp, pm + max(pm) * 0.04, sprintf("%.3f", pm), cex = 0.75)

  v_col_for <- function(v) switch(as.character(v),
                                   PASS = "darkgreen", WARN = "darkorange",
                                   FAIL = "red", "grey60")
  # Helper: per-class verdict + percentile table (single method).
  one_method_panel <- function(per_class_list, class_names, predicted_class,
                               title_main) {
    plot.new(); title(main = title_main)
    text(0.5, 0.94, "PASS=inside  WARN=edge  FAIL=outside",
         cex = 0.72, col = "grey40")
    text(0.18, 0.80, "class",   adj = 0, cex = 0.9, font = 2)
    text(0.62, 0.80, "verdict", adj = 0, cex = 0.9, font = 2)
    text(0.92, 0.80, "pct",     adj = 1, cex = 0.9, font = 2)
    K_local <- length(class_names)
    y_top  <- 0.65
    y_step <- 0.55 / max(1, K_local)
    for (i in seq_len(K_local)) {
      cn <- class_names[i]
      e  <- per_class_list[[cn]]
      yi <- y_top - (i - 1) * y_step
      v  <- if (!is.null(e$verdict) && !is.na(e$verdict)) e$verdict else "--"
      p  <- if (!is.null(e$p_value) && !is.na(e$p_value)) e$p_value else NA_real_
      pct_label <- if (is.na(p)) "--"
                   else if (v == "FAIL") ">100"
                   else sprintf("%.1f", 100 * (1 - p))
      star <- if (cn == predicted_class) " *" else ""
      text(0.18, yi, paste0(cn, star), adj = 0, cex = 0.95)
      text(0.62, yi, v, adj = 0, cex = 1.0, font = 2, col = v_col_for(v))
      text(0.92, yi, pct_label, adj = 1, cex = 0.95, col = "grey20")
    }
    text(0.5, 0.06, "* = predicted class", cex = 0.7, col = "grey40")
  }

  # (1,2) k-NN per-class verdict
  one_method_panel(x$latent_nn, x$class_names, x$predicted_class,
                   sprintf("k-NN k=%d  (Sun 2022)",
                           if (!is.null(x$knn.k)) x$knn.k else 10L))

  # (1,3) Mahalanobis per-class verdict
  one_method_panel(x$latent_mahal, x$class_names, x$predicted_class,
                   "Mahalanobis  (Lee 2018)")

  # ---- Bottom-row scatter helper ----
  # X = Latent Z1, Y = per-sim LOO distance in the chosen metric. Per-class
  # color; obs marked once per class at (Z1_obs, obs-to-class-k distance).
  scatter_panel <- function(sc, so, per_class_list, ylab_text, title_main) {
    sim_loo <- rep(NA_real_, nrow(sc))
    for (k in seq_len(K)) {
      cn <- x$class_names[k]
      rows <- which(x$latent_class_id == k)
      nd <- per_class_list[[cn]]$null_distribution
      if (length(rows) >= 2L && length(nd) == length(rows))
        sim_loo[rows] <- nd
    }
    obs_d <- vapply(per_class_list, function(e)
      if (is.null(e$d_obs)) NA_real_ else e$d_obs, numeric(1))
    cls_v <- cols[x$class_names[x$latent_class_id]]
    y_all <- c(sim_loo, obs_d)
    y_all <- y_all[is.finite(y_all) & y_all > 0]
    use_log_y <- length(y_all) > 0 &&
                 (max(y_all) / min(y_all)) > 100
    plot(sc[, 1], sim_loo,
         col = adjustcolor(cls_v, 0.5), pch = 16, cex = 0.5,
         log = if (use_log_y) "y" else "",
         xlab = "Latent Z1 (z-scored)",
         ylab = ylab_text,
         main = title_main)
    for (k in seq_len(K)) {
      cn <- x$class_names[k]
      if (!is.na(obs_d[k]))
        points(so[1, 1], obs_d[k],
               col = cols[cn], pch = 8, cex = 2.0, lwd = 2.0)
    }
    legend("topright", legend = c(x$class_names, "obs (per class)"),
           col = c(cols, "black"), pch = c(rep(16, K), 8),
           pt.cex = c(rep(1, K), 1.6), cex = 0.65, bty = "n")
  }

  # (2,1) k-NN scatter
  if (!is.null(x$latent_scores) && ncol(x$latent_scores) >= 2L) {
    scatter_panel(
      x$latent_scores, x$latent_score_obs, x$latent_nn,
      sprintf("LOO k-NN distance (full %d-D latent)", x$n_latent_dims),
      sprintf("Latent: Z1 vs LOO k-NN (k=%d, %d dims)",
              if (!is.null(x$knn.k)) x$knn.k else 10L, x$n_latent_dims))
  } else {
    plot.new(); title(main = "k-NN scatter -- unavailable")
  }

  # (2,2) Mahalanobis scatter
  if (!is.null(x$latent_scores) && ncol(x$latent_scores) >= 2L) {
    scatter_panel(
      x$latent_scores, x$latent_score_obs, x$latent_mahal,
      "LOO Mahalanobis distance (pooled Sigma)",
      "Latent: Z1 vs LOO Mahalanobis")
  } else {
    plot.new(); title(main = "Mahalanobis scatter -- unavailable")
  }

  # (2,3) Prediction summary text panel
  plot.new(); title(main = "Prediction summary")
  text(0.5, 0.92, sprintf("Predicted class : %s", x$predicted_class),
       cex = 1.0, font = 2)
  text(0.5, 0.82, sprintf("P(class | obs) : %.3f", x$predicted_prob), cex = 0.9)
  text(0.5, 0.74, sprintf("Margin (top1-top2) : %.3f", x$margin), cex = 0.9)
  if (!is.null(x$log10_bf_runner_up) && is.finite(x$log10_bf_runner_up)) {
    bf_lab <- if (!is.null(x$bf_runner_up_label)) x$bf_runner_up_label else "--"
    text(0.5, 0.65, sprintf("log10 BF (vs runner-up) : %.2f  [%s]",
                            x$log10_bf_runner_up, bf_lab),
         cex = 0.9, col = "grey20")
  }
  # Per-class BF list (predicted vs each other)
  if (!is.null(x$bayes_factor)) {
    text(0.5, 0.55, "Bayes factor (pred vs class) :",
         cex = 0.8, col = "grey40")
    y_b <- 0.48
    for (cn in x$class_names) {
      if (cn == x$predicted_class) next
      lbf <- x$log10_bayes_factor[[cn]]
      bf  <- x$bayes_factor[[cn]]
      bf_str <- if (!is.finite(bf)) "Inf" else
                if (bf >= 100) sprintf("%.1f", bf) else
                sprintf("%.2f", bf)
      lbf_str <- if (!is.finite(lbf)) "Inf" else sprintf("%.2f", lbf)
      text(0.5, y_b,
           sprintf("vs %s : BF=%s  log10=%s", cn, bf_str, lbf_str),
           cex = 0.8, col = "grey20")
      y_b <- y_b - 0.07
    }
  }
  text(0.5, 0.20, sprintf("Latent dims : %d", x$n_latent_dims),
       cex = 0.8, col = "grey40")
  text(0.5, 0.12,
       sprintf("Obs in predicted class : %s",
               ifelse(x$obs_in_predicted_class, "YES", "NO")), cex = 0.85)
  text(0.5, 0.05,
       sprintf("Obs in any class : %s",
               ifelse(x$obs_in_any_class, "YES", "NO")), cex = 0.85)
}


#' @export
print.OOD_posttrain_classify <- function(x, ...) {
  .ood.posttrain.classify.print.summary(x)
  invisible(x)
}

#' @export
plot.OOD_posttrain_classify <- function(x, ...) {
  .ood.posttrain.classify.plot(x)
  invisible(NULL)
}

#' @export
summary.OOD_posttrain_classify <- function(object, ...) {
  out <- data.frame(
    class     = object$class_names,
    p_class   = unname(object$prediction$prob_mean[object$class_names]),
    d_obs     = vapply(object$latent_nn, function(e) e$d_obs, numeric(1)),
    p_value   = vapply(object$latent_nn, function(e) e$p_value, numeric(1)),
    pass      = vapply(object$latent_nn,
                       function(e) isTRUE(e$pass), logical(1)),
    predicted = object$class_names == object$predicted_class,
    stringsAsFactors = FALSE
  )
  cat(sprintf("Verdict: %s | predicted=%s (P=%.3f, margin=%.3f)\n\n",
              toupper(object$verdict),
              object$predicted_class,
              object$predicted_prob,
              object$margin))
  print(out, row.names = FALSE)
  invisible(out)
}


# ============================================================================
# OOD.attribute.classify — per-bin attribution explaining classifier OOD
#
# Forensic follow-up to OOD.posttrain.classify: when obs lies outside every
# class's latent cloud, attribution tells us *which* SFS bins drive the
# extrapolation under three target functions per class:
#   (a) centroid — distance from obs to class-k centroid in latent space
#   (b) nearest  — distance from obs to its nearest class-k sim in latent
#   (c) logit    — class-k logit (pre-softmax) at obs
#
# Method: integrated gradients (Sundararajan et al. 2017) in the classifier's
# z-scored augmented input space (raw stats + log1p|stats|). Per-target
# baselines:
#   centroid -> class-k mean (per class)
#   nearest  -> the actual nearest class-k sim (per class)
#   logit    -> global mean across all sims
#
# Aggregation: per-augmented-feature attributions for (raw_j, log1p_j) are
# summed to give a single per-raw-bin attribution, in z-scored input units.
# ============================================================================


# --- Internal: forward through everything except output_layer; return
#     penultimate latent and (optionally) logits (= output_layer applied to
#     latent). Differentiable; intended for autograd. ---
.ood.classify.fwd.both <- function(model, X, output_layer) {
  orig_output <- model$output_layer
  model$output_layer <- torch::nn_identity()
  on.exit(model$output_layer <- orig_output, add = TRUE)
  z <- model(X)                # penultimate (linear), differentiable
  logits <- output_layer(z)    # apply trained head separately
  list(z = z, logits = logits)
}


# --- Internal: integrated-gradients attribution for one (target, class).
#     target_fn(z, logits, k) returns a scalar tensor scoring "OOD-ness"
#     of x given class k. Returns numeric vector length n_aug_features. ---
.ood.classify.ig <- function(model, output_layer, x_obs_t, x_base_t,
                             target_fn, k, n_steps, device) {
  alphas <- (seq_len(n_steps) - 0.5) / n_steps        # midpoint rule
  acc <- torch::torch_zeros_like(x_obs_t)
  delta <- x_obs_t - x_base_t
  for (a in alphas) {
    x <- x_base_t + a * delta
    x$requires_grad_(TRUE)
    out <- .ood.classify.fwd.both(model, x, output_layer)
    loss <- target_fn(out$z, out$logits, k)
    g <- torch::autograd_grad(loss, x, retain_graph = FALSE,
                               create_graph = FALSE)[[1]]
    acc <- acc + g$detach()
  }
  attr_aug <- as.numeric((delta * (acc / n_steps))$cpu())
  attr_aug
}


#' Per-bin Attribution for Classifier OOD (Integrated Gradients)
#'
#' Forensic complement to \code{\link{OOD.posttrain.classify}}: when obs
#' is flagged as OOD (lies outside every class's latent cloud), this
#' function attributes the OOD-ness to individual summary statistics
#' (e.g. SFS bins) using integrated gradients through the trained
#' classifier.
#'
#' Three target functions per class \eqn{k}:
#' \itemize{
#'   \item \strong{centroid}: \eqn{\|z(\mathbf{s}) - \bar z_k\|^2} —
#'         distance to class-\eqn{k} centroid in latent space.
#'         Baseline = class-\eqn{k} mean SFS.
#'   \item \strong{nearest}: \eqn{\|z(\mathbf{s}) - z(\mathbf{s}_k^*)\|^2}
#'         — distance to obs's nearest class-\eqn{k} sim in latent space.
#'         Baseline = that nearest sim.
#'   \item \strong{logit}: classifier's class-\eqn{k} logit (pre-softmax).
#'         Positive attribution = bin pushes prediction toward class
#'         \eqn{k}. Baseline = global mean SFS across all sims.
#' }
#'
#' Attribution is computed in the classifier's z-scored augmented input
#' space (raw stats + log1p|stats|), then the (raw, log1p) pair for each
#' original bin is summed to give a per-bin attribution.
#'
#' @param classifier Object from \code{\link{tune.nn.classify}}.
#' @param observed Numeric vector or 1-row data.frame/matrix.
#' @param reftable Reference table containing classifier features +
#'   label column.
#' @param model_col Label column name.
#' @param classes Character vector — subset of class names to attribute.
#'   NULL = all.
#' @param targets Character vector — subset of \code{c("centroid",
#'   "nearest", "logit")}. Default all three.
#' @param n_steps Integration steps for integrated gradients (default 50).
#' @param subsample Per-class subsample for nearest-sim search and
#'   centroid computation (default NULL = use all sims).
#' @param seed Random seed.
#' @param plot Logical — draw the heatmap.
#' @param verbose Logical — print progress.
#'
#' @return A list of class
#'   \code{c("OOD_attribution_classify", "OOD_diagnostic")} with:
#' \describe{
#'   \item{class_names, n_classes, n_bins}{metadata}
#'   \item{bin_names}{character vector of raw stat names (length n_bins)}
#'   \item{obs}{numeric obs values (length n_bins)}
#'   \item{class_means}{K x n_bins matrix of per-class mean raw stats}
#'   \item{predicted_class}{best class from \code{nn.predict.classify}}
#'   \item{attribution}{named list with one matrix per target — each
#'     \code{K_used x n_bins}, attribution in z-scored input units}
#'   \item{nearest_idx}{(if \code{nearest} requested) named integer of
#'     nearest sim index per class (1-based, into the kept subsample)}
#'   \item{n_steps}{IG steps used}
#' }
#'
#' @section Plot:
#' Three side-by-side heatmaps (one per target), classes on Y, bins on X,
#' diverging blue-white-red color = signed attribution. Predicted class
#' labelled with an asterisk.
#'
#' @seealso \code{\link{OOD.posttrain.classify}} for the OOD test that
#'   motivates running attribution.
#'
#' @export
OOD.attribute.classify <- function(classifier, observed, reftable,
                                    model_col = "model",
                                    classes = NULL,
                                    targets = c("centroid", "nearest", "logit"),
                                    n_steps = 50L,
                                    subsample = NULL,
                                    seed = 42,
                                    plot = TRUE,
                                    verbose = TRUE) {
  if (!inherits(classifier, "model_selection_classifier"))
    stop("classifier must come from tune.nn.classify()")
  if (!model_col %in% colnames(reftable))
    stop(sprintf("OOD.attribute.classify: model_col '%s' not in reftable",
                 model_col))
  if (!requireNamespace("torch", quietly = TRUE))
    stop("OOD.attribute.classify requires the 'torch' package.")

  targets <- match.arg(targets, c("centroid", "nearest", "logit"),
                       several.ok = TRUE)

  feat_cols   <- classifier$feature.cols
  feat_mu     <- classifier$feat_mu
  feat_sd     <- classifier$feat_sd
  K_all       <- classifier$n_classes
  class_names_all <- classifier$class_names
  device      <- if (!is.null(classifier$device)) classifier$device else "cpu"
  dev         <- torch::torch_device(device)
  model       <- classifier$model
  output_layer <- model$output_layer  # capture trained head once
  n_bins      <- length(feat_cols)

  if (is.null(classes)) {
    use_classes <- class_names_all
  } else {
    miss_cls <- setdiff(classes, class_names_all)
    if (length(miss_cls) > 0L)
      stop("classes not in classifier: ", paste(miss_cls, collapse = ", "))
    use_classes <- classes
  }
  K_used <- length(use_classes)

  # --- Resolve labels and (optionally) subsample ---
  raw_labels <- as.character(reftable[[model_col]])
  miss_lbl <- setdiff(unique(raw_labels), class_names_all)
  if (length(miss_lbl) > 0L)
    stop("reftable contains labels not seen by classifier: ",
         paste(miss_lbl, collapse = ", "))
  model_id_int <- match(raw_labels, class_names_all)

  set.seed(seed)
  if (!is.null(subsample)) {
    counts <- table(model_id_int)
    if (any(counts > subsample)) {
      keep_idx <- unlist(lapply(seq_len(K_all), function(k) {
        r <- which(model_id_int == k)
        if (length(r) <= subsample) r else sample(r, subsample)
      }))
      reftable <- reftable[keep_idx, , drop = FALSE]
      model_id_int <- model_id_int[keep_idx]
    }
  }

  # --- Resolve observed ---
  if (is.data.frame(observed) || is.matrix(observed)) {
    if (nrow(observed) > 1L) stop("observed must be a single row")
    obs_named <- setNames(as.numeric(observed[1, ]), colnames(observed))
  } else {
    if (is.null(names(observed)))
      stop("observed must be a named numeric vector")
    obs_named <- setNames(as.numeric(observed), names(observed))
  }
  miss_feat <- setdiff(feat_cols, names(obs_named))
  if (length(miss_feat) > 0L)
    stop("observed missing feature columns: ",
         paste(miss_feat, collapse = ", "))
  obs_raw <- obs_named[feat_cols]

  # --- Build augmented z-scored sim matrix ---
  S_raw <- as.matrix(reftable[, feat_cols, drop = FALSE])
  bad <- apply(S_raw, 1, function(x) any(!is.finite(x)))
  if (any(bad)) {
    S_raw <- S_raw[!bad, , drop = FALSE]
    model_id_int <- model_id_int[!bad]
  }
  S_aug <- cbind(S_raw, log1p(abs(S_raw)))
  S_z   <- sweep(sweep(S_aug, 2, feat_mu, "-"), 2, feat_sd, "/")

  obs_aug <- c(obs_raw, log1p(abs(obs_raw)))
  obs_z   <- (as.numeric(obs_aug) - feat_mu) / feat_sd
  obs_z_mat <- matrix(obs_z, nrow = 1L)

  # --- Class means in raw space (for plotting + centroid baseline) ---
  class_means_raw <- t(vapply(seq_len(K_all), function(k) {
    rows <- which(model_id_int == k)
    if (length(rows) == 0L) rep(NA_real_, n_bins)
    else colMeans(S_raw[rows, , drop = FALSE])
  }, numeric(n_bins)))
  rownames(class_means_raw) <- class_names_all

  # --- Class mean and global mean in aug-z space (baselines) ---
  class_means_z <- t(vapply(seq_len(K_all), function(k) {
    rows <- which(model_id_int == k)
    if (length(rows) == 0L) rep(NA_real_, ncol(S_z))
    else colMeans(S_z[rows, , drop = FALSE])
  }, numeric(ncol(S_z))))
  rownames(class_means_z) <- class_names_all
  global_mean_z <- colMeans(S_z)

  # --- Latent embeddings of all sims + obs (no grad), used for nearest ---
  if (verbose)
    cat(sprintf("PipeMaster:: extracting latent embeddings (%d sims)\n",
                nrow(S_z)))
  Z_sim <- .torch.penultimate(model, S_z, device = dev)
  Z_obs <- .torch.penultimate(model, obs_z_mat, device = dev)

  # --- Per-class centroid in latent + nearest sim per class ---
  z_centroids <- t(vapply(seq_len(K_all), function(k) {
    rows <- which(model_id_int == k)
    if (length(rows) == 0L) rep(NA_real_, ncol(Z_sim))
    else colMeans(Z_sim[rows, , drop = FALSE])
  }, numeric(ncol(Z_sim))))
  rownames(z_centroids) <- class_names_all

  nearest_idx <- integer(K_all)
  names(nearest_idx) <- class_names_all
  if ("nearest" %in% targets) {
    if (verbose) cat("PipeMaster:: locating nearest sim per class\n")
    for (k in seq_len(K_all)) {
      rows <- which(model_id_int == k)
      if (length(rows) == 0L) { nearest_idx[k] <- NA_integer_; next }
      Z_k <- Z_sim[rows, , drop = FALSE]
      nn <- RANN::nn2(Z_k, Z_obs, k = 1L)
      nearest_idx[k] <- rows[nn$nn.idx[1, 1]]
    }
  }

  # --- Predicted class via classifier ---
  obs_df <- as.data.frame(matrix(obs_raw, nrow = 1L,
                                 dimnames = list(NULL, feat_cols)))
  prediction <- nn.predict.classify(classifier, obs_df,
                                    mc_samples = 0L, verbose = FALSE)

  # --- Set up attribution loop ---
  obs_t <- torch::torch_tensor(obs_z_mat, dtype = torch::torch_float(),
                                device = dev)

  attribution <- list()
  for (tg in targets) attribution[[tg]] <-
    matrix(NA_real_, nrow = K_used, ncol = n_bins,
           dimnames = list(use_classes, feat_cols))

  if (verbose)
    cat(sprintf("PipeMaster:: integrated gradients (%d steps, %d classes, %d targets)\n",
                n_steps, K_used, length(targets)))

  for (i in seq_along(use_classes)) {
    cn <- use_classes[i]
    k_global <- match(cn, class_names_all)

    if ("centroid" %in% targets) {
      base_z <- class_means_z[k_global, ]
      base_t <- torch::torch_tensor(matrix(base_z, nrow = 1L),
                                     dtype = torch::torch_float(),
                                     device = dev)
      z_c <- torch::torch_tensor(matrix(z_centroids[k_global, ], nrow = 1L),
                                  dtype = torch::torch_float(),
                                  device = dev)
      tgt_fn <- function(z, logits, k) (z - z_c)$pow(2)$sum()
      attr_aug <- .ood.classify.ig(model, output_layer, obs_t, base_t,
                                    tgt_fn, k_global, n_steps, dev)
      attribution$centroid[i, ] <- attr_aug[seq_len(n_bins)] +
                                    attr_aug[(n_bins + 1L):(2L * n_bins)]
    }

    if ("nearest" %in% targets) {
      idx <- nearest_idx[k_global]
      if (!is.na(idx)) {
        base_z <- S_z[idx, ]
        base_t <- torch::torch_tensor(matrix(base_z, nrow = 1L),
                                       dtype = torch::torch_float(),
                                       device = dev)
        z_n <- torch::torch_tensor(matrix(Z_sim[idx, ], nrow = 1L),
                                    dtype = torch::torch_float(),
                                    device = dev)
        tgt_fn <- function(z, logits, k) (z - z_n)$pow(2)$sum()
        attr_aug <- .ood.classify.ig(model, output_layer, obs_t, base_t,
                                      tgt_fn, k_global, n_steps, dev)
        attribution$nearest[i, ] <- attr_aug[seq_len(n_bins)] +
                                     attr_aug[(n_bins + 1L):(2L * n_bins)]
      }
    }

    if ("logit" %in% targets) {
      base_t <- torch::torch_tensor(matrix(global_mean_z, nrow = 1L),
                                     dtype = torch::torch_float(),
                                     device = dev)
      kk <- k_global
      tgt_fn <- function(z, logits, k) logits[1, kk]
      attr_aug <- .ood.classify.ig(model, output_layer, obs_t, base_t,
                                    tgt_fn, k_global, n_steps, dev)
      attribution$logit[i, ] <- attr_aug[seq_len(n_bins)] +
                                 attr_aug[(n_bins + 1L):(2L * n_bins)]
    }

    if (verbose)
      cat(sprintf("  [%d/%d] %s done\n", i, K_used, cn))
  }

  results <- list(
    class_names     = use_classes,
    n_classes       = K_used,
    n_bins          = n_bins,
    bin_names       = feat_cols,
    obs             = obs_raw,
    class_means     = class_means_raw[use_classes, , drop = FALSE],
    predicted_class = prediction$best_model,
    prediction      = prediction,
    attribution     = attribution,
    nearest_idx     = if ("nearest" %in% targets) nearest_idx else NULL,
    n_steps         = as.integer(n_steps),
    targets         = targets
  )
  class(results) <- c("OOD_attribution_classify", "OOD_diagnostic")

  if (verbose) .ood.attribute.classify.print.summary(results)
  if (plot)    .ood.attribute.classify.plot(results)

  invisible(results)
}


# --- Internal: print summary listing top-3 bins per target, signed ---
.ood.attribute.classify.print.summary <- function(x) {
  cat("\nOOD.attribute.classify summary\n")
  cat(sprintf("  Predicted class : %s\n", x$predicted_class))
  cat(sprintf("  IG steps        : %d\n", x$n_steps))
  cat(sprintf("  Bins            : %d (%s)\n", x$n_bins,
              paste(utils::head(x$bin_names, 3), collapse = ", ")))
  cat("\n  Top-3 |attribution| per (target, class):\n")
  for (tg in x$targets) {
    cat(sprintf("\n  %s:\n", toupper(tg)))
    M <- x$attribution[[tg]]
    for (cn in rownames(M)) {
      a <- M[cn, ]
      if (all(is.na(a))) next
      ord <- order(abs(a), decreasing = TRUE)[seq_len(min(3L, length(a)))]
      mark <- if (cn == x$predicted_class) " <- predicted" else ""
      txt <- paste(sprintf("%s=%+.3f", x$bin_names[ord], a[ord]),
                   collapse = "  ")
      cat(sprintf("    %-15s %s%s\n", cn, txt, mark))
    }
  }
  cat("\n")
}


# --- Internal: 1x3 heatmap layout (one heatmap per target) ---
.ood.attribute.classify.plot <- function(x) {
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  ntg <- length(x$targets)
  par(mfrow = c(1, ntg), mar = c(5, 6, 3, 5))

  # Diverging palette (blue-white-red)
  ramp <- grDevices::colorRampPalette(c("#2166ac", "#67a9cf", "#f7f7f7",
                                         "#ef8a62", "#b2182b"))
  n_col <- 101L
  pal <- ramp(n_col)

  for (tg in x$targets) {
    M <- x$attribution[[tg]]
    if (all(is.na(M))) {
      plot.new(); title(main = sprintf("%s -- unavailable", tg)); next
    }
    rng <- max(abs(M), na.rm = TRUE)
    if (!is.finite(rng) || rng == 0) rng <- 1
    breaks <- seq(-rng, rng, length.out = n_col + 1L)

    image(seq_len(x$n_bins), seq_len(x$n_classes),
          t(M), col = pal, breaks = breaks,
          xlab = "", ylab = "", axes = FALSE,
          main = sprintf("%s attribution\n(units: z-scored input)", tg))
    axis(1, at = seq_len(x$n_bins), labels = x$bin_names,
         las = 2, cex.axis = 0.7)
    cls_lab <- ifelse(x$class_names == x$predicted_class,
                      paste0(x$class_names, "*"),
                      x$class_names)
    axis(2, at = seq_len(x$n_classes), labels = cls_lab,
         las = 1, cex.axis = 0.8)
    box()

    # color legend
    usr <- par("usr")
    leg_x <- usr[2] + diff(usr[1:2]) * 0.04
    leg_w <- diff(usr[1:2]) * 0.03
    leg_y <- seq(usr[3], usr[4], length.out = n_col + 1L)
    rect(leg_x, leg_y[-(n_col + 1L)],
         leg_x + leg_w, leg_y[-1L],
         col = pal, border = NA, xpd = NA)
    text(leg_x + leg_w * 1.5, leg_y[1L],
         sprintf("%+.2f", -rng), pos = 4, cex = 0.7, xpd = NA)
    text(leg_x + leg_w * 1.5, leg_y[n_col + 1L],
         sprintf("%+.2f", rng), pos = 4, cex = 0.7, xpd = NA)
    text(leg_x + leg_w * 1.5, mean(leg_y),
         "0", pos = 4, cex = 0.7, xpd = NA)
  }
}


#' @export
print.OOD_attribution_classify <- function(x, ...) {
  .ood.attribute.classify.print.summary(x)
  invisible(x)
}

#' @export
plot.OOD_attribution_classify <- function(x, ...) {
  .ood.attribute.classify.plot(x)
  invisible(NULL)
}

#' @export
summary.OOD_attribution_classify <- function(object, ...) {
  out <- do.call(rbind, lapply(object$targets, function(tg) {
    M <- object$attribution[[tg]]
    if (all(is.na(M))) return(NULL)
    apply(M, 1, function(a) {
      ord <- order(abs(a), decreasing = TRUE)[1L]
      sprintf("%s=%+.3f", object$bin_names[ord], a[ord])
    })
  }))
  if (!is.null(out)) {
    rownames(out) <- toupper(object$targets[seq_len(nrow(out))])
    cat(sprintf("Predicted class: %s | targets: %s | bins: %d | IG steps: %d\n\n",
                object$predicted_class,
                paste(object$targets, collapse = ", "),
                object$n_bins, object$n_steps))
    cat("Top |attribution| bin per (target, class):\n")
    print(out, quote = FALSE)
  }
  invisible(out)
}

