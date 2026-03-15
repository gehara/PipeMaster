# ============================================================================
# ABC Rejection — standalone ABC rejection algorithm
#
# abc.rejection() — accepts closest simulations to observed data
# ============================================================================

#' ABC Rejection Algorithm
#'
#' Performs approximate Bayesian computation (ABC) by rejection sampling.
#' Accepts the closest simulations to the observed data based on Euclidean
#' distance in summary statistic space, optionally normalized by standard
#' deviation, median absolute deviation, or Mahalanobis distance.
#'
#' @param reftable data.frame or character path — reference table from
#'   \code{sim.sumstat()} or \code{sim.sumstats()} containing parameter
#'   and summary statistic columns. If a file path, read with
#'   \code{read.table(header = TRUE)}.
#' @param observed numeric vector or 1-row data.frame — observed summary
#'   statistics.
#' @param param.cols character vector — names of parameter columns.
#' @param tol numeric — tolerance (acceptance fraction, default 0.01).
#' @param distance character — distance metric: \code{"sd"} (default),
#'   \code{"mad"}, or \code{"mahalanobis"}.
#' @param prior data.frame or NULL — prior samples for plotting (optional).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return An object of class \code{"nn.posterior"} with:
#' \describe{
#'   \item{point_estimate}{weighted mean of accepted samples (weight = 1/distance)}
#'   \item{abc_rejection}{matrix of accepted parameter samples}
#'   \item{abc_distance}{distance metric used}
#'   \item{prior}{prior samples (if provided)}
#'   \item{param_names}{parameter column names}
#' }
#'
#' @export
abc.rejection <- function(reftable, observed,
                          param.cols,
                          tol = 0.01,
                          distance = c("sd", "mad", "mahalanobis"),
                          prior = NULL,
                          verbose = TRUE) {

  distance <- match.arg(distance)

  # --- Load from file if path ---
  if (is.character(reftable) && length(reftable) == 1L) {
    if (verbose) cat(sprintf("PipeMaster:: Reading reftable from %s\n", reftable))
    reftable <- read.table(reftable, header = TRUE)
  }

  # --- Split params and stats ---
  nuisance  <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance))
  n_samples <- nrow(reftable)
  n_accept  <- max(1L, floor(n_samples * tol))
  n_stats   <- length(stat_cols)

  theta_orig <- as.matrix(reftable[, param.cols, drop = FALSE])
  pred_orig  <- as.matrix(reftable[, stat_cols, drop = FALSE])

  # --- Prepare observed ---
  if (is.data.frame(observed)) {
    if (all(stat_cols %in% colnames(observed)))
      observed <- as.numeric(observed[1, stat_cols])
    else
      observed <- as.numeric(observed[1, ])
  }
  if (length(observed) != n_stats)
    stop(sprintf("observed has %d elements but reftable has %d stat columns.",
                 length(observed), n_stats))

  if (verbose) {
    cat(sprintf("PipeMaster:: abc.rejection: %s samples, tol=%.4f → keep %s, distance=%s\n",
                format(n_samples, big.mark = ","), tol,
                format(n_accept, big.mark = ","), distance))
  }

  # --- Compute distances ---
  if (distance == "mahalanobis") {
    Sigma_pred <- cov(pred_orig)
    eig <- eigen(Sigma_pred, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) < 1e-8 * max(eig) || min(eig) < 0) {
      ridge <- max(eig) * 1e-6
      Sigma_pred <- Sigma_pred + diag(ridge, n_stats)
    }
    Sigma_inv <- chol2inv(chol(Sigma_pred))
    resid <- t(t(pred_orig) - observed)
    dists <- sqrt(rowSums(resid %*% Sigma_inv * resid))
  } else if (distance == "mad") {
    stat_scale <- apply(pred_orig, 2, function(x) {
      m <- median(x, na.rm = TRUE)
      mad_val <- median(abs(x - m), na.rm = TRUE)
      if (mad_val == 0) sd(x, na.rm = TRUE) else mad_val
    })
    stat_scale[stat_scale == 0 | !is.finite(stat_scale)] <- 1
    resid <- t((t(pred_orig) - observed) / stat_scale)
    dists <- sqrt(rowSums(resid^2))
  } else {
    # SD-normalized (default)
    stat_scale <- apply(pred_orig, 2, sd, na.rm = TRUE)
    stat_scale[stat_scale == 0 | !is.finite(stat_scale)] <- 1
    resid <- t((t(pred_orig) - observed) / stat_scale)
    dists <- sqrt(rowSums(resid^2))
  }

  # --- Keep top tol fraction ---
  accept_idx <- order(dists)[seq_len(n_accept)]
  accepted   <- theta_orig[accept_idx, , drop = FALSE]
  acc_dists  <- dists[accept_idx]

  if (verbose) {
    cat(sprintf("PipeMaster:: Accepted %s samples (tol=%.4f)\n",
                format(n_accept, big.mark = ","), tol))
    cat(sprintf("PipeMaster:: Distance range: [%.4f, %.4f], median=%.4f\n",
                min(acc_dists), max(acc_dists), median(acc_dists)))
  }

  # --- Point estimate: weighted mean of accepted (weight = 1/distance) ---
  w <- 1 / (acc_dists + 1e-12)
  w <- w / sum(w)
  point_est <- colSums(accepted * w)
  names(point_est) <- param.cols

  if (verbose) {
    est_str <- paste(sprintf("%s=%.1f", param.cols, point_est), collapse = " ")
    cat(sprintf("PipeMaster:: Point estimate: %s\n", est_str))
  }

  # --- Prior samples ---
  prior_samples <- NULL
  if (!is.null(prior)) {
    pcols <- intersect(param.cols, colnames(prior))
    pcols <- setdiff(pcols, nuisance)
    if (length(pcols) > 0) {
      prior_samples <- as.matrix(prior[, pcols, drop = FALSE])
      colnames(prior_samples) <- pcols
    }
  }

  # --- Return nn.posterior object ---
  result <- list(
    point_estimate = point_est,
    conformal      = NULL,
    bootstrap      = NULL,
    mc_dropout     = NULL,
    abc_rejection  = accepted,
    abc_distance   = distance,
    prior          = prior_samples,
    param_names    = param.cols
  )
  class(result) <- "nn.posterior"
  result
}
