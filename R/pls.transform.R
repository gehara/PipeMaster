#' PLS Dimensionality Reduction for ABC
#'
#' Fits a Partial Least Squares (PLS) regression of parameters on summary
#' statistics and returns the projection. PLS components maximize covariance
#' between stats and parameters, keeping the most informative dimensions
#' for inference. Implements NIPALS algorithm (no external dependencies).
#'
#' @param stats Numeric matrix (n x p) — summary statistics (predictors).
#' @param params Numeric matrix (n x q) — parameters (response). Used to
#'   find the projection directions. Not needed for \code{pls.project()}.
#' @param n.comp Integer — number of PLS components to retain (default 20).
#' @param scale Logical — center and scale stats before PLS (default TRUE).
#' @param max.rows Integer — maximum rows for NIPALS fitting (default 10000).
#'   When \code{nrow(stats) > max.rows}, a random subsample is used for
#'   fitting. Centering/scaling parameters are computed from the full data.
#'   This avoids excessive memory use on large reftables.
#'
#' @return A list with class \code{"pls.fit"}:
#' \describe{
#'   \item{projection}{matrix (p x n.comp) — projection from stat space to PLS space}
#'   \item{center}{numeric (p) — column means used for centering}
#'   \item{scale}{numeric (p) — column SDs used for scaling (1 if scale=FALSE)}
#'   \item{n.comp}{integer — number of components retained}
#'   \item{scores}{matrix (n x n.comp) — training scores (stats projected)}
#'   \item{var.explained}{numeric (n.comp) — fraction of X variance per component}
#' }
#'
#' @details
#' Uses the NIPALS (Nonlinear Iterative Partial Least Squares) algorithm.
#' For each component: find direction w in X-space that maximizes
#' cov(Xw, Y), project X onto w, deflate X and Y, repeat.
#'
#' Columns with zero variance are dropped before fitting.
#'
#' @examples
#' \dontrun{
#' # Fit PLS on reftable
#' pls_fit <- pls.fit(stats = as.matrix(reftable[, stat_cols]),
#'                    params = as.matrix(reftable[, param_cols]),
#'                    n.comp = 20)
#'
#' # Project reftable and observed into PLS space
#' stats_pls <- pls.project(pls_fit, as.matrix(reftable[, stat_cols]))
#' obs_pls   <- pls.project(pls_fit, matrix(obs_raw, nrow = 1))
#'
#' # Use in ABC
#' reftable_pls <- cbind(reftable[, param_cols], as.data.frame(stats_pls))
#' abc_result <- abc.rejection(reftable_pls, obs_pls, param.cols = param_cols)
#' }
#'
#' @export
pls.fit <- function(stats, params, n.comp = 20L, scale = TRUE, max.rows = 10000L) {

  stats  <- as.matrix(stats)
  params <- as.matrix(params)
  n <- nrow(stats)
  p <- ncol(stats)
  q <- ncol(params)

  if (nrow(params) != n)
    stop("stats and params must have the same number of rows")

  # Center and scale from full data
  col_means <- colMeans(stats, na.rm = TRUE)
  if (scale) {
    col_sds <- apply(stats, 2, sd, na.rm = TRUE)
    col_sds[col_sds == 0 | !is.finite(col_sds)] <- 1
  } else {
    col_sds <- rep(1, p)
  }

  # Subsample for NIPALS fitting if too many rows
  if (n > max.rows) {
    sub_idx <- sample(n, max.rows)
    stats_fit  <- stats[sub_idx, , drop = FALSE]
    params_fit <- params[sub_idx, , drop = FALSE]
    n_fit <- max.rows
    cat(sprintf("PipeMaster:: pls.fit: subsampling %d/%d rows for fitting\n", max.rows, n))
  } else {
    stats_fit  <- stats
    params_fit <- params
    n_fit <- n
  }

  # Center and scale using full-data parameters
  X <- t((t(stats_fit) - col_means) / col_sds)
  X[!is.finite(X)] <- 0

  # Drop zero-variance columns
  x_var <- colSums(X^2)
  keep <- x_var > 1e-12
  if (sum(keep) < n.comp) {
    n.comp <- sum(keep)
    if (n.comp == 0) stop("All stat columns have zero variance")
  }
  X <- X[, keep, drop = FALSE]

  # Center params
  Y_means <- colMeans(params_fit, na.rm = TRUE)
  Y <- t(t(params_fit) - Y_means)
  Y[!is.finite(Y)] <- 0

  # Clamp n.comp
  n.comp <- min(n.comp, ncol(X), n_fit - 1L)

  # NIPALS PLS2
  W <- matrix(0, nrow = ncol(X), ncol = n.comp)  # weights
  P <- matrix(0, nrow = ncol(X), ncol = n.comp)  # X loadings
  T_mat <- matrix(0, nrow = n_fit, ncol = n.comp)  # X scores
  C <- matrix(0, nrow = q, ncol = n.comp)  # Y loadings (per-component, per-param)
  total_var <- sum(X^2)

  for (a in seq_len(n.comp)) {
    # Initial u = first column of Y (or largest singular vector)
    u <- Y[, 1, drop = TRUE]

    for (iter in 1:100) {
      # X weight: w = X'u / u'u
      w <- crossprod(X, u)
      w <- w / sqrt(sum(w^2))

      # X score: t = Xw
      t_new <- X %*% w

      # Y loading: c = Y't / t't
      cc <- crossprod(Y, t_new) / sum(t_new^2)

      # Y score: u = Yc / c'c
      u_new <- Y %*% cc / sum(cc^2)

      # Check convergence
      if (sum((u_new - u)^2) < 1e-10 * sum(u_new^2 + 1e-20)) break
      u <- u_new
    }

    # X loading: p = X't / t't
    p_load <- crossprod(X, t_new) / sum(t_new^2)

    # Store
    W[, a] <- w
    P[, a] <- p_load
    T_mat[, a] <- t_new
    C[, a] <- as.numeric(cc)

    # Deflate X and Y
    X <- X - t_new %*% t(p_load)
    Y <- Y - t_new %*% t(cc)
  }

  # Build projection matrix: maps centered/scaled stats to PLS scores
  # R = W (P'W)^{-1}, so scores = X_centered_scaled %*% R
  PtW <- crossprod(P, W)
  R <- W %*% solve(PtW)

  # Expand R back to full column space (re-insert zeros for dropped cols)
  R_full <- matrix(0, nrow = sum(keep) + sum(!keep), ncol = n.comp)
  R_full[keep, ] <- R

  # Variance explained per component
  var_comp <- colSums(T_mat^2)
  var_explained <- var_comp / total_var

  # Y-loadings per component, with param-name rownames so callers can map
  # each PLS component back to the params it predicts. In normalized
  # (NIPALS-internal) units — sign and relative magnitude across params at
  # a given component are interpretable; absolute scale is not.
  Y_names <- colnames(params)
  if (is.null(Y_names)) Y_names <- paste0("Y", seq_len(q))
  rownames(C) <- Y_names
  colnames(C) <- paste0("PLS", seq_len(n.comp))

  result <- list(
    projection    = R_full,
    center        = col_means,
    scale         = col_sds,
    n.comp        = n.comp,
    scores        = T_mat,
    var.explained = var_explained,
    Y_loadings    = C,
    Y_names       = Y_names,
    Y_means       = Y_means
  )
  class(result) <- "pls.fit"

  cat(sprintf("PipeMaster:: pls.fit: %d stats -> %d PLS components (%.1f%% variance)\n",
              p, n.comp, 100 * sum(var_explained)))

  result
}


#' Project new data into PLS component space
#'
#' @param pls.fit A \code{pls.fit} object from \code{pls.fit()}.
#' @param newdata Numeric matrix or vector — new summary statistics to project.
#'   Must have the same columns as the training stats.
#'
#' @return Numeric matrix (n x n.comp) of PLS scores, or named vector if
#'   input was a single observation.
#'
#' @export
pls.project <- function(pls.fit, newdata) {

  single <- is.null(dim(newdata))
  if (single) newdata <- matrix(newdata, nrow = 1)
  newdata <- as.matrix(newdata)

  # Center and scale using training parameters
  X <- t((t(newdata) - pls.fit$center) / pls.fit$scale)
  X[!is.finite(X)] <- 0

  # Project
  scores <- X %*% pls.fit$projection
  colnames(scores) <- paste0("PLS", seq_len(pls.fit$n.comp))

  if (single) return(as.numeric(scores))
  scores
}
