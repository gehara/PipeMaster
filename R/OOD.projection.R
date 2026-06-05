# ============================================================================
# OOD Projection Forensics + Per-Outlier Sim Distributions
#
# OOD.projection.diagnose() — class-dispatched on OOD_pretrain / OOD_posttrain.
#   For OOD_pretrain : basis in {"pca", "pls"}.
#   For OOD_posttrain: basis in {"pca", "pls", "nn_latent"}; default nn_latent.
#
# OOD.outliers() — small-multiples plot of the sim distribution for each
#   outlier stat with the observed value overlaid.
# ============================================================================


#' Diagnose why the "obs projected (outliers neutralized)" arrow lands where
#' it does in the chosen projection basis.
#'
#' Runs three complementary checks on the chosen basis:
#' \enumerate{
#'   \item Full-component NN distance: does the projection bring obs closer
#'     to the sim cloud when all retained components are considered, or only
#'     in the first two? Cumulative distance vs number of components shown.
#'   \item Per-component obs extremity: absolute z-score of obs (and
#'     projected obs) on each component (score / sim component sd). Reveals
#'     extremity along components beyond the first two that are not visible
#'     in scatter plots.
#'   \item Top outlier-stat loadings on components 1 and 2 (PCA / PLS only):
#'     which specific outlier stats are pulling obs along the first two
#'     components, driving the visible arrow direction when neutralized.
#'     Skipped for the NN-latent basis (the mapping is nonlinear; loadings
#'     do not exist).
#' }
#'
#' Class-dispatched on the input:
#' \itemize{
#'   \item \code{OOD_pretrain} accepts \code{basis} in \code{"pca"} or
#'     \code{"pls"} (default \code{"pca"}).
#'   \item \code{OOD_posttrain} accepts \code{"pca"}, \code{"pls"}, or
#'     \code{"nn_latent"} (default \code{"nn_latent"} since that is the
#'     post-tier-native basis).
#' }
#'
#' @param ood_result output of \code{OOD.pretrain()} or \code{OOD.posttrain()}.
#' @param basis character — projection basis to diagnose. See class dispatch.
#' @param pdf_file character or NULL — write the plot to this PDF.
#' @param top_n integer — how many top stats to show in the loading plots.
#' @param pdf_width,pdf_height numeric — PDF dimensions if \code{pdf_file} given.
#'
#' @return Invisibly a list with:
#' \describe{
#'   \item{nn_obs_full}{Full-component NN distance from obs to nearest sim}
#'   \item{nn_proj_full}{Full-component NN distance from projected obs to
#'     nearest sim}
#'   \item{cumulative_distance}{matrix [n_comp x 2] — cumulative distance
#'     to nearest sim as components are added one by one}
#'   \item{obs_pc_z, proj_pc_z}{|score| / sim component sd per component}
#'   \item{loadings_pc1, loadings_pc2}{outlier-stat loadings on components
#'     1 and 2, sorted by |loading|. NULL for NN-latent basis.}
#'   \item{basis}{the basis that was used}
#' }
#' @export
OOD.projection.diagnose <- function(ood_result, basis = NULL,
                                    pdf_file = NULL, top_n = 15,
                                    pdf_width = 12, pdf_height = 10) {
  if (!inherits(ood_result, "OOD_diagnostic"))
    stop("ood_result must be an OOD_pretrain or OOD_posttrain object.")

  is_post <- inherits(ood_result, "OOD_posttrain")
  valid_bases <- if (is_post) c("nn_latent", "pls", "pca")
                 else c("pls", "pca")
  if (is.null(basis)) basis <- valid_bases[1]
  basis <- match.arg(basis, valid_bases)

  # Resolve the basis result and metadata
  br <- .ood.proj.resolve.basis(ood_result, basis)
  if (is.null(br$scores))
    stop(sprintf("ood_result$%s has no scores — basis unavailable.", br$slot))
  if (is.null(br$score_obs_projected))
    stop("No projected obs in result — no outlier stats, nothing to diagnose.")

  scores     <- br$scores
  score_obs  <- as.numeric(br$score_obs[1, ])
  score_proj <- as.numeric(br$score_obs_projected[1, ])
  n_comp     <- ncol(scores)

  # --- (A) Full-component NN distance: obs vs projected ---
  nn_obs_full  <- min(sqrt(rowSums(sweep(scores, 2, score_obs,  "-")^2)))
  nn_proj_full <- min(sqrt(rowSums(sweep(scores, 2, score_proj, "-")^2)))

  cumdist <- matrix(NA_real_, nrow = n_comp, ncol = 2,
                    dimnames = list(NULL, c("obs", "projected")))
  for (k in seq_len(n_comp)) {
    cumdist[k, 1] <- min(sqrt(rowSums(
      sweep(scores[, 1:k, drop = FALSE], 2, score_obs[1:k],  "-")^2)))
    cumdist[k, 2] <- min(sqrt(rowSums(
      sweep(scores[, 1:k, drop = FALSE], 2, score_proj[1:k], "-")^2)))
  }

  # --- (B) per-component obs extremity ---
  comp_sd <- apply(scores, 2, sd)
  obs_pc_z  <- abs(score_obs)  / comp_sd
  proj_pc_z <- abs(score_proj) / comp_sd

  # --- (C) loadings: stat (PCA/PLS) + param (PLS only) ---
  load_pc1 <- NULL; load_pc2 <- NULL
  param_load_pc1 <- NULL; param_load_pc2 <- NULL
  if (basis %in% c("pca", "pls")) {
    loads <- .ood.proj.loadings(ood_result, basis, br, top_n)
    load_pc1 <- loads$pc1; load_pc2 <- loads$pc2
  }
  if (basis == "pls") {
    pload <- .ood.proj.param.loadings(br, top_n)
    param_load_pc1 <- pload$pc1
    param_load_pc2 <- pload$pc2
  }

  # --- plotting ---
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  # PLS basis: 2x3 with stat + param loadings. Other bases: 2x2.
  if (basis == "pls") {
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  } else {
    par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  }

  comp_label <- switch(basis, pca = "PC", pls = "PLS", nn_latent = "dim")

  # Panel: cumulative NN distance vs # components
  plot(seq_len(n_comp), cumdist[, 1], type = "b", col = "red", pch = 16,
       ylim = range(cumdist, na.rm = TRUE),
       xlab = sprintf("number of %ss included", comp_label),
       ylab = "NN distance to nearest sim",
       main = sprintf("Full-%s NN: obs=%.3f  projected=%.3f  (%s basis)",
                      comp_label, nn_obs_full, nn_proj_full, basis))
  lines(seq_len(n_comp), cumdist[, 2], type = "b", col = "orange", pch = 2)
  legend("bottomright",
         legend = c("obs", "obs projected (outliers neutralized)"),
         col = c("red", "orange"), pch = c(16, 2), lty = 1, cex = 0.8)

  # PLS-only: stat loadings on PLS1 (top row col 2)
  if (basis == "pls") {
    .ood.proj.plot.loadings(load_pc1, comp_label, 1, kind = "stat")
    .ood.proj.plot.loadings(param_load_pc1, comp_label, 1, kind = "param")
  }

  # Panel: per-component obs z-score (bottom row col 1 for PLS, top right for others)
  bar_mat <- rbind(obs = obs_pc_z, projected = proj_pc_z)
  barplot(bar_mat, beside = TRUE,
          names.arg = paste0(comp_label, seq_len(n_comp)),
          col = c("red", "orange"), border = NA,
          main = sprintf("Per-%s obs extremity (|score| / sim sd)  [n=%d]",
                         comp_label, n_comp),
          ylab = "z-score", las = 1, cex.names = 0.8)
  abline(h = c(2, 3), lty = 2, col = "grey60")
  legend("topright", legend = c("obs", "projected"),
         fill = c("red", "orange"), cex = 0.8, bty = "n")

  if (basis == "pls") {
    # bottom row cols 2-3: stat loadings PLS2 + param loadings PLS2
    .ood.proj.plot.loadings(load_pc2, comp_label, 2, kind = "stat")
    .ood.proj.plot.loadings(param_load_pc2, comp_label, 2, kind = "param")
  } else if (basis == "pca") {
    # bottom row cols 1-2: stat loadings PC1, PC2
    .ood.proj.plot.loadings(load_pc1, comp_label, 1, kind = "stat")
    .ood.proj.plot.loadings(load_pc2, comp_label, 2, kind = "stat")
  } else {
    # NN-latent: no loadings exist
    plot.new(); title(main = "Loadings: not defined for NN-latent")
    plot.new(); title(main = "(nonlinear basis — use OOD.outliers())")
  }

  invisible(list(
    nn_obs_full = nn_obs_full, nn_proj_full = nn_proj_full,
    cumulative_distance = cumdist,
    obs_pc_z = obs_pc_z, proj_pc_z = proj_pc_z,
    loadings_pc1 = load_pc1, loadings_pc2 = load_pc2,
    param_loadings_pc1 = param_load_pc1,
    param_loadings_pc2 = param_load_pc2,
    basis = basis
  ))
}


# ---------------------------------------------------------------------------
# Internal: resolve which slot to read for the chosen basis
# ---------------------------------------------------------------------------

.ood.proj.resolve.basis <- function(ood_result, basis) {
  slot <- switch(basis,
                 pca       = "pca_all",
                 pls       = "pls_all",
                 nn_latent = "nn_latent")
  br <- ood_result[[slot]]
  if (is.null(br))
    stop(sprintf("ood_result$%s is NULL — basis '%s' unavailable.",
                 slot, basis))
  list(
    slot                = slot,
    scores              = br$scores,
    score_obs           = br$score_obs,
    score_obs_projected = br$score_obs_projected,
    keep_mask           = br$keep_mask,
    pca_fit             = br$pca_fit,
    pls_fit             = br$pls_fit
  )
}


# ---------------------------------------------------------------------------
# Internal: outlier-stat loadings on components 1 and 2 for PCA / PLS
# ---------------------------------------------------------------------------

.ood.proj.loadings <- function(ood_result, basis, br, top_n) {
  outlier_stats <- ood_result$percentiles$stat[
    ood_result$percentiles$reason == "outlier"]

  # Recover the stat names that correspond to the rows of the loading matrix.
  # For PCA: pca_fit$rotation has rownames = stat names used.
  # For PLS: pls_fit$projection has nrow = ncol(input stats); we use the
  # keep_mask + ood_result$.context$sim_cols.
  if (basis == "pca") {
    rotation <- br$pca_fit$rotation
    stat_names_basis <- rownames(rotation)
  } else {
    proj_mat <- br$pls_fit$projection
    sim_cols <- ood_result$.context$sim_cols
    stat_names_basis <- sim_cols[br$keep_mask]
    rotation <- proj_mat
    rownames(rotation) <- stat_names_basis
  }

  outlier_in_basis <- intersect(outlier_stats, stat_names_basis)
  if (length(outlier_in_basis) == 0) {
    return(list(pc1 = data.frame(stat = character(0), loading = numeric(0)),
                pc2 = data.frame(stat = character(0), loading = numeric(0))))
  }

  get_load_df <- function(comp_idx) {
    loads <- rotation[outlier_in_basis, comp_idx]
    df <- data.frame(stat = names(loads), loading = as.numeric(loads),
                     stringsAsFactors = FALSE)
    df <- df[order(-abs(df$loading)), ]
    head(df, top_n)
  }
  list(pc1 = get_load_df(1), pc2 = get_load_df(2))
}

.ood.proj.plot.loadings <- function(load_df, comp_label, comp_idx,
                                    kind = c("stat", "param")) {
  kind <- match.arg(kind)
  if (is.null(load_df) || nrow(load_df) == 0) {
    plot.new()
    title(main = sprintf("No %s loadings on %s%d", kind, comp_label, comp_idx))
    return(invisible(NULL))
  }
  par(mar = c(4, 9, 3, 1))
  load_plot <- load_df[rev(seq_len(nrow(load_df))), ]
  label_col <- if (!is.null(load_plot$stat)) load_plot$stat else load_plot$param
  bar_col <- if (kind == "stat")
    ifelse(load_plot$loading > 0, "firebrick", "royalblue")
  else
    ifelse(load_plot$loading > 0, "darkgreen", "darkorange3")
  title_prefix <- if (kind == "stat") "outlier-stat" else "param"
  barplot(load_plot$loading,
          names.arg = label_col, horiz = TRUE, las = 1,
          col = bar_col, border = NA, cex.names = 0.7,
          main = sprintf("Top %d %s loadings on %s%d",
                         nrow(load_df), title_prefix, comp_label, comp_idx),
          xlab = sprintf("%s%d loading", comp_label, comp_idx))
  abline(v = 0, col = "grey60")
  par(mar = c(4, 4, 3, 1))
}


# Param Y-loadings on each PLS component (requires extended pls.fit).
# Returns sorted top-N data.frames for components 1 and 2.
.ood.proj.param.loadings <- function(br, top_n) {
  pls_fit <- br$pls_fit
  if (is.null(pls_fit) || is.null(pls_fit$Y_loadings)) {
    return(list(pc1 = NULL, pc2 = NULL))
  }
  Y <- pls_fit$Y_loadings
  param_names <- pls_fit$Y_names
  if (is.null(param_names)) param_names <- rownames(Y)
  if (is.null(param_names))
    param_names <- paste0("Y", seq_len(nrow(Y)))

  get_load_df <- function(comp_idx) {
    if (comp_idx > ncol(Y)) return(NULL)
    loads <- Y[, comp_idx]
    df <- data.frame(param = param_names, loading = as.numeric(loads),
                     stringsAsFactors = FALSE)
    df <- df[order(-abs(df$loading)), ]
    head(df, top_n)
  }
  list(pc1 = get_load_df(1), pc2 = get_load_df(2))
}


# ---------------------------------------------------------------------------
# OOD.priors.bestfit: prior-bound pressure on best-fit ABC posterior
# ---------------------------------------------------------------------------

#' Prior Boundary Pressure on Best-Fit (ABC Posterior) Samples
#'
#' Given the parameter values of the simulations that fit obs best (the
#' rejection-ABC posterior approximation), reports where those values sit
#' within the prior range and flags boundary pressure. The pressed
#' direction is the actionable answer to "which prior bound is the data
#' pushing against?"
#'
#' Two input modes (mutually exclusive):
#' \enumerate{
#'   \item \code{ood_result =} (in-pipeline): the function performs K-NN
#'     on the cached PLS / PCA scores from \code{OOD.pretrain()}. Mechanically
#'     equivalent to rejection ABC with PLS reduction and tolerance
#'     \code{K_frac}, but reuses the projection that was already fit.
#'   \item \code{posterior =} (external): the function consumes any
#'     posterior object with an \code{$abc} matrix of accepted parameter
#'     samples — e.g., from \code{nn.predict(method = "ABC_NN_regression")},
#'     or any future ABC/NPE method that returns class \code{"posterior"}.
#'     The \code{reftable} is required to extract empirical prior bounds.
#' }
#'
#' Prior bounds are inferred from the empirical \code{min}/\code{max} of
#' each param column (assumes uniform priors — the PipeMaster default).
#'
#' Boundary pressure is flagged when the median of the accepted samples
#' sits within \code{edge_threshold * (prior_hi - prior_lo)} of either
#' bound.
#'
#' Pairs naturally with \code{OOD.projection.diagnose(basis = "pls")} (the
#' Y-loadings panel says which param a stat-direction wants) and with
#' \code{OOD.priors.outliers()} (per-outlier-stat marginal voting).
#' Agreement across diagnostics = strong signal to widen.
#'
#' @param ood_result output of \code{OOD.pretrain()} (or \code{OOD.posttrain()}).
#'   Mutually exclusive with \code{posterior}.
#' @param posterior object of class \code{"posterior"} (e.g., from
#'   \code{nn.predict(method = "ABC_NN_regression")}). Mutually exclusive
#'   with \code{ood_result}.
#' @param reftable data.frame — required when \code{posterior} is given,
#'   for extracting empirical prior bounds.
#' @param K_frac numeric — fraction of sims to keep as nearest neighbors
#'   when \code{ood_result} is given (default 0.025). Ignored in
#'   \code{posterior} mode.
#' @param edge_threshold numeric — fraction of prior range considered
#'   "near boundary" (default 0.05).
#' @param basis character — projection basis for the in-pipeline K-NN.
#'   \code{"pls"} (default) is param-aligned; \code{"pca"} is
#'   variance-aligned. Ignored in \code{posterior} mode.
#' @param plot logical — produce a per-param grid showing the prior range
#'   with accepted samples overlaid as a boxplot (default TRUE).
#' @param pdf_file character or NULL — write the plot to this PDF.
#' @param pdf_width,pdf_height numeric — PDF dimensions if \code{pdf_file}.
#'
#' @return A list with:
#' \describe{
#'   \item{table}{data.frame: \code{param, prior_lo, prior_hi,
#'     neighbor_lo, neighbor_q25, neighbor_med, neighbor_q75, neighbor_hi,
#'     edge_pressure, med_dist_lo, med_dist_hi}}
#'   \item{method}{character — short description of how \code{accepted}
#'     was obtained}
#'   \item{K}{integer — number of accepted samples used}
#' }
#' @export
OOD.priors.bestfit <- function(ood_result = NULL,
                                posterior = NULL,
                                reftable = NULL,
                                K_frac = 0.025,
                                edge_threshold = 0.05,
                                basis = c("pls", "pca"),
                                plot = TRUE,
                                pdf_file = NULL,
                                pdf_width = 12, pdf_height = 9) {
  if (is.null(ood_result) && is.null(posterior))
    stop("Provide either ood_result or posterior.")
  if (!is.null(ood_result) && !is.null(posterior))
    stop("Provide only one of ood_result or posterior, not both.")

  if (!is.null(posterior)) {
    # Mode B: external posterior
    if (!inherits(posterior, "posterior"))
      stop("posterior must be of class 'posterior' (e.g., from nn.predict(method = 'ABC_NN_regression')).")
    if (is.null(posterior$abc))
      stop("posterior must contain $abc (matrix of accepted samples).")
    if (is.null(reftable))
      stop("reftable required when using posterior= (for prior bounds).")
    accepted <- as.matrix(posterior$abc)
    param_names <- if (!is.null(posterior$param_names))
                     posterior$param_names else colnames(accepted)
    if (is.null(param_names))
      stop("Cannot determine param names from posterior.")
    missing_p <- setdiff(param_names, colnames(reftable))
    if (length(missing_p) > 0)
      stop(sprintf("Param columns missing from reftable: %s",
                   paste(missing_p, collapse = ", ")))
    prior_full <- as.matrix(reftable[, param_names, drop = FALSE])
    method_label <- sprintf("user-supplied posterior (%d accepted)",
                            nrow(accepted))
    K <- nrow(accepted)
    n_sims <- nrow(prior_full)

  } else {
    # Mode A: OOD-integrated K-NN on cached projection
    if (!inherits(ood_result, "OOD_diagnostic"))
      stop("ood_result must be an OOD_pretrain or OOD_posttrain object.")
    basis <- match.arg(basis)
    P_sim <- ood_result$.context$P_sim
    if (is.null(P_sim))
      stop("ood_result must have been called with param.cols.")
    br <- if (basis == "pls" && !is.null(ood_result$pls_all))
            list(scores = ood_result$pls_all$scores,
                 score_obs = ood_result$pls_all$score_obs,
                 method = "PLS")
          else
            list(scores = ood_result$pca_all$scores,
                 score_obs = ood_result$pca_all$score_obs,
                 method = "PCA")
    if (is.null(br$scores))
      stop(sprintf("Basis '%s' has no scores in ood_result.", basis))
    scores    <- br$scores
    score_obs <- as.numeric(br$score_obs[1, ])
    n_sims <- nrow(scores)
    K <- min(max(1L, ceiling(K_frac * n_sims)), n_sims)
    d <- sqrt(rowSums(sweep(scores, 2, score_obs, "-")^2))
    nn_idx <- order(d)[seq_len(K)]
    accepted   <- P_sim[nn_idx, , drop = FALSE]
    prior_full <- P_sim
    param_names <- colnames(P_sim)
    method_label <- sprintf("internal K-NN on cached %s scores (K=%d / %d)",
                            br$method, K, n_sims)
  }

  # --- Common downstream: boundary diagnostic ---
  prior_lo <- apply(prior_full, 2, min)
  prior_hi <- apply(prior_full, 2, max)
  prior_range <- prior_hi - prior_lo

  q <- apply(accepted, 2, function(x)
    quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1), names = FALSE))
  med <- q[3, ]
  med_dist_lo <- (med - prior_lo) / pmax(prior_range, .Machine$double.eps)
  med_dist_hi <- (prior_hi - med) / pmax(prior_range, .Machine$double.eps)

  edge_pressure <- ifelse(med_dist_lo < edge_threshold &
                          med_dist_hi < edge_threshold, "BOTH",
                  ifelse(med_dist_lo < edge_threshold, "LO",
                  ifelse(med_dist_hi < edge_threshold, "HI", "ok")))

  table <- data.frame(
    param = param_names,
    prior_lo = prior_lo, prior_hi = prior_hi,
    neighbor_lo  = q[1, ], neighbor_q25 = q[2, ],
    neighbor_med = q[3, ], neighbor_q75 = q[4, ],
    neighbor_hi  = q[5, ],
    edge_pressure = edge_pressure,
    med_dist_lo = round(med_dist_lo, 4),
    med_dist_hi = round(med_dist_hi, 4),
    stringsAsFactors = FALSE
  )
  rownames(table) <- NULL

  cat(sprintf("PipeMaster:: OOD.priors.bestfit (%s, edge=%.0f%%)\n",
              method_label, 100 * edge_threshold))
  pressed <- table[table$edge_pressure != "ok", ]
  if (nrow(pressed) > 0) {
    cat(sprintf("  Boundary pressure on %d/%d params:\n",
                nrow(pressed), nrow(table)))
    for (r in seq_len(nrow(pressed)))
      cat(sprintf("    %-18s [%s] med=%.4g  prior=[%.4g, %.4g]\n",
                  pressed$param[r], pressed$edge_pressure[r],
                  pressed$neighbor_med[r],
                  pressed$prior_lo[r], pressed$prior_hi[r]))
  } else {
    cat("  No boundary pressure detected.\n")
  }

  if (plot) {
    .ood.priors.bestfit.plot(table, accepted, method_label, K,
                              edge_threshold,
                              pdf_file = pdf_file,
                              pdf_width = pdf_width,
                              pdf_height = pdf_height)
  }

  invisible(list(table = table, method = method_label, K = K))
}


# Per-param grid: prior range bar + accepted-sample boxplot, colored by pressure.
.ood.priors.bestfit.plot <- function(out, accepted, method_label, K,
                                     edge_threshold,
                                     pdf_file = NULL,
                                     pdf_width = 12, pdf_height = 9) {
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  n_params <- nrow(out)
  n_col <- min(4, max(1, ceiling(sqrt(n_params))))
  n_row <- ceiling(n_params / n_col)
  par(mfrow = c(n_row, n_col), mar = c(3, 4, 3, 1),
      oma = c(0, 0, 2, 0))

  pressure_col <- function(p) switch(p,
    LO = "royalblue", HI = "firebrick",
    BOTH = "purple", ok = "grey60")

  for (i in seq_len(n_params)) {
    p_lo <- out$prior_lo[i]; p_hi <- out$prior_hi[i]
    nb_vals <- accepted[, i]
    pr <- out$edge_pressure[i]
    col <- pressure_col(pr)

    plot.new()
    plot.window(xlim = c(p_lo, p_hi),
                ylim = c(0.3, 1.3), xaxs = "i")
    axis(1, cex.axis = 0.8)
    title(main = sprintf("%s [%s]", out$param[i], pr),
          col.main = col, cex.main = 1.0)

    # Prior range bar
    rect(p_lo, 1.05, p_hi, 1.20, col = "grey85", border = "grey60")
    # Edge zones (within edge_threshold of bounds), translucent
    edge_w <- (p_hi - p_lo) * edge_threshold
    rect(p_lo, 1.05, p_lo + edge_w, 1.20,
         col = rgb(0.25, 0.41, 0.88, 0.25), border = NA)
    rect(p_hi - edge_w, 1.05, p_hi, 1.20,
         col = rgb(0.70, 0.13, 0.13, 0.25), border = NA)
    # Accepted-sample boxplot below
    boxplot(nb_vals, horizontal = TRUE, add = TRUE,
            at = 0.6, boxwex = 0.4, axes = FALSE,
            col = col, border = "grey20", outline = TRUE,
            outcex = 0.4)
    points(out$neighbor_med[i], 1.125, pch = 18, col = col, cex = 1.2)
  }

  mtext(sprintf("OOD.priors.bestfit  (%s, edge=%.0f%%)",
                method_label, 100 * edge_threshold),
        outer = TRUE, cex = 1.1)
}


# ---------------------------------------------------------------------------
# OOD.priors.outliers: per-outlier-stat marginal prior-change voting
# ---------------------------------------------------------------------------

#' Suggest Prior Changes from Per-Outlier-Stat Marginal Correlations
#'
#' For each outlier statistic (obs strictly outside sim support), finds the
#' top-K most-correlated parameters in the reftable and translates each
#' \code{(corr sign x outlier direction)} pair into a vote to widen that
#' parameter's HI or LO bound. Aggregates votes into a per-parameter
#' consensus table — the actionable answer to "if I want to keep these
#' outlier stats, which priors should I widen and in which direction?"
#'
#' Complementary to \code{OOD.priors.bestfit()}, which uses joint K-NN-of-obs
#' (rejection ABC) boundary pressure (one signal averaged over all stats).
#' This function is per-outlier-stat marginal — different signal, different
#' failure mode. Use both.
#'
#' Caveat: marginal correlations miss param interactions. Strong when one
#' param dominantly drives a stat, weak when two params jointly determine
#' it (e.g., bottleneck depth depends on Ne0/Ne1 ratio). Treat as
#' first-order guidance, not gospel.
#'
#' @param ood_result output of \code{OOD.pretrain()} (or \code{OOD.posttrain()})
#'   that was called with \code{param.cols} so \code{P_sim} was stashed.
#' @param top_k integer — number of top-correlated params per outlier stat
#'   that get a vote (default 3).
#' @param corr_threshold numeric — minimum \code{|corr|} for a param to vote
#'   (default 0.2). Params with weaker signal than this are ignored even if
#'   in the top-k.
#' @param ambiguous_threshold numeric — when \code{|net| / total <
#'   ambiguous_threshold}, the consensus is reported as \code{"ambiguous"}
#'   rather than HI/LO (default 0.3).
#' @param plot logical — produce a per-param visual (default TRUE):
#'   side-by-side diverging "net" bars (HI − LO) and "total" bars.
#' @param pdf_file character or NULL — write the plot to this PDF.
#' @param pdf_width,pdf_height numeric — PDF dimensions if \code{pdf_file}.
#'
#' @return A list with:
#' \describe{
#'   \item{consensus}{data.frame, one row per param implicated by any
#'     outlier: \code{param, n_HI, n_LO, evidence_HI, evidence_LO,
#'     consensus, net, total}, sorted by \code{total} descending}
#'   \item{per_outlier}{data.frame, one row per (outlier, param) vote:
#'     \code{outlier_stat, obs_dir, param, corr, vote}}
#' }
#'
#' @export
OOD.priors.outliers <- function(ood_result, top_k = 3L,
                               corr_threshold = 0.2,
                               ambiguous_threshold = 0.3,
                               plot = TRUE,
                               pdf_file = NULL,
                               pdf_width = 12, pdf_height = 8) {
  if (!inherits(ood_result, "OOD_diagnostic"))
    stop("ood_result must be an OOD_pretrain or OOD_posttrain object.")

  ctx <- ood_result$.context
  P_sim <- ctx$P_sim
  if (is.null(P_sim))
    stop("OOD.priors.outliers requires param.cols to have been passed to OOD.pretrain.")
  S_sim    <- ctx$S_sim
  observed <- ctx$observed
  sim_cols <- ctx$sim_cols
  outlier_mask <- ctx$outlier_mask

  if (is.null(outlier_mask) || !any(outlier_mask)) {
    cat("PipeMaster:: OOD.priors.outliers: no outlier stats — nothing to vote on.\n")
    return(invisible(list(
      consensus = data.frame(),
      per_outlier = data.frame())))
  }

  outlier_idx <- which(outlier_mask)
  n_params <- ncol(P_sim)
  param_names <- colnames(P_sim)
  if (is.null(param_names)) param_names <- paste0("p", seq_len(n_params))

  # Per-outlier voting
  votes <- list()
  for (k in outlier_idx) {
    s_vals <- S_sim[, k]
    sim_min <- min(s_vals); sim_max <- max(s_vals)
    obs_v <- observed[k]
    direction <- if (obs_v > sim_max) 1L
                 else if (obs_v < sim_min) -1L
                 else next  # shouldn't happen given outlier_mask
    obs_dir_lab <- if (direction > 0) "above_max" else "below_min"

    # Correlations of this stat with each param across the reftable
    corrs <- suppressWarnings(
      vapply(seq_len(n_params),
             function(p) cor(P_sim[, p], s_vals, use = "complete.obs"),
             numeric(1)))
    corrs[!is.finite(corrs)] <- 0

    keep <- which(abs(corrs) >= corr_threshold)
    if (length(keep) == 0) next
    keep <- keep[order(-abs(corrs[keep]))][seq_len(min(top_k, length(keep)))]

    for (p in keep) {
      eff <- sign(corrs[p]) * direction
      vote <- if (eff > 0) "widen_HI" else "widen_LO"
      votes[[length(votes) + 1L]] <- data.frame(
        outlier_stat = sim_cols[k],
        obs_dir      = obs_dir_lab,
        param        = param_names[p],
        corr         = round(corrs[p], 4),
        vote         = vote,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(votes) == 0) {
    cat(sprintf("PipeMaster:: OOD.priors.outliers: no params with |corr| >= %.2f for any outlier.\n",
                corr_threshold))
    return(invisible(list(
      consensus = data.frame(),
      per_outlier = data.frame())))
  }
  per_outlier <- do.call(rbind, votes)

  # Aggregate per param
  agg <- lapply(unique(per_outlier$param), function(p) {
    d <- per_outlier[per_outlier$param == p, ]
    n_HI <- sum(d$vote == "widen_HI")
    n_LO <- sum(d$vote == "widen_LO")
    e_HI <- sum(abs(d$corr[d$vote == "widen_HI"]))
    e_LO <- sum(abs(d$corr[d$vote == "widen_LO"]))
    net   <- n_HI - n_LO
    total <- n_HI + n_LO
    cons <- if (total == 0) "no_signal"
            else if (abs(net) / total < ambiguous_threshold) "ambiguous"
            else if (net > 0) "widen_HI" else "widen_LO"
    data.frame(
      param = p, n_HI = n_HI, n_LO = n_LO,
      evidence_HI = round(e_HI, 4), evidence_LO = round(e_LO, 4),
      consensus = cons, net = net, total = total,
      stringsAsFactors = FALSE
    )
  })
  consensus <- do.call(rbind, agg)
  consensus <- consensus[order(-consensus$total, -abs(consensus$net)), ]
  rownames(consensus) <- NULL

  # Stdout summary
  n_outliers <- length(unique(per_outlier$outlier_stat))
  cat(sprintf("PipeMaster:: OOD.priors.outliers (top_k=%d, |corr|>=%.2f, ambig<%.0f%%)\n",
              top_k, corr_threshold, 100 * ambiguous_threshold))
  cat(sprintf("  %d outlier stats voted on %d params (%d total votes)\n",
              n_outliers, nrow(consensus), nrow(per_outlier)))
  for (r in seq_len(nrow(consensus))) {
    cat(sprintf("    %-18s [%-10s] HI=%d LO=%d  (evidence HI=%.2f, LO=%.2f)\n",
                consensus$param[r], consensus$consensus[r],
                consensus$n_HI[r], consensus$n_LO[r],
                consensus$evidence_HI[r], consensus$evidence_LO[r]))
  }

  if (plot) {
    .ood.priors.outliers.plot(consensus, ambiguous_threshold,
                               pdf_file = pdf_file,
                               pdf_width = pdf_width,
                               pdf_height = pdf_height)
  }

  invisible(list(consensus = consensus, per_outlier = per_outlier))
}


# Side-by-side panels: diverging net (HI - LO) bars + total bars.
.ood.priors.outliers.plot <- function(consensus, ambiguous_threshold,
                                     pdf_file = NULL,
                                     pdf_width = 12, pdf_height = 8) {
  if (nrow(consensus) == 0) {
    plot.new(); title(main = "No prior-change suggestions")
    return(invisible(NULL))
  }
  if (!is.null(pdf_file)) {
    pdf(pdf_file, width = pdf_width, height = pdf_height)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4, 9, 3, 1), oma = c(0, 0, 0, 0))

  # Plot from bottom (largest total) to top (smallest), so reverse for
  # horizontal barplot which fills bottom-to-top
  d <- consensus[rev(seq_len(nrow(consensus))), ]
  cons_col <- function(c) switch(c,
    widen_HI  = "firebrick",
    widen_LO  = "royalblue",
    ambiguous = "grey55",
    no_signal = "grey80")
  bar_cols <- vapply(d$consensus, cons_col, character(1))

  # Panel 1: net (HI - LO), diverging
  net_max <- max(abs(d$net), 1)
  barplot(d$net, names.arg = d$param, horiz = TRUE, las = 1,
          col = bar_cols, border = NA, cex.names = 0.7,
          xlim = c(-net_max, net_max),
          xlab = "net votes  (n_HI - n_LO)",
          main = sprintf("Prior-change net votes  (ambig if |net|/total < %.0f%%)",
                         100 * ambiguous_threshold))
  abline(v = 0, col = "grey30")

  # Panel 2: total votes (single direction)
  barplot(d$total, names.arg = d$param, horiz = TRUE, las = 1,
          col = bar_cols, border = NA, cex.names = 0.7,
          xlab = "total votes  (n_HI + n_LO)",
          main = "Total outliers implicating this param")
  legend("topright", bty = "n", cex = 0.8,
         legend = c("widen HI", "widen LO", "ambiguous", "no signal"),
         fill = c("firebrick", "royalblue", "grey55", "grey80"))
}


# ---------------------------------------------------------------------------
# OOD.outliers: per-outlier-stat sim-distribution histograms
# ---------------------------------------------------------------------------

#' Plot simulated distributions for OOD outlier statistics
#'
#' Small-multiples diagnostic: one histogram per outlier statistic showing the
#' simulation distribution with the observed value overlaid as a vertical line.
#' Complements the meta-histogram in \code{OOD.pretrain()} /
#' \code{OOD.posttrain()} by letting the user inspect each problematic stat
#' in detail.
#'
#' Outliers are sorted by extremity (\code{max(percentile, 1 - percentile)})
#' and truncated to \code{max_stats}.
#'
#' @param ood_result output of \code{OOD.pretrain()} or \code{OOD.posttrain()}.
#' @param observed numeric vector or 1-row data.frame of observed stats, same
#'   as passed to \code{OOD.pretrain()} / \code{OOD.posttrain()}.
#' @param reftable data.frame of simulated stats (the reference table).
#' @param max_stats integer — maximum number of outlier stats to plot
#'   (default 24). Most extreme first.
#' @param ncol integer — number of histograms per row in the grid (default 4).
#' @param pdf_file character or NULL — if provided, writes the plot to this
#'   PDF file. Otherwise plots to the active device.
#' @param pdf_width,pdf_height numeric — PDF dimensions if \code{pdf_file} given.
#'
#' @return Invisibly returns the outlier data.frame (sorted by extremity, truncated).
#' @export
OOD.outliers <- function(ood_result, observed, reftable,
                         max_stats = 24, ncol = 4,
                         pdf_file = NULL,
                         pdf_width = 12, pdf_height = 9) {
  if (!inherits(ood_result, "OOD_diagnostic"))
    stop("ood_result must be an OOD_pretrain or OOD_posttrain object.")

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
