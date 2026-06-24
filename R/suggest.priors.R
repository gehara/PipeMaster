#' Suggest prior bounds from observed summary statistics
#'
#' Derives point estimates for current Ne, ancestral Ne, divergence time,
#' and migration from classical population-genetic moment estimators, then
#' writes wide prior bounds around those estimates into the supplied Model
#' object. Bounds are \code{est / expand_factor} to \code{est * expand_factor}.
#'
#' Estimators used:
#' \itemize{
#'   \item \strong{Ne0.pop\eqn{p}}: average of \eqn{Ne_\pi = \pi/(4\mu)} and
#'     \eqn{Ne_W = \theta_W/(4\mu)} using per-site values
#'     (\code{s_mean_pi_p / locus_length} etc.).
#'   \item \strong{Ne1.pop\eqn{p}}: same bounds as Ne0.pop\eqn{p}.
#'     The function does NOT call a demographic syndrome (expansion / bottleneck
#'     / constant) from Tajima's D -- \eqn{|D| < 0.5} is too noisy to call,
#'     and even \eqn{|D| > 0.5} can be selection, structure, or other
#'     departures. The raw TajD value is reported in the suggestions object
#'     for the user to inspect.
#'   \item \strong{Ne.anc_i_j}: \code{mean(Ne_i, Ne_j)}.
#'   \item \strong{t.Ne1.pop\eqn{p}}: heuristic \code{[1, expand_factor * Ne0]}.
#'   \item \strong{join_i_j}: \code{T = -2 * Ne_avg * log(1 - Fst_i_j)}
#'     (pure-isolation reading). Used as an upper bound when migration is
#'     also in the model.
#'   \item \strong{mig_i_j}: \code{Nm = (1 - Fst_i_j) / (4 * Fst_i_j)}
#'     (drift-migration equilibrium reading). Used as an upper bound when
#'     divergence is also in the model.
#' }
#'
#' Moment estimators are time-averaged: under expansion (Ne0 >> Ne1),
#' \eqn{\pi / (4\mu)} sits closer to Ne1 than Ne0. Treat the suggested
#' Ne bounds as a starting point -- if the classifier or posterior pushes
#' against the upper bound, widen.
#'
#' When \code{model$flags$ej} and \code{model$flags$m} are both present,
#' the Fst-derived T and Nm cannot be disentangled from Fst alone. The
#' function reports both readings, uses the ratio
#' \code{fixed_i_j / (fixed_i_j + shared_i_j)} as a soft tiebreaker
#' (> 0.5 -> ancient/isolation tilt, < 0.1 -> recent/migration tilt),
#' and emits a per-pair warning in the printed summary.
#'
#' This function does NOT modify \code{model$conds}, \code{model$tree},
#' \code{model$labels}, or \code{model$sum_anc_ne}. Distribution column
#' (column 6 of each flag matrix) is left unchanged.
#'
#' @param model A PipeMaster Model object (class "Model").
#' @param obs Single-row data.frame or named numeric vector of observed
#'   summary statistics, in the format returned by \code{observed.sumstats()}
#'   (columns prefixed \code{s_mean_}).
#' @param mu Per-site per-generation mutation rate.
#' @param locus_length Locus length in bp used to convert per-locus pi /
#'   thetaW to per-site values. If NULL, taken as \code{mean(model$loci[,2])}.
#' @param expand_factor Multiplicative half-width for prior bounds.
#'   \code{lo = est / expand_factor}, \code{hi = est * expand_factor}.
#'   Default 10.
#' @param verbose Print the suggestion summary to the console. Default TRUE.
#' @return The model object with \code{model$flags} lo/hi columns rewritten,
#'   and a \code{model$prior_suggestions} list attached recording the point
#'   estimates, per-pair Fst readings, and any warnings.
#' @export
suggest.priors <- function(model, obs, mu,
                           locus_length = NULL,
                           expand_factor = 10,
                           verbose = TRUE) {

  if (!inherits(model, "Model"))
    warning("model does not inherit from class 'Model'; proceeding anyway.")
  if (!is.numeric(mu) || length(mu) != 1 || !is.finite(mu) || mu <= 0)
    stop("mu must be a single positive finite number")
  if (!is.numeric(expand_factor) || expand_factor <= 1)
    stop("expand_factor must be > 1")

  ## --- coerce obs to a 1-row data.frame ---
  if (is.data.frame(obs)) {
    if (nrow(obs) > 1) obs <- obs[1, , drop = FALSE]
  } else if (is.matrix(obs)) {
    obs <- as.data.frame(obs[1, , drop = FALSE], stringsAsFactors = FALSE)
  } else {
    obs <- as.data.frame(as.list(obs), stringsAsFactors = FALSE)
  }

  ## --- locus length ---
  if (is.null(locus_length)) {
    bp <- suppressWarnings(as.numeric(model$loci[, 2]))
    bp <- bp[is.finite(bp) & bp > 0]
    if (length(bp) == 0)
      stop("Could not infer locus_length from model$loci; pass it explicitly.")
    locus_length <- mean(bp)
  }
  if (!is.finite(locus_length) || locus_length <= 0)
    stop("locus_length must be a positive finite number")

  pull <- function(name) {
    v <- obs[[name]]
    if (is.null(v)) return(NA_real_)
    suppressWarnings(as.numeric(v))
  }

  npops <- as.integer(model$I[1, 3])
  nsam_pop <- suppressWarnings(as.integer(model$I[1, 3 + seq_len(npops)]))
  nsam_pop[!is.finite(nsam_pop) | nsam_pop <= 0] <- 2L

  mig_scale <- if (is.null(model$mig_scale)) "Nm" else model$mig_scale
  if (!mig_scale %in% c("Nm", "m"))
    warning("Unknown mig_scale '", mig_scale, "'; treating as 'Nm'.")

  ## --- per-pop Ne0 from pi and thetaW ---
  ne_est       <- rep(NA_real_, npops)
  ne_est_pi    <- rep(NA_real_, npops)
  ne_est_w     <- rep(NA_real_, npops)
  hd_obs       <- rep(NA_real_, npops)
  for (p in seq_len(npops)) {
    pi_loc      <- pull(sprintf("s_mean_pi_%d", p))
    thetaW_loc  <- pull(sprintf("s_mean_thetaW_%d", p))
    hd_obs[p]   <- pull(sprintf("s_mean_Hd_%d", p))

    pi_per_site     <- pi_loc / locus_length
    thetaW_per_site <- thetaW_loc / locus_length
    ne_est_pi[p] <- pi_per_site     / (4 * mu)
    ne_est_w[p]  <- thetaW_per_site / (4 * mu)

    ests <- c(ne_est_pi[p], ne_est_w[p])
    ests <- ests[is.finite(ests) & ests > 0]
    if (length(ests) > 0) ne_est[p] <- mean(ests)
  }

  ## --- per-pop Tajima's D (reported only, no syndrome call) ---
  tajd_obs <- rep(NA_real_, npops)
  for (p in seq_len(npops)) {
    tajd_obs[p] <- pull(sprintf("s_mean_tajd_%d", p))
  }

  ## --- pairwise Fst-driven T_iso and Nm_eq ---
  pair_estimates <- list()
  if (npops >= 2) {
    for (i in 1:(npops - 1)) for (j in (i + 1):npops) {
      fst <- pull(sprintf("s_mean_Fst_%d_%d", i, j))
      if (!is.finite(fst) || fst <= 0 || fst >= 1) next

      fixed  <- pull(sprintf("s_mean_fixed_%d_%d",  i, j))
      shared <- pull(sprintf("s_mean_shared_%d_%d", i, j))

      ne_ij <- mean(c(ne_est[i], ne_est[j]), na.rm = TRUE)
      if (!is.finite(ne_ij)) next

      T_iso <- -2 * ne_ij * log(1 - fst)
      Nm_eq <- (1 - fst) / (4 * fst)

      r_fs <- NA_real_
      tilt <- "balanced"
      T_tilt  <- 1.0
      Nm_tilt <- 1.0
      if (is.finite(fixed) && is.finite(shared) && (fixed + shared) > 0) {
        r_fs <- fixed / (fixed + shared)
        if (r_fs > 0.5) {
          Nm_tilt <- 0.5
          tilt <- sprintf("ancient/isolation (r=%.2f, Nm hi shrunk 0.5x)", r_fs)
        } else if (r_fs < 0.1) {
          T_tilt <- 0.5
          tilt <- sprintf("recent/migration (r=%.2f, T hi shrunk 0.5x)", r_fs)
        } else {
          tilt <- sprintf("genuinely ambiguous (r=%.2f)", r_fs)
        }
      }

      pair_estimates[[paste0(i, "_", j)]] <- list(
        i = i, j = j, fst = fst,
        fixed = fixed, shared = shared, r_fs = r_fs,
        T_iso = T_iso, Nm_eq = Nm_eq, ne_avg = ne_ij,
        T_tilt = T_tilt, Nm_tilt = Nm_tilt, tilt = tilt
      )
    }
  }

  ## --- helper to rewrite a flag matrix row by name ---
  changed <- character(0)
  rewrite <- function(mat, key, lo, hi) {
    if (is.null(mat) || nrow(mat) == 0) return(mat)
    idx <- which(mat[, 1] == key)
    if (length(idx) == 0) return(mat)
    if (is.finite(lo) && lo >= 0) mat[idx, 4] <- format(lo, scientific = FALSE)
    if (is.finite(hi) && hi >= 0) mat[idx, 5] <- format(hi, scientific = FALSE)
    changed <<- c(changed, key)
    mat
  }

  ## --- Ne0 ---
  for (p in seq_len(npops)) {
    if (!is.finite(ne_est[p])) next
    floor_ne <- max(nsam_pop[p] * 2, 50)
    lo <- max(floor_ne, ne_est[p] / expand_factor)
    hi <- ne_est[p] * expand_factor
    model$flags$n <- rewrite(model$flags$n, sprintf("Ne0.pop%d", p), lo, hi)
  }

  ## --- Ne1 (same bounds as Ne0; no syndrome assumption) ---
  if (!is.null(model$flags$en) && !is.null(model$flags$en$size)) {
    for (p in seq_len(npops)) {
      if (!is.finite(ne_est[p])) next
      floor_ne <- max(nsam_pop[p] * 2, 50)
      lo <- max(floor_ne, ne_est[p] / expand_factor)
      hi <- ne_est[p] * expand_factor
      model$flags$en$size <- rewrite(model$flags$en$size,
                                     sprintf("Ne1.pop%d", p), lo, hi)
    }
  }

  ## --- t.Ne1 (heuristic) ---
  if (!is.null(model$flags$en) && !is.null(model$flags$en$time)) {
    for (p in seq_len(npops)) {
      if (!is.finite(ne_est[p])) next
      lo <- 1
      hi <- ne_est[p] * expand_factor
      model$flags$en$time <- rewrite(model$flags$en$time,
                                     sprintf("t.Ne1.pop%d", p), lo, hi)
    }
  }

  ## --- Ne.anc_i_j (average daughters' Ne) ---
  if (!is.null(model$flags$en) && !is.null(model$flags$en$size)) {
    anc_rows <- grep("^Ne\\.anc_", model$flags$en$size[, 1])
    for (k in anc_rows) {
      nm   <- model$flags$en$size[k, 1]
      parts <- strsplit(sub("^Ne\\.anc_", "", nm), "_")[[1]]
      if (length(parts) != 2) next
      i <- suppressWarnings(as.integer(parts[1]))
      j <- suppressWarnings(as.integer(parts[2]))
      if (!is.finite(i) || !is.finite(j) || i > npops || j > npops) next
      ne_anc <- mean(c(ne_est[i], ne_est[j]), na.rm = TRUE)
      if (!is.finite(ne_anc)) next
      floor_ne <- max(nsam_pop[i] * 2, nsam_pop[j] * 2, 50)
      lo <- max(floor_ne, ne_anc / expand_factor)
      hi <- ne_anc * expand_factor
      model$flags$en$size[k, 4] <- format(lo, scientific = FALSE)
      model$flags$en$size[k, 5] <- format(hi, scientific = FALSE)
      changed <- c(changed, nm)
    }
  }

  ## --- Fst-derived join times and migrants ---
  warnings_list <- character(0)
  has_ej <- !is.null(model$flags$ej) && nrow(model$flags$ej) > 0
  has_m  <- !is.null(model$flags$m)  && nrow(model$flags$m)  > 0
  ambiguous <- has_ej && has_m

  if (length(pair_estimates) > 0) {
    for (key in names(pair_estimates)) {
      pe <- pair_estimates[[key]]
      i <- pe$i; j <- pe$j

      ## join time: search ej rows for "i j" or "j i"
      if (has_ej) {
        for (k in seq_len(nrow(model$flags$ej))) {
          parts <- strsplit(model$flags$ej[k, 3], " ")[[1]]
          src <- suppressWarnings(as.integer(parts[1]))
          dst <- suppressWarnings(as.integer(parts[2]))
          if (!is.finite(src) || !is.finite(dst)) next
          if ((src == i && dst == j) || (src == j && dst == i)) {
            lo <- 1
            hi <- pe$T_iso * expand_factor * pe$T_tilt
            model$flags$ej[k, 4] <- format(lo, scientific = FALSE)
            model$flags$ej[k, 5] <- format(hi, scientific = FALSE)
            changed <- c(changed, model$flags$ej[k, 1])

            ## sync t.Ne.anc time row, if present
            if (!is.null(model$flags$en$time)) {
              ant <- which(grepl(sprintf("^t\\.Ne\\.anc_%d_%d$|^t\\.Ne\\.anc_%d_%d$",
                                          i, j, j, i),
                                  model$flags$en$time[, 1]))
              for (a in ant) {
                model$flags$en$time[a, 4] <- format(lo, scientific = FALSE)
                model$flags$en$time[a, 5] <- format(hi, scientific = FALSE)
                changed <- c(changed, model$flags$en$time[a, 1])
              }
            }
          }
        }
      }

      ## migration: search m rows for "i j" or "j i", set both directions
      if (has_m) {
        for (k in seq_len(nrow(model$flags$m))) {
          parts <- strsplit(model$flags$m[k, 3], " ")[[1]]
          src <- suppressWarnings(as.integer(parts[1]))
          dst <- suppressWarnings(as.integer(parts[2]))
          if (!is.finite(src) || !is.finite(dst)) next
          if ((src == i && dst == j) || (src == j && dst == i)) {
            lo <- 0
            hi <- pe$Nm_eq * expand_factor * pe$Nm_tilt
            if (mig_scale == "m") {
              ## convert Nm hi to per-gen rate using receiving pop Ne
              ne_recv <- ne_est[dst]
              if (is.finite(ne_recv) && ne_recv > 0) hi <- hi / ne_recv
            }
            model$flags$m[k, 4] <- format(lo, scientific = FALSE)
            model$flags$m[k, 5] <- format(hi, scientific = FALSE)
            changed <- c(changed, model$flags$m[k, 1])
          }
        }
      }

      if (ambiguous) {
        msg <- sprintf(
          paste0("Pair (%d,%d) Fst=%.3f: isolation reading T=%.0f gen | ",
                 "migration reading Nm=%.3f\n",
                 "    tilt: %s\n",
                 "    WARNING: a single Fst value cannot separate T from Nm; ",
                 "both priors set wide."),
          i, j, pe$fst, pe$T_iso, pe$Nm_eq, pe$tilt)
        warnings_list <- c(warnings_list, msg)
      }
    }
  }

  ## --- inconsistency checks ---
  for (p in seq_len(npops)) {
    if (!is.finite(ne_est_pi[p]) || !is.finite(ne_est_w[p])) next
    avg <- mean(c(ne_est_pi[p], ne_est_w[p]))
    rel <- abs(ne_est_pi[p] - ne_est_w[p]) / avg
    if (rel > 0.5) {
      warnings_list <- c(warnings_list, sprintf(
        "pop %d: |Ne_pi - Ne_W| / mean = %.2f (>0.5) -> data shows non-stationary signal; widen Ne0 prior or trust TajD-tilted Ne1.",
        p, rel))
    }
  }

  ## --- assemble suggestions object ---
  suggestions <- list(
    mu = mu, locus_length = locus_length, expand_factor = expand_factor,
    npops = npops, mig_scale = mig_scale,
    ne_est = ne_est, ne_est_pi = ne_est_pi, ne_est_W = ne_est_w,
    tajd_obs = tajd_obs, hd_obs = hd_obs,
    pair_estimates = pair_estimates,
    ambiguous_T_Nm = ambiguous,
    changed_parameters = unique(changed),
    warnings = warnings_list
  )
  model$prior_suggestions <- suggestions

  if (verbose) {
    cat(sprintf("PipeMaster:: suggest.priors  mu=%g  locus_length=%g  expand_factor=%g  npops=%d  mig_scale=%s\n",
                mu, locus_length, expand_factor, npops, mig_scale))
    cat("  NOTE: these are STARTING-POINT priors from classical moment estimators.\n")
    cat("        Moment estimators (pi/(4mu), thetaW/(4mu), Fst-based T and Nm)\n")
    cat("        assume stationarity/equilibrium and are time-averaged. Under\n")
    cat("        expansion, bottleneck, structure, or selection, they can be\n")
    cat("        biased by an order of magnitude. Do not trust them blindly:\n")
    cat("          - run OOD.pretrain() to check obs lies inside the prior support;\n")
    cat("          - inspect whether the posterior pushes against any prior bound;\n")
    cat("          - widen expand_factor or set bounds manually if it does.\n")
    for (p in seq_len(npops)) {
      cat(sprintf("  pop %d: Ne_pi=%s  Ne_W=%s  Ne_est=%s  Hd=%s  TajD=%s\n",
                  p,
                  format(ne_est_pi[p], digits = 4, scientific = FALSE),
                  format(ne_est_w[p],  digits = 4, scientific = FALSE),
                  format(ne_est[p],    digits = 4, scientific = FALSE),
                  format(hd_obs[p],    digits = 3, scientific = FALSE),
                  format(tajd_obs[p],  digits = 3, scientific = FALSE)))
    }
    if (length(pair_estimates) > 0) {
      cat("  Pairwise (Fst -> T_iso & Nm_eq):\n")
      for (key in names(pair_estimates)) {
        pe <- pair_estimates[[key]]
        cat(sprintf("    %s: Fst=%.4f  T_iso=%.0f gen  Nm_eq=%.3f  r_fixed=%s  tilt=%s\n",
                    key, pe$fst, pe$T_iso, pe$Nm_eq,
                    if (is.finite(pe$r_fs)) sprintf("%.2f", pe$r_fs) else "NA",
                    pe$tilt))
      }
    }
    if (length(unique(changed)) > 0)
      cat(sprintf("  Updated %d prior rows: %s\n",
                  length(unique(changed)),
                  paste(unique(changed), collapse = ", ")))
    if (length(warnings_list) > 0) {
      cat("  WARNINGS:\n")
      for (w in warnings_list) cat("    - ", w, "\n", sep = "")
    }
  }

  invisible(model)
}
