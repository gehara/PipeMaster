# ============================================================================
# Bayesian Posterior Estimation via Neural Emulator MCMC
#
# emulator.MCMC()          — Python-accelerated multi-chain MH-MCMC with
#                                  emulator as likelihood surrogate
# .estimate.residual.variance() — per-stat residual variance from validation set
# .hpd.interval()               — highest posterior density interval
# .ess()                         — effective sample size via autocorrelation
# .rhat()                        — Gelman-Rubin convergence diagnostic
# ============================================================================

# ============================================================================
# Internal: estimate per-stat residual variance from validation data
# ============================================================================

.estimate.residual.variance <- function(emulator) {
  data      <- emulator$data
  em_model  <- emulator$best_model
  target_mu <- data$target_mu
  target_sd <- data$target_sd

  # Predict on validation set (already Z-scored)
  Y_val_hat_z <- as.matrix(predict(em_model, data$X_val, verbose = 0L))

  # Inverse-transform both to original scale
  Y_val_orig     <- .inv.transform.emulator(data$Y_val, target_mu, target_sd)
  Y_val_hat_orig <- .inv.transform.emulator(Y_val_hat_z, target_mu, target_sd)

  # Per-stat MSE = residual variance
  sigma2 <- colMeans((Y_val_orig - Y_val_hat_orig)^2)
  names(sigma2) <- data$stat_cols
  sigma2
}

# ============================================================================
# Internal: highest posterior density interval
# ============================================================================

.hpd.interval <- function(samples, prob = 0.95) {
  n <- length(samples)
  samples <- sort(samples)
  n_in <- ceiling(prob * n)
  widths <- samples[(n_in + 1):n] - samples[1:(n - n_in)]
  best <- which.min(widths)
  c(lower = samples[best], upper = samples[best + n_in])
}

# ============================================================================
# Internal: effective sample size via autocorrelation
# ============================================================================

.ess <- function(x) {
  n <- length(x)
  if (n < 10L) return(n)
  a <- acf(x, lag.max = min(n - 1L, 500L), plot = FALSE)$acf[-1]
  # Truncate at first negative autocorrelation (Geyer's initial positive rule)
  neg_idx <- which(a < 0)
  if (length(neg_idx) > 0L) a <- a[seq_len(neg_idx[1L] - 1L)]
  if (length(a) == 0L) return(n)
  ess <- n / (1 + 2 * sum(a))
  max(1, ess)
}

# ============================================================================
# Internal: Gelman-Rubin R-hat convergence diagnostic
# ============================================================================

.rhat <- function(chain_list) {
  # chain_list: list of numeric vectors (one per chain), all same length
  m <- length(chain_list)
  if (m < 2L) return(NA_real_)
  n <- length(chain_list[[1L]])
  if (n < 2L) return(NA_real_)

  chain_means <- vapply(chain_list, mean, numeric(1))
  chain_vars  <- vapply(chain_list, var, numeric(1))

  grand_mean <- mean(chain_means)
  B <- n * var(chain_means)   # between-chain variance
  W <- mean(chain_vars)       # within-chain variance

  if (W == 0) return(NA_real_)
  var_hat <- ((n - 1) * W + B) / n
  sqrt(var_hat / W)
}

# ============================================================================
# Main: Bayesian posterior via batched multi-chain Metropolis-Hastings MCMC
# ============================================================================

#' Bayesian Posterior Estimation via Emulator-Assisted MCMC
#'
#' Estimates the full posterior distribution p(theta|S_obs) using Metropolis-Hastings
#' MCMC with the trained neural emulator as a likelihood surrogate. The emulator
#' defines a Gaussian approximate likelihood where the residual variance is estimated
#' from the validation set.
#'
#' The MCMC loop runs entirely in Python (via \code{reticulate::py_run_string()})
#' to eliminate per-iteration R-to-Python overhead. The emulator model is called
#' directly as \code{model(tensor, training=False)} in TensorFlow/numpy, giving
#' ~20-50x speedup over an R-side \code{predict()} loop. Multiple chains are
#' evaluated in a single batched model call per iteration.
#'
#' All sampling is done in log-space (demographic parameters are strictly positive),
#' with a symmetric Gaussian proposal. Adaptive Metropolis tuning during burn-in
#' targets ~23\% acceptance rate.
#'
#' @param emulator list -- output from \code{train.emulator()}.
#' @param observed numeric vector or 1-row data.frame -- observed summary statistics.
#' @param model PipeMaster model object -- for extracting prior bounds via
#'   \code{get.prior.table()}. Either \code{model} or \code{prior.bounds}
#'   must be provided.
#' @param prior.bounds data.frame with columns \code{Parameter}, \code{prior.1},
#'   \code{prior.2} -- alternative to \code{model} for specifying parameter bounds.
#' @param n_samples integer -- total MCMC samples to retain across all chains
#'   (default 10000).
#' @param burnin integer -- burn-in samples to discard per chain (default 2000).
#' @param thin integer -- thinning interval (default 1).
#' @param n_chains integer -- number of parallel MCMC chains (default 10). All
#'   chains are evaluated in a single batched predict call per iteration, so
#'   increasing \code{n_chains} adds minimal wall-time cost while giving more
#'   independent samples and enabling convergence diagnostics (R-hat).
#' @param proposal_sd numeric vector or NULL -- proposal standard deviations per
#'   parameter in LOG space. If NULL, auto-calibrated from prior range.
#' @param sigma2 numeric vector or NULL -- residual variance per stat. If NULL,
#'   estimated from validation data via \code{.estimate.residual.variance()}.
#' @param start numeric vector, matrix, or NULL -- starting point(s) for MCMC
#'   (original scale). If a vector, all chains start from the same point (with
#'   jitter). If a matrix (n_chains x n_params), each row is a chain's start.
#'   If NULL, uses diverse optima from \code{emulator.optimize()}.
#' @param device character -- device for MCMC forward passes: \code{"auto"}
#'   (default) uses the model as-is (same device it was trained on),
#'   \code{"cpu"} clones model weights to CPU before MCMC (faster for small
#'   batches when model was trained on GPU), \code{"gpu"} keeps model on GPU.
#' @param adaptive logical -- if TRUE, use adaptive Metropolis (tune proposal_sd
#'   during burn-in to target ~23\% acceptance rate). Default TRUE.
#' @param verbose logical -- print progress. Default TRUE.
#'
#' @return A list of class \code{"emulator_posterior"} with:
#' \describe{
#'   \item{samples}{matrix (n_samples x n_params) of posterior samples (original scale),
#'     combined from all chains}
#'   \item{log_posterior}{numeric vector of log-posterior values per sample}
#'   \item{acceptance_rate}{overall acceptance rate (mean across chains)}
#'   \item{chain_acceptance}{per-chain acceptance rates}
#'   \item{sigma2}{residual variances used}
#'   \item{prior_bounds}{bounds used (data.frame)}
#'   \item{param_cols}{parameter names}
#'   \item{point_estimate}{posterior mean}
#'   \item{credible_intervals}{95\% HPD intervals per parameter}
#'   \item{ess}{effective sample size per parameter}
#'   \item{n_chains}{number of chains}
#'   \item{chain_samples}{array (samples_per_chain x n_params x n_chains) of
#'     per-chain samples for diagnostics}
#'   \item{rhat}{Gelman-Rubin R-hat per parameter (if n_chains >= 2)}
#' }
#'
#' @export
emulator.MCMC <- function(emulator, observed,
                               model = NULL, prior.bounds = NULL,
                               n_samples = 10000L, burnin = 2000L,
                               thin = 1L, n_chains = 10L,
                               proposal_sd = NULL,
                               sigma2 = NULL, start = NULL,
                               device = "auto",
                               adaptive = TRUE, verbose = TRUE) {

  if (!requireNamespace("keras", quietly = TRUE))
    stop("emulator.MCMC() requires the 'keras' R package.")

  # --- Extract bounds (same pattern as emulator.optimize) ---
  if (!is.null(model)) {
    ptab <- get.prior.table(model)
    idx  <- match(emulator$param_cols, ptab$Parameter)
    if (any(is.na(idx)))
      stop("Parameters not found in model: ",
           paste(emulator$param_cols[is.na(idx)], collapse = ", "))
    bounds <- data.frame(
      Parameter = ptab$Parameter[idx],
      prior.1   = as.numeric(ptab$prior.1[idx]),
      prior.2   = as.numeric(ptab$prior.2[idx]),
      stringsAsFactors = FALSE
    )
  } else if (!is.null(prior.bounds)) {
    bounds <- prior.bounds
  } else {
    stop("Either 'model' or 'prior.bounds' must be provided.")
  }

  # --- Extract emulator components ---
  data       <- emulator$data
  param_cols <- emulator$param_cols
  stat_cols  <- data$stat_cols
  feat_mu    <- data$feat_mu
  feat_sd    <- data$feat_sd
  target_mu  <- data$target_mu
  target_sd  <- data$target_sd
  n_params   <- length(param_cols)
  n_stats    <- length(stat_cols)
  em_model   <- emulator$best_model

  # --- Device selection: optionally clone model to CPU ---
  device <- match.arg(device, c("auto", "cpu", "gpu"))
  if (device == "cpu") {
    if (!requireNamespace("reticulate", quietly = TRUE))
      stop("emulator.MCMC() requires the 'reticulate' R package.")
    if (verbose) cat("PipeMaster:: Cloning model weights to CPU for MCMC\n")
    reticulate::py_run_string("
import tensorflow as tf
def _clone_model_to_cpu(model):
    with tf.device('/CPU:0'):
        cpu_model = tf.keras.models.clone_model(model)
        cpu_model.set_weights(model.get_weights())
    return cpu_model
", convert = FALSE)
    em_model <- reticulate::py$`_clone_model_to_cpu`(em_model)
  }

  lo <- as.numeric(bounds$prior.1)
  hi <- as.numeric(bounds$prior.2)
  bounds_log_lo <- log(lo)   # (n_params,)
  bounds_log_hi <- log(hi)   # (n_params,)

  # --- Prepare observed stats ---
  if (is.data.frame(observed)) {
    if (!is.null(stat_cols) && all(stat_cols %in% colnames(observed)))
      observed <- as.numeric(observed[1, stat_cols])
    else
      observed <- as.numeric(observed[1, ])
  }
  if (length(observed) != n_stats)
    stop(sprintf("observed has %d elements but emulator expects %d stats.",
                 length(observed), n_stats))

  # --- Estimate residual variance if not provided ---
  if (is.null(sigma2)) {
    if (verbose) cat("PipeMaster:: Estimating residual variance from validation data\n")
    sigma2 <- .estimate.residual.variance(emulator)
  }
  # Floor sigma2 to avoid division by near-zero
  sigma2[sigma2 < 1e-12] <- 1e-12

  # --- Starting points (n_chains x n_params in log-space) ---
  n_chains <- as.integer(n_chains)
  if (is.null(start)) {
    # Use diverse optima from emulator.optimize
    if (verbose) cat("PipeMaster:: Getting starting points via emulator.optimize()\n")
    n_opt_starts <- max(50L, n_chains * 5L)
    opt <- emulator.optimize(emulator, observed,
                             prior.bounds = bounds,
                             n_starts = n_opt_starts,
                             n_steps = 1000L,
                             verbose = verbose)
    # Pick top n_chains optima (diverse starts)
    top_idx <- order(opt$losses)[seq_len(min(n_chains, length(opt$losses)))]
    starts_orig <- opt$optima[top_idx, , drop = FALSE]
    if (nrow(starts_orig) < n_chains) {
      # Replicate with jitter if not enough optima
      extra <- n_chains - nrow(starts_orig)
      base <- starts_orig[rep(1L, extra), , drop = FALSE]
      jitter_sd <- (log(hi) - log(lo)) / 40
      base_log <- log(base) + matrix(rnorm(extra * n_params), nrow = extra) *
                  rep(jitter_sd, each = extra)
      base_log <- pmax(pmin(base_log, rep(log(hi), each = extra)),
                       rep(log(lo), each = extra))
      starts_orig <- rbind(starts_orig, exp(base_log))
    }
    theta_current_log <- log(starts_orig)
  } else if (is.matrix(start)) {
    # User provided matrix (n_chains x n_params)
    if (nrow(start) != n_chains)
      stop(sprintf("start matrix has %d rows but n_chains=%d", nrow(start), n_chains))
    theta_current_log <- log(start)
  } else {
    # User provided single vector — jitter for multiple chains
    start_log <- log(as.numeric(start))
    jitter_sd <- (log(hi) - log(lo)) / 40
    theta_current_log <- matrix(start_log, nrow = n_chains, ncol = n_params,
                                byrow = TRUE)
    if (n_chains > 1L) {
      noise <- matrix(rnorm((n_chains - 1L) * n_params), nrow = n_chains - 1L)
      theta_current_log[-1L, ] <- theta_current_log[-1L, , drop = FALSE] +
                                  t(t(noise) * jitter_sd)
      # Clamp to bounds
      theta_current_log[-1L, ] <- pmax(
        pmin(theta_current_log[-1L, , drop = FALSE],
             rep(log(hi), each = n_chains - 1L)),
        rep(log(lo), each = n_chains - 1L))
    }
  }

  # --- Auto-calibrate proposal_sd ---
  if (is.null(proposal_sd)) {
    proposal_sd <- (log(hi) - log(lo)) / 20
  }
  if (length(proposal_sd) == 1L)
    proposal_sd <- rep(proposal_sd, n_params)

  # --- MCMC setup ---
  n_samples       <- as.integer(n_samples)
  burnin          <- as.integer(burnin)
  thin            <- as.integer(thin)
  samples_per_chain <- ceiling(n_samples / n_chains)
  iter_per_chain  <- burnin + samples_per_chain * thin

  if (verbose) {
    cat("PipeMaster:: emulator.MCMC (Python-accelerated multi-chain MH-MCMC)\n")
    cat(sprintf("PipeMaster:: %d chains x %d samples/chain = %d total, %d burn-in, thin=%d\n",
                n_chains, samples_per_chain, samples_per_chain * n_chains,
                burnin, thin))
    cat(sprintf("PipeMaster:: %d iterations/chain (%d total model calls)\n",
                iter_per_chain, iter_per_chain))
    cat(sprintf("PipeMaster:: %d params, %d stats\n", n_params, n_stats))
    for (j in seq_len(n_params))
      cat(sprintf("  %s: [%.1f, %.1f], proposal_sd=%.4f\n",
                  param_cols[j], lo[j], hi[j], proposal_sd[j]))
    cat(sprintf("PipeMaster:: Residual variance (sigma2) range: [%.2f, %.2f]\n",
                min(sigma2), max(sigma2)))
  }

  # --- Define Python MCMC functions (once per call) ---
  if (!requireNamespace("reticulate", quietly = TRUE))
    stop("emulator.MCMC() requires the 'reticulate' R package.")

  reticulate::py_run_string("
import numpy as np
import time

def _log_posterior_batch(model, theta_log, bounds_lo, bounds_hi,
                         feat_mu, feat_sd, target_mu, target_sd,
                         observed, sigma2):
    '''Vectorized log-posterior for (n_chains, n_params) batch.'''
    # Bounds check (vectorized)
    in_bounds = np.all((theta_log >= bounds_lo) & (theta_log <= bounds_hi), axis=1)
    lp = np.full(theta_log.shape[0], -np.inf)
    if not np.any(in_bounds):
        return lp
    idx = np.where(in_bounds)[0]
    theta_in = theta_log[idx]

    # Z-score to emulator input space
    theta_z = (theta_in - feat_mu) / feat_sd

    # Model predict — direct Python call, no R round-trip
    import tensorflow as tf
    pred_z = model(tf.constant(theta_z, dtype=tf.float32), training=False).numpy()

    # Inverse transform: un-Z-score then inverse signed-log1p
    pred_t = pred_z * target_sd + target_mu
    pred_orig = np.sign(pred_t) * np.expm1(np.abs(pred_t))

    # Gaussian log-likelihood (per-chain)
    resid = pred_orig - observed
    ll = -0.5 * np.sum(resid**2 / sigma2, axis=1)
    lp[idx] = ll
    return lp


def mcmc_batched(model, theta_start_log, bounds_lo, bounds_hi,
                 feat_mu, feat_sd, target_mu, target_sd,
                 observed, sigma2, proposal_sd,
                 n_iter, burnin, thin, n_chains, adaptive, verbose):
    '''Full batched multi-chain MH-MCMC loop in Python.

    Returns dict with chain_samples, chain_lp, accepts, proposal_sd.
    '''
    n_params = theta_start_log.shape[1]
    samples_per_chain = int(np.ceil((n_iter - burnin) / thin))

    # Ensure all inputs are numpy float64
    theta_current = theta_start_log.astype(np.float64).copy()
    bounds_lo = bounds_lo.astype(np.float64).ravel()
    bounds_hi = bounds_hi.astype(np.float64).ravel()
    feat_mu = feat_mu.astype(np.float64).ravel()
    feat_sd = feat_sd.astype(np.float64).ravel()
    target_mu = target_mu.astype(np.float64).ravel()
    target_sd = target_sd.astype(np.float64).ravel()
    observed = observed.astype(np.float64).ravel()
    sigma2 = sigma2.astype(np.float64).ravel()
    proposal_sd = proposal_sd.astype(np.float64).ravel().copy()

    # Storage
    chain_samples = np.empty((samples_per_chain, n_params, n_chains), dtype=np.float64)
    chain_lp = np.empty((samples_per_chain, n_chains), dtype=np.float64)
    accepts = np.zeros(n_chains, dtype=np.int64)
    adapt_accepts = np.zeros(n_chains, dtype=np.int64)

    # Evaluate starting points
    lp_current = _log_posterior_batch(model, theta_current,
                                      bounds_lo, bounds_hi,
                                      feat_mu, feat_sd,
                                      target_mu, target_sd,
                                      observed, sigma2)

    log_interval = max(100, n_iter // 40)
    t0 = time.time()

    if verbose:
        print('PipeMaster:: Burn-in phase', flush=True)

    for i in range(1, n_iter + 1):
        # Propose
        noise = np.random.randn(n_chains, n_params) * proposal_sd
        theta_proposed = theta_current + noise

        # Evaluate
        lp_proposed = _log_posterior_batch(model, theta_proposed,
                                           bounds_lo, bounds_hi,
                                           feat_mu, feat_sd,
                                           target_mu, target_sd,
                                           observed, sigma2)

        # Accept/reject (vectorized)
        log_alpha = lp_proposed - lp_current
        log_u = np.log(np.random.rand(n_chains))
        accept = np.isfinite(log_alpha) & (log_u < log_alpha)

        if np.any(accept):
            theta_current[accept] = theta_proposed[accept]
            lp_current[accept] = lp_proposed[accept]
            accepts[accept] += 1
            adapt_accepts[accept] += 1

        # Adaptive Metropolis during burn-in
        if adaptive and i <= burnin and i % 100 == 0:
            adapt_rate = np.sum(adapt_accepts) / (100.0 * n_chains)
            if adapt_rate > 0.30:
                proposal_sd *= 1.1
            elif adapt_rate < 0.15:
                proposal_sd *= 0.9
            adapt_accepts[:] = 0

        # Transition message
        if verbose and i == burnin + 1:
            print('PipeMaster:: Sampling phase', flush=True)

        # Store (after burn-in, respecting thinning)
        if i > burnin and (i - burnin) % thin == 0:
            idx = (i - burnin) // thin - 1
            if idx < samples_per_chain:
                chain_samples[idx, :, :] = np.exp(theta_current).T
                chain_lp[idx, :] = lp_current

        # Progress
        if verbose and i % log_interval == 0:
            elapsed = time.time() - t0
            rate = i / elapsed if elapsed > 0 else float('inf')
            eta = (n_iter - i) / rate if rate > 0 else 0
            phase = 'burn-in' if i <= burnin else 'sampling'
            mean_acc = 100.0 * np.sum(accepts) / (i * n_chains)
            finite_lp = lp_current[np.isfinite(lp_current)]
            mean_lp = np.mean(finite_lp) if len(finite_lp) > 0 else float('-inf')
            print(f'  [{i}/{n_iter}] {phase} | acc={mean_acc:.1f}% | mean logP={mean_lp:.1f} | {rate:.0f} it/s | ETA {eta:.0f}s',
                  flush=True)

    elapsed = time.time() - t0
    if verbose:
        rate = n_iter / elapsed if elapsed > 0 else float('inf')
        print(f'PipeMaster:: MCMC done in {elapsed:.1f}s ({rate:.0f} it/s, {n_chains} chains x {n_iter} it)',
              flush=True)

    return {
        'chain_samples': chain_samples,
        'chain_lp': chain_lp,
        'accepts': accepts,
        'proposal_sd': proposal_sd
    }
", convert = FALSE)

  # --- Run MCMC in Python ---
  np <- reticulate::import("numpy", convert = FALSE)

  py_result <- reticulate::py$mcmc_batched(
    model       = em_model,
    theta_start_log = np$array(theta_current_log),
    bounds_lo   = np$array(bounds_log_lo),
    bounds_hi   = np$array(bounds_log_hi),
    feat_mu     = np$array(feat_mu),
    feat_sd     = np$array(feat_sd),
    target_mu   = np$array(target_mu),
    target_sd   = np$array(target_sd),
    observed    = np$array(observed),
    sigma2      = np$array(sigma2),
    proposal_sd = np$array(proposal_sd),
    n_iter      = as.integer(iter_per_chain),
    burnin      = burnin,
    thin        = thin,
    n_chains    = n_chains,
    adaptive    = adaptive,
    verbose     = verbose
  )

  # Extract results back to R
  chain_samples <- reticulate::py_to_r(py_result$chain_samples)
  chain_lp      <- reticulate::py_to_r(py_result$chain_lp)
  accepts       <- as.integer(reticulate::py_to_r(py_result$accepts))
  proposal_sd   <- as.numeric(reticulate::py_to_r(py_result$proposal_sd))

  # --- Combine chains ---
  # Interleave: chain 1 sample 1, chain 2 sample 1, ..., chain N sample 1,
  #             chain 1 sample 2, ...
  total_from_chains <- samples_per_chain * n_chains
  all_samples <- matrix(NA_real_, nrow = total_from_chains, ncol = n_params)
  all_lp      <- numeric(total_from_chains)
  for (c in seq_len(n_chains)) {
    rows <- seq(c, total_from_chains, by = n_chains)
    all_samples[rows, ] <- chain_samples[, , c]
    all_lp[rows]        <- chain_lp[, c]
  }
  # Truncate to exactly n_samples
  if (total_from_chains > n_samples) {
    all_samples <- all_samples[seq_len(n_samples), , drop = FALSE]
    all_lp      <- all_lp[seq_len(n_samples)]
  }
  colnames(all_samples) <- param_cols

  # --- Diagnostics ---
  # Per-chain acceptance rates
  chain_acc <- accepts / iter_per_chain
  acceptance_rate <- mean(chain_acc)

  # Posterior mean
  point_estimate <- colMeans(all_samples)
  names(point_estimate) <- param_cols

  # 95% HPD intervals
  ci <- t(apply(all_samples, 2, .hpd.interval, prob = 0.95))
  colnames(ci) <- c("lower", "upper")
  rownames(ci) <- param_cols

  # ESS: compute per-chain then sum (interleaved samples break autocorrelation)
  ess_vals <- rep(0, n_params)
  for (j in seq_len(n_params)) {
    for (cc in seq_len(n_chains))
      ess_vals[j] <- ess_vals[j] + .ess(chain_samples[, j, cc])
  }
  names(ess_vals) <- param_cols

  # Gelman-Rubin R-hat (per-param, across chains)
  rhat_vals <- rep(NA_real_, n_params)
  names(rhat_vals) <- param_cols
  if (n_chains >= 2L) {
    for (j in seq_len(n_params)) {
      chain_list <- lapply(seq_len(n_chains), function(c) chain_samples[, j, c])
      rhat_vals[j] <- .rhat(chain_list)
    }
  }

  if (verbose) {
    cat(sprintf("PipeMaster:: Acceptance rate: %.1f%% (mean across %d chains)\n",
                100 * acceptance_rate, n_chains))
    # Per-chain summary
    cat(sprintf("PipeMaster:: Per-chain acceptance: %s\n",
                paste(sprintf("%.0f%%", 100 * chain_acc), collapse = " ")))
    if (acceptance_rate < 0.15 || acceptance_rate > 0.35)
      cat("PipeMaster:: WARNING: acceptance rate outside [15%, 35%] range\n")

    cat("PipeMaster:: Posterior summary:\n")
    hdr <- sprintf("  %-15s %12s %15s %8s", "Param", "Mean", "95% HPD", "ESS")
    if (n_chains >= 2L) hdr <- paste0(hdr, sprintf(" %6s", "R-hat"))
    cat(hdr, "\n")
    for (j in seq_len(n_params)) {
      line <- sprintf("  %-15s %12.1f [%7.1f, %7.1f] %8.0f",
                      param_cols[j], point_estimate[j],
                      ci[j, 1], ci[j, 2], ess_vals[j])
      if (n_chains >= 2L)
        line <- paste0(line, sprintf(" %6.3f", rhat_vals[j]))
      cat(line, "\n")
    }

    low_ess <- ess_vals[ess_vals < 200]
    if (length(low_ess) > 0L)
      cat(sprintf("PipeMaster:: WARNING: low ESS for: %s\n",
                  paste(names(low_ess), collapse = ", ")))
    if (n_chains >= 2L) {
      bad_rhat <- rhat_vals[!is.na(rhat_vals) & rhat_vals > 1.1]
      if (length(bad_rhat) > 0L)
        cat(sprintf("PipeMaster:: WARNING: R-hat > 1.1 for: %s (chains may not have converged)\n",
                    paste(names(bad_rhat), collapse = ", ")))
    }
  }

  result <- list(
    samples            = all_samples,
    log_posterior      = all_lp,
    acceptance_rate    = acceptance_rate,
    chain_acceptance   = chain_acc,
    sigma2             = sigma2,
    prior_bounds       = bounds,
    param_cols         = param_cols,
    point_estimate     = point_estimate,
    credible_intervals = ci,
    ess                = ess_vals,
    n_chains           = n_chains,
    chain_samples      = chain_samples,
    rhat               = rhat_vals
  )
  class(result) <- "emulator_posterior"
  result
}

# ============================================================================
# S3: summary.emulator_posterior
# ============================================================================

#' @export
summary.emulator_posterior <- function(object, ...) {
  samples    <- object$samples
  param_cols <- object$param_cols
  n_params   <- length(param_cols)

  tab <- data.frame(
    Parameter = param_cols,
    Mean      = colMeans(samples),
    Median    = apply(samples, 2, median),
    SD        = apply(samples, 2, sd),
    HPD_lower = object$credible_intervals[, "lower"],
    HPD_upper = object$credible_intervals[, "upper"],
    ESS       = object$ess,
    stringsAsFactors = FALSE
  )
  if (!is.null(object$rhat) && any(!is.na(object$rhat)))
    tab$Rhat <- object$rhat
  rownames(tab) <- param_cols

  result <- list(
    table            = tab,
    acceptance_rate  = object$acceptance_rate,
    chain_acceptance = object$chain_acceptance,
    n_samples        = nrow(samples),
    n_chains         = object$n_chains
  )
  class(result) <- "summary.emulator_posterior"
  result
}

# ============================================================================
# S3: print.summary.emulator_posterior
# ============================================================================

#' @export
print.summary.emulator_posterior <- function(x, ...) {
  cat("Emulator MCMC Posterior Summary\n")
  cat(sprintf("  Chains: %d, Samples: %d, Acceptance rate: %.1f%%\n\n",
              x$n_chains, x$n_samples, 100 * x$acceptance_rate))

  tab <- x$table
  has_rhat <- "Rhat" %in% colnames(tab)
  hdr <- sprintf("  %-15s %10s %10s %10s %20s %8s",
                 "Parameter", "Mean", "Median", "SD", "95% HPD", "ESS")
  if (has_rhat) hdr <- paste0(hdr, sprintf(" %6s", "R-hat"))
  cat(hdr, "\n")
  cat(paste(rep("-", if (has_rhat) 88 else 80), collapse = ""), "\n")
  for (i in seq_len(nrow(tab))) {
    line <- sprintf("  %-15s %10.1f %10.1f %10.1f [%8.1f, %8.1f] %8.0f",
                    tab$Parameter[i], tab$Mean[i], tab$Median[i], tab$SD[i],
                    tab$HPD_lower[i], tab$HPD_upper[i], tab$ESS[i])
    if (has_rhat)
      line <- paste0(line, sprintf(" %6.3f", tab$Rhat[i]))
    cat(line, "\n")
  }
  cat("\n")
  invisible(x)
}

# ============================================================================
# S3: print.emulator_posterior
# ============================================================================

#' @export
print.emulator_posterior <- function(x, ...) {
  cat(sprintf("Emulator MCMC Posterior (%d chains, %d samples, acceptance=%.1f%%)\n",
              x$n_chains, nrow(x$samples), 100 * x$acceptance_rate))
  cat(sprintf("Parameters: %s\n", paste(x$param_cols, collapse = ", ")))
  est_str <- paste(sprintf("%s=%.1f", x$param_cols, x$point_estimate),
                   collapse = " ")
  cat(sprintf("Posterior mean: %s\n", est_str))
  invisible(x)
}

# ============================================================================
# S3: plot.emulator_posterior
# ============================================================================

#' @export
plot.emulator_posterior <- function(x, true_values = NULL, ...) {
  samples    <- x$samples
  param_cols <- x$param_cols
  n_params   <- length(param_cols)
  bounds     <- x$prior_bounds
  n_chains   <- x$n_chains
  chain_samp <- x$chain_samples   # (samples_per_chain, n_params, n_chains)

  lo <- as.numeric(bounds$prior.1)
  hi <- as.numeric(bounds$prior.2)

  # Determine layout: traces + histograms, optionally + pairwise
  has_pairs <- n_params <= 6 && n_params >= 2
  if (has_pairs) {
    n_pairs <- n_params * (n_params - 1) / 2
    n_pair_rows <- ceiling(n_pairs / n_params)
    total_rows <- 2 + n_pair_rows
  } else {
    total_rows <- 2
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(total_rows, n_params), mar = c(4, 4, 2.5, 1))

  # Chain colors
  if (n_chains <= 8L) {
    chain_cols <- c("steelblue", "firebrick", "darkgreen", "darkorange",
                    "purple", "brown", "deeppink", "darkgoldenrod")[seq_len(n_chains)]
  } else {
    chain_cols <- rainbow(n_chains, s = 0.7, v = 0.8)
  }
  chain_cols_alpha <- adjustcolor(chain_cols, alpha.f = 0.5)

  # Row 1: Trace plots (per-chain overlaid)
  for (j in seq_len(n_params)) {
    spc <- dim(chain_samp)[1]  # samples_per_chain
    # Compute y-range across all chains
    yrange <- range(chain_samp[, j, ], na.rm = TRUE)
    if (!is.null(true_values) && param_cols[j] %in% names(true_values))
      yrange <- range(yrange, true_values[param_cols[j]])
    plot(seq_len(spc), chain_samp[, j, 1], type = "l",
         col = chain_cols_alpha[1],
         xlab = "Iteration", ylab = param_cols[j],
         main = paste("Trace:", param_cols[j]),
         ylim = yrange)
    if (n_chains > 1L)
      for (c in 2:n_chains)
        lines(seq_len(spc), chain_samp[, j, c], col = chain_cols_alpha[c])
    # Posterior mean
    abline(h = x$point_estimate[j], col = "black", lwd = 2)
    # True value
    if (!is.null(true_values) && param_cols[j] %in% names(true_values))
      abline(h = true_values[param_cols[j]], col = "red", lwd = 2, lty = 2)
  }

  # Row 2: Marginal posterior histograms (combined samples)
  for (j in seq_len(n_params)) {
    hist(samples[, j], breaks = 50, freq = FALSE,
         col = "lightblue", border = "white",
         xlab = param_cols[j], ylab = "Density",
         main = paste("Posterior:", param_cols[j]))
    # Prior range shading
    usr <- par("usr")
    rect(lo[j], usr[3], hi[j], usr[4],
         col = adjustcolor("grey", alpha.f = 0.15), border = NA)
    # Re-draw histogram on top
    hist(samples[, j], breaks = 50, freq = FALSE,
         col = adjustcolor("steelblue", alpha.f = 0.5),
         border = "white", add = TRUE)
    # Density curve
    d <- density(samples[, j])
    lines(d, col = "darkblue", lwd = 2)
    # Posterior mean
    abline(v = x$point_estimate[j], col = "steelblue", lwd = 2)
    # True value
    if (!is.null(true_values) && param_cols[j] %in% names(true_values))
      abline(v = true_values[param_cols[j]], col = "red", lwd = 2, lty = 2)
    # HPD interval
    abline(v = x$credible_intervals[j, ], col = "darkblue", lwd = 1, lty = 3)
  }

  # Row 3+: Pairwise scatter plots
  if (has_pairs) {
    pair_idx <- 0
    for (i in 1:(n_params - 1)) {
      for (k in (i + 1):n_params) {
        pair_idx <- pair_idx + 1
        n_plot <- min(nrow(samples), 2000L)
        plot_idx <- if (nrow(samples) > n_plot)
          sample(nrow(samples), n_plot) else seq_len(nrow(samples))
        plot(samples[plot_idx, i], samples[plot_idx, k],
             pch = 16, cex = 0.3,
             col = adjustcolor("steelblue", alpha.f = 0.3),
             xlab = param_cols[i], ylab = param_cols[k],
             main = paste(param_cols[i], "vs", param_cols[k]))
        if (!is.null(true_values)) {
          if (param_cols[i] %in% names(true_values) &&
              param_cols[k] %in% names(true_values))
            points(true_values[param_cols[i]], true_values[param_cols[k]],
                   pch = 4, cex = 2, col = "red", lwd = 2)
        }
      }
    }
    remaining <- n_params * n_pair_rows - pair_idx
    if (remaining > 0)
      for (r in seq_len(remaining)) plot.new()
  }

  invisible(x)
}
