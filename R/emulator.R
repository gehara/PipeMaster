# ============================================================================
# Neural Emulator: Forward model θ→S with gradient-based optimization
#
# train.emulator()    — train a NN that maps parameters → summary statistics
# emulator.optimize() — gradient-descent optimization of θ given observed S
# save/load functions  — persist emulator results to disk
# ============================================================================

#' Train a Neural Emulator (Forward Model θ→S)
#'
#' Trains a neural network that maps demographic parameters to summary statistics
#' (the reverse direction of \code{tune.nn()}, which maps S→θ). Uses the same
#' Hyperband search infrastructure as \code{tune.nn()} with a dense ResNet
#' architecture, but with swapped input/output: parameters are features and
#' summary statistics are targets.
#'
#' The trained emulator can be used with \code{emulator.optimize()} to find
#' parameter values that best reproduce observed summary statistics via
#' gradient-based optimization.
#'
#' @param reftable data.frame — output of \code{sim.sumstat()} containing both
#'   parameter columns and statistic columns.
#' @param param.cols character vector — names of demographic parameter columns.
#'   These become the \emph{input} features of the emulator.
#' @param search.space named list — overrides default HP ranges. NULL uses
#'   emulator-specific defaults (see \code{.emulator.search.space()}).
#' @param exclude.cols character vector — additional column names to exclude from
#'   targets (e.g., nuisance columns). Default NULL.
#' @param val.frac numeric — validation fraction (default 0.1).
#' @param max_epochs integer — maximum epochs for Hyperband (default 500).
#' @param eta numeric — Hyperband halving factor (default 3).
#' @param n_searches integer — number of independent Hyperband searches (default 5).
#' @param top_k integer — number of top models to keep (default 3).
#' @param cores integer — max concurrent search workers (default 5).
#' @param gpus integer — number of GPUs (default 0, CPU-only).
#' @param gpu.threshold integer — max searches per GPU (default 4).
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list of class \code{"emulator_result"} with:
#' \describe{
#'   \item{best_hp}{named list of best hyperparameters}
#'   \item{best_val_loss}{best validation loss achieved}
#'   \item{all_results}{data.frame of all evaluated configs}
#'   \item{best_model}{trained keras model (best config)}
#'   \item{models}{list of top-K keras models, sorted by val_loss}
#'   \item{models_val_loss}{numeric vector of validation losses for each model}
#'   \item{models_metrics}{list of per-model metrics (R², MPE per stat column)}
#'   \item{data}{preprocessed data list (transforms, splits)}
#'   \item{param_cols}{character vector of parameter column names (input)}
#'   \item{stat_cols}{character vector of statistic column names (output)}
#'   \item{type}{\code{"emulator"}}
#' }
#'
#' @export
train.emulator <- function(reftable, param.cols,
                           search.space = NULL,
                           exclude.cols = NULL,
                           val.frac = 0.1,
                           max_epochs = 500, eta = 3,
                           n_searches = 5L, top_k = 3L,
                           cores = 5L, gpus = 0L,
                           gpu.threshold = 4L, greedy = TRUE,
                           seed = 42, verbose = TRUE) {

  # Delegate to tune.nn with type = "emulator"
  result <- tune.nn(reftable, param.cols = param.cols,
                    type = "emulator",
                    search_space = search.space,
                    exclude.cols = exclude.cols,
                    val.frac = val.frac,
                    max_epochs = max_epochs, eta = eta,
                    n_searches = n_searches, top_k = top_k,
                    cores = cores, gpus = gpus,
                    gpu.threshold = gpu.threshold, greedy = greedy,
                    seed = seed, verbose = verbose)

  # Augment with emulator-specific metadata
  result$param_cols <- param.cols
  result$stat_cols  <- result$data$stat_cols
  structure(result, class = "emulator_result")
}

# ============================================================================
# Internal: prepare emulator data (params → stats)
# ============================================================================

.prep.emulator.data <- function(reftable, param.cols, exclude.cols, val.frac) {
  nuisance <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance, exclude.cols))

  params <- as.matrix(reftable[, param.cols, drop = FALSE])
  stats  <- as.matrix(reftable[, stat_cols, drop = FALSE])

  # Remove rows with non-finite values
  bad <- apply(params, 1, function(x) any(!is.finite(x) | x <= 0)) |
         apply(stats, 1, function(x) any(!is.finite(x)))
  if (any(bad)) {
    params <- params[!bad, , drop = FALSE]
    stats  <- stats[!bad, , drop = FALSE]
  }

  # Split train/val
  n <- nrow(params)
  idx <- sample(n)
  n_val <- floor(val.frac * n)
  n_train <- n - n_val
  tr <- idx[1:n_train]; va <- idx[(n_train + 1):n]

  # Features: log(params) → Z-score (demographic params are always positive)
  feat_log <- log(params)
  feat_mu <- colMeans(feat_log[tr, , drop = FALSE])
  feat_sd <- apply(feat_log[tr, , drop = FALSE], 2, sd)
  feat_sd[feat_sd == 0] <- 1
  zscore_feat <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore_feat(feat_log[tr, , drop = FALSE])
  X_val   <- zscore_feat(feat_log[va, , drop = FALSE])

  # Targets: sign(x) * log1p(abs(x)) → Z-score (handles negative stats)
  signed_log1p <- function(x) sign(x) * log1p(abs(x))
  stats_t <- signed_log1p(stats)

  t_mu <- colMeans(stats_t[tr, , drop = FALSE])
  t_sd <- apply(stats_t[tr, , drop = FALSE], 2, sd)
  t_sd[t_sd == 0] <- 1
  zscore_tar <- function(Y) t((t(Y) - t_mu) / t_sd)

  list(
    X_train    = X_train,
    X_val      = X_val,
    Y_train    = zscore_tar(stats_t[tr, , drop = FALSE]),
    Y_val      = zscore_tar(stats_t[va, , drop = FALSE]),
    n_features = ncol(X_train),
    feat_mu    = feat_mu,
    feat_sd    = feat_sd,
    target_mu  = t_mu,
    target_sd  = t_sd,
    stat_cols  = stat_cols
  )
}

# ============================================================================
# Internal: default search space for emulator
# ============================================================================

.emulator.search.space <- function() {
  list(
    units_1       = c(64, 128, 256),
    n_resblocks_1 = 1:3,
    units_2       = c(64, 128, 256),
    n_resblocks_2 = 0:2,
    units_3       = c(32, 64, 128),
    learning_rate = c(1e-4, 1e-2),
    l2_reg        = c(1e-5, 1e-3),
    dropout       = c(0, 0.3),
    use_dropout   = c(TRUE, FALSE),
    batch_size    = c(256, 512, 1024),
    huber_delta   = c(0.5, 1.0, 1.5)
  )
}

# ============================================================================
# Internal: inverse transform for emulator targets
#
# Emulator targets use signed log1p:  sign(x) * log1p(abs(x))
# Inverse: sign(y) * expm1(abs(y))
# ============================================================================

.inv.transform.emulator <- function(Y_z, target_mu, target_sd) {
  # Un-Z-score
  Y_t <- t(t(Y_z) * target_sd + target_mu)
  # Inverse signed log1p
  sign(Y_t) * expm1(abs(Y_t))
}

# ============================================================================
# Gradient-based optimization: find θ that minimizes ||f(θ) - S_obs||²
# ============================================================================

#' Optimize Parameters via Neural Emulator Gradient Descent
#'
#' Given a trained neural emulator (from \code{train.emulator()}) and observed
#' summary statistics, finds the parameter values θ that minimize the squared
#' difference between the emulator's predicted statistics f(θ) and the observed
#' statistics. Uses TensorFlow GradientTape for backpropagation with projected
#' gradient descent (clipping to prior bounds) and multi-start to avoid local
#' optima.
#'
#' @param emulator list — output from \code{train.emulator()}.
#' @param observed numeric vector or 1-row data.frame — observed summary statistics.
#' @param model PipeMaster model object — used to extract prior bounds via
#'   \code{get.prior.table()}. Either \code{model} or \code{prior.bounds}
#'   must be provided.
#' @param prior.bounds data.frame with columns \code{Parameter}, \code{prior.1},
#'   \code{prior.2} — alternative to \code{model} for specifying parameter bounds.
#' @param n_starts integer — number of random starting points (default 100).
#' @param n_steps integer — number of gradient descent steps per start (default 2000).
#' @param learning_rate numeric — gradient descent learning rate (default 0.01).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{optima}{matrix (n_starts x n_params) of optimized parameter values (original scale)}
#'   \item{losses}{numeric vector of final MSE losses per start}
#'   \item{point_estimate}{mean of top-10\% optima (lowest loss)}
#'   \item{prior_bounds}{the bounds used (data.frame)}
#' }
#'
#' @export
emulator.optimize <- function(emulator, observed,
                              model = NULL, prior.bounds = NULL,
                              n_starts = 100L, n_steps = 2000L,
                              learning_rate = 0.01,
                              verbose = TRUE) {

  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("emulator.optimize() requires the 'keras' and 'tensorflow' R packages.")

  tf <- tensorflow::tf

  # --- Extract bounds ---
  if (!is.null(model)) {
    ptab <- get.prior.table(model)
    # Match to emulator param_cols
    idx <- match(emulator$param_cols, ptab$Parameter)
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
    stop("Either 'model' or 'prior.bounds' must be provided")
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

  # Use best model by default
  em_model <- emulator$best_model

  # --- Prepare observed stats ---
  if (is.data.frame(observed)) {
    if (!is.null(stat_cols) && all(stat_cols %in% colnames(observed)))
      observed <- as.numeric(observed[1, stat_cols])
    else
      observed <- as.numeric(observed[1, ])
  }

  # Transform observed to same space as training targets:
  # signed_log1p → Z-score
  obs_t <- sign(observed) * log1p(abs(observed))
  obs_z <- (obs_t - target_mu) / target_sd
  S_obs <- tf$constant(matrix(obs_z, nrow = 1), dtype = "float32")

  # --- Compute bounds in Z-scored log-space ---
  lo <- as.numeric(bounds$prior.1)
  hi <- as.numeric(bounds$prior.2)
  lo_z <- (log(lo) - feat_mu) / feat_sd
  hi_z <- (log(hi) - feat_mu) / feat_sd

  # Ensure lo_z < hi_z
  bounds_lo <- pmin(lo_z, hi_z)
  bounds_hi <- pmax(lo_z, hi_z)
  bounds_lo_tf <- tf$constant(matrix(bounds_lo, nrow = 1), dtype = "float32")
  bounds_hi_tf <- tf$constant(matrix(bounds_hi, nrow = 1), dtype = "float32")

  if (verbose) {
    cat("PipeMaster:: emulator.optimize\n")
    cat(sprintf("PipeMaster:: %d starts x %d steps, lr=%.4f\n",
                n_starts, n_steps, learning_rate))
    cat(sprintf("PipeMaster:: %d params, %d stats\n", n_params, length(stat_cols)))
    for (i in seq_len(n_params))
      cat(sprintf("  %s: [%.1f, %.1f]\n", param_cols[i], lo[i], hi[i]))
  }

  # --- Generate random starting points (uniform in log-space, Z-scored) ---
  set.seed(42L)
  starts_log <- matrix(NA_real_, nrow = n_starts, ncol = n_params)
  for (j in seq_len(n_params))
    starts_log[, j] <- runif(n_starts, log(lo[j]), log(hi[j]))

  # Z-score the starting points
  starts_z <- t((t(starts_log) - feat_mu) / feat_sd)

  # --- Batched gradient descent ---
  # All n_starts run in parallel as a single (n_starts, n_params) tensor
  theta <- tf$Variable(starts_z, dtype = "float32")
  S_obs_batch <- tf$tile(S_obs, c(as.integer(n_starts), 1L))
  lr_tf <- tf$constant(learning_rate, dtype = "float32")

  `%as%` <- reticulate::`%as%`

  if (verbose) cat("PipeMaster:: Optimizing")

  for (step in seq_len(n_steps)) {
    with(tf$GradientTape() %as% tape, {
      pred <- em_model(theta, training = FALSE)           # (n_starts, n_stats)
      per_start_loss <- tf$reduce_mean((pred - S_obs_batch)^2, axis = 1L)
      total_loss <- tf$reduce_sum(per_start_loss)
    })
    grad <- tape$gradient(total_loss, theta)
    theta$assign_sub(lr_tf * grad)
    # Project to bounds (broadcasting: bounds are (1, n_params))
    theta$assign(tf$clip_by_value(theta, bounds_lo_tf, bounds_hi_tf))

    if (verbose && step %% 500L == 0L)
      cat(sprintf(" [step %d, mean_loss=%.4f]", step,
                  as.numeric(tf$reduce_mean(per_start_loss)$numpy())))
  }
  if (verbose) cat(" done\n")

  # --- Extract results ---
  optima_z <- as.matrix(theta$numpy())
  losses   <- as.numeric(per_start_loss$numpy())

  # Inverse transform: un-Z-score → exp → original scale
  optima_log <- t(t(optima_z) * feat_sd + feat_mu)
  optima     <- exp(optima_log)
  colnames(optima) <- param_cols

  # --- Point estimate: mean of top-10% ---
  n_top <- max(1L, floor(n_starts * 0.1))
  top_idx <- order(losses)[1:n_top]
  point_estimate <- colMeans(optima[top_idx, , drop = FALSE])
  names(point_estimate) <- param_cols

  if (verbose) {
    cat(sprintf("PipeMaster:: Top-%d mean loss: %.6f\n", n_top,
                mean(losses[top_idx])))
    est_str <- paste(sprintf("%s=%.1f", param_cols, point_estimate), collapse = " ")
    cat(sprintf("PipeMaster:: Point estimate: %s\n", est_str))
  }

  list(
    optima         = optima,
    losses         = losses,
    point_estimate = point_estimate,
    prior_bounds   = bounds
  )
}

# ============================================================================
# ABC Rejection Posterior via Neural Emulator
# ============================================================================

#' ABC Rejection Posterior via Neural Emulator
#'
#' Performs ABC (Approximate Bayesian Computation) rejection sampling using a
#' trained neural emulator as a fast simulator surrogate. Samples parameter
#' vectors from the prior, predicts summary statistics through the emulator
#' ensemble, and retains the closest \code{tol} fraction as approximate
#' posterior samples.
#'
#' Returns an \code{nn.posterior} object so that existing \code{summary()},
#' \code{density()}, and \code{plot()} methods work automatically.
#'
#' @param emulator list — output from \code{train.emulator()}.
#' @param observed numeric vector or 1-row data.frame — observed summary statistics.
#' @param model PipeMaster model object — used to extract prior bounds via
#'   \code{get.prior.table()}. Either \code{model} or \code{prior.bounds}
#'   must be provided.
#' @param prior.bounds data.frame with columns \code{Parameter}, \code{prior.1},
#'   \code{prior.2} — alternative to \code{model} for specifying parameter bounds.
#' @param reftable optional data.frame — if provided, extracts prior samples for
#'   the \code{nn.posterior} object (used by \code{plot()} with \code{show_prior}).
#' @param n_samples integer — number of parameter vectors to draw from the prior
#'   (default 1e6).
#' @param tol numeric — fraction of closest samples to retain as posterior
#'   (default 0.01, i.e. top 1\%).
#' @param distance character — distance metric: \code{"euclidean"} (default) or
#'   \code{"mse"}.
#' @param verbose logical — print progress (default TRUE).
#'
#' @return An object of class \code{"nn.posterior"} with:
#' \describe{
#'   \item{point_estimate}{weighted mean of accepted samples (1/distance weights)}
#'   \item{abc_rejection}{matrix (n_accepted x n_params) of posterior samples}
#'   \item{prior}{prior samples from reftable (if provided)}
#'   \item{param_names}{character vector of parameter names}
#' }
#'
#' @export
emulator.predict <- function(emulator, observed,
                             model = NULL, prior.bounds = NULL,
                             reftable = NULL,
                             n_samples = 1e6, tol = 0.01,
                             distance = c("euclidean", "mse"),
                             verbose = TRUE) {

  if (!requireNamespace("keras", quietly = TRUE))
    stop("emulator.predict() requires the 'keras' R package.")

  distance <- match.arg(distance)
  n_samples <- as.integer(n_samples)
  n_accept  <- max(1L, floor(n_samples * tol))

  # --- Extract bounds (same logic as emulator.optimize) ---
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

  lo <- as.numeric(bounds$prior.1)
  hi <- as.numeric(bounds$prior.2)

  if (verbose) {
    cat("PipeMaster:: emulator.predict (ABC rejection)\n")
    cat(sprintf("PipeMaster:: %s samples, tol=%.4f → keep %s, distance=%s\n",
                format(n_samples, big.mark = ","), tol,
                format(n_accept, big.mark = ","), distance))
    cat(sprintf("PipeMaster:: %d params, %d stats\n", n_params, n_stats))
    for (i in seq_len(n_params))
      cat(sprintf("  %s: [%.1f, %.1f]\n", param_cols[i], lo[i], hi[i]))
  }

  # --- Sample θ from prior (uniform in original scale) ---
  theta_orig <- matrix(NA_real_, nrow = n_samples, ncol = n_params)
  for (j in seq_len(n_params))
    theta_orig[, j] <- runif(n_samples, lo[j], hi[j])
  colnames(theta_orig) <- param_cols

  # --- Transform to emulator input space: log → Z-score ---
  theta_z <- t((t(log(theta_orig)) - feat_mu) / feat_sd)

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

  # --- Ensemble prediction in batches ---
  models    <- emulator$models
  val_losses <- emulator$models_val_loss
  if (is.null(models) || length(models) == 0L) {
    models     <- list(emulator$best_model)
    val_losses <- emulator$best_val_loss
  }
  # Inverse val_loss weights (normalized)
  weights <- 1 / val_losses
  weights <- weights / sum(weights)
  n_models <- length(models)

  batch_size <- 50000L
  n_batches  <- ceiling(n_samples / batch_size)

  # Allocate prediction matrix (Z-scored target space)
  pred_z <- matrix(0, nrow = n_samples, ncol = n_stats)

  if (verbose) cat("PipeMaster:: Predicting")
  for (b in seq_len(n_batches)) {
    i_start <- (b - 1L) * batch_size + 1L
    i_end   <- min(b * batch_size, n_samples)
    X_batch <- theta_z[i_start:i_end, , drop = FALSE]

    # Weighted ensemble average
    batch_pred <- matrix(0, nrow = i_end - i_start + 1L, ncol = n_stats)
    for (m in seq_len(n_models)) {
      p <- as.matrix(predict(models[[m]], X_batch, verbose = 0L))
      batch_pred <- batch_pred + weights[m] * p
    }
    pred_z[i_start:i_end, ] <- batch_pred

    if (verbose && (b %% 5L == 0L || b == n_batches))
      cat(sprintf(" [%d/%d]", b, n_batches))
  }
  if (verbose) cat(" done\n")

  # --- Inverse transform predictions to original scale ---
  pred_orig <- .inv.transform.emulator(pred_z, target_mu, target_sd)

  # --- Compute distances (in original scale, MAD-normalized) ---
  # Normalize each stat by MAD to give equal weight
  stat_mad <- apply(pred_orig, 2, function(x) {
    m <- median(x, na.rm = TRUE)
    mad_val <- median(abs(x - m), na.rm = TRUE)
    if (mad_val == 0) sd(x, na.rm = TRUE) else mad_val
  })
  stat_mad[stat_mad == 0 | !is.finite(stat_mad)] <- 1

  # Normalized residuals
  resid <- t((t(pred_orig) - observed) / stat_mad)

  if (distance == "euclidean") {
    dists <- sqrt(rowSums(resid^2))
  } else {
    dists <- rowMeans(resid^2)
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
  names(point_est) <- param_cols

  if (verbose) {
    est_str <- paste(sprintf("%s=%.1f", param_cols, point_est), collapse = " ")
    cat(sprintf("PipeMaster:: Point estimate: %s\n", est_str))
  }

  # --- Prior samples from reftable ---
  prior_samples <- NULL
  if (!is.null(reftable)) {
    nuisance <- c("mean.rate", "sd.rate")
    pcols <- intersect(param_cols, colnames(reftable))
    pcols <- setdiff(pcols, nuisance)
    if (length(pcols) > 0) {
      prior_samples <- as.matrix(reftable[, pcols, drop = FALSE])
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
    prior          = prior_samples,
    param_names    = param_cols
  )
  class(result) <- "nn.posterior"
  result
}

# ============================================================================
# Save / Load emulator results
# ============================================================================

#' Save Emulator Result to Disk
#'
#' Saves a \code{train.emulator()} result to a directory. Keras models are saved
#' in TensorFlow SavedModel format; metadata is saved as RDS.
#'
#' @param emulator list — output from \code{train.emulator()}.
#' @param path character — directory path to save to.
#'
#' @export
save.emulator.result <- function(emulator, path) {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("save.emulator.result() requires the 'keras' package.")

  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  # Save best_model
  model_dir <- file.path(path, "best_model")
  keras::save_model_tf(emulator$best_model, model_dir)

  # Save all top-K models
  models <- emulator$models
  if (!is.null(models) && length(models) > 0L) {
    models_dir <- file.path(path, "models")
    dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_along(models))
      keras::save_model_tf(models[[i]],
                           file.path(models_dir, sprintf("model_%03d", i)))
  }

  # Save metadata as RDS (strip keras objects)
  result_no_model <- emulator
  result_no_model$best_model <- NULL
  result_no_model$models <- NULL
  class(result_no_model) <- "list"
  saveRDS(result_no_model, file.path(path, "emulator_result.rds"))

  cat(sprintf("PipeMaster:: Saved emulator result to %s (%d models)\n",
              path, length(models)))
}

#' Load Emulator Result from Disk
#'
#' Loads a \code{train.emulator()} result previously saved with
#' \code{save.emulator.result()}.
#'
#' @param path character — directory path where files were saved.
#'
#' @return A list identical to the output of \code{train.emulator()}.
#'
#' @export
load.emulator.result <- function(path) {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("load.emulator.result() requires the 'keras' package.")

  rds_file  <- file.path(path, "emulator_result.rds")
  model_dir <- file.path(path, "best_model")

  if (!file.exists(rds_file))
    stop("emulator_result.rds not found in: ", path)
  if (!dir.exists(model_dir))
    stop("best_model/ directory not found in: ", path)

  result <- readRDS(rds_file)
  result$best_model <- keras::load_model_tf(model_dir)

  # Load top-K models
  models_dir <- file.path(path, "models")
  if (dir.exists(models_dir)) {
    model_subdirs <- sort(list.dirs(models_dir, recursive = FALSE,
                                    full.names = TRUE))
    if (length(model_subdirs) > 0L) {
      result$models <- lapply(model_subdirs, function(d) {
        tryCatch(keras::load_model_tf(d), error = function(e) {
          warning("Could not load model from ", d, ": ",
                  conditionMessage(e), call. = FALSE)
          NULL
        })
      })
      keep <- !vapply(result$models, is.null, logical(1))
      result$models <- result$models[keep]
      if (!is.null(result$models_val_loss))
        result$models_val_loss <- result$models_val_loss[keep]
      cat(sprintf("PipeMaster:: Loaded emulator result from %s (%d models)\n",
                  path, length(result$models)))
    } else {
      result$models <- list(result$best_model)
      cat(sprintf("PipeMaster:: Loaded emulator result from %s\n", path))
    }
  } else {
    result$models <- list(result$best_model)
    if (is.null(result$models_val_loss))
      result$models_val_loss <- result$best_val_loss
    cat(sprintf("PipeMaster:: Loaded emulator result from %s\n", path))
  }

  structure(result, class = "emulator_result")
}
