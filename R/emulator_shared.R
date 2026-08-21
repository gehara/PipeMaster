# ============================================================================
# Emulator data preparation and hyperparameter search space
#
# WHY THIS FILE EXISTS (2026-08-21)
#   These two functions used to live in R/emulator.R, which is gitignored
#   (.gitignore:60) and therefore absent from the GitHub repository. But
#   R/tune.nn.R -- which IS tracked -- calls them:
#     .emulator.search.space()  at R/tune.nn.R:139 and :1800
#     .prep.emulator.data()     at R/tune.nn.R:197
#   so every install from GitHub produced a package whose tune.nn() failed on
#   first use with "could not find function '.emulator.search.space'".
#   tune.nn.classify() was broken by the same route: torch_classifier.R:515
#   calls .default.search.space(), whose list literal at tune.nn.R:139
#   eagerly evaluates .emulator.search.space().
#
#   It went unnoticed because lagarto and segovia both install from a working
#   directory where R/emulator.R is present via Dropbox, so R CMD INSTALL picks
#   it up. Only a git-based install (tarrega) exposes it.
#
#   Both functions are pure base R with no keras/reticulate dependency, so they
#   are moved here rather than un-ignoring R/emulator.R -- that keeps the keras
#   removal in 50d83f4 coherent.
# ============================================================================

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
