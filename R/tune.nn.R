#' Hyperband Tuner for Neural Network Architectures
#'
#' Pure R implementation of the Hyperband algorithm for tuning keras neural
#' networks used in demographic parameter estimation.
#' Supports ResNet (summary statistics), 1D CNN (single-pop SFS), and
#' 2D CNN (joint SFS) architectures.
#'
#' @param reftable data.frame — output of \code{sim.sumstat()} or \code{sim.sfs()}
#'   containing both parameter columns and statistic columns.
#' @param param.cols character vector — names of parameter columns (targets to predict).
#' @param type character — architecture type: \code{"sumstat"} (ResNet),
#'   \code{"sfs1d"} (1D CNN), or \code{"sfs2d"} (2D CNN).
#' @param sfs.dims integer vector — for \code{type = "sfs2d"} only: \code{c(dim1, dim2)}
#'   of the joint SFS matrix.
#' @param max_epochs integer — maximum epochs for Hyperband (default 500).
#' @param eta numeric — Hyperband halving factor (default 3).
#' @param hyperband_iterations integer — number of full Hyperband sweeps (default 1).
#' @param search_space named list — overrides default HP ranges. NULL uses
#'   architecture-specific defaults.
#' @param val.frac numeric — validation fraction (default 0.1).
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{best_hp}{named list of best hyperparameters}
#'   \item{best_val_loss}{best validation loss achieved}
#'   \item{all_results}{data.frame of all evaluated configs (hp_string, val_loss, bracket, round)}
#'   \item{best_model}{trained keras model (best config, retrained to max_epochs)}
#' }
#'
#' @export
tune.nn <- function(reftable, param.cols,
                    type = c("sumstat", "sfs1d", "sfs2d"),
                    sfs.dims = NULL,
                    max_epochs = 500, eta = 3, hyperband_iterations = 1,
                    search_space = NULL,
                    val.frac = 0.1, seed = 42, verbose = TRUE) {

  # --- Dependency check ---
  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("tune.nn() requires the 'keras' and 'tensorflow' R packages.\n",
         "Install with: install.packages(c('keras', 'tensorflow'))\n",
         "Then run: keras::install_keras()")

  type <- match.arg(type)

  if (type == "sfs2d" && (is.null(sfs.dims) || length(sfs.dims) != 2))
    stop("sfs.dims must be c(dim1, dim2) for type='sfs2d'")

  # --- Search space ---
  ss <- if (is.null(search_space)) .default.search.space(type) else search_space

  # --- Prepare data ---
  if (verbose) cat(sprintf("PipeMaster:: tune.nn \u2014 Hyperband (%s)\n",
                           switch(type, sumstat = "ResNet", sfs1d = "1D CNN", sfs2d = "2D CNN")))

  data <- .prep.data(reftable, param.cols, type, sfs.dims, val.frac, seed)

  n_feat <- if (type == "sumstat") ncol(data$X_train) else data$n_features
  n_targ <- ncol(data$Y_train)

  if (verbose) cat(sprintf("PipeMaster:: %d features, %d targets | %d train, %d val\n",
                           n_feat, n_targ, nrow(data$X_train), nrow(data$X_val)))

  # --- Run Hyperband ---
  hb <- .hyperband(
    search_space = ss,
    data         = data,
    type         = type,
    sfs.dims     = sfs.dims,
    max_epochs   = max_epochs,
    eta          = eta,
    n_iter       = hyperband_iterations,
    seed         = seed,
    verbose      = verbose
  )

  # --- Retrain best config (warm-start from saved weights) ---
  if (verbose)
    cat(sprintf("\nPipeMaster:: Best config: %s\n", .hp.to.string(hb$best_hp, type)))

  tensorflow::tf$random$set_seed(as.integer(seed))
  best_model <- .build.nn(hb$best_hp, data, type, sfs.dims)

  # Load weights from best Hyperband model
  weights_loaded <- tryCatch({
    keras::load_model_weights_tf(best_model, hb$best_weights_path)
    TRUE
  }, error = function(e) {
    if (verbose) cat("PipeMaster:: [warn] Could not load saved weights, training from scratch\n")
    FALSE
  })

  # Continue training if there are remaining epochs
  start_epoch <- if (weights_loaded) as.integer(hb$best_epochs) else 0L
  final_val_loss <- hb$best_val_loss

  if (start_epoch < as.integer(max_epochs)) {
    if (verbose)
      cat(sprintf("PipeMaster:: Retraining best for %d epochs (warm-start from epoch %d)...\n",
                  max_epochs, start_epoch))

    bs <- as.integer(hb$best_hp$batch_size)
    retrain_history <- best_model |> keras::fit(
      x = data$X_train, y = data$Y_train,
      validation_data = list(data$X_val, data$Y_val),
      epochs        = as.integer(max_epochs),
      initial_epoch = as.integer(start_epoch),
      batch_size    = bs,
      callbacks  = list(
        keras::callback_early_stopping(monitor = "val_loss", patience = 30L,
                                       restore_best_weights = TRUE),
        keras::callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,
                                             factor = 0.5, min_lr = 1e-6, verbose = 0L)
      ),
      verbose = 0L
    )

    retrain_vl <- retrain_history$metrics$val_loss
    if (is.null(retrain_vl)) retrain_vl <- retrain_history$history$val_loss
    retrain_vl <- unlist(retrain_vl)
    if (length(retrain_vl) > 0 && any(is.finite(retrain_vl)))
      final_val_loss <- min(final_val_loss, min(retrain_vl[is.finite(retrain_vl)]))
  } else {
    if (verbose) cat("PipeMaster:: Best model already trained to max_epochs, skipping retrain\n")
  }

  # Clean up temp weights
  unlink(hb$best_weights_path, recursive = TRUE)
  if (verbose) cat(sprintf("PipeMaster:: Final val_loss: %.6f\n", final_val_loss))

  list(
    best_hp       = hb$best_hp,
    best_val_loss = final_val_loss,
    all_results   = hb$all_results,
    best_model    = best_model,
    data          = data,
    type          = type,
    sfs.dims      = sfs.dims
  )
}

# ============================================================================
# Internal: default search spaces
# ============================================================================

.default.search.space <- function(type) {
  switch(type,
    sumstat = list(
      units_1       = c(128, 256, 512),
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
    ),
    sfs1d = list(
      n_blocks      = 2:5,
      base_filters  = c(16, 32, 64),
      kernel_start  = c(3, 5, 7, 9),
      use_residual  = c(TRUE, FALSE),
      n_dense       = 1:2,
      dense_units   = c(32, 64, 128, 256),
      dropout       = c(0, 0.4),
      l2_reg        = c(1e-5, 1e-3),
      learning_rate = c(1e-4, 1e-2),
      batch_size    = c(64, 128, 256, 512),
      loss          = c("huber", "mse")
    ),
    sfs2d = list(
      n_blocks      = 2:5,
      base_filters  = c(16, 32, 64),
      kernel_start  = c(3, 5, 7),
      use_residual  = c(TRUE, FALSE),
      n_dense       = 1:2,
      dense_units   = c(64, 128, 256),
      dropout       = c(0, 0.4),
      l2_reg        = c(1e-5, 1e-3),
      learning_rate = c(1e-4, 1e-2),
      batch_size    = c(64, 128, 256),
      loss          = c("huber", "mse")
    )
  )
}

# ============================================================================
# Internal: sample one HP configuration
# ============================================================================

.sample.config <- function(search_space) {
  hp <- list()
  log_uniform_params <- c("learning_rate", "l2_reg")

  for (nm in names(search_space)) {
    vals <- search_space[[nm]]

    if (nm %in% log_uniform_params && length(vals) == 2 && is.numeric(vals)) {
      # Log-uniform: 2-element numeric range
      hp[[nm]] <- 10^runif(1, log10(vals[1]), log10(vals[2]))
    } else if (nm == "dropout" && length(vals) == 2 && is.numeric(vals)) {
      # Uniform: 2-element numeric range
      hp[[nm]] <- runif(1, vals[1], vals[2])
    } else {
      # Discrete: sample from values
      hp[[nm]] <- sample(vals, 1)
    }
  }
  hp
}

# ============================================================================
# Internal: HP config to string
# ============================================================================

.hp.to.string <- function(hp, type) {
  if (type == "sumstat") {
    sprintf("u1=%d rb1=%d u2=%d rb2=%d u3=%d lr=%.1e l2=%.1e do=%.2f bs=%d hd=%.1f",
            hp$units_1, hp$n_resblocks_1, hp$units_2, hp$n_resblocks_2,
            hp$units_3, hp$learning_rate, hp$l2_reg, hp$dropout,
            hp$batch_size, hp$huber_delta)
  } else {
    sprintf("blk=%d filt=%d ks=%d res=%s dns=%dx%d do=%.2f l2=%.1e lr=%.1e bs=%d loss=%s",
            hp$n_blocks, hp$base_filters, hp$kernel_start,
            if (isTRUE(hp$use_residual)) "Y" else "N",
            hp$n_dense, hp$dense_units, hp$dropout, hp$l2_reg,
            hp$learning_rate, hp$batch_size, hp$loss)
  }
}

# ============================================================================
# Internal: build and compile a keras model from an HP config
# ============================================================================

.build.nn <- function(hp, data, type, sfs.dims) {
  switch(type,
    sumstat = .build.resnet(hp, data),
    sfs1d   = .build.cnn1d(hp, data),
    sfs2d   = .build.cnn2d(hp, data, sfs.dims)
  )
}

# --- ResNet for summary statistics ---
.build.resnet <- function(hp, data) {
  n_features <- ncol(data$X_train)
  n_targets  <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  # Residual block helper
  res_block <- function(x, units) {
    skip <- x
    x <- x |>
      keras::layer_dense(units = as.integer(units), activation = "relu",
                         kernel_regularizer = l2) |>
      keras::layer_batch_normalization() |>
      keras::layer_dense(units = as.integer(units), activation = "linear",
                         kernel_regularizer = l2) |>
      keras::layer_batch_normalization()
    x <- keras::layer_add(list(x, skip))
    x <- keras::layer_activation(x, activation = "relu")
    x
  }

  inp <- keras::layer_input(shape = n_features)

  # First dense + residual group
  x <- inp |>
    keras::layer_dense(units = as.integer(hp$units_1), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  for (i in seq_len(hp$n_resblocks_1))
    x <- res_block(x, hp$units_1)

  # Middle transition + optional dropout
  x <- x |>
    keras::layer_dense(units = as.integer(hp$units_2), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  # Second residual group
  if (hp$n_resblocks_2 > 0)
    for (i in seq_len(hp$n_resblocks_2))
      x <- res_block(x, hp$units_2)

  # Final dense + optional dropout
  x <- x |>
    keras::layer_dense(units = as.integer(hp$units_3), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  out <- x |> keras::layer_dense(units = n_targets, activation = "linear")

  model <- keras::keras_model(inp, out)
  model |> keras::compile(
    loss      = keras::loss_huber(delta = hp$huber_delta),
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate)
  )
  model
}

# --- 1D CNN for single-population SFS ---
.build.cnn1d <- function(hp, data) {
  n_bins    <- data$n_bins
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)
  input <- keras::layer_input(shape = c(n_bins, 1L))
  x <- input

  for (b in seq_len(hp$n_blocks)) {
    filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))
    ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))

    x <- x |>
      keras::layer_conv_1d(filters = filters, kernel_size = ks,
                           padding = "same", kernel_regularizer = l2) |>
      keras::layer_batch_normalization() |>
      keras::layer_activation("relu")

    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {
      skip <- x
      x <- x |>
        keras::layer_conv_1d(filters = filters, kernel_size = ks,
                             padding = "same", kernel_regularizer = l2) |>
        keras::layer_batch_normalization() |>
        keras::layer_activation("relu") |>
        keras::layer_conv_1d(filters = filters, kernel_size = ks,
                             padding = "same", kernel_regularizer = l2) |>
        keras::layer_batch_normalization()
      x <- keras::layer_add(list(x, skip)) |> keras::layer_activation("relu")
    }
  }

  x <- x |> keras::layer_global_average_pooling_1d()

  for (d in seq_len(hp$n_dense)) {
    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
    x <- x |>
      keras::layer_dense(units = units, activation = "relu",
                         kernel_regularizer = l2) |>
      keras::layer_batch_normalization() |>
      keras::layer_dropout(rate = hp$dropout)
  }

  output <- x |> keras::layer_dense(units = n_targets, activation = "linear")
  model <- keras::keras_model(input, output)

  loss_fn <- if (identical(hp$loss, "huber")) keras::loss_huber(delta = 1.0) else "mse"
  model |> keras::compile(
    loss      = loss_fn,
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate),
    metrics   = list("mae")
  )
  model
}

# --- 2D CNN for joint SFS ---
.build.cnn2d <- function(hp, data, sfs.dims) {
  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)
  input <- keras::layer_input(shape = c(dim1, dim2, 1L))
  x <- input

  for (b in seq_len(hp$n_blocks)) {
    filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))
    ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))

    x <- x |>
      keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                           padding = "same", kernel_regularizer = l2) |>
      keras::layer_batch_normalization() |>
      keras::layer_activation("relu")

    # Max pooling on first and middle blocks
    if (b == 1 || b == (hp$n_blocks %/% 2 + 1))
      x <- keras::layer_max_pooling_2d(x, pool_size = c(2L, 2L))

    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {
      skip <- x
      x <- x |>
        keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                             padding = "same", kernel_regularizer = l2) |>
        keras::layer_batch_normalization() |>
        keras::layer_activation("relu") |>
        keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                             padding = "same", kernel_regularizer = l2) |>
        keras::layer_batch_normalization()
      x <- keras::layer_add(list(x, skip)) |> keras::layer_activation("relu")
    }
  }

  x <- x |> keras::layer_global_average_pooling_2d()

  for (d in seq_len(hp$n_dense)) {
    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
    x <- x |>
      keras::layer_dense(units = units, activation = "relu",
                         kernel_regularizer = l2) |>
      keras::layer_batch_normalization() |>
      keras::layer_dropout(rate = hp$dropout)
  }

  output <- x |> keras::layer_dense(units = n_targets, activation = "linear")
  model <- keras::keras_model(input, output)

  loss_fn <- if (identical(hp$loss, "huber")) keras::loss_huber(delta = 1.0) else "mse"
  model |> keras::compile(
    loss      = loss_fn,
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate),
    metrics   = list("mae")
  )
  model
}

# ============================================================================
# Internal: data preparation
# ============================================================================

.prep.data <- function(reftable, param.cols, type, sfs.dims, val.frac, seed) {
  set.seed(seed)

  if (type == "sumstat") {
    .prep.sumstat(reftable, param.cols, val.frac)
  } else if (type == "sfs1d") {
    .prep.sfs1d(reftable, param.cols, val.frac)
  } else {
    .prep.sfs2d(reftable, param.cols, sfs.dims, val.frac)
  }
}

# --- sumstat: augment features with log1p(abs(x)), Z-score, log-targets ---
.prep.sumstat <- function(reftable, param.cols, val.frac) {
  nuisance <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance))
  target_cols <- setdiff(param.cols, nuisance)

  targets <- as.matrix(reftable[, target_cols])
  features_raw <- as.matrix(reftable[, stat_cols])

  # Feature augmentation
  features <- cbind(features_raw, log1p(abs(features_raw)))

  # Remove bad rows
  bad <- apply(features, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { features <- features[!bad, ]; targets <- targets[!bad, ] }

  # Split
  n <- nrow(features)
  idx <- sample(n)
  n_val <- floor(val.frac * n)
  n_train <- n - n_val

  tr <- idx[1:n_train]; va <- idx[(n_train + 1):n]

  # Z-score features
  feat_mu <- colMeans(features[tr, ])
  feat_sd <- apply(features[tr, ], 2, sd); feat_sd[feat_sd == 0] <- 1
  zscore_feat <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore_feat(features[tr, ])
  X_val   <- zscore_feat(features[va, ])

  # Log-transform + Z-score targets
  Y_tr_log <- log(targets[tr, ]); Y_va_log <- log(targets[va, ])
  t_mu <- colMeans(Y_tr_log); t_sd <- apply(Y_tr_log, 2, sd); t_sd[t_sd == 0] <- 1
  zscore_tar <- function(Y) t((t(Y) - t_mu) / t_sd)

  list(
    X_train = X_train, X_val = X_val,
    Y_train = zscore_tar(Y_tr_log), Y_val = zscore_tar(Y_va_log),
    n_features = ncol(X_train),
    feat_mu = feat_mu, feat_sd = feat_sd,
    target_mu = t_mu, target_sd = t_sd,
    stat_cols = stat_cols
  )
}

# --- sfs1d: proportions + Z-score, reshape to (n, bins, 1) ---
.prep.sfs1d <- function(reftable, param.cols, val.frac) {
  sfs_cols <- grep("^sfs_", colnames(reftable), value = TRUE)
  n_bins <- length(sfs_cols)

  targets <- as.matrix(reftable[, param.cols])
  sfs_raw <- as.matrix(reftable[, sfs_cols])

  # Log1p transform
  sfs_log <- log1p(sfs_raw)

  bad <- apply(sfs_log, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { sfs_log <- sfs_log[!bad, ]; targets <- targets[!bad, ] }

  n <- nrow(sfs_log)
  idx <- sample(n)
  n_val <- floor(val.frac * n)
  n_train <- n - n_val
  tr <- idx[1:n_train]; va <- idx[(n_train + 1):n]

  feat_mu <- colMeans(sfs_log[tr, ])
  feat_sd <- apply(sfs_log[tr, ], 2, sd); feat_sd[feat_sd == 0] <- 1
  zscore <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore(sfs_log[tr, ]); X_val <- zscore(sfs_log[va, ])

  # Reshape to (n, bins, 1)
  dim(X_train) <- c(nrow(X_train), n_bins, 1L)
  dim(X_val)   <- c(nrow(X_val), n_bins, 1L)

  # Log + Z-score targets
  Y_tr_log <- log(targets[tr, ]); Y_va_log <- log(targets[va, ])
  t_mu <- colMeans(Y_tr_log); t_sd <- apply(Y_tr_log, 2, sd); t_sd[t_sd == 0] <- 1
  zscore_tar <- function(Y) t((t(Y) - t_mu) / t_sd)

  list(
    X_train = X_train, X_val = X_val,
    Y_train = zscore_tar(Y_tr_log), Y_val = zscore_tar(Y_va_log),
    n_bins = n_bins, n_features = n_bins,
    feat_mu = feat_mu, feat_sd = feat_sd,
    target_mu = t_mu, target_sd = t_sd,
    stat_cols = sfs_cols
  )
}

# --- sfs2d: log1p + Z-score, reshape to (n, dim1, dim2, 1) ---
.prep.sfs2d <- function(reftable, param.cols, sfs.dims, val.frac) {
  sfs_cols <- grep("^sfs_", colnames(reftable), value = TRUE)
  n_sfs <- length(sfs_cols)
  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]
  if (n_sfs != dim1 * dim2)
    stop(sprintf("Number of SFS columns (%d) doesn't match sfs.dims %dx%d = %d",
                 n_sfs, dim1, dim2, dim1 * dim2))

  targets <- as.matrix(reftable[, param.cols])
  sfs_raw <- as.matrix(reftable[, sfs_cols])

  sfs_log <- log1p(sfs_raw)
  bad <- apply(sfs_log, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { sfs_log <- sfs_log[!bad, ]; targets <- targets[!bad, ] }

  n <- nrow(sfs_log)
  idx <- sample(n)
  n_val <- floor(val.frac * n)
  n_train <- n - n_val
  tr <- idx[1:n_train]; va <- idx[(n_train + 1):n]

  feat_mu <- colMeans(sfs_log[tr, ])
  feat_sd <- apply(sfs_log[tr, ], 2, sd); feat_sd[feat_sd == 0] <- 1
  zscore <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore(sfs_log[tr, ]); X_val <- zscore(sfs_log[va, ])

  # Reshape to (n, dim1, dim2, 1)
  X_train <- array(X_train, dim = c(nrow(X_train), dim1, dim2, 1L))
  X_val   <- array(X_val,   dim = c(nrow(X_val),   dim1, dim2, 1L))

  # Log + Z-score targets
  Y_tr_log <- log(targets[tr, ]); Y_va_log <- log(targets[va, ])
  t_mu <- colMeans(Y_tr_log); t_sd <- apply(Y_tr_log, 2, sd); t_sd[t_sd == 0] <- 1
  zscore_tar <- function(Y) t((t(Y) - t_mu) / t_sd)

  list(
    X_train = X_train, X_val = X_val,
    Y_train = zscore_tar(Y_tr_log), Y_val = zscore_tar(Y_va_log),
    n_bins = n_sfs, n_features = n_sfs,
    feat_mu = feat_mu, feat_sd = feat_sd,
    target_mu = t_mu, target_sd = t_sd,
    stat_cols = sfs_cols
  )
}

# ============================================================================
# Internal: Hyperband algorithm
# ============================================================================

.hyperband <- function(search_space, data, type, sfs.dims,
                       max_epochs, eta, n_iter, seed, verbose) {

  s_max <- floor(log(max_epochs) / log(eta))

  if (verbose) cat(sprintf("PipeMaster:: max_epochs=%d, eta=%d, s_max=%d, %d brackets\n\n",
                           max_epochs, eta, s_max, s_max + 1))

  all_results <- data.frame(
    hp_string = character(), val_loss = numeric(),
    bracket = integer(), round = integer(),
    stringsAsFactors = FALSE
  )
  global_best_loss   <- Inf
  global_best_hp     <- NULL
  global_best_epochs <- 0L
  best_weights_path  <- tempfile("hb_best_weights_")

  for (iter in seq_len(n_iter)) {
    if (verbose && n_iter > 1) cat(sprintf("=== Iteration %d ===\n", iter))

    for (s in s_max:0) {
      n <- ceiling((s_max + 1) / (s + 1)) * as.integer(eta^s)
      r <- max_epochs / eta^s

      if (verbose) cat(sprintf("  Bracket %d | %d configs \u00d7 %d epochs\n",
                               s, n, round(r)))

      # Sample configs and build models
      set.seed(seed + iter * 1000 + s)
      configs <- lapply(seq_len(n), function(i) .sample.config(search_space))

      tensorflow::tf$random$set_seed(as.integer(seed + iter * 1000 + s))
      models <- lapply(seq_along(configs), function(k) {
        tryCatch(.build.nn(configs[[k]], data, type, sfs.dims),
                 error = function(e) {
                   if (verbose) cat(sprintf("    [warn] config %d failed to build: %s\n",
                                            k, conditionMessage(e)))
                   NULL
                 })
      })

      # Track how many epochs each model has been trained
      prev_epochs <- rep(0L, n)

      for (i in 0:s) {
        r_i <- round(r * eta^i)
        n_i <- max(1, floor(n / eta^i))
        n_keep <- max(1, ceiling(n_i / eta))

        # Train each surviving model
        val_losses <- rep(Inf, length(models))

        for (j in seq_along(models)) {
          if (is.null(models[[j]])) next

          tryCatch({
            history <- models[[j]] |> keras::fit(
              x = data$X_train, y = data$Y_train,
              validation_data = list(data$X_val, data$Y_val),
              epochs        = as.integer(r_i),
              initial_epoch = as.integer(prev_epochs[j]),
              batch_size    = as.integer(configs[[j]]$batch_size),
              callbacks     = list(
                keras::callback_early_stopping(monitor = "val_loss",
                                               patience = 10L,
                                               restore_best_weights = TRUE),
                keras::callback_reduce_lr_on_plateau(monitor = "val_loss",
                                                      patience = 5L,
                                                      factor = 0.5,
                                                      min_lr = 1e-6,
                                                      verbose = 0L)
              ),
              verbose = 0L
            )
            vl <- history$metrics$val_loss
            if (is.null(vl)) vl <- history$history$val_loss
            val_losses[j] <- min(unlist(vl))
          }, error = function(e) {
            if (verbose) cat(sprintf("    [warn] config %d train error: %s\n",
                                     j, conditionMessage(e)))
            val_losses[j] <<- Inf
          })
        }

        prev_epochs[seq_along(models)] <- r_i

        # Record results
        for (j in seq_along(models)) {
          if (!is.null(models[[j]]) && is.finite(val_losses[j])) {
            all_results <- rbind(all_results, data.frame(
              hp_string = .hp.to.string(configs[[j]], type),
              val_loss  = val_losses[j],
              bracket   = s,
              round     = i,
              stringsAsFactors = FALSE
            ))
          }
        }

        best_round_loss <- min(val_losses[is.finite(val_losses)])
        if (verbose)
          cat(sprintf("    Round %d: %d configs, %d ep \u2192 best val_loss=%.4f",
                      i, sum(!sapply(models, is.null)), r_i, best_round_loss))

        # Update global best — save weights so we can resume after k_clear_session
        best_idx <- which.min(val_losses)
        if (val_losses[best_idx] < global_best_loss) {
          global_best_loss   <- val_losses[best_idx]
          global_best_hp     <- configs[[best_idx]]
          global_best_epochs <- as.integer(r_i)
          tryCatch(
            keras::save_model_weights_tf(models[[best_idx]], best_weights_path),
            error = function(e) NULL
          )
        }

        # Prune: keep top n_keep
        if (i < s) {
          ranking <- order(val_losses)
          keep <- ranking[1:min(n_keep, length(ranking))]
          discard <- setdiff(seq_along(models), keep)

          if (verbose) cat(sprintf(" | pruning to %d\n", length(keep)))

          # Compact lists (subset first, then clear session)
          models      <- models[keep]
          configs     <- configs[keep]
          prev_epochs <- prev_epochs[keep]
          val_losses  <- val_losses[keep]
        } else {
          if (verbose) cat("\n")
        }
      }

      # Clear remaining models from this bracket
      rm(models)
      gc()
      tryCatch(keras::k_clear_session(), error = function(e) NULL)
    }
  }

  list(
    best_hp           = global_best_hp,
    best_val_loss     = global_best_loss,
    best_epochs       = global_best_epochs,
    best_weights_path = best_weights_path,
    all_results       = all_results
  )
}

# ============================================================================
# nn.predict — Posterior estimation via conformal prediction and bootstrap
# ============================================================================

#' Posterior Estimation via Conformal Prediction and Bootstrap
#'
#' Estimates posterior distributions for observed data using a trained neural
#' network from \code{tune.nn()}, with conformal prediction and/or bootstrap
#' methods for uncertainty quantification.
#'
#' @param tune.result list — output from \code{tune.nn()}.
#' @param observed named numeric vector or 1-row data.frame of observed summary
#'   statistics (or SFS bins).
#' @param reftable data.frame — original reference table (needed for bootstrap
#'   retraining and conformal calibration).
#' @param param.cols character vector — parameter column names.
#' @param type character — architecture type: \code{"sumstat"}, \code{"sfs1d"},
#'   or \code{"sfs2d"}. If NULL, uses the type stored in tune.result.
#' @param sfs.dims integer vector — for 2D CNN only. If NULL, uses tune.result.
#' @param method character — \code{"conformal"}, \code{"bootstrap"}, or both.
#' @param n_boot integer — number of bootstrap replicates (default 50).
#' @param n_ensemble integer — number of ensemble models for conformal (default 10).
#' @param cal.frac numeric — fraction of reftable used as calibration set (default 0.1).
#' @param max_epochs integer — max training epochs for conformal/bootstrap models (default 1000).
#' @param cores integer — number of parallel Rscript workers for bootstrap (default 1).
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{point_estimate}{named numeric vector — inverse-transformed point
#'     prediction from best model}
#'   \item{conformal}{matrix of posterior samples (n_samples x n_params), or NULL}
#'   \item{bootstrap}{matrix of posterior samples (n_boot x n_params), or NULL}
#'   \item{param_names}{character vector of parameter column names}
#' }
#'
#' @export
nn.predict <- function(tune.result, observed, reftable = NULL, param.cols = NULL,
                       type = NULL, sfs.dims = NULL,
                       method = c("conformal", "bootstrap"),
                       n_boot = 50, n_ensemble = 10, cal.frac = 0.1,
                       max_epochs = 1000, cores = 1, seed = 42, verbose = TRUE) {

  # --- Dependency check ---
  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("nn.predict() requires the 'keras' and 'tensorflow' R packages.")

  method <- match.arg(method, c("conformal", "bootstrap"), several.ok = TRUE)

  # Extract from tune.result
  best_hp    <- tune.result$best_hp
  best_model <- tune.result$best_model
  data       <- tune.result$data

  if (is.null(best_hp) || is.null(data))
    stop("tune.result must be output from tune.nn() with $best_hp and $data")

  if (is.null(type)) type <- tune.result$type
  if (is.null(type)) stop("type must be specified or present in tune.result")
  type <- match.arg(type, c("sumstat", "sfs1d", "sfs2d"))

  if (is.null(sfs.dims)) sfs.dims <- tune.result$sfs.dims

  if (("conformal" %in% method || "bootstrap" %in% method) && is.null(reftable))
    stop("reftable is required for conformal and bootstrap methods")
  if (!is.null(reftable) && is.null(param.cols))
    stop("param.cols is required when reftable is provided")

  param_names <- if (!is.null(param.cols)) param.cols else colnames(data$Y_train)

  # Coerce observed to numeric vector, reordering to match reftable stat column order
  if (is.data.frame(observed)) {
    stat_cols <- data$stat_cols
    if (!is.null(stat_cols) && all(stat_cols %in% colnames(observed))) {
      observed <- as.numeric(observed[1, stat_cols])
    } else {
      observed <- as.numeric(observed[1, ])
    }
  }

  # --- Header ---
  method_str <- paste(
    ifelse("conformal" %in% method, "Conformal", ""),
    ifelse(length(method) == 2, "+", ""),
    ifelse("bootstrap" %in% method, "Bootstrap", ""),
    sep = " "
  )
  method_str <- trimws(gsub("\\s+", " ", method_str))

  if (verbose) {
    cat(sprintf("PipeMaster:: nn.predict \u2014 %s\n", method_str))
  }

  # --- Prepare observed data ---
  X_obs <- .prep.observed(observed, data, type, sfs.dims)

  if (verbose) {
    n_obs_feat <- if (type == "sumstat") length(observed) else length(observed)
    cat(sprintf("PipeMaster:: Observed: %d %s\n", n_obs_feat,
                if (type == "sumstat") "summary statistics" else "SFS bins"))
  }

  # --- Point estimate from best model ---
  Y_z <- predict(best_model, X_obs, verbose = 0L)
  point_est <- as.numeric(.inv.transform(Y_z, data$target_mu, data$target_sd))
  names(point_est) <- param_names

  if (verbose) {
    est_str <- paste(sprintf("%s=%.0f", param_names, point_est), collapse = " ")
    cat(sprintf("PipeMaster:: Point estimate: %s\n", est_str))
  }

  # --- Conformal prediction ---
  conf_samples <- NULL
  if ("conformal" %in% method) {
    conf_samples <- .run.conformal(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      n_ensemble, cal.frac, max_epochs, seed, verbose
    )
    colnames(conf_samples) <- param_names
  }

  # --- Bootstrap ---
  boot_samples <- NULL
  if ("bootstrap" %in% method) {
    boot_samples <- .run.bootstrap(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      n_boot, max_epochs, cores, seed, verbose
    )
    colnames(boot_samples) <- param_names
  }

  if (verbose) cat("PipeMaster:: Done.\n")

  list(
    point_estimate = point_est,
    conformal      = conf_samples,
    bootstrap      = boot_samples,
    param_names    = param_names
  )
}

# ============================================================================
# Internal: prepare observed data for prediction
# ============================================================================

.prep.observed <- function(observed, data, type, sfs.dims) {
  feat_mu <- data$feat_mu
  feat_sd <- data$feat_sd

  if (type == "sumstat") {
    obs_aug <- c(observed, log1p(abs(observed)))
    X_obs <- matrix((obs_aug - feat_mu) / feat_sd, nrow = 1)
    X_obs[!is.finite(X_obs)] <- 0
  } else if (type == "sfs1d") {
    obs_log <- log1p(observed)
    obs_z <- (obs_log - feat_mu) / feat_sd
    obs_z[!is.finite(obs_z)] <- 0
    X_obs <- array(obs_z, dim = c(1L, length(obs_z), 1L))
  } else {
    obs_log <- log1p(observed)
    obs_z <- (obs_log - feat_mu) / feat_sd
    obs_z[!is.finite(obs_z)] <- 0
    X_obs <- array(obs_z, dim = c(1L, sfs.dims[1], sfs.dims[2], 1L))
  }
  X_obs
}

# ============================================================================
# Internal: normalize observed with arbitrary normalization params
# ============================================================================

.prep.observed.with <- function(observed, feat_mu, feat_sd, type, sfs.dims) {
  if (type == "sumstat") {
    obs_aug <- c(observed, log1p(abs(observed)))
    X_obs <- matrix((obs_aug - feat_mu) / feat_sd, nrow = 1)
    X_obs[!is.finite(X_obs)] <- 0
  } else if (type == "sfs1d") {
    obs_log <- log1p(observed)
    obs_z <- (obs_log - feat_mu) / feat_sd
    obs_z[!is.finite(obs_z)] <- 0
    X_obs <- array(obs_z, dim = c(1L, length(obs_z), 1L))
  } else {
    obs_log <- log1p(observed)
    obs_z <- (obs_log - feat_mu) / feat_sd
    obs_z[!is.finite(obs_z)] <- 0
    X_obs <- array(obs_z, dim = c(1L, sfs.dims[1], sfs.dims[2], 1L))
  }
  X_obs
}

# ============================================================================
# Internal: inverse transform predictions from normalized log-space
# ============================================================================

.inv.transform <- function(Y_z, target_mu, target_sd) {
  exp(t(t(Y_z) * target_sd + target_mu))
}

# ============================================================================
# Internal: write standalone model-builder R script for parallel workers
# ============================================================================

.write.builder.script <- function(filepath, type) {
  # Emit self-contained build_nn() that only needs keras/tensorflow loaded.
  # Only writes the architecture(s) actually needed.

  lines <- c(
    '# Auto-generated model builder for nn.predict parallel workers',
    'build_nn <- function(hp, data, type, sfs.dims) {',
    '  switch(type,',
    '    sumstat = build_resnet(hp, data),',
    '    sfs1d   = build_cnn1d(hp, data),',
    '    sfs2d   = build_cnn2d(hp, data, sfs.dims)',
    '  )',
    '}'
  )

  if (type == "sumstat" || type == "all") {
    lines <- c(lines, '',
    'build_resnet <- function(hp, data) {',
    '  n_features <- ncol(data$X_train)',
    '  n_targets  <- ncol(data$Y_train)',
    '  l2 <- regularizer_l2(hp$l2_reg)',
    '  res_block <- function(x, units) {',
    '    skip <- x',
    '    x <- x |>',
    '      layer_dense(units = as.integer(units), activation = "relu",',
    '                  kernel_regularizer = l2) |>',
    '      layer_batch_normalization() |>',
    '      layer_dense(units = as.integer(units), activation = "linear",',
    '                  kernel_regularizer = l2) |>',
    '      layer_batch_normalization()',
    '    x <- layer_add(list(x, skip))',
    '    layer_activation(x, activation = "relu")',
    '  }',
    '  inp <- layer_input(shape = n_features)',
    '  x <- inp |>',
    '    layer_dense(units = as.integer(hp$units_1), activation = "relu",',
    '                kernel_regularizer = l2) |>',
    '    layer_batch_normalization()',
    '  for (i in seq_len(hp$n_resblocks_1)) x <- res_block(x, hp$units_1)',
    '  x <- x |>',
    '    layer_dense(units = as.integer(hp$units_2), activation = "relu",',
    '                kernel_regularizer = l2) |>',
    '    layer_batch_normalization()',
    '  if (isTRUE(hp$use_dropout) && hp$dropout > 0)',
    '    x <- layer_dropout(x, rate = hp$dropout)',
    '  if (hp$n_resblocks_2 > 0)',
    '    for (i in seq_len(hp$n_resblocks_2)) x <- res_block(x, hp$units_2)',
    '  x <- x |>',
    '    layer_dense(units = as.integer(hp$units_3), activation = "relu",',
    '                kernel_regularizer = l2) |>',
    '    layer_batch_normalization()',
    '  if (isTRUE(hp$use_dropout) && hp$dropout > 0)',
    '    x <- layer_dropout(x, rate = hp$dropout)',
    '  out <- x |> layer_dense(units = n_targets, activation = "linear")',
    '  model <- keras_model(inp, out)',
    '  model |> compile(',
    '    loss = loss_huber(delta = hp$huber_delta),',
    '    optimizer = optimizer_adam(learning_rate = hp$learning_rate)',
    '  )',
    '  model',
    '}')
  }

  if (type == "sfs1d" || type == "all") {
    lines <- c(lines, '',
    'build_cnn1d <- function(hp, data) {',
    '  n_bins <- data$n_bins',
    '  n_targets <- ncol(data$Y_train)',
    '  l2 <- regularizer_l2(hp$l2_reg)',
    '  input <- layer_input(shape = c(n_bins, 1L))',
    '  x <- input',
    '  for (b in seq_len(hp$n_blocks)) {',
    '    filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))',
    '    ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))',
    '    x <- x |>',
    '      layer_conv_1d(filters = filters, kernel_size = ks,',
    '                    padding = "same", kernel_regularizer = l2) |>',
    '      layer_batch_normalization() |>',
    '      layer_activation("relu")',
    '    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {',
    '      skip <- x',
    '      x <- x |>',
    '        layer_conv_1d(filters = filters, kernel_size = ks,',
    '                      padding = "same", kernel_regularizer = l2) |>',
    '        layer_batch_normalization() |>',
    '        layer_activation("relu") |>',
    '        layer_conv_1d(filters = filters, kernel_size = ks,',
    '                      padding = "same", kernel_regularizer = l2) |>',
    '        layer_batch_normalization()',
    '      x <- layer_add(list(x, skip)) |> layer_activation("relu")',
    '    }',
    '  }',
    '  x <- x |> layer_global_average_pooling_1d()',
    '  for (d in seq_len(hp$n_dense)) {',
    '    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))',
    '    x <- x |>',
    '      layer_dense(units = units, activation = "relu",',
    '                  kernel_regularizer = l2) |>',
    '      layer_batch_normalization() |>',
    '      layer_dropout(rate = hp$dropout)',
    '  }',
    '  output <- x |> layer_dense(units = n_targets, activation = "linear")',
    '  model <- keras_model(input, output)',
    '  loss_fn <- if (identical(hp$loss, "huber")) loss_huber(delta = 1.0) else "mse"',
    '  model |> compile(loss = loss_fn,',
    '    optimizer = optimizer_adam(learning_rate = hp$learning_rate),',
    '    metrics = list("mae"))',
    '  model',
    '}')
  }

  if (type == "sfs2d" || type == "all") {
    lines <- c(lines, '',
    'build_cnn2d <- function(hp, data, sfs.dims) {',
    '  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]',
    '  n_targets <- ncol(data$Y_train)',
    '  l2 <- regularizer_l2(hp$l2_reg)',
    '  input <- layer_input(shape = c(dim1, dim2, 1L))',
    '  x <- input',
    '  for (b in seq_len(hp$n_blocks)) {',
    '    filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))',
    '    ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))',
    '    x <- x |>',
    '      layer_conv_2d(filters = filters, kernel_size = c(ks, ks),',
    '                    padding = "same", kernel_regularizer = l2) |>',
    '      layer_batch_normalization() |>',
    '      layer_activation("relu")',
    '    if (b == 1 || b == (hp$n_blocks %/% 2 + 1))',
    '      x <- layer_max_pooling_2d(x, pool_size = c(2L, 2L))',
    '    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {',
    '      skip <- x',
    '      x <- x |>',
    '        layer_conv_2d(filters = filters, kernel_size = c(ks, ks),',
    '                      padding = "same", kernel_regularizer = l2) |>',
    '        layer_batch_normalization() |>',
    '        layer_activation("relu") |>',
    '        layer_conv_2d(filters = filters, kernel_size = c(ks, ks),',
    '                      padding = "same", kernel_regularizer = l2) |>',
    '        layer_batch_normalization()',
    '      x <- layer_add(list(x, skip)) |> layer_activation("relu")',
    '    }',
    '  }',
    '  x <- x |> layer_global_average_pooling_2d()',
    '  for (d in seq_len(hp$n_dense)) {',
    '    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))',
    '    x <- x |>',
    '      layer_dense(units = units, activation = "relu",',
    '                  kernel_regularizer = l2) |>',
    '      layer_batch_normalization() |>',
    '      layer_dropout(rate = hp$dropout)',
    '  }',
    '  output <- x |> layer_dense(units = n_targets, activation = "linear")',
    '  model <- keras_model(input, output)',
    '  loss_fn <- if (identical(hp$loss, "huber")) loss_huber(delta = 1.0) else "mse"',
    '  model |> compile(loss = loss_fn,',
    '    optimizer = optimizer_adam(learning_rate = hp$learning_rate),',
    '    metrics = list("mae"))',
    '  model',
    '}')
  }

  writeLines(lines, filepath)
}

# ============================================================================
# Internal: prepare reftable split for conformal/bootstrap training
# ============================================================================

.prep.reftable.split <- function(reftable, param.cols, row_idx, type, sfs.dims) {
  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)

  targets <- as.matrix(reftable[row_idx, target_cols, drop = FALSE])

  if (type == "sumstat") {
    stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance))
    features_raw <- as.matrix(reftable[row_idx, stat_cols])
    features <- cbind(features_raw, log1p(abs(features_raw)))

    bad <- apply(features, 1, function(x) any(!is.finite(x)))
    if (any(bad)) { features <- features[!bad, ]; targets <- targets[!bad, , drop = FALSE] }
  } else {
    stat_cols <- grep("^sfs_", colnames(reftable), value = TRUE)
    sfs_raw <- as.matrix(reftable[row_idx, stat_cols])
    features <- log1p(sfs_raw)

    bad <- apply(features, 1, function(x) any(!is.finite(x)))
    if (any(bad)) { features <- features[!bad, ]; targets <- targets[!bad, , drop = FALSE] }
  }

  list(features = features, targets = targets)
}

# ============================================================================
# Internal: conformal prediction (sequential ensemble)
# ============================================================================

.run.conformal <- function(reftable, param.cols, observed, best_hp, type, sfs.dims,
                           n_ensemble, cal.frac, max_epochs, seed, verbose) {

  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params <- length(target_cols)

  n_total <- nrow(reftable)
  n_cal   <- floor(cal.frac * n_total)
  n_train <- n_total - n_cal

  if (verbose)
    cat(sprintf("\nPipeMaster:: Conformal prediction (%d ensemble \u00d7 %d cal samples)...\n",
                n_ensemble, n_cal))

  all_conf_samples <- list()

  for (ens_i in seq_len(n_ensemble)) {
    # Split reftable into train / calibration
    set.seed(seed + ens_i * 100)
    idx <- sample(n_total)
    tr_idx  <- idx[1:n_train]
    cal_idx <- idx[(n_train + 1):n_total]

    # Prepare train and calibration splits
    tr_split  <- .prep.reftable.split(reftable, param.cols, tr_idx, type, sfs.dims)
    cal_split <- .prep.reftable.split(reftable, param.cols, cal_idx, type, sfs.dims)

    feat_tr  <- tr_split$features
    targ_tr  <- tr_split$targets
    feat_cal <- cal_split$features
    params_cal <- cal_split$targets

    # Z-score features using train stats
    f_mu <- colMeans(feat_tr)
    f_sd <- apply(feat_tr, 2, sd); f_sd[f_sd == 0] <- 1
    X_tr  <- t((t(feat_tr)  - f_mu) / f_sd)
    X_cal <- t((t(feat_cal) - f_mu) / f_sd)

    # Log + Z-score targets using train stats
    Y_tr_log <- log(targ_tr)
    t_mu <- colMeans(Y_tr_log)
    t_sd <- apply(Y_tr_log, 2, sd); t_sd[t_sd == 0] <- 1
    Y_tr <- t((t(Y_tr_log) - t_mu) / t_sd)

    # Reshape for CNN architectures
    if (type == "sfs1d") {
      n_bins <- ncol(X_tr)
      dim(X_tr)  <- c(nrow(X_tr), n_bins, 1L)
      dim(X_cal) <- c(nrow(X_cal), n_bins, 1L)
    } else if (type == "sfs2d") {
      X_tr  <- array(X_tr,  dim = c(nrow(X_tr),  sfs.dims[1], sfs.dims[2], 1L))
      X_cal <- array(X_cal, dim = c(nrow(X_cal), sfs.dims[1], sfs.dims[2], 1L))
    }

    # Split train into train/val for early stopping
    n_tr_rows <- nrow(X_tr)
    n_va <- max(1L, floor(0.1 * n_tr_rows))
    va_rows <- seq_len(n_va)
    tr_rows <- seq(n_va + 1L, n_tr_rows)

    if (type == "sumstat") {
      X_va_s <- X_tr[va_rows, , drop = FALSE]
      X_tr_s <- X_tr[tr_rows, , drop = FALSE]
    } else if (type == "sfs1d") {
      X_va_s <- X_tr[va_rows, , , drop = FALSE]
      X_tr_s <- X_tr[tr_rows, , , drop = FALSE]
    } else {
      X_va_s <- X_tr[va_rows, , , , drop = FALSE]
      X_tr_s <- X_tr[tr_rows, , , , drop = FALSE]
    }
    Y_va_s <- Y_tr[va_rows, , drop = FALSE]
    Y_tr_s <- Y_tr[tr_rows, , drop = FALSE]

    # Build data list for .build.nn
    ens_data <- list(
      X_train = X_tr_s, X_val = X_va_s,
      Y_train = Y_tr_s, Y_val = Y_va_s,
      n_features = ncol(feat_tr),
      n_bins = if (type != "sumstat") ncol(feat_cal) else NULL,
      feat_mu = f_mu, feat_sd = f_sd,
      target_mu = t_mu, target_sd = t_sd
    )

    # Build and train model
    tensorflow::tf$random$set_seed(as.integer(seed + ens_i * 1000))
    model <- .build.nn(best_hp, ens_data, type, sfs.dims)

    bs <- as.integer(best_hp$batch_size)
    history <- model |> keras::fit(
      x = ens_data$X_train, y = ens_data$Y_train,
      validation_data = list(ens_data$X_val, ens_data$Y_val),
      epochs     = as.integer(max_epochs),
      batch_size = bs,
      callbacks  = list(
        keras::callback_early_stopping(monitor = "val_loss", patience = 30L,
                                       restore_best_weights = TRUE),
        keras::callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,
                                             factor = 0.5, min_lr = 1e-6, verbose = 0L)
      ),
      verbose = 0L
    )

    # Count epochs trained
    vl <- history$metrics$val_loss
    if (is.null(vl)) vl <- history$history$val_loss
    n_ep <- length(unlist(vl))

    # Predict on calibration set -> inverse transform
    cal_pred_z <- predict(model, X_cal, verbose = 0L)
    cal_pred   <- .inv.transform(cal_pred_z, t_mu, t_sd)

    # Residuals: true params - predicted params (original scale)
    cal_resid <- params_cal - cal_pred

    # Predict on observed (normalized with this split's params)
    X_obs <- .prep.observed.with(observed, f_mu, f_sd, type, sfs.dims)
    obs_pred_z <- predict(model, X_obs, verbose = 0L)
    obs_est <- as.numeric(.inv.transform(obs_pred_z, t_mu, t_sd))

    # Conformal samples: obs_est[j] + cal_resid[, j]
    ens_samples <- matrix(NA, nrow = nrow(cal_resid), ncol = n_params)
    for (j in seq_len(n_params))
      ens_samples[, j] <- obs_est[j] + cal_resid[, j]

    all_conf_samples[[ens_i]] <- ens_samples

    if (verbose)
      cat(sprintf("  Ensemble %d/%d \u2014 trained %d ep, %d conformal samples\n",
                  ens_i, n_ensemble, n_ep, nrow(ens_samples)))

    rm(model)
    gc()
    tryCatch(keras::k_clear_session(), error = function(e) NULL)
  }

  # Stack all ensemble samples
  do.call(rbind, all_conf_samples)
}

# ============================================================================
# Internal: bootstrap posterior estimation
# ============================================================================

.run.bootstrap <- function(reftable, param.cols, observed, best_hp, type, sfs.dims,
                           n_boot, max_epochs, cores, seed, verbose) {

  if (cores <= 1) {
    return(.run.bootstrap.sequential(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      n_boot, max_epochs, seed, verbose
    ))
  }

  # Parallel bootstrap via Rscript workers
  .run.bootstrap.parallel(
    reftable, param.cols, observed, best_hp, type, sfs.dims,
    n_boot, max_epochs, cores, seed, verbose
  )
}

# ============================================================================
# Internal: sequential bootstrap (cores = 1)
# ============================================================================

.run.bootstrap.sequential <- function(reftable, param.cols, observed, best_hp,
                                      type, sfs.dims, n_boot, max_epochs,
                                      seed, verbose) {

  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params <- length(target_cols)

  if (verbose)
    cat(sprintf("\nPipeMaster:: Bootstrap (%d replicates, sequential)...\n", n_boot))

  boot_matrix <- matrix(NA, nrow = n_boot, ncol = n_params)

  # Prepare full data once
  all_rows <- seq_len(nrow(reftable))
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims)
  features_all <- full_split$features
  targets_all  <- full_split$targets
  n_rows <- nrow(features_all)

  for (b in seq_len(n_boot)) {
    set.seed(seed + b)

    # Bootstrap sample with replacement
    boot_idx <- sample(n_rows, replace = TRUE)
    feat_boot <- features_all[boot_idx, ]
    targ_boot <- targets_all[boot_idx, , drop = FALSE]

    # Z-score features
    f_mu <- colMeans(feat_boot)
    f_sd <- apply(feat_boot, 2, sd); f_sd[f_sd == 0] <- 1
    X_boot <- t((t(feat_boot) - f_mu) / f_sd)

    # Log + Z-score targets
    Y_boot_log <- log(targ_boot)
    t_mu <- colMeans(Y_boot_log)
    t_sd <- apply(Y_boot_log, 2, sd); t_sd[t_sd == 0] <- 1
    Y_boot <- t((t(Y_boot_log) - t_mu) / t_sd)

    # Reshape for CNN
    if (type == "sfs1d") {
      n_bins <- ncol(X_boot)
      dim(X_boot) <- c(nrow(X_boot), n_bins, 1L)
    } else if (type == "sfs2d") {
      X_boot <- array(X_boot, dim = c(nrow(X_boot), sfs.dims[1], sfs.dims[2], 1L))
    }

    # Split off validation for early stopping
    n_b <- nrow(X_boot)
    n_va <- max(1L, floor(0.1 * n_b))
    va_rows <- seq_len(n_va)
    tr_rows <- seq(n_va + 1L, n_b)

    if (type == "sumstat") {
      X_va_b <- X_boot[va_rows, , drop = FALSE]
      X_tr_b <- X_boot[tr_rows, , drop = FALSE]
    } else if (type == "sfs1d") {
      X_va_b <- X_boot[va_rows, , , drop = FALSE]
      X_tr_b <- X_boot[tr_rows, , , drop = FALSE]
    } else {
      X_va_b <- X_boot[va_rows, , , , drop = FALSE]
      X_tr_b <- X_boot[tr_rows, , , , drop = FALSE]
    }
    Y_va_b <- Y_boot[va_rows, , drop = FALSE]
    Y_tr_b <- Y_boot[tr_rows, , drop = FALSE]

    boot_data <- list(
      X_train = X_tr_b, X_val = X_va_b,
      Y_train = Y_tr_b, Y_val = Y_va_b,
      n_features = ncol(feat_boot),
      n_bins = if (type != "sumstat") ncol(feat_boot) else NULL,
      feat_mu = f_mu, feat_sd = f_sd,
      target_mu = t_mu, target_sd = t_sd
    )

    tensorflow::tf$random$set_seed(as.integer(seed + b))
    model <- .build.nn(best_hp, boot_data, type, sfs.dims)

    bs <- as.integer(best_hp$batch_size)
    model |> keras::fit(
      x = boot_data$X_train, y = boot_data$Y_train,
      validation_data = list(boot_data$X_val, boot_data$Y_val),
      epochs     = as.integer(max_epochs),
      batch_size = bs,
      callbacks  = list(
        keras::callback_early_stopping(monitor = "val_loss", patience = 30L,
                                       restore_best_weights = TRUE),
        keras::callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,
                                             factor = 0.5, min_lr = 1e-6, verbose = 0L)
      ),
      verbose = 0L
    )

    # Predict on observed (normalized with this bootstrap's params)
    X_obs <- .prep.observed.with(observed, f_mu, f_sd, type, sfs.dims)
    obs_pred_z <- predict(model, X_obs, verbose = 0L)
    boot_matrix[b, ] <- as.numeric(.inv.transform(obs_pred_z, t_mu, t_sd))

    if (verbose && (b %% 10 == 0 || b == n_boot))
      cat(sprintf("  Progress: %d/%d\n", b, n_boot))

    rm(model)
    gc()
    tryCatch(keras::k_clear_session(), error = function(e) NULL)
  }

  boot_matrix
}

# ============================================================================
# Internal: parallel bootstrap via Rscript workers
# ============================================================================

.run.bootstrap.parallel <- function(reftable, param.cols, observed, best_hp,
                                    type, sfs.dims, n_boot, max_epochs,
                                    cores, seed, verbose) {

  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params <- length(target_cols)

  if (verbose)
    cat(sprintf("\nPipeMaster:: Bootstrap (%d replicates, %d cores)...\n", n_boot, cores))

  # Create temp working directory
  work_dir <- tempfile("nn_boot_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)

  # Prepare full data and save shared .RData
  all_rows <- seq_len(nrow(reftable))
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims)
  features_all <- full_split$features
  targets_all  <- full_split$targets
  n_rows <- nrow(features_all)

  # Save everything the worker needs
  shared_file <- file.path(work_dir, "shared_data.RData")
  save(features_all, targets_all, observed, best_hp,
       type, sfs.dims, n_rows, n_params, max_epochs, seed, target_cols,
       file = shared_file)

  # Save model-building functions so workers are self-contained (no PipeMaster dep)
  builder_script <- file.path(work_dir, "_build_model.R")
  .write.builder.script(builder_script, type)

  # Generate worker script
  worker_script <- file.path(work_dir, "_nn_worker.R")
  writeLines(c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    'Sys.setenv(TF_NUM_INTRAOP_THREADS = "1",',
    '           TF_NUM_INTEROP_THREADS = "1",',
    '           OMP_NUM_THREADS = "1")',
    '',
    'load("shared_data.RData")',
    '',
    'out_file <- file.path("results", sprintf("boot_%04d.csv", task_id))',
    'if (file.exists(out_file)) { cat("skip\\n"); q("no") }',
    '',
    'suppressPackageStartupMessages({',
    '  library(keras)',
    '  library(tensorflow)',
    '})',
    'source("_build_model.R")',
    '',
    'set.seed(seed + task_id)',
    'boot_idx <- sample(n_rows, replace = TRUE)',
    'feat_boot <- features_all[boot_idx, ]',
    'targ_boot <- targets_all[boot_idx, , drop = FALSE]',
    '',
    'f_mu <- colMeans(feat_boot)',
    'f_sd <- apply(feat_boot, 2, sd); f_sd[f_sd == 0] <- 1',
    'X_boot <- t((t(feat_boot) - f_mu) / f_sd)',
    '',
    'Y_boot_log <- log(targ_boot)',
    't_mu <- colMeans(Y_boot_log)',
    't_sd <- apply(Y_boot_log, 2, sd); t_sd[t_sd == 0] <- 1',
    'Y_boot <- t((t(Y_boot_log) - t_mu) / t_sd)',
    '',
    'if (type == "sfs1d") {',
    '  n_bins <- ncol(X_boot)',
    '  dim(X_boot) <- c(nrow(X_boot), n_bins, 1L)',
    '} else if (type == "sfs2d") {',
    '  X_boot <- array(X_boot, dim = c(nrow(X_boot), sfs.dims[1], sfs.dims[2], 1L))',
    '}',
    '',
    'n_b <- nrow(X_boot)',
    'n_va <- max(1L, floor(0.1 * n_b))',
    'va_rows <- seq_len(n_va)',
    'tr_rows <- seq(n_va + 1L, n_b)',
    '',
    'if (type == "sumstat") {',
    '  X_va_b <- X_boot[va_rows, , drop = FALSE]',
    '  X_tr_b <- X_boot[tr_rows, , drop = FALSE]',
    '} else if (type == "sfs1d") {',
    '  X_va_b <- X_boot[va_rows, , , drop = FALSE]',
    '  X_tr_b <- X_boot[tr_rows, , , drop = FALSE]',
    '} else {',
    '  X_va_b <- X_boot[va_rows, , , , drop = FALSE]',
    '  X_tr_b <- X_boot[tr_rows, , , , drop = FALSE]',
    '}',
    'Y_va_b <- Y_boot[va_rows, , drop = FALSE]',
    'Y_tr_b <- Y_boot[tr_rows, , drop = FALSE]',
    '',
    'boot_data <- list(',
    '  X_train = X_tr_b, X_val = X_va_b,',
    '  Y_train = Y_tr_b, Y_val = Y_va_b,',
    '  n_features = ncol(feat_boot),',
    '  n_bins = if (type != "sumstat") ncol(feat_boot) else NULL,',
    '  feat_mu = f_mu, feat_sd = f_sd,',
    '  target_mu = t_mu, target_sd = t_sd',
    ')',
    '',
    'tf$random$set_seed(as.integer(seed + task_id))',
    'model <- build_nn(best_hp, boot_data, type, sfs.dims)',
    '',
    'bs <- as.integer(best_hp$batch_size)',
    'model |> fit(',
    '  x = boot_data$X_train, y = boot_data$Y_train,',
    '  validation_data = list(boot_data$X_val, boot_data$Y_val),',
    '  epochs = as.integer(max_epochs), batch_size = bs,',
    '  callbacks = list(',
    '    callback_early_stopping(monitor = "val_loss", patience = 30L,',
    '                            restore_best_weights = TRUE),',
    '    callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,',
    '                                  factor = 0.5, min_lr = 1e-6, verbose = 0L)',
    '  ),',
    '  verbose = 0L',
    ')',
    '',
    '# Normalize observed and predict',
    'if (type == "sumstat") {',
    '  obs_aug <- c(observed, log1p(abs(observed)))',
    '  X_obs <- matrix((obs_aug - f_mu) / f_sd, nrow = 1)',
    '  X_obs[!is.finite(X_obs)] <- 0',
    '} else if (type == "sfs1d") {',
    '  obs_z <- (log1p(observed) - f_mu) / f_sd',
    '  obs_z[!is.finite(obs_z)] <- 0',
    '  X_obs <- array(obs_z, dim = c(1L, length(obs_z), 1L))',
    '} else {',
    '  obs_z <- (log1p(observed) - f_mu) / f_sd',
    '  obs_z[!is.finite(obs_z)] <- 0',
    '  X_obs <- array(obs_z, dim = c(1L, sfs.dims[1], sfs.dims[2], 1L))',
    '}',
    '',
    'obs_pred_z <- predict(model, X_obs, verbose = 0L)',
    'est <- as.numeric(exp(t(t(obs_pred_z) * t_sd + t_mu)))',
    '',
    'df <- as.data.frame(matrix(est, nrow = 1))',
    'colnames(df) <- target_cols',
    'write.csv(df, out_file, row.names = FALSE)',
    'cat(sprintf("  boot %d done\\n", task_id))',
    'k_clear_session()'
  ), worker_script)

  # Launch workers in batches
  if (verbose)
    cat(sprintf("  Launching %d Rscript workers, %d concurrent...\n", n_boot, cores))

  old_wd <- getwd()
  setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  # Process in batches of `cores`
  batch_starts <- seq(1, n_boot, by = cores)

  for (batch_start in batch_starts) {
    batch_end <- min(batch_start + cores - 1, n_boot)
    batch_ids <- batch_start:batch_end

    # Launch batch: background processes via system()
    for (id in batch_ids) {
      out_csv <- file.path("results", sprintf("boot_%04d.csv", id))
      if (file.exists(out_csv)) next
      cmd <- sprintf("Rscript %s %d > /dev/null 2>&1 &",
                     shQuote(worker_script), id)
      system(cmd, wait = FALSE)
    }

    # Poll until all output CSVs exist or timeout
    expected <- file.path("results", sprintf("boot_%04d.csv", batch_ids))
    timeout <- max_epochs * 10  # generous timeout in seconds
    elapsed <- 0
    while (elapsed < timeout) {
      if (all(file.exists(expected))) break
      Sys.sleep(2)
      elapsed <- elapsed + 2
    }

    if (verbose)
      cat(sprintf("  Progress: %d/%d\n", batch_end, n_boot))
  }

  # Collect results
  boot_matrix <- matrix(NA, nrow = n_boot, ncol = n_params)
  for (b in seq_len(n_boot)) {
    csv_file <- file.path(results_dir, sprintf("boot_%04d.csv", b))
    if (file.exists(csv_file)) {
      row <- read.csv(csv_file)
      boot_matrix[b, ] <- as.numeric(row[1, ])
    }
  }

  actual_failures <- sum(is.na(boot_matrix[, 1]))
  if (verbose)
    cat(sprintf("  Failures: %d/%d\n", actual_failures, n_boot))

  # Clean up temp files
  unlink(work_dir, recursive = TRUE)

  boot_matrix
}
