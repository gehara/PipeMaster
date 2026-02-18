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

  # --- Retrain best config for max_epochs ---
  if (verbose) {
    cat(sprintf("\nPipeMaster:: Best config: %s\n", .hp.to.string(hb$best_hp, type)))
    cat(sprintf("PipeMaster:: Retraining best for %d epochs...\n", max_epochs))
  }

  tensorflow::tf$random$set_seed(as.integer(seed))
  best_model <- .build.nn(hb$best_hp, data, type, sfs.dims)

  bs <- as.integer(hb$best_hp$batch_size)
  retrain_history <- best_model %>% keras::fit(
    x = data$X_train, y = data$Y_train,
    validation_data = list(data$X_val, data$Y_val),
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

  final_val_loss <- min(unlist(retrain_history$history$val_loss))
  if (verbose) cat(sprintf("PipeMaster:: Final val_loss: %.6f\n", final_val_loss))

  list(
    best_hp       = hb$best_hp,
    best_val_loss = final_val_loss,
    all_results   = hb$all_results,
    best_model    = best_model
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
    x <- x %>%
      keras::layer_dense(units = as.integer(units), activation = "relu",
                         kernel_regularizer = l2) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_dense(units = as.integer(units), activation = "linear",
                         kernel_regularizer = l2) %>%
      keras::layer_batch_normalization()
    x <- keras::layer_add(list(x, skip))
    x <- keras::layer_activation(x, activation = "relu")
    x
  }

  inp <- keras::layer_input(shape = n_features)

  # First dense + residual group
  x <- inp %>%
    keras::layer_dense(units = as.integer(hp$units_1), activation = "relu",
                       kernel_regularizer = l2) %>%
    keras::layer_batch_normalization()
  for (i in seq_len(hp$n_resblocks_1))
    x <- res_block(x, hp$units_1)

  # Middle transition + optional dropout
  x <- x %>%
    keras::layer_dense(units = as.integer(hp$units_2), activation = "relu",
                       kernel_regularizer = l2) %>%
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  # Second residual group
  if (hp$n_resblocks_2 > 0)
    for (i in seq_len(hp$n_resblocks_2))
      x <- res_block(x, hp$units_2)

  # Final dense + optional dropout
  x <- x %>%
    keras::layer_dense(units = as.integer(hp$units_3), activation = "relu",
                       kernel_regularizer = l2) %>%
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  out <- x %>% keras::layer_dense(units = n_targets, activation = "linear")

  model <- keras::keras_model(inp, out)
  model %>% keras::compile(
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

    x <- x %>%
      keras::layer_conv_1d(filters = filters, kernel_size = ks,
                           padding = "same", kernel_regularizer = l2) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_activation("relu")

    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {
      skip <- x
      x <- x %>%
        keras::layer_conv_1d(filters = filters, kernel_size = ks,
                             padding = "same", kernel_regularizer = l2) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_activation("relu") %>%
        keras::layer_conv_1d(filters = filters, kernel_size = ks,
                             padding = "same", kernel_regularizer = l2) %>%
        keras::layer_batch_normalization()
      x <- keras::layer_add(list(x, skip)) %>% keras::layer_activation("relu")
    }
  }

  x <- x %>% keras::layer_global_average_pooling_1d()

  for (d in seq_len(hp$n_dense)) {
    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
    x <- x %>%
      keras::layer_dense(units = units, activation = "relu",
                         kernel_regularizer = l2) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_dropout(rate = hp$dropout)
  }

  output <- x %>% keras::layer_dense(units = n_targets, activation = "linear")
  model <- keras::keras_model(input, output)

  loss_fn <- if (identical(hp$loss, "huber")) keras::loss_huber(delta = 1.0) else "mse"
  model %>% keras::compile(
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

    x <- x %>%
      keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                           padding = "same", kernel_regularizer = l2) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_activation("relu")

    # Max pooling on first and middle blocks
    if (b == 1 || b == (hp$n_blocks %/% 2 + 1))
      x <- keras::layer_max_pooling_2d(x, pool_size = c(2L, 2L))

    if (isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks) {
      skip <- x
      x <- x %>%
        keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                             padding = "same", kernel_regularizer = l2) %>%
        keras::layer_batch_normalization() %>%
        keras::layer_activation("relu") %>%
        keras::layer_conv_2d(filters = filters, kernel_size = c(ks, ks),
                             padding = "same", kernel_regularizer = l2) %>%
        keras::layer_batch_normalization()
      x <- keras::layer_add(list(x, skip)) %>% keras::layer_activation("relu")
    }
  }

  x <- x %>% keras::layer_global_average_pooling_2d()

  for (d in seq_len(hp$n_dense)) {
    units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
    x <- x %>%
      keras::layer_dense(units = units, activation = "relu",
                         kernel_regularizer = l2) %>%
      keras::layer_batch_normalization() %>%
      keras::layer_dropout(rate = hp$dropout)
  }

  output <- x %>% keras::layer_dense(units = n_targets, activation = "linear")
  model <- keras::keras_model(input, output)

  loss_fn <- if (identical(hp$loss, "huber")) keras::loss_huber(delta = 1.0) else "mse"
  model %>% keras::compile(
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
    target_mu = t_mu, target_sd = t_sd
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
    target_mu = t_mu, target_sd = t_sd
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
    target_mu = t_mu, target_sd = t_sd
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
  global_best_loss <- Inf
  global_best_hp   <- NULL

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
      models <- lapply(configs, function(hp) {
        tryCatch(.build.nn(hp, data, type, sfs.dims),
                 error = function(e) NULL)
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
            history <- models[[j]] %>% keras::fit(
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
            val_losses[j] <- min(unlist(history$history$val_loss))
          }, error = function(e) {
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

        # Update global best
        best_idx <- which.min(val_losses)
        if (val_losses[best_idx] < global_best_loss) {
          global_best_loss <- val_losses[best_idx]
          global_best_hp   <- configs[[best_idx]]
        }

        # Prune: keep top n_keep
        if (i < s) {
          ranking <- order(val_losses)
          keep <- ranking[1:min(n_keep, length(ranking))]
          discard <- setdiff(seq_along(models), keep)

          if (verbose) cat(sprintf(" | pruning to %d\n", length(keep)))

          # Free memory for eliminated models
          for (d in discard) {
            if (!is.null(models[[d]])) {
              tryCatch(keras::k_clear_session(), error = function(e) NULL)
              models[[d]] <- NULL
            }
          }

          # Compact lists
          models      <- models[keep]
          configs     <- configs[keep]
          prev_epochs <- prev_epochs[keep]
          val_losses  <- val_losses[keep]
        } else {
          if (verbose) cat("\n")
        }
      }

      # Clear remaining models from this bracket
      for (m in models) {
        if (!is.null(m))
          tryCatch(keras::k_clear_session(), error = function(e) NULL)
      }
    }
  }

  list(
    best_hp      = global_best_hp,
    best_val_loss = global_best_loss,
    all_results  = all_results
  )
}
