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
#' @param exclude.cols character vector — additional column names to exclude from
#'   features (e.g., other parameter columns not in \code{param.cols}). Required
#'   when estimating a single parameter from a reftable that contains multiple
#'   parameter columns, to prevent other parameters from leaking into the feature
#'   set. Default NULL (only \code{param.cols} and nuisance columns are excluded).
#' @param val.frac numeric — validation fraction (default 0.1).
#' @param cores integer — number of parallel Rscript workers for Hyperband rounds
#'   (default 1, sequential). Values > 1 train multiple configs simultaneously.
#'   Each worker spawns a separate R process that loads TensorFlow and a full
#'   copy of the training data, so RAM usage scales linearly with \code{cores}
#'   (~1.5 GB per worker). When \code{gpus > 0}, workers are assigned to GPUs
#'   via round-robin: multiple workers can share a GPU (they allocate VRAM
#'   incrementally), but too many workers per GPU may cause GPU out-of-memory
#'   errors.
#' @param gpus integer — number of GPUs to round-robin across workers (default 0,
#'   CPU-only). When \code{gpus = 0}, all workers run on CPU with GPUs disabled.
#'   Ignored when \code{cores = 1}.
#' @param gpu.threshold integer — maximum number of workers per GPU before
#'   switching to CPU-only for a round (default 4). When \code{gpus > 0},
#'   Hyperband rounds with more than \code{gpu.threshold * gpus} configs run
#'   on CPU (all \code{cores} workers, GPUs disabled), while rounds with fewer
#'   configs use the GPU. This allows early brackets (many configs, few epochs)
#'   to exploit CPU parallelism and later brackets (few configs, many epochs) to
#'   use GPU throughput. Ignored when \code{gpus = 0}.
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
                    exclude.cols = NULL,
                    val.frac = 0.1, cores = 1L, gpus = 0L,
                    gpu.threshold = 4L,
                    seed = 42, verbose = TRUE) {

  # --- Dependency check ---
  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("tune.nn() requires the 'keras' and 'tensorflow' R packages.\n",
         "Install with: install.packages(c('keras', 'tensorflow'))\n",
         "Then run: keras::install_keras()")

  # --- Memory guard for parallel mode ---
  if (cores > 1L) {
    avail_gb <- tryCatch({
      mem_info <- system("free -b 2>/dev/null", intern = TRUE)
      if (length(mem_info) >= 2) {
        fields <- as.numeric(strsplit(trimws(mem_info[2]), "\\s+")[[1]])
        fields[7] / 1e9  # "available" column (7th field)
      } else {
        NA_real_
      }
    }, error = function(e) NA_real_)

    if (!is.na(avail_gb)) {
      est_per_worker <- 1.5  # ~1.5 GB per TF worker (runtime + data copy)
      est_total <- cores * est_per_worker
      if (est_total > avail_gb * 0.85) {
        warning(sprintf(
          paste0("cores=%d workers may exceed available RAM ",
                 "(%.1f GB free, ~%.0f GB estimated). ",
                 "This can cause swapping and severe slowdowns. ",
                 "Reduce cores if you experience memory issues."),
          cores, avail_gb, est_total),
          call. = FALSE)
      }
    }

  }

  # --- Enable GPU memory growth (prevent TF from grabbing all VRAM) ---
  # Only in sequential mode; parallel workers set CUDA_VISIBLE_DEVICES per-process
  if (cores <= 1L) {
    tryCatch({
      tf_gpus <- tensorflow::tf$config$list_physical_devices("GPU")
      for (gpu in tf_gpus)
        tensorflow::tf$config$experimental$set_memory_growth(gpu, TRUE)
    }, error = function(e) NULL)
  }

  type <- match.arg(type)

  if (type == "sfs2d" && (is.null(sfs.dims) || length(sfs.dims) != 2))
    stop("sfs.dims must be c(dim1, dim2) for type='sfs2d'")

  # --- Search space ---
  ss <- if (is.null(search_space)) .default.search.space(type) else search_space

  # --- Prepare data ---
  if (verbose) cat(sprintf("PipeMaster:: tune.nn \u2014 Hyperband (%s)\n",
                           switch(type, sumstat = "ResNet", sfs1d = "1D CNN", sfs2d = "2D CNN")))

  data <- .prep.data(reftable, param.cols, type, sfs.dims, exclude.cols, val.frac, seed)

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
    cores        = cores,
    gpus         = gpus,
    gpu.threshold = gpu.threshold,
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
    keras::load_model_weights_tf(best_model, file.path(hb$best_weights_path, "ckpt"))
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
    sfs.dims      = sfs.dims,
    exclude.cols  = exclude.cols
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

.build.nn <- function(hp, data, type, sfs.dims, mc_dropout = FALSE) {
  switch(type,
    sumstat = .build.resnet(hp, data, mc_dropout = mc_dropout),
    sfs1d   = .build.cnn1d(hp, data, mc_dropout = mc_dropout),
    sfs2d   = .build.cnn2d(hp, data, sfs.dims, mc_dropout = mc_dropout)
  )
}

# --- ResNet for summary statistics ---
.build.resnet <- function(hp, data, mc_dropout = FALSE) {
  n_features <- ncol(data$X_train)
  n_targets  <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  # Dropout layer: standard or always-on (for MC dropout inference)
  .drop <- if (mc_dropout) {
    function(x, rate) {
      r <- rate
      keras::layer_lambda(x, f = function(x) tensorflow::tf$nn$dropout(x, rate = r))
    }
  } else {
    function(x, rate) keras::layer_dropout(x, rate = rate)
  }

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
    x <- .drop(x, hp$dropout)

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
    x <- .drop(x, hp$dropout)

  out <- x |> keras::layer_dense(units = n_targets, activation = "linear")

  model <- keras::keras_model(inp, out)
  model |> keras::compile(
    loss      = keras::loss_huber(delta = hp$huber_delta),
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate)
  )
  model
}

# --- 1D CNN for single-population SFS ---
.build.cnn1d <- function(hp, data, mc_dropout = FALSE) {
  n_bins    <- data$n_bins
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  .drop <- if (mc_dropout && hp$dropout > 0) {
    function(x, rate) {
      r <- rate
      keras::layer_lambda(x, f = function(x) tensorflow::tf$nn$dropout(x, rate = r))
    }
  } else {
    function(x, rate) keras::layer_dropout(x, rate = rate)
  }

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
      .drop(hp$dropout)
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
.build.cnn2d <- function(hp, data, sfs.dims, mc_dropout = FALSE) {
  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  .drop <- if (mc_dropout && hp$dropout > 0) {
    function(x, rate) {
      r <- rate
      keras::layer_lambda(x, f = function(x) tensorflow::tf$nn$dropout(x, rate = r))
    }
  } else {
    function(x, rate) keras::layer_dropout(x, rate = rate)
  }

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
      .drop(hp$dropout)
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

.prep.data <- function(reftable, param.cols, type, sfs.dims, exclude.cols, val.frac, seed) {
  set.seed(seed)

  if (type == "sumstat") {
    .prep.sumstat(reftable, param.cols, exclude.cols, val.frac)
  } else if (type == "sfs1d") {
    .prep.sfs1d(reftable, param.cols, val.frac)
  } else {
    .prep.sfs2d(reftable, param.cols, sfs.dims, val.frac)
  }
}

# --- sumstat: augment features with log1p(abs(x)), Z-score, log-targets ---
.prep.sumstat <- function(reftable, param.cols, exclude.cols, val.frac) {
  nuisance <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance, exclude.cols))
  target_cols <- setdiff(param.cols, nuisance)

  targets <- as.matrix(reftable[, target_cols, drop = FALSE])
  features_raw <- as.matrix(reftable[, stat_cols, drop = FALSE])

  # Feature augmentation
  features <- cbind(features_raw, log1p(abs(features_raw)))

  # Remove bad rows
  bad <- apply(features, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { features <- features[!bad, , drop = FALSE]; targets <- targets[!bad, , drop = FALSE] }

  # Split
  n <- nrow(features)
  idx <- sample(n)
  n_val <- floor(val.frac * n)
  n_train <- n - n_val

  tr <- idx[1:n_train]; va <- idx[(n_train + 1):n]

  # Z-score features
  feat_mu <- colMeans(features[tr, , drop = FALSE])
  feat_sd <- apply(features[tr, , drop = FALSE], 2, sd); feat_sd[feat_sd == 0] <- 1
  zscore_feat <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore_feat(features[tr, , drop = FALSE])
  X_val   <- zscore_feat(features[va, , drop = FALSE])

  # Log-transform + Z-score targets
  Y_tr_log <- log(targets[tr, , drop = FALSE]); Y_va_log <- log(targets[va, , drop = FALSE])
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

  targets <- as.matrix(reftable[, param.cols, drop = FALSE])
  sfs_raw <- as.matrix(reftable[, sfs_cols])

  # Log1p transform
  sfs_log <- log1p(sfs_raw)

  bad <- apply(sfs_log, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { sfs_log <- sfs_log[!bad, ]; targets <- targets[!bad, , drop = FALSE] }

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
  Y_tr_log <- log(targets[tr, , drop = FALSE]); Y_va_log <- log(targets[va, , drop = FALSE])
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

  targets <- as.matrix(reftable[, param.cols, drop = FALSE])
  sfs_raw <- as.matrix(reftable[, sfs_cols])

  sfs_log <- log1p(sfs_raw)
  bad <- apply(sfs_log, 1, function(x) any(!is.finite(x)))
  if (any(bad)) { sfs_log <- sfs_log[!bad, ]; targets <- targets[!bad, , drop = FALSE] }

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
  Y_tr_log <- log(targets[tr, , drop = FALSE]); Y_va_log <- log(targets[va, , drop = FALSE])
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
                       max_epochs, eta, n_iter, cores, gpus,
                       gpu.threshold, seed, verbose) {

  s_max <- min(floor(log(max_epochs) / log(eta)), 3L)  # cap at 3 (max 27 configs/bracket)

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

      # Sample configs
      set.seed(seed + iter * 1000 + s)
      configs <- lapply(seq_len(n), function(i) .sample.config(search_space))

      # Track how many epochs each model has been trained
      prev_epochs <- rep(0L, n)
      # Temp dir for per-config weight files (avoids GPU OOM from keeping all models)
      weight_dir <- tempfile("hb_weights_")
      dir.create(weight_dir, recursive = TRUE)

      for (i in 0:s) {
        r_i <- round(r * eta^i)
        n_i <- max(1, floor(n / eta^i))
        n_keep <- max(1, ceiling(n_i / eta))

        val_losses <- rep(Inf, length(configs))

        if (cores > 1L && length(configs) > 1L) {
          # PARALLEL: choose CPU vs GPU based on config count per round
          if (gpus > 0L && length(configs) <= gpus * gpu.threshold) {
            round_cores <- min(cores, gpus * gpu.threshold)
            round_gpus  <- gpus
          } else {
            round_cores <- cores
            round_gpus  <- 0L
          }

          val_losses <- .hyperband.round.parallel(
            configs, data, type, sfs.dims, r_i, prev_epochs,
            weight_dir, seed, iter, s, round_cores, round_gpus, verbose)

          # Update global best from parallel results
          for (j in seq_along(configs)) {
            if (val_losses[j] < global_best_loss) {
              global_best_loss   <- val_losses[j]
              global_best_hp     <- configs[[j]]
              global_best_epochs <- as.integer(r_i)
              wpath <- file.path(weight_dir, sprintf("cfg_%d", j))
              if (dir.exists(wpath)) {
                unlink(best_weights_path, recursive = TRUE)
                dir.create(best_weights_path, recursive = TRUE, showWarnings = FALSE)
                file.copy(list.files(wpath, full.names = TRUE),
                          best_weights_path, recursive = TRUE)
              }
            }
          }
        } else {
          # SEQUENTIAL: train each config one at a time
          for (j in seq_along(configs)) {
            tryCatch({
              tensorflow::tf$random$set_seed(as.integer(seed + iter * 1000 + s + j))
              model <- .build.nn(configs[[j]], data, type, sfs.dims)

              # Load saved weights from previous round (if any)
              wpath <- file.path(weight_dir, sprintf("cfg_%d", j))
              wfile <- file.path(wpath, "ckpt")
              if (prev_epochs[j] > 0 && dir.exists(wpath)) {
                keras::load_model_weights_tf(model, wfile)
              }

              history <- model |> keras::fit(
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

              # Save weights for this config (overwrite previous round's)
              dir.create(wpath, recursive = TRUE, showWarnings = FALSE)
              keras::save_model_weights_tf(model, wfile)

              # Update global best
              if (val_losses[j] < global_best_loss) {
                global_best_loss   <- val_losses[j]
                global_best_hp     <- configs[[j]]
                global_best_epochs <- as.integer(r_i)
                tryCatch({
                  dir.create(best_weights_path, recursive = TRUE, showWarnings = FALSE)
                  keras::save_model_weights_tf(model, file.path(best_weights_path, "ckpt"))
                },
                  error = function(e) NULL
                )
              }

              # Free GPU memory immediately
              rm(model); gc()
              tryCatch(keras::k_clear_session(), error = function(e) NULL)

            }, error = function(e) {
              if (verbose) cat(sprintf("    [warn] config %d error: %s\n",
                                       j, conditionMessage(e)))
              val_losses[j] <<- Inf
              tryCatch({ rm(model); gc(); keras::k_clear_session() },
                       error = function(e2) NULL)
            })
          }
        }

        prev_epochs[seq_along(configs)] <- r_i

        # Record results
        for (j in seq_along(configs)) {
          if (is.finite(val_losses[j])) {
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
                      i, length(configs), r_i, best_round_loss))

        # Prune: keep top n_keep
        if (i < s) {
          ranking <- order(val_losses)
          keep <- ranking[1:min(n_keep, length(ranking))]

          if (verbose) cat(sprintf(" | pruning to %d\n", length(keep)))

          # Remove weight files for discarded configs
          discard <- setdiff(seq_along(configs), keep)
          for (d in discard) {
            wpath <- file.path(weight_dir, sprintf("cfg_%d", d))
            unlink(wpath, recursive = TRUE)
          }

          # Renumber surviving configs/weights 1..n_keep
          new_configs     <- configs[keep]
          new_prev_epochs <- prev_epochs[keep]
          for (k in seq_along(keep)) {
            old_path <- file.path(weight_dir, sprintf("cfg_%d", keep[k]))
            new_path <- file.path(weight_dir, sprintf("cfg_new_%d", k))
            if (dir.exists(old_path)) file.rename(old_path, new_path)
          }
          for (k in seq_along(keep)) {
            old_path <- file.path(weight_dir, sprintf("cfg_new_%d", k))
            new_path <- file.path(weight_dir, sprintf("cfg_%d", k))
            if (dir.exists(old_path)) file.rename(old_path, new_path)
          }
          configs     <- new_configs
          prev_epochs <- new_prev_epochs
        } else {
          if (verbose) cat("\n")
        }
      }

      # Clean up bracket weight files
      unlink(weight_dir, recursive = TRUE)
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
# Internal: parallel training of one Hyperband round via Rscript worker pool
# ============================================================================

.hyperband.round.parallel <- function(configs, data, type, sfs.dims, r_i,
                                      prev_epochs, weight_dir, seed, iter, s,
                                      cores, gpus, verbose) {
  n_configs <- length(configs)
  val_losses <- rep(Inf, n_configs)

  # Create temp working directory

  work_dir <- tempfile("hb_round_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)
  weights_work <- file.path(work_dir, "weights")
  dir.create(weights_work)

  # Save shared training data
  X_train <- data$X_train
  X_val   <- data$X_val
  Y_train <- data$Y_train
  Y_val   <- data$Y_val
  n_features <- if (type == "sumstat") ncol(data$X_train) else data$n_features
  n_bins     <- data$n_bins
  shared_file <- file.path(work_dir, "shared_data.RData")
  save(X_train, X_val, Y_train, Y_val,
       type, sfs.dims, n_features, n_bins,
       file = shared_file)

  # Save round-specific task data (configs, epochs, seeds)
  round_configs <- configs
  round_r_i     <- as.integer(r_i)
  round_prev_epochs <- as.integer(prev_epochs)
  round_seeds <- vapply(seq_len(n_configs), function(j) {
    as.integer(seed + iter * 1000 + s + j)
  }, integer(1))
  round_file <- file.path(work_dir, "round_tasks.RData")
  save(round_configs, round_r_i, round_prev_epochs, round_seeds,
       file = round_file)

  # Copy existing weights from weight_dir/cfg_* → work_dir/weights/cfg_*
  for (j in seq_len(n_configs)) {
    src <- file.path(weight_dir, sprintf("cfg_%d", j))
    if (prev_epochs[j] > 0 && dir.exists(src)) {
      dst <- file.path(weights_work, sprintf("cfg_%d", j))
      dir.create(dst, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
    }
  }

  # Write scripts
  .write.builder.script(file.path(work_dir, "_build_model.R"), type)
  .write.hyperband.worker.script(file.path(work_dir, "_hb_worker.R"))

  # Build task list for the pool
  tasks <- lapply(seq_len(n_configs), function(j) {
    list(
      script = "_hb_worker.R",
      id     = j,
      result = sprintf("results/hb_%04d.csv", j),
      prefix = "hb"
    )
  })

  if (verbose)
    cat(sprintf("    [parallel] %d configs on %d cores%s\n",
                n_configs, cores,
                if (gpus > 0) sprintf(", %d GPUs", gpus) else ""))

  # Launch pool
  pool_result <- .launch.rscript.pool(tasks, cores, work_dir,
                                       timeout_per_task = r_i * 10,
                                       gpus = gpus, verbose = verbose,
                                       max_retries = 1L)

  # Collect val_losses from result CSVs
  for (j in seq_len(n_configs)) {
    csv_file <- file.path(results_dir, sprintf("hb_%04d.csv", j))
    if (file.exists(csv_file)) {
      row <- tryCatch(read.csv(csv_file), error = function(e) NULL)
      if (!is.null(row) && "val_loss" %in% names(row))
        val_losses[j] <- as.numeric(row$val_loss[1])
    }
  }

  # Copy updated weights back: work_dir/weights/cfg_* → weight_dir/cfg_*
  for (j in seq_len(n_configs)) {
    src <- file.path(weights_work, sprintf("cfg_%d", j))
    if (dir.exists(src)) {
      dst <- file.path(weight_dir, sprintf("cfg_%d", j))
      unlink(dst, recursive = TRUE)
      dir.create(dst, recursive = TRUE, showWarnings = FALSE)
      file.copy(list.files(src, full.names = TRUE), dst, recursive = TRUE)
    }
  }

  # Clean up
  unlink(work_dir, recursive = TRUE)

  val_losses
}

# ============================================================================
# save / load tune.nn results (keras model needs special serialization)
# ============================================================================

#' Save tune.nn Result to Disk
#'
#' Saves the output of \code{tune.nn()} so it can be loaded on another machine.
#' The keras model is serialized separately using \code{keras::save_model_tf()}.
#'
#' @param tune.result list — output from \code{tune.nn()}.
#' @param path character — directory path where files will be saved (created if needed).
#'
#' @export
save.tune.result <- function(tune.result, path) {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("save.tune.result() requires the 'keras' package.")

  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  # Save keras model
  model_dir <- file.path(path, "best_model")
  keras::save_model_tf(tune.result$best_model, model_dir)

  # Save everything else as RData
  result_no_model <- tune.result
  result_no_model$best_model <- NULL
  saveRDS(result_no_model, file.path(path, "tune_result.rds"))

  cat(sprintf("PipeMaster:: Saved tune.nn result to %s\n", path))
}

#' Load tune.nn Result from Disk
#'
#' Loads a tune.nn() result previously saved with \code{save.tune.result()}.
#'
#' @param path character — directory path where files were saved.
#'
#' @return A list identical to the output of \code{tune.nn()}.
#'
#' @export
load.tune.result <- function(path) {
  if (!requireNamespace("keras", quietly = TRUE))
    stop("load.tune.result() requires the 'keras' package.")

  rds_file  <- file.path(path, "tune_result.rds")
  model_dir <- file.path(path, "best_model")

  if (!file.exists(rds_file))
    stop("tune_result.rds not found in: ", path)
  if (!dir.exists(model_dir))
    stop("best_model/ directory not found in: ", path)

  result <- readRDS(rds_file)
  result$best_model <- keras::load_model_tf(model_dir)

  cat(sprintf("PipeMaster:: Loaded tune.nn result from %s\n", path))
  result
}

# ============================================================================
# nn.predict — Posterior estimation via conformal prediction and bootstrap
# ============================================================================

#' Posterior Estimation via Neural Network Uncertainty Quantification
#'
#' Estimates posterior distributions for observed data using a trained neural
#' network from \code{tune.nn()}, with multiple methods for uncertainty
#' quantification: conformal prediction, bootstrap, MC dropout, and quantile
#' regression.
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
#' @param method character — one or more of \code{"conformal"}, \code{"bootstrap"},
#'   \code{"mc_dropout"}, \code{"quantile"}.
#' @param n_boot integer — number of bootstrap replicates (default 20).
#' @param n_ensemble integer — number of ensemble models for conformal (default 1).
#' @param cal.frac numeric — fraction of reftable used as calibration set (default 0.1).
#' @param n_mc integer — number of MC dropout forward passes (default 1000).
#'   Only used when \code{method} includes \code{"mc_dropout"}.
#' @param mc_dropout_rate numeric — dropout rate for MC dropout (default NULL =
#'   use tuned rate). When set (e.g. 0.1), the model is retrained with this
#'   dropout rate before running forward passes. Useful when the tuned dropout
#'   is too low to produce meaningful posterior spread.
#' @param n_quantiles integer — number of evenly-spaced quantile bins for quantile
#'   regression (default 19). Higher values (e.g. 50, 100) produce smoother
#'   posteriors. Ignored if \code{q_probs} is provided explicitly.
#' @param q_probs numeric vector — custom quantile probabilities. If NULL
#'   (default), generated automatically from \code{n_quantiles}. Only used when
#'   \code{method} includes \code{"quantile"}.
#' @param max_epochs integer — max training epochs for conformal/bootstrap/quantile
#'   models (default 1000).
#' @param cores integer — number of parallel Rscript workers (default 1).
#' @param gpus integer — number of GPUs for parallel workers (default 0 = CPU-only).
#'   When \code{gpus > 0}, workers are assigned GPUs round-robin via
#'   \code{CUDA_VISIBLE_DEVICES}. Each worker sets
#'   \code{TF_FORCE_GPU_ALLOW_GROWTH=true} so multiple workers sharing a GPU
#'   don't OOM. \code{cores} and \code{gpus} are independent — \code{cores}
#'   controls max concurrent workers, \code{gpus} controls GPU assignment.
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return An object of class \code{"nn.posterior"} (a list) with:
#' \describe{
#'   \item{point_estimate}{named numeric vector — inverse-transformed point
#'     prediction from best model}
#'   \item{conformal}{matrix of posterior samples (n_samples x n_params), or NULL}
#'   \item{bootstrap}{matrix of posterior samples (n_boot x n_params), or NULL}
#'   \item{mc_dropout}{matrix of posterior samples (n_mc x n_params), or NULL}
#'   \item{quantile}{matrix of quantile values (n_quantiles x n_params), or NULL}
#'   \item{q_probs}{numeric vector of quantile probabilities used, or NULL}
#'   \item{param_names}{character vector of parameter column names}
#' }
#' Use \code{summary()} to get a table of median, mean, and quantiles.
#'
#' @export
nn.predict <- function(tune.result, observed, reftable = NULL, param.cols = NULL,
                       type = NULL, sfs.dims = NULL,
                       method = c("conformal", "bootstrap", "mc_dropout", "quantile"),
                       n_boot = 20, n_ensemble = 1, cal.frac = 0.1,
                       n_mc = 1000L, mc_dropout_rate = NULL,
                       n_quantiles = 19L, q_probs = NULL,
                       max_epochs = 1000, cores = 1, gpus = 0, seed = 42,
                       verbose = TRUE) {

  # --- Dependency check ---
  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("nn.predict() requires the 'keras' and 'tensorflow' R packages.")

  method <- match.arg(method, c("conformal", "bootstrap", "mc_dropout", "quantile"),
                      several.ok = TRUE)

  # Generate q_probs from n_quantiles if not provided
  if (is.null(q_probs)) {
    q_probs <- seq(0.01, 0.99, length.out = as.integer(n_quantiles))
  }

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
  exclude.cols <- tune.result$exclude.cols

  needs_reftable <- any(c("conformal", "bootstrap", "quantile") %in% method)
  if (needs_reftable && is.null(reftable))
    stop("reftable is required for conformal, bootstrap, and quantile methods")
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
  method_labels <- c(conformal = "Conformal", bootstrap = "Bootstrap",
                     mc_dropout = "MC Dropout", quantile = "Quantile Regression")
  method_str <- paste(method_labels[method], collapse = " + ")

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

  # --- Uncertainty quantification ---
  conf_samples   <- NULL
  boot_samples   <- NULL
  mcdrop_samples <- NULL
  quant_matrix   <- NULL
  quant_probs    <- NULL
  do_conformal  <- "conformal" %in% method
  do_bootstrap  <- "bootstrap" %in% method
  do_mc_dropout <- "mc_dropout" %in% method
  do_quantile   <- "quantile" %in% method

  # --- MC Dropout ---
  if (do_mc_dropout) {
    tf <- tensorflow::tf
    mc_hp <- best_hp

    # Override dropout rate if requested
    needs_retrain <- FALSE
    if (!is.null(mc_dropout_rate)) {
      mc_hp <- best_hp
      mc_hp$dropout     <- mc_dropout_rate
      mc_hp$use_dropout <- TRUE
      needs_retrain <- (mc_dropout_rate != best_hp$dropout)
    }

    # Check for dropout
    has_dropout <- (isTRUE(mc_hp$use_dropout) && mc_hp$dropout > 0) ||
                   (!is.null(mc_hp$dropout) && mc_hp$dropout > 0)
    if (!has_dropout)
      warning("MC Dropout requested but dropout rate is 0. ",
              "Predictions will have no stochastic variation. ",
              "Set mc_dropout_rate (e.g. 0.1) to override.", call. = FALSE)

    if (verbose)
      cat(sprintf("\nPipeMaster:: MC Dropout (%d forward passes, dropout=%.3f%s)...\n",
                  n_mc, mc_hp$dropout,
                  if (needs_retrain) ", retraining" else ""))

    # Build MC model: same architecture but with always-on dropout (lambda layers)
    # so predict() keeps BN in inference mode while dropout stays active
    mc_model <- .build.nn(mc_hp, data, type, sfs.dims, mc_dropout = TRUE)

    if (needs_retrain) {
      # Retrain with the new dropout rate using tune.nn's data splits
      tensorflow::tf$random$set_seed(as.integer(seed + 555L))
      bs <- as.integer(mc_hp$batch_size)

      if (verbose) cat("PipeMaster:: Retraining model with new dropout rate...\n")

      mc_model |> keras::fit(
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

      if (verbose) cat("PipeMaster:: Retrain done. Running forward passes...\n")
    } else {
      # Same dropout rate — just copy weights from best model
      mc_model$set_weights(best_model$get_weights())
    }

    n_params_mc <- length(param_names)
    mc_raw <- matrix(NA_real_, nrow = as.integer(n_mc), ncol = n_params_mc)
    for (i in seq_len(n_mc)) {
      pred_z <- predict(mc_model, X_obs, verbose = 0L)
      mc_raw[i, ] <- as.numeric(.inv.transform(pred_z, data$target_mu, data$target_sd))
    }

    rm(mc_model); gc()

    mcdrop_samples <- mc_raw
    colnames(mcdrop_samples) <- param_names

    if (verbose)
      cat(sprintf("PipeMaster:: MC Dropout done — %d samples\n", n_mc))
  }

  # --- Quantile Regression (trains 1 model) ---
  if (do_quantile) {
    quant_probs <- q_probs
    quant_matrix <- .run.quantile.regression(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      q_probs, max_epochs, seed, verbose, data, best_model,
      exclude.cols = exclude.cols
    )
    colnames(quant_matrix) <- param_names
  }

  if (cores > 1 && (do_conformal || do_bootstrap)) {
    # Unified parallel dispatch via priority pool
    pool_out <- .run.parallel.pool(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      do_conformal, do_bootstrap, n_ensemble, n_boot,
      cal.frac, max_epochs, cores, gpus, seed, verbose,
      point_est = point_est, exclude.cols = exclude.cols
    )
    if (do_conformal && !is.null(pool_out$conformal)) {
      conf_samples <- pool_out$conformal
      colnames(conf_samples) <- param_names
    }
    if (do_bootstrap && !is.null(pool_out$bootstrap)) {
      boot_samples <- pool_out$bootstrap
      colnames(boot_samples) <- param_names
    }
  } else {
    # Sequential fallback (cores <= 1)
    if (do_conformal) {
      conf_samples <- .run.conformal.sequential(
        reftable, param.cols, observed, best_hp, type, sfs.dims,
        n_ensemble, cal.frac, max_epochs, seed, verbose,
        point_est = point_est, exclude.cols = exclude.cols
      )
      colnames(conf_samples) <- param_names
    }
    if (do_bootstrap) {
      boot_samples <- .run.bootstrap.sequential(
        reftable, param.cols, observed, best_hp, type, sfs.dims,
        n_boot, max_epochs, seed, verbose,
        exclude.cols = exclude.cols
      )
      colnames(boot_samples) <- param_names
    }
  }

  # Clip to prior range — use reftable parameter columns as bounds
  if (!is.null(reftable) && !is.null(param_names)) {
    for (j in seq_along(param_names)) {
      p <- param_names[j]
      if (p %in% colnames(reftable)) {
        lo <- min(reftable[[p]])
        hi <- max(reftable[[p]])
        if (!is.null(conf_samples))   conf_samples[, j]   <- pmax(lo, pmin(hi, conf_samples[, j]))
        if (!is.null(boot_samples))   boot_samples[, j]   <- pmax(lo, pmin(hi, boot_samples[, j]))
        if (!is.null(mcdrop_samples)) mcdrop_samples[, j] <- pmax(lo, pmin(hi, mcdrop_samples[, j]))
        if (!is.null(quant_matrix))   quant_matrix[, j]   <- pmax(lo, pmin(hi, quant_matrix[, j]))
      }
    }
  } else {
    # Fallback: at minimum clip negatives
    if (!is.null(conf_samples))   conf_samples[conf_samples < 0]     <- 0
    if (!is.null(boot_samples))   boot_samples[boot_samples < 0]     <- 0
    if (!is.null(mcdrop_samples)) mcdrop_samples[mcdrop_samples < 0] <- 0
    if (!is.null(quant_matrix))   quant_matrix[quant_matrix < 0]     <- 0
  }

  if (verbose) cat("PipeMaster:: Done.\n")

  # Store prior samples from reftable
  prior_samples <- NULL
  if (!is.null(reftable) && !is.null(param_names)) {
    nuisance <- c("mean.rate", "sd.rate")
    pcols <- setdiff(param_names, nuisance)
    prior_samples <- as.matrix(reftable[, pcols, drop = FALSE])
    colnames(prior_samples) <- pcols
  }

  result <- list(
    point_estimate = point_est,
    conformal      = conf_samples,
    bootstrap      = boot_samples,
    mc_dropout     = mcdrop_samples,
    quantile       = quant_matrix,
    q_probs        = quant_probs,
    prior          = prior_samples,
    param_names    = param_names
  )
  class(result) <- "nn.posterior"
  result
}

#' @export
summary.nn.posterior <- function(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  param_names <- object$param_names
  n_params <- length(param_names)
  q_labels <- paste0(formatC(probs * 100, format = "fg"), "%")

  .summarize <- function(mat) {
    tbl <- matrix(NA, nrow = n_params, ncol = length(probs) + 2)
    colnames(tbl) <- c("Mean", "Median", q_labels)
    rownames(tbl) <- param_names
    for (j in seq_len(n_params)) {
      vals <- mat[, j]
      vals <- vals[is.finite(vals)]
      tbl[j, "Mean"]   <- mean(vals)
      tbl[j, "Median"] <- median(vals)
      tbl[j, -(1:2)]   <- quantile(vals, probs = probs)
    }
    tbl
  }

  out <- list(point_estimate = object$point_estimate)
  if (!is.null(object$conformal))
    out$conformal <- .summarize(object$conformal)
  if (!is.null(object$bootstrap))
    out$bootstrap <- .summarize(object$bootstrap)
  if (!is.null(object$mc_dropout))
    out$mc_dropout <- .summarize(object$mc_dropout)
  if (!is.null(object$quantile)) {
    # Quantile matrix is (n_quantiles x n_params) — report directly
    tbl <- t(object$quantile)
    colnames(tbl) <- paste0("Q", formatC(object$q_probs * 100, format = "fg"), "%")
    rownames(tbl) <- param_names
    out$quantile <- tbl
    out$q_probs  <- object$q_probs
  }
  class(out) <- "summary.nn.posterior"
  out
}

#' @export
print.summary.nn.posterior <- function(x, digits = 2, ...) {
  cat("Point estimate:\n")
  print(round(x$point_estimate, digits))
  if (!is.null(x$conformal)) {
    cat("\nConformal posterior:\n")
    print(round(x$conformal, digits))
  }
  if (!is.null(x$bootstrap)) {
    cat("\nBootstrap posterior:\n")
    print(round(x$bootstrap, digits))
  }
  if (!is.null(x$mc_dropout)) {
    cat("\nMC Dropout posterior:\n")
    print(round(x$mc_dropout, digits))
  }
  if (!is.null(x$quantile)) {
    cat("\nQuantile Regression:\n")
    print(round(x$quantile, digits))
  }
  invisible(x)
}

#' @export
print.nn.posterior <- function(x, ...) {
  methods <- c()
  if (!is.null(x$conformal))  methods <- c(methods, sprintf("conformal (%d samples)", nrow(x$conformal)))
  if (!is.null(x$bootstrap))  methods <- c(methods, sprintf("bootstrap (%d samples)", nrow(x$bootstrap)))
  if (!is.null(x$mc_dropout)) methods <- c(methods, sprintf("mc_dropout (%d samples)", nrow(x$mc_dropout)))
  if (!is.null(x$quantile))   methods <- c(methods, sprintf("quantile (%d quantiles)", nrow(x$quantile)))
  cat(sprintf("nn.posterior object — %s\n", paste(methods, collapse = " + ")))
  cat("Point estimate:\n")
  print(round(x$point_estimate, 2))
  cat("\nUse summary() for posterior quantiles, density() for density objects.\n")
  invisible(x)
}

#' @export
density.nn.posterior <- function(x, method = NULL, ...) {
  param_names <- x$param_names
  sample_methods <- c("conformal", "bootstrap", "mc_dropout")

  if (is.null(method)) {
    for (m in sample_methods) {
      if (!is.null(x[[m]])) { method <- m; break }
    }
    if (is.null(method)) stop("No posterior samples available.")
  }
  method <- match.arg(method, sample_methods)
  mat <- x[[method]]
  if (is.null(mat)) stop(sprintf("No %s samples available.", method))

  densities <- lapply(seq_along(param_names), function(j) {
    vals <- mat[, j]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 2) density(vals, ...) else NULL
  })
  names(densities) <- param_names
  densities
}

#' @export
plot.nn.posterior <- function(x, method = NULL, col = "red", lwd = 2,
                             show_point_est = TRUE, show_prior = FALSE, ...) {
  param_names <- x$param_names
  n_params <- length(param_names)

  all_methods <- c("prior", "conformal", "bootstrap", "mc_dropout", "quantile")

  # Pick methods
  if (is.null(method)) {
    methods <- c()
    for (m in all_methods) {
      if (!is.null(x[[m]])) methods <- c(methods, m)
    }
  } else {
    methods <- match.arg(method, all_methods, several.ok = TRUE)
  }
  if (show_prior && !"prior" %in% methods && !is.null(x$prior))
    methods <- c("prior", methods)

  post_methods <- setdiff(methods, "prior")
  if (length(post_methods) == 0 && !"prior" %in% methods)
    stop("No posterior samples available.")

  # Colors for methods
  method_cols <- c(prior = "grey50", conformal = "red", bootstrap = "blue",
                   mc_dropout = "darkgreen", quantile = "purple")
  if (length(col) == 1 && length(post_methods) == 1)
    method_cols[post_methods] <- col

  # Separate sample-based methods from quantile method
  sample_methods <- intersect(methods, c("prior", "conformal", "bootstrap", "mc_dropout"))
  has_quantile <- "quantile" %in% methods

  par(mfrow = c(1, n_params))
  for (j in seq_len(n_params)) {
    pname <- param_names[j]

    # Use prior range as valid bounds for density estimation and x-axis
    prior_range <- if (!is.null(x$prior)) range(x$prior[, j], na.rm = TRUE) else NULL

    # --- Pre-compute all density objects ---
    densities   <- list()
    density_lty <- list()

    for (m in sample_methods) {
      mat <- x[[m]]
      if (is.null(mat)) next
      vals <- mat[, j]
      vals <- vals[is.finite(vals)]
      if (length(vals) <= 2) next
      if (!is.null(prior_range)) {
        d <- density(vals, from = prior_range[1], to = prior_range[2])
      } else {
        d <- density(vals)
      }
      densities[[m]]   <- d
      density_lty[[m]] <- if (m == "prior") 3 else 1
    }

    if (has_quantile && !is.null(x$quantile)) {
      qvals  <- x$quantile[, j]
      qprobs <- x$q_probs
      pseudo <- approx(qprobs, qvals, xout = seq(qprobs[1], qprobs[length(qprobs)],
                                                   length.out = 10000),
                        rule = 2)$y
      if (!is.null(prior_range)) {
        d_q <- density(pseudo, n = 1024, adjust = 1.5,
                       from = prior_range[1], to = prior_range[2])
      } else {
        d_q <- density(pseudo, n = 1024, adjust = 1.5)
      }
      densities[["quantile"]]   <- d_q
      density_lty[["quantile"]] <- 1
    }

    # --- Derive common y-range across all densities ---
    combined_ylim <- c(0, max(vapply(densities, function(d) max(d$y), numeric(1))))

    # --- Plot all curves with shared axes ---
    first <- TRUE
    for (m in names(densities)) {
      d <- densities[[m]]
      if (first) {
        plot(d, main = pname, xlab = pname, ylab = "Density",
             col = method_cols[m], lwd = lwd, lty = density_lty[[m]],
             xlim = prior_range, ylim = combined_ylim, ...)
        first <- FALSE
      } else {
        lines(d, col = method_cols[m], lwd = lwd, lty = density_lty[[m]])
      }
    }

    if (show_point_est)
      abline(v = x$point_estimate[j], lty = 2, lwd = lwd)

    # Legend
    legend_labels <- c()
    legend_cols   <- c()
    legend_lty    <- c()
    for (m in sample_methods) {
      if (!is.null(densities[[m]])) {
        legend_labels <- c(legend_labels, m)
        legend_cols   <- c(legend_cols, method_cols[m])
        legend_lty    <- c(legend_lty, if (m == "prior") 3 else 1)
      }
    }
    if (has_quantile && !is.null(densities[["quantile"]])) {
      legend_labels <- c(legend_labels, "quantile")
      legend_cols   <- c(legend_cols, method_cols["quantile"])
      legend_lty    <- c(legend_lty, 1)
    }
    if (show_point_est) {
      legend_labels <- c(legend_labels, "Point est.")
      legend_cols   <- c(legend_cols, "black")
      legend_lty    <- c(legend_lty, 2)
    }

    legend("topright", legend = legend_labels, col = legend_cols,
           lty = legend_lty, lwd = lwd, cex = 0.8)
  }
}

# ============================================================================
# Internal: prepare observed data for prediction
# ============================================================================

.prep.observed <- function(observed, data, type, sfs.dims) {
  feat_mu <- data$feat_mu
  feat_sd <- data$feat_sd

  if (type == "sumstat") {
    obs_aug <- c(observed, log1p(abs(observed)))
    if (length(obs_aug) != length(feat_mu))
      stop(sprintf("Dimension mismatch: observed has %d features (augmented) but model expects %d.\n  This usually means the tune.nn() model was trained with extra columns as features.\n  Use exclude.cols in tune.nn() to exclude non-target parameter columns.",
                    length(obs_aug), length(feat_mu)))
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
    if (length(obs_aug) != length(feat_mu))
      stop(sprintf("Dimension mismatch: observed has %d features (augmented) but model expects %d.\n  This usually means the tune.nn() model was trained with extra columns as features.\n  Use exclude.cols in tune.nn() to exclude non-target parameter columns.",
                    length(obs_aug), length(feat_mu)))
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
# Internal: build quantile regression model (same architecture, different head)
# ============================================================================

.build.nn.qr <- function(hp, data, type, sfs.dims, n_quantiles, n_params, q_probs) {
  tf <- tensorflow::tf

  # Total outputs: one per (param, quantile) pair
  n_out <- as.integer(n_params * n_quantiles)

  # Build the backbone — reuse same architecture but with expanded output layer
  switch(type,
    sumstat = .build.resnet.qr(hp, data, n_out, n_params, n_quantiles, q_probs),
    sfs1d   = .build.cnn.qr(hp, data, n_out, n_params, n_quantiles, q_probs, "1d"),
    sfs2d   = .build.cnn2d.qr(hp, data, n_out, n_params, n_quantiles, q_probs, sfs.dims)
  )
}

.build.resnet.qr <- function(hp, data, n_out, n_params, n_quantiles, q_probs) {
  tf <- tensorflow::tf
  n_features <- ncol(data$X_train)
  l2 <- keras::regularizer_l2(hp$l2_reg)

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
    keras::layer_activation(x, activation = "relu")
  }

  inp <- keras::layer_input(shape = n_features)
  x <- inp |>
    keras::layer_dense(units = as.integer(hp$units_1), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  for (i in seq_len(hp$n_resblocks_1)) x <- res_block(x, hp$units_1)

  x <- x |>
    keras::layer_dense(units = as.integer(hp$units_2), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  if (hp$n_resblocks_2 > 0)
    for (i in seq_len(hp$n_resblocks_2)) x <- res_block(x, hp$units_2)

  x <- x |>
    keras::layer_dense(units = as.integer(hp$units_3), activation = "relu",
                       kernel_regularizer = l2) |>
    keras::layer_batch_normalization()
  if (isTRUE(hp$use_dropout) && hp$dropout > 0)
    x <- keras::layer_dropout(x, rate = hp$dropout)

  out <- x |> keras::layer_dense(units = n_out, activation = "linear")
  model <- keras::keras_model(inp, out)

  pinball <- .make.pinball.loss(n_params, n_quantiles, q_probs)
  model |> keras::compile(
    loss      = pinball,
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate)
  )
  model
}

.build.cnn.qr <- function(hp, data, n_out, n_params, n_quantiles, q_probs, dim_type) {
  tf <- tensorflow::tf
  n_bins    <- data$n_bins
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

  output <- x |> keras::layer_dense(units = n_out, activation = "linear")
  model <- keras::keras_model(input, output)

  pinball <- .make.pinball.loss(n_params, n_quantiles, q_probs)
  model |> keras::compile(
    loss      = pinball,
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate),
    metrics   = list("mae")
  )
  model
}

.build.cnn2d.qr <- function(hp, data, n_out, n_params, n_quantiles, q_probs, sfs.dims) {
  tf <- tensorflow::tf
  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]
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

  output <- x |> keras::layer_dense(units = n_out, activation = "linear")
  model <- keras::keras_model(input, output)

  pinball <- .make.pinball.loss(n_params, n_quantiles, q_probs)
  model |> keras::compile(
    loss      = pinball,
    optimizer = keras::optimizer_adam(learning_rate = hp$learning_rate),
    metrics   = list("mae")
  )
  model
}

# ============================================================================
# Internal: pinball (quantile) loss function
# ============================================================================

.make.pinball.loss <- function(n_params, n_quantiles, q_probs) {
  tf <- tensorflow::tf
  # q_probs is a numeric vector of length n_quantiles
  # output shape: (batch, n_params * n_quantiles)
  # y_true shape: (batch, n_params) — replicated targets

  q_tensor <- tf$constant(q_probs, dtype = "float32")
  np <- as.integer(n_params)
  nq <- as.integer(n_quantiles)

  function(y_true, y_pred) {
    # y_true: (batch, n_params) — expand to (batch, n_params, n_quantiles)
    y_true_exp <- tf$expand_dims(y_true, axis = -1L)
    y_true_rep <- tf$`repeat`(y_true_exp, repeats = nq, axis = -1L)

    # y_pred: (batch, n_params * n_quantiles) — reshape to (batch, n_params, n_quantiles)
    y_pred_3d <- tf$reshape(y_pred, c(-1L, np, nq))

    # Pinball loss
    err <- y_true_rep - y_pred_3d
    loss <- tf$maximum(q_tensor * err, (q_tensor - 1.0) * err)
    tf$reduce_mean(loss)
  }
}

# ============================================================================
# Internal: quantile regression — train and predict
# ============================================================================

.run.quantile.regression <- function(reftable, param.cols, observed, best_hp,
                                     type, sfs.dims, q_probs, max_epochs,
                                     seed, verbose, data, best_model = NULL,
                                     exclude.cols = NULL) {

  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params    <- length(target_cols)
  n_quantiles <- length(q_probs)

  warm_start <- !is.null(best_model)

  if (verbose)
    cat(sprintf("\nPipeMaster:: Quantile Regression (%d quantiles, warm_start=%s)...\n",
                n_quantiles, warm_start))

  # Prepare features and targets from reftable
  all_rows   <- seq_len(nrow(reftable))
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims,
                                      exclude.cols = exclude.cols)
  features_all <- full_split$features
  targets_all  <- full_split$targets
  n_total <- nrow(features_all)

  # Split into train/val
  set.seed(seed + 999L)
  idx <- sample(n_total)
  n_val   <- max(1L, floor(0.1 * n_total))
  n_train <- n_total - n_val
  tr_idx <- idx[1:n_train]
  va_idx <- idx[(n_train + 1):n_total]

  feat_tr <- features_all[tr_idx, ]
  feat_va <- features_all[va_idx, ]
  targ_tr <- targets_all[tr_idx, , drop = FALSE]
  targ_va <- targets_all[va_idx, , drop = FALSE]

  if (warm_start) {
    # Use same normalization as best_model (from tune.nn)
    f_mu <- data$feat_mu; f_sd <- data$feat_sd
    t_mu <- data$target_mu; t_sd <- data$target_sd
  } else {
    # Compute normalization from scratch
    f_mu <- colMeans(feat_tr)
    f_sd <- apply(feat_tr, 2, sd); f_sd[f_sd == 0] <- 1
    t_mu <- colMeans(log(targ_tr))
    t_sd <- apply(log(targ_tr), 2, sd); t_sd[t_sd == 0] <- 1
  }

  # Z-score features
  X_tr <- t((t(feat_tr) - f_mu) / f_sd)
  X_va <- t((t(feat_va) - f_mu) / f_sd)

  # Log + Z-score targets
  Y_tr <- t((t(log(targ_tr)) - t_mu) / t_sd)
  Y_va <- t((t(log(targ_va)) - t_mu) / t_sd)

  # Reshape for CNN architectures
  if (type == "sfs1d") {
    n_bins <- ncol(X_tr)
    dim(X_tr) <- c(nrow(X_tr), n_bins, 1L)
    dim(X_va) <- c(nrow(X_va), n_bins, 1L)
  } else if (type == "sfs2d") {
    X_tr <- array(X_tr, dim = c(nrow(X_tr), sfs.dims[1], sfs.dims[2], 1L))
    X_va <- array(X_va, dim = c(nrow(X_va), sfs.dims[1], sfs.dims[2], 1L))
  }

  qr_data <- list(
    X_train = X_tr, X_val = X_va,
    Y_train = Y_tr, Y_val = Y_va,
    n_features = ncol(features_all),
    n_bins = if (type != "sumstat") ncol(features_all) else NULL,
    feat_mu = f_mu, feat_sd = f_sd,
    target_mu = t_mu, target_sd = t_sd
  )

  # Build quantile regression model
  tensorflow::tf$random$set_seed(as.integer(seed + 777L))
  model <- .build.nn.qr(best_hp, qr_data, type, sfs.dims,
                         n_quantiles, n_params, q_probs)

  # --- Warm-start: copy backbone weights from best_model ---
  if (warm_start) {
    orig_w <- best_model$get_weights()
    qr_w   <- model$get_weights()

    # All weights match except the last 2 (output Dense kernel + bias)
    # which differ in shape: (units_3, n_params) vs (units_3, n_params*n_quantiles)
    n_backbone <- length(qr_w) - 2L

    if (length(orig_w) - 2L == n_backbone) {
      for (i in seq_len(n_backbone)) {
        qr_w[[i]] <- orig_w[[i]]
      }
      model$set_weights(qr_w)
      if (verbose)
        cat("PipeMaster:: Warm-start: copied backbone weights from best model\n")
    } else {
      warning("Warm-start: backbone weight count mismatch (",
              length(orig_w) - 2L, " vs ", n_backbone,
              "). Training from scratch.")
      warm_start <- FALSE
    }
  }

  pinball <- .make.pinball.loss(n_params, n_quantiles, q_probs)
  bs <- as.integer(best_hp$batch_size)

  if (warm_start) {
    # Phase 1: freeze backbone, train only output head at full LR
    n_layers <- length(model$layers)
    for (i in seq_len(n_layers - 1L)) {
      layer <- model$layers[[i]]
      layer$trainable <- FALSE
    }
    model |> keras::compile(
      loss      = pinball,
      optimizer = keras::optimizer_adam(learning_rate = best_hp$learning_rate)
    )
    if (verbose) cat("PipeMaster:: Phase 1: training output head (backbone frozen)...\n")
    model |> keras::fit(
      x = qr_data$X_train, y = qr_data$Y_train,
      validation_data = list(qr_data$X_val, qr_data$Y_val),
      epochs     = 100L,
      batch_size = bs,
      callbacks  = list(
        keras::callback_early_stopping(monitor = "val_loss", patience = 20L,
                                       restore_best_weights = TRUE)
      ),
      verbose = 0L
    )

    # Phase 2: unfreeze all, fine-tune at lower LR
    for (i in seq_len(n_layers - 1L)) {
      layer <- model$layers[[i]]
      layer$trainable <- TRUE
    }
    model |> keras::compile(
      loss      = pinball,
      optimizer = keras::optimizer_adam(learning_rate = best_hp$learning_rate / 10)
    )
    if (verbose) cat("PipeMaster:: Phase 2: fine-tuning full model...\n")
  } else {
    model |> keras::compile(
      loss      = pinball,
      optimizer = keras::optimizer_adam(learning_rate = best_hp$learning_rate)
    )
  }

  history <- model |> keras::fit(
    x = qr_data$X_train, y = qr_data$Y_train,
    validation_data = list(qr_data$X_val, qr_data$Y_val),
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

  vl <- history$metrics$val_loss
  if (is.null(vl)) vl <- history$history$val_loss
  n_ep <- length(unlist(vl))
  final_vl <- min(unlist(vl))

  if (verbose)
    cat(sprintf("PipeMaster:: QR model trained %d epochs (val_loss=%.6f)\n", n_ep, final_vl))

  # Predict on observed
  X_obs_qr <- .prep.observed.with(observed, f_mu, f_sd, type, sfs.dims)
  pred_z <- predict(model, X_obs_qr, verbose = 0L)   # (1, n_params * n_quantiles)
  pred_z <- as.numeric(pred_z)

  # Reshape to (n_params, n_quantiles) — TF reshape is row-major, so use byrow=TRUE
  pred_mat_z <- matrix(pred_z, nrow = n_params, ncol = n_quantiles, byrow = TRUE)

  # Inverse transform each value: exp(z * t_sd + t_mu)
  quant_mat <- matrix(NA_real_, nrow = n_quantiles, ncol = n_params)
  for (j in seq_len(n_params)) {
    quant_mat[, j] <- exp(pred_mat_z[j, ] * t_sd[j] + t_mu[j])
  }

  # Enforce monotonicity: sort quantiles per parameter
  for (j in seq_len(n_params)) {
    quant_mat[, j] <- sort(quant_mat[, j])
  }

  rm(model); gc()
  tryCatch(keras::k_clear_session(), error = function(e) NULL)

  if (verbose)
    cat("PipeMaster:: Quantile Regression done\n")

  quant_mat
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
# Internal: write standalone Hyperband worker Rscript
# ============================================================================

.write.hyperband.worker.script <- function(filepath) {
  writeLines(c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    '# Threading env (GPU env set externally by pool launcher)',
    'Sys.setenv(TF_NUM_INTRAOP_THREADS = "1",',
    '           TF_NUM_INTEROP_THREADS = "1",',
    '           OMP_NUM_THREADS = "1")',
    '',
    'load("shared_data.RData")',
    'load("round_tasks.RData")',
    '',
    'out_file <- file.path("results", sprintf("hb_%04d.csv", task_id))',
    'if (file.exists(out_file)) { cat("skip\\n"); q("no") }',
    '',
    'suppressPackageStartupMessages({',
    '  library(keras)',
    '  library(tensorflow)',
    '})',
    'source("_build_model.R")',
    '',
    '# Build data list for builder',
    'hb_data <- list(',
    '  X_train = X_train, X_val = X_val,',
    '  Y_train = Y_train, Y_val = Y_val,',
    '  n_features = n_features,',
    '  n_bins = n_bins',
    ')',
    '',
    'cfg <- round_configs[[task_id]]',
    'prev_ep <- round_prev_epochs[task_id]',
    'target_ep <- round_r_i',
    'worker_seed <- round_seeds[task_id]',
    '',
    'tf$random$set_seed(as.integer(worker_seed))',
    'model <- build_nn(cfg, hb_data, type, sfs.dims)',
    '',
    '# Load weights from previous round (if any)',
    'wpath <- file.path("weights", sprintf("cfg_%d", task_id))',
    'wfile <- file.path(wpath, "ckpt")',
    'if (prev_ep > 0 && dir.exists(wpath)) {',
    '  load_model_weights_tf(model, wfile)',
    '}',
    '',
    'history <- model |> fit(',
    '  x = hb_data$X_train, y = hb_data$Y_train,',
    '  validation_data = list(hb_data$X_val, hb_data$Y_val),',
    '  epochs        = as.integer(target_ep),',
    '  initial_epoch = as.integer(prev_ep),',
    '  batch_size    = as.integer(cfg$batch_size),',
    '  callbacks     = list(',
    '    callback_early_stopping(monitor = "val_loss",',
    '                            patience = 10L,',
    '                            restore_best_weights = TRUE),',
    '    callback_reduce_lr_on_plateau(monitor = "val_loss",',
    '                                  patience = 5L,',
    '                                  factor = 0.5,',
    '                                  min_lr = 1e-6,',
    '                                  verbose = 0L)',
    '  ),',
    '  verbose = 0L',
    ')',
    '',
    'vl <- history$metrics$val_loss',
    'if (is.null(vl)) vl <- history$history$val_loss',
    'best_vl <- min(unlist(vl))',
    '',
    '# Save weights',
    'dir.create(wpath, recursive = TRUE, showWarnings = FALSE)',
    'save_model_weights_tf(model, wfile)',
    '',
    '# Write result',
    'write.csv(data.frame(task_id = task_id, val_loss = best_vl),',
    '          out_file, row.names = FALSE)',
    'cat(sprintf("  hb %d done (val_loss=%.4f)\\n", task_id, best_vl))',
    'k_clear_session()'
  ), filepath)
}

# ============================================================================
# Internal: write standalone conformal worker Rscript
# ============================================================================

.write.conformal.worker.script <- function(filepath) {
  lines <- c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    '# Threading env (GPU env set externally by pool launcher)',
    'Sys.setenv(TF_NUM_INTRAOP_THREADS = "1",',
    '           TF_NUM_INTEROP_THREADS = "1",',
    '           OMP_NUM_THREADS = "1")',
    '',
    'load("shared_data.RData")',
    '',
    'out_file <- file.path("results", sprintf("conf_%04d.csv", task_id))',
    'if (file.exists(out_file)) { cat("skip\\n"); q("no") }',
    '',
    'suppressPackageStartupMessages({',
    '  library(keras)',
    '  library(tensorflow)',
    '})',
    'source("_build_model.R")',
    '',
    '# Split into train / calibration',
    'n_cal   <- floor(cal_frac * n_rows)',
    'n_train <- n_rows - n_cal',
    '',
    'set.seed(seed + task_id * 100)',
    'idx     <- sample(n_rows)',
    'tr_idx  <- idx[1:n_train]',
    'cal_idx <- idx[(n_train + 1):n_rows]',
    '',
    'feat_tr    <- features_all[tr_idx, ]',
    'targ_tr    <- targets_all[tr_idx, , drop = FALSE]',
    'feat_cal   <- features_all[cal_idx, ]',
    'params_cal <- targets_all[cal_idx, , drop = FALSE]',
    '',
    '# Z-score features using train stats',
    'f_mu  <- colMeans(feat_tr)',
    'f_sd  <- apply(feat_tr, 2, sd); f_sd[f_sd == 0] <- 1',
    'X_tr  <- t((t(feat_tr)  - f_mu) / f_sd)',
    'X_cal <- t((t(feat_cal) - f_mu) / f_sd)',
    '',
    '# Log + Z-score targets using train stats',
    'Y_tr_log <- log(targ_tr)',
    't_mu <- colMeans(Y_tr_log)',
    't_sd <- apply(Y_tr_log, 2, sd); t_sd[t_sd == 0] <- 1',
    'Y_tr <- t((t(Y_tr_log) - t_mu) / t_sd)',
    '',
    '# Reshape for CNN if needed',
    'if (type == "sfs1d") {',
    '  n_bins <- ncol(X_tr)',
    '  dim(X_tr)  <- c(nrow(X_tr), n_bins, 1L)',
    '  dim(X_cal) <- c(nrow(X_cal), n_bins, 1L)',
    '} else if (type == "sfs2d") {',
    '  X_tr  <- array(X_tr,  dim = c(nrow(X_tr),  sfs.dims[1], sfs.dims[2], 1L))',
    '  X_cal <- array(X_cal, dim = c(nrow(X_cal), sfs.dims[1], sfs.dims[2], 1L))',
    '}',
    '',
    '# Split train into train/val for early stopping',
    'n_tr_rows <- nrow(X_tr)',
    'n_va <- max(1L, floor(0.1 * n_tr_rows))',
    'va_rows <- seq_len(n_va)',
    'tr_rows <- seq(n_va + 1L, n_tr_rows)',
    '',
    'if (type == "sumstat") {',
    '  X_va_s <- X_tr[va_rows, , drop = FALSE]',
    '  X_tr_s <- X_tr[tr_rows, , drop = FALSE]',
    '} else if (type == "sfs1d") {',
    '  X_va_s <- X_tr[va_rows, , , drop = FALSE]',
    '  X_tr_s <- X_tr[tr_rows, , , drop = FALSE]',
    '} else {',
    '  X_va_s <- X_tr[va_rows, , , , drop = FALSE]',
    '  X_tr_s <- X_tr[tr_rows, , , , drop = FALSE]',
    '}',
    'Y_va_s <- Y_tr[va_rows, , drop = FALSE]',
    'Y_tr_s <- Y_tr[tr_rows, , drop = FALSE]',
    '',
    'ens_data <- list(',
    '  X_train = X_tr_s, X_val = X_va_s,',
    '  Y_train = Y_tr_s, Y_val = Y_va_s,',
    '  n_features = ncol(feat_tr),',
    '  n_bins = if (type != "sumstat") ncol(feat_tr) else NULL,',
    '  feat_mu = f_mu, feat_sd = f_sd,',
    '  target_mu = t_mu, target_sd = t_sd',
    ')',
    '',
    'tf$random$set_seed(as.integer(seed + task_id * 1000))',
    'model <- build_nn(best_hp, ens_data, type, sfs.dims)',
    '',
    'bs <- as.integer(best_hp$batch_size)',
    'model |> fit(',
    '  x = ens_data$X_train, y = ens_data$Y_train,',
    '  validation_data = list(ens_data$X_val, ens_data$Y_val),',
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
    '# Predict on calibration set -> inverse transform',
    'cal_pred_z <- predict(model, X_cal, verbose = 0L)',
    'cal_pred   <- exp(t(t(cal_pred_z) * t_sd + t_mu))',
    '',
    '# Residuals: true params - predicted params',
    'cal_resid <- params_cal - cal_pred',
    '',
    '# Center on best model point estimate (avoids multi-peaked posteriors)',
    'if (!is.null(point_est)) {',
    '  center <- point_est',
    '} else {',
    '  if (type == "sumstat") {',
    '    obs_aug <- c(observed, log1p(abs(observed)))',
    '    X_obs <- matrix((obs_aug - f_mu) / f_sd, nrow = 1)',
    '    X_obs[!is.finite(X_obs)] <- 0',
    '  } else if (type == "sfs1d") {',
    '    obs_z <- (log1p(observed) - f_mu) / f_sd',
    '    obs_z[!is.finite(obs_z)] <- 0',
    '    X_obs <- array(obs_z, dim = c(1L, length(obs_z), 1L))',
    '  } else {',
    '    obs_z <- (log1p(observed) - f_mu) / f_sd',
    '    obs_z[!is.finite(obs_z)] <- 0',
    '    X_obs <- array(obs_z, dim = c(1L, sfs.dims[1], sfs.dims[2], 1L))',
    '  }',
    '  obs_pred_z <- predict(model, X_obs, verbose = 0L)',
    '  center <- as.numeric(exp(t(t(obs_pred_z) * t_sd + t_mu)))',
    '}',
    '',
    '# Conformal samples: center[j] + cal_resid[, j]',
    'conf_samples <- matrix(NA, nrow = nrow(cal_resid), ncol = n_params)',
    'for (j in seq_len(n_params))',
    '  conf_samples[, j] <- center[j] + cal_resid[, j]',
    '',
    'df <- as.data.frame(conf_samples)',
    'colnames(df) <- target_cols',
    'write.csv(df, out_file, row.names = FALSE)',
    'cat(sprintf("  conf %d done (%d samples)\\n", task_id, nrow(conf_samples)))',
    'k_clear_session()'
  )

  writeLines(lines, filepath)
}

# ============================================================================
# Internal: prepare reftable split for conformal/bootstrap training
# ============================================================================

.prep.reftable.split <- function(reftable, param.cols, row_idx, type, sfs.dims,
                                  exclude.cols = NULL) {
  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)

  targets <- as.matrix(reftable[row_idx, target_cols, drop = FALSE])

  if (type == "sumstat") {
    stat_cols <- setdiff(colnames(reftable), c(param.cols, nuisance, exclude.cols))
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
# Internal: conformal prediction (sequential ensemble, cores <= 1 fallback)
# ============================================================================

.run.conformal.sequential <- function(reftable, param.cols, observed, best_hp, type, sfs.dims,
                           n_ensemble, cal.frac, max_epochs, seed, verbose,
                           point_est = NULL, exclude.cols = NULL) {

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
    tr_split  <- .prep.reftable.split(reftable, param.cols, tr_idx, type, sfs.dims,
                                       exclude.cols = exclude.cols)
    cal_split <- .prep.reftable.split(reftable, param.cols, cal_idx, type, sfs.dims,
                                       exclude.cols = exclude.cols)

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

    # Center conformal samples on best model's point estimate (if provided)
    # to avoid multi-peaked posteriors from varying per-ensemble obs_est
    if (!is.null(point_est)) {
      center <- point_est
    } else {
      X_obs <- .prep.observed.with(observed, f_mu, f_sd, type, sfs.dims)
      obs_pred_z <- predict(model, X_obs, verbose = 0L)
      center <- as.numeric(.inv.transform(obs_pred_z, t_mu, t_sd))
    }

    # Conformal samples: center[j] + cal_resid[, j]
    ens_samples <- matrix(NA, nrow = nrow(cal_resid), ncol = n_params)
    for (j in seq_len(n_params))
      ens_samples[, j] <- center[j] + cal_resid[, j]

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
                           n_boot, max_epochs, cores, gpus = 0, seed, verbose,
                           exclude.cols = NULL) {

  if (cores <= 1) {
    return(.run.bootstrap.sequential(
      reftable, param.cols, observed, best_hp, type, sfs.dims,
      n_boot, max_epochs, seed, verbose,
      exclude.cols = exclude.cols
    ))
  }

  # Parallel bootstrap via Rscript workers
  .run.bootstrap.parallel(
    reftable, param.cols, observed, best_hp, type, sfs.dims,
    n_boot, max_epochs, cores, gpus = gpus, seed, verbose,
    exclude.cols = exclude.cols
  )
}

# ============================================================================
# Internal: sequential bootstrap (cores = 1)
# ============================================================================

.run.bootstrap.sequential <- function(reftable, param.cols, observed, best_hp,
                                      type, sfs.dims, n_boot, max_epochs,
                                      seed, verbose, exclude.cols = NULL) {

  nuisance <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params <- length(target_cols)

  if (verbose)
    cat(sprintf("\nPipeMaster:: Bootstrap (%d replicates, sequential)...\n", n_boot))

  boot_matrix <- matrix(NA, nrow = n_boot, ncol = n_params)

  # Prepare full data once
  all_rows <- seq_len(nrow(reftable))
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims,
                                      exclude.cols = exclude.cols)
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
# Internal: write standalone bootstrap worker Rscript
# ============================================================================

.write.bootstrap.worker.script <- function(filepath) {
  writeLines(c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    '# Threading env (GPU env set externally by pool launcher)',
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
  ), filepath)
}

# ============================================================================
# Internal: continuous worker pool with priority scheduling
# ============================================================================

.launch.rscript.pool <- function(tasks, cores, work_dir,
                                  timeout_per_task, gpus, verbose,
                                  max_retries = 2L) {
  old_wd <- getwd()
  setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  logs_dir <- file.path(work_dir, "logs")
  dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)
  sentinels_dir <- file.path(work_dir, "sentinels")
  dir.create(sentinels_dir, showWarnings = FALSE, recursive = TRUE)

  n_tasks <- length(tasks)
  pending   <- seq_len(n_tasks)
  active    <- list()
  completed <- integer(0)
  failed    <- integer(0)
  gpu_counter <- 0L
  retry_count <- integer(n_tasks)  # starts at 0 for each task

  launch_one <- function(task_idx) {
    task <- tasks[[task_idx]]
    attempt <- retry_count[task_idx] + 1L
    log_file <- file.path("logs", sprintf("%s_%04d_a%d.log",
                                          task$prefix, task$id, attempt))

    # Sentinel file paths (include attempt to avoid stale files)
    done_file <- file.path(sentinels_dir,
                           sprintf("%s_%04d_a%d.done", task$prefix, task$id, attempt))
    pid_file  <- file.path(sentinels_dir,
                           sprintf("%s_%04d_a%d.pid", task$prefix, task$id, attempt))

    # GPU environment
    if (gpus > 0) {
      gpu_id <- gpu_counter %% gpus
      gpu_counter <<- gpu_counter + 1L
      gpu_env <- sprintf("CUDA_VISIBLE_DEVICES=%d TF_FORCE_GPU_ALLOW_GROWTH=true",
                         gpu_id)
    } else {
      gpu_env <- "CUDA_VISIBLE_DEVICES=-1"
    }

    # Shell wrapper: run Rscript, capture exit code in .done file, record PID
    cmd <- sprintf(
      "{ env %s Rscript %s %d > %s 2>&1; echo $? > %s; } & echo $! > %s",
      gpu_env, shQuote(task$script), task$id,
      shQuote(log_file), shQuote(done_file), shQuote(pid_file))
    system(cmd, wait = FALSE)

    list(task_idx = task_idx, start_time = Sys.time(), log_file = log_file,
         done_file = done_file, pid_file = pid_file, attempt = attempt)
  }

  # Helper: read last N lines from a log file for diagnostics
  .log_tail <- function(log_rel_path, n = 20L) {
    log_path <- file.path(work_dir, log_rel_path)
    if (!file.exists(log_path)) return(NULL)
    log_lines <- readLines(log_path, warn = FALSE)
    n_lines <- length(log_lines)
    if (n_lines == 0L) return(NULL)
    log_lines[max(1L, n_lines - n + 1L):n_lines]
  }

  # Fill initial slots (conformal tasks come first = priority)
  while (length(active) < cores && length(pending) > 0) {
    tidx <- pending[1]
    pending <- pending[-1]
    active[[length(active) + 1]] <- launch_one(tidx)
  }

  if (verbose)
    cat(sprintf("  [pool] Launched %d / %d tasks (%d pending)\n",
                length(active), n_tasks, length(pending)))

  # Poll loop — check every second, 4-way detection
  while (length(active) > 0) {
    Sys.sleep(1)

    still_active <- list()
    for (a in active) {
      task <- tasks[[a$task_idx]]
      elapsed_sec <- as.numeric(difftime(Sys.time(), a$start_time, units = "secs"))

      if (file.exists(task$result)) {
        # --- SUCCESS: result CSV exists ---
        completed <- c(completed, a$task_idx)
        if (verbose)
          cat(sprintf("  [pool] %s %d done (%d/%d, %.0fs)\n",
                      task$prefix, task$id,
                      length(completed), n_tasks, elapsed_sec))

      } else if (file.exists(a$done_file)) {
        # --- CRASH: process exited but no result CSV ---
        exit_code <- tryCatch(
          as.integer(trimws(readLines(a$done_file, n = 1L, warn = FALSE))),
          error = function(e) NA_integer_)
        tail_lines <- .log_tail(a$log_file)

        if (retry_count[a$task_idx] < max_retries) {
          retry_count[a$task_idx] <- retry_count[a$task_idx] + 1L
          if (verbose)
            cat(sprintf("  [pool] %s %d CRASHED (exit %s, %.0fs) — retry %d/%d\n",
                        task$prefix, task$id,
                        if (is.na(exit_code)) "?" else as.character(exit_code),
                        elapsed_sec,
                        retry_count[a$task_idx], max_retries))
          # Re-queue at front of pending so it gets the next free slot
          pending <- c(a$task_idx, pending)
        } else {
          failed <- c(failed, a$task_idx)
          tail_msg <- if (!is.null(tail_lines))
            paste(tail_lines, collapse = "\n") else "(no log)"
          warning(sprintf(
            "Worker %s %d crashed (exit %s) after %d attempt(s), %.0fs. Last log:\n%s",
            task$prefix, task$id,
            if (is.na(exit_code)) "?" else as.character(exit_code),
            max_retries + 1L, elapsed_sec, tail_msg),
            call. = FALSE)
        }

      } else if (elapsed_sec > timeout_per_task) {
        # --- TIMEOUT: kill process and fail/retry ---
        # Try to kill the child process
        if (file.exists(a$pid_file)) {
          pid <- tryCatch(
            as.integer(trimws(readLines(a$pid_file, n = 1L, warn = FALSE))),
            error = function(e) NA_integer_)
          if (!is.na(pid)) tryCatch(tools::pskill(pid), error = function(e) NULL)
        }

        if (retry_count[a$task_idx] < max_retries) {
          retry_count[a$task_idx] <- retry_count[a$task_idx] + 1L
          if (verbose)
            cat(sprintf("  [pool] %s %d TIMEOUT (%.0fs) — retry %d/%d\n",
                        task$prefix, task$id, elapsed_sec,
                        retry_count[a$task_idx], max_retries))
          pending <- c(a$task_idx, pending)
        } else {
          failed <- c(failed, a$task_idx)
          tail_lines <- .log_tail(a$log_file)
          tail_msg <- if (!is.null(tail_lines))
            paste(tail_lines, collapse = "\n") else "(no log)"
          warning(sprintf(
            "Worker %s %d timed out after %d attempt(s), %.0fs. Last log:\n%s",
            task$prefix, task$id, max_retries + 1L, elapsed_sec, tail_msg),
            call. = FALSE)
        }

      } else {
        # --- RUNNING: still within timeout ---
        still_active[[length(still_active) + 1]] <- a
      }
    }
    active <- still_active

    # Backfill freed slots from pending queue
    while (length(active) < cores && length(pending) > 0) {
      tidx <- pending[1]
      pending <- pending[-1]
      active[[length(active) + 1]] <- launch_one(tidx)
    }
  }

  list(completed = completed, failed = failed)
}

# ============================================================================
# Internal: unified parallel pool for conformal + bootstrap
# ============================================================================

.run.parallel.pool <- function(reftable, param.cols, observed, best_hp,
                                type, sfs.dims, do_conformal, do_bootstrap,
                                n_ensemble, n_boot, cal.frac, max_epochs,
                                cores, gpus, seed, verbose,
                                point_est = NULL, exclude.cols = NULL) {

  nuisance    <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params    <- length(target_cols)

  # Create temp working directory
  work_dir    <- tempfile("nn_pool_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)

  # Prepare and save shared data
  all_rows   <- seq_len(nrow(reftable))
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims,
                                      exclude.cols = exclude.cols)
  features_all <- full_split$features
  targets_all  <- full_split$targets
  n_rows     <- nrow(features_all)
  cal_frac   <- cal.frac

  shared_file <- file.path(work_dir, "shared_data.RData")
  save(features_all, targets_all, observed, best_hp,
       type, sfs.dims, n_rows, n_params, max_epochs, seed, target_cols,
       cal_frac, point_est,
       file = shared_file)

  # Write model builder script
  .write.builder.script(file.path(work_dir, "_build_model.R"), type)

  # Build task list: conformal first (priority), then bootstrap
  tasks <- list()

  if (do_conformal) {
    .write.conformal.worker.script(file.path(work_dir, "_conf_worker.R"))
    for (i in seq_len(n_ensemble)) {
      tasks[[length(tasks) + 1]] <- list(
        script = "_conf_worker.R", id = i,
        result = sprintf("results/conf_%04d.csv", i),
        prefix = "conf"
      )
    }
  }

  if (do_bootstrap) {
    .write.bootstrap.worker.script(file.path(work_dir, "_boot_worker.R"))
    for (b in seq_len(n_boot)) {
      tasks[[length(tasks) + 1]] <- list(
        script = "_boot_worker.R", id = b,
        result = sprintf("results/boot_%04d.csv", b),
        prefix = "boot"
      )
    }
  }

  n_conf <- if (do_conformal) n_ensemble else 0L
  n_boot_tasks <- if (do_bootstrap) n_boot else 0L

  if (verbose)
    cat(sprintf("\nPipeMaster:: Unified pool: %d conformal + %d bootstrap tasks, %d cores%s\n",
                n_conf, n_boot_tasks, cores,
                if (gpus > 0) sprintf(", %d GPUs", gpus) else ""))

  # Launch pool
  pool_result <- .launch.rscript.pool(tasks, cores, work_dir,
                                       timeout_per_task = max_epochs * 10,
                                       gpus = gpus, verbose = verbose)

  # Collect conformal results
  result <- list(conformal = NULL, bootstrap = NULL)

  if (do_conformal) {
    conf_list <- list()
    for (i in seq_len(n_ensemble)) {
      csv_file <- file.path(results_dir, sprintf("conf_%04d.csv", i))
      if (file.exists(csv_file)) {
        conf_list[[length(conf_list) + 1]] <- as.matrix(read.csv(csv_file))
      }
    }
    if (length(conf_list) > 0)
      result$conformal <- do.call(rbind, conf_list)

    if (verbose)
      cat(sprintf("PipeMaster:: Conformal: %d samples from %d/%d ensemble models\n",
                  if (!is.null(result$conformal)) nrow(result$conformal) else 0L,
                  length(conf_list), n_ensemble))
  }

  # Collect bootstrap results
  if (do_bootstrap) {
    boot_matrix <- matrix(NA, nrow = n_boot, ncol = n_params)
    for (b in seq_len(n_boot)) {
      csv_file <- file.path(results_dir, sprintf("boot_%04d.csv", b))
      if (file.exists(csv_file)) {
        row <- read.csv(csv_file)
        boot_matrix[b, ] <- as.numeric(row[1, ])
      }
    }
    result$bootstrap <- boot_matrix

    n_fail <- sum(is.na(boot_matrix[, 1]))
    if (verbose)
      cat(sprintf("PipeMaster:: Bootstrap: %d/%d successful\n",
                  n_boot - n_fail, n_boot))
  }

  # Clean up temp files
  unlink(work_dir, recursive = TRUE)

  result
}

# ============================================================================
# Internal: parallel bootstrap via Rscript workers
# ============================================================================

.run.bootstrap.parallel <- function(reftable, param.cols, observed, best_hp,
                                    type, sfs.dims, n_boot, max_epochs,
                                    cores, gpus = 0, seed, verbose,
                                    exclude.cols = NULL) {

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
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims,
                                      exclude.cols = exclude.cols)
  features_all <- full_split$features
  targets_all  <- full_split$targets
  n_rows <- nrow(features_all)

  shared_file <- file.path(work_dir, "shared_data.RData")
  save(features_all, targets_all, observed, best_hp,
       type, sfs.dims, n_rows, n_params, max_epochs, seed, target_cols,
       file = shared_file)

  # Write model builder and worker scripts
  .write.builder.script(file.path(work_dir, "_build_model.R"), type)
  .write.bootstrap.worker.script(file.path(work_dir, "_boot_worker.R"))

  # Build task list
  tasks <- list()
  for (b in seq_len(n_boot)) {
    tasks[[b]] <- list(
      script = "_boot_worker.R", id = b,
      result = sprintf("results/boot_%04d.csv", b),
      prefix = "boot"
    )
  }

  # Launch via continuous pool
  pool_result <- .launch.rscript.pool(tasks, cores, work_dir,
                                       timeout_per_task = max_epochs * 10,
                                       gpus = gpus, verbose = verbose)

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
