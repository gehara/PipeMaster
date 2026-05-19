#' Hyperband Tuner for Neural Network Architectures
#'
#' Pure R implementation of the Hyperband algorithm for tuning neural
#' networks used in demographic parameter estimation.
#' Supports ResNet (summary statistics), 1D CNN (single-pop SFS), and
#' 2D CNN (joint SFS) architectures. Two backends available:
#' \code{"torch"} (pure R, no Python dependency) and \code{"keras"}
#' (TensorFlow via reticulate, faster XLA-compiled training).
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
#' @param search_space named list — overrides default HP ranges. NULL uses
#'   architecture-specific defaults.
#' @param exclude.cols character vector — additional column names to exclude from
#'   features (e.g., other parameter columns not in \code{param.cols}). Required
#'   when estimating a single parameter from a reftable that contains multiple
#'   parameter columns, to prevent other parameters from leaking into the feature
#'   set. Default NULL (only \code{param.cols} and nuisance columns are excluded).
#' @param val.frac numeric — validation fraction (default 0.1).
#' @param n_searches integer — number of independent Hyperband searches to run
#'   (default 1). Each search explores a different population of random HP
#'   configurations (unique seed). Higher values increase the chance of finding a
#'   good architecture. With \code{cores > 1}, searches run in parallel via
#'   separate Rscript workers. With \code{cores = 1}, searches run sequentially
#'   in the main process.
#' @param cores integer — maximum number of concurrent search workers
#'   (default 1, sequential). Each worker spawns a separate R process with its
#'   own TensorFlow runtime and a full copy of the training data, so RAM usage
#'   scales linearly (~1.5-2 GB per worker). Ignored when \code{n_searches = 1}.
#' @param gpus integer — number of GPUs to distribute searches across
#'   (default 0, CPU-only). Searches are assigned to GPUs round-robin, up to
#'   \code{gpu.threshold} per GPU. Remaining searches run CPU-only. When
#'   \code{n_searches = 1}, the single search uses all available GPUs directly.
#' @param gpu.threshold integer — maximum searches per GPU (default 4).
#'   Total GPU searches = \code{min(n_searches, gpu.threshold * gpus)}. Excess
#'   searches run on CPU. Example: \code{n_searches=10, gpus=1, gpu.threshold=4}
#'   assigns 4 searches to the GPU and 6 to CPU. Too many searches per GPU may
#'   cause GPU out-of-memory errors. Ignored when \code{gpus = 0}.
#' @param top_k integer — number of top models to keep from parallel searches
#'   (default 1). When \code{top_k > 1} and \code{n_searches > 1}, the best K
#'   models (by validation loss) are retained. At prediction time,
#'   \code{nn.predict()} uses a proximity-weighted ensemble of these models for
#'   more stable point estimates. Capped at \code{n_searches}.
#' @param backend character — \code{"torch"} (default), \code{"keras"}, or
#'   \code{"auto"}. Torch is a pure R backend (no Python dependency).
#'   Keras requires TensorFlow via reticulate. \code{"auto"} uses torch
#'   if available, else keras.
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list with:
#' \describe{
#'   \item{best_hp}{named list of best hyperparameters}
#'   \item{best_val_loss}{best validation loss achieved}
#'   \item{all_results}{data.frame of all evaluated configs (hp_string, val_loss, bracket, round)}
#'   \item{best_model}{trained model (keras or torch, depending on backend)}
#'   \item{models}{list of top-K models, sorted by val_loss (length 1 when \code{top_k = 1})}
#'   \item{models_val_loss}{numeric vector of validation losses for each model in \code{models}}
#'   \item{models_metrics}{list of per-model metrics (each element: \code{list(r2, mpe)} with named vectors in original parameter scale)}
#'   \item{backend}{character — \code{"torch"} or \code{"keras"}}
#' }
#'
#' @export
tune.nn <- function(reftable, param.cols,
                    type = c("sumstat", "sfs1d", "sfs2d", "emulator"),
                    sfs.dims = NULL,
                    max_epochs = 500, eta = 3,
                    search_space = NULL,
                    exclude.cols = NULL,
                    val.frac = 0.1,
                    n_searches = 1L, cores = 1L, gpus = 0L,
                    gpu.threshold = 4L, greedy = TRUE,
                    top_k = 1L,
                    backend = c("torch", "keras", "auto"),
                    seed = 42, verbose = TRUE) {

  # --- Resolve backend ---
  backend <- match.arg(backend)
  if (backend == "auto") {
    backend <- if (requireNamespace("torch", quietly = TRUE) &&
                   torch::torch_is_installed()) "torch" else "keras"
  }
  if (verbose) cat(sprintf("PipeMaster:: tune.nn backend: %s\n", backend))

  # --- Dispatch to torch backend ---
  if (backend == "torch") {
    res <- .tune.nn.torch(reftable = reftable, param.cols = param.cols,
                          type = type, sfs.dims = sfs.dims,
                          max_epochs = max_epochs, eta = eta,
                          search_space = search_space,
                          exclude.cols = exclude.cols,
                          val.frac = val.frac,
                          n_searches = n_searches, cores = cores,
                          gpus = gpus, gpu.threshold = gpu.threshold,
                          greedy = greedy,
                          top_k = top_k,
                          seed = seed, verbose = verbose)
    return(res)
  }

  # --- Keras dependency check ---
  if (!requireNamespace("keras", quietly = TRUE) ||
      !requireNamespace("tensorflow", quietly = TRUE))
    stop("tune.nn(backend='keras') requires the 'keras' and 'tensorflow' R packages.\n",
         "Install with: install.packages(c('keras', 'tensorflow'))\n",
         "Then run: keras::install_keras()")

  # --- Validate top_k ---
  top_k <- as.integer(top_k)
  if (top_k < 1L) stop("top_k must be >= 1")
  if (n_searches > 1L && top_k > n_searches)
    top_k <- n_searches

  n_concurrent <- min(as.integer(cores), as.integer(n_searches))

  # --- Memory guard ---
  if (n_concurrent > 1L) {
    avail_gb <- tryCatch({
      mem_info <- system("free -b 2>/dev/null", intern = TRUE)
      if (length(mem_info) >= 2) {
        fields <- as.numeric(strsplit(trimws(mem_info[2]), "\\s+")[[1]])
        fields[7] / 1e9
      } else NA_real_
    }, error = function(e) NA_real_)

    if (!is.na(avail_gb)) {
      est_per_worker <- 1.5
      est_total <- n_concurrent * est_per_worker
      if (est_total > avail_gb * 0.85) {
        warning(sprintf(
          paste0("cores=%d concurrent searches may exceed available RAM ",
                 "(%.1f GB free, ~%.0f GB estimated). ",
                 "Each search loads TensorFlow + a full data copy (~1.5 GB). ",
                 "Reduce cores or n_searches if you experience memory issues."),
          n_concurrent, avail_gb, est_total), call. = FALSE)
      }
    }
  }

  # --- GPU OOM warning ---
  if (gpus > 0L && n_searches > 1L) {
    n_gpu_searches <- min(n_searches, gpu.threshold * gpus)
    per_gpu <- ceiling(n_gpu_searches / gpus)
    if (per_gpu > 2L) {
      warning(sprintf(
        paste0("gpu.threshold=%d with %d GPU(s) puts up to %d searches per GPU. ",
               "This may cause GPU out-of-memory. Reduce gpu.threshold ",
               "if you experience OOM errors."),
        gpu.threshold, gpus, per_gpu), call. = FALSE)
    }
  }

  # --- GPU memory growth (single-search mode only) ---
  if (n_searches <= 1L) {
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
                           switch(type, sumstat = "ResNet", sfs1d = "1D CNN",
                                  sfs2d = "2D CNN", emulator = "Emulator ResNet")))

  data <- .prep.data(reftable, param.cols, type, sfs.dims, exclude.cols, val.frac, seed)

  n_feat <- if (type %in% c("sumstat", "emulator")) ncol(data$X_train) else data$n_features
  n_targ <- ncol(data$Y_train)

  if (verbose) cat(sprintf("PipeMaster:: %d features, %d targets | %d train, %d val\n",
                           n_feat, n_targ, nrow(data$X_train), nrow(data$X_val)))

  # --- Run Hyperband ---
  if (n_searches > 1L) {
    result <- .parallel.hyperband.search(
      data = data, search_space = ss, type = type, sfs.dims = sfs.dims,
      max_epochs = max_epochs, eta = eta,
      n_searches = n_searches, cores = cores, gpus = gpus,
      gpu.threshold = gpu.threshold, greedy = greedy, top_k = top_k,
      seed = seed, verbose = verbose)

    return(list(
      best_hp         = result$best_hp,
      best_val_loss   = result$best_val_loss,
      all_results     = result$all_results,
      best_model      = result$best_model,
      models          = result$models,
      models_val_loss = result$models_val_loss,
      models_metrics  = result$models_metrics,
      data            = data,
      type            = type,
      sfs.dims        = sfs.dims,
      exclude.cols    = exclude.cols,
      backend         = "keras"
    ))
  }

  # --- Single search (sequential path) ---
  hb <- .hyperband(
    search_space = ss,
    data         = data,
    type         = type,
    sfs.dims     = sfs.dims,
    max_epochs   = max_epochs,
    eta          = eta,
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

  # Compute R² and MPE
  models_metrics <- list(.compute.model.metrics(best_model, data, type))
  if (verbose) .print.model.metrics(models_metrics, final_val_loss, verbose)

  list(
    best_hp         = hb$best_hp,
    best_val_loss   = final_val_loss,
    all_results     = hb$all_results,
    best_model      = best_model,
    models          = list(best_model),
    models_val_loss = final_val_loss,
    models_metrics  = models_metrics,
    data            = data,
    type            = type,
    sfs.dims        = sfs.dims,
    exclude.cols    = exclude.cols,
    backend         = "keras"
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
    ),
    emulator = .emulator.search.space()
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
  if (type %in% c("sumstat", "emulator")) {
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
    sumstat  = .build.resnet(hp, data),
    sfs1d    = .build.cnn1d(hp, data),
    sfs2d    = .build.cnn2d(hp, data, sfs.dims),
    emulator = .build.resnet(hp, data)
  )
}

# --- ResNet for summary statistics ---
.build.resnet <- function(hp, data) {
  n_features <- ncol(data$X_train)
  n_targets  <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  .drop <- function(x, rate) keras::layer_dropout(x, rate = rate)

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
.build.cnn1d <- function(hp, data) {
  n_bins    <- data$n_bins
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  .drop <- function(x, rate) keras::layer_dropout(x, rate = rate)

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
.build.cnn2d <- function(hp, data, sfs.dims) {
  dim1 <- sfs.dims[1]; dim2 <- sfs.dims[2]
  n_targets <- ncol(data$Y_train)

  l2 <- keras::regularizer_l2(hp$l2_reg)

  .drop <- function(x, rate) keras::layer_dropout(x, rate = rate)

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
  } else if (type == "emulator") {
    .prep.emulator.data(reftable, param.cols, exclude.cols, val.frac)
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
                       max_epochs, eta, seed, verbose) {

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

  # Precompute total Hyperband budget (config-epochs) for ETA estimation
  total_budget <- 0
  for (ss in s_max:0) {
    nn <- ceiling((s_max + 1) / (ss + 1)) * as.integer(eta^ss)
    rr <- max_epochs / eta^ss
    for (ii in 0:ss) {
      r_ii   <- round(rr * eta^ii)
      prev_r <- if (ii > 0) round(rr * eta^(ii - 1)) else 0L
      n_ii   <- max(1L, floor(nn / eta^ii))
      total_budget <- total_budget + n_ii * (r_ii - prev_r)
    }
  }
  consumed_budget <- 0
  global_t0 <- proc.time()[3]

  for (s in s_max:0) {
    bracket_t0 <- proc.time()[3]
    n <- ceiling((s_max + 1) / (s + 1)) * as.integer(eta^s)
    r <- max_epochs / eta^s

    if (verbose) cat(sprintf("  Bracket %d | %d configs \u00d7 %d epochs\n",
                             s, n, round(r)))

    # Sample configs
    set.seed(seed + s)
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
      cfg_times  <- rep(NA_real_, length(configs))

      for (j in seq_along(configs)) {
        tryCatch({
          cfg_t0 <- proc.time()[3]
          tensorflow::tf$random$set_seed(as.integer(seed + s + j))
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
          cfg_times[j]  <- proc.time()[3] - cfg_t0

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
            if (verbose) cat(sprintf("    \u2605 new best: val_loss=%.4f (bracket %d, round %d)\n",
                                      val_losses[j], s, i))
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
        consumed_budget <- consumed_budget + (r_i - prev_epochs[j])
      }

      # Log per-config val_loss and timing
      if (verbose) {
        items <- character(length(configs))
        for (k in seq_along(configs)) {
          vl_str <- if (is.finite(val_losses[k])) sprintf("%.4f", val_losses[k]) else "  NA "
          t_str  <- if (!is.na(cfg_times[k])) sprintf("%.1fs", cfg_times[k]) else "?"
          items[k] <- sprintf("cfg %d: %s (%s)", k, vl_str, t_str)
        }
        for (start in seq(1, length(items), by = 4)) {
          end <- min(start + 3, length(items))
          cat("    ", paste(items[start:end], collapse = "   "), "\n")
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

        if (verbose && length(configs) > 1) {
          sorted_losses <- val_losses[ranking]
          kept_str <- paste(sprintf("%.4f", sorted_losses[seq_len(length(keep))]),
                            collapse = " ")
          if (length(ranking) > length(keep)) {
            pruned_str <- paste(sprintf("%.4f",
                                sorted_losses[(length(keep) + 1):length(ranking)]),
                                collapse = " ")
            cat(sprintf("    Ranking: %s | %s\n", kept_str, pruned_str))
            pad <- nchar(sprintf("    Ranking: %s ", kept_str))
            cat(sprintf("%s^-- cutoff (top %d kept)\n", strrep(" ", pad), length(keep)))
          }
        }

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

    if (verbose) {
      bracket_elapsed <- (proc.time()[3] - bracket_t0) / 60
      global_elapsed  <- (proc.time()[3] - global_t0) / 60
      if (consumed_budget > 0 && consumed_budget < total_budget) {
        eta_min <- global_elapsed * (total_budget - consumed_budget) / consumed_budget
        cat(sprintf("    Bracket %d done in %.1f min | elapsed %.1f min, ~%.0f min remaining\n\n",
                    s, bracket_elapsed, global_elapsed, eta_min))
      } else {
        cat(sprintf("    Bracket %d done in %.1f min | elapsed %.1f min\n\n",
                    s, bracket_elapsed, global_elapsed))
      }
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
# Internal: generate diverse (eta, max_epochs) configs for parallel searches
# ============================================================================

.generate.search.configs <- function(n_searches, base_max_epochs, base_eta) {
  base_max_epochs <- as.integer(base_max_epochs)
  base_eta        <- as.integer(base_eta)

  if (n_searches == 1L)
    return(data.frame(eta = base_eta, max_epochs = base_max_epochs))

  # Pool of diverse (eta, max_epochs) combos
  pool_eta <- c(2L, 3L, 4L)
  pool_me  <- unique(as.integer(round(
    base_max_epochs * c(0.6, 1.0, 2.0)
  )))
  pool <- expand.grid(eta = pool_eta, max_epochs = pool_me,
                      KEEP.OUT.ATTRS = FALSE)

  # Put the user's base config first
  base_idx <- which(pool$eta == base_eta & pool$max_epochs == base_max_epochs)
  if (length(base_idx) > 0L) {
    pool <- rbind(pool[base_idx[1L], ], pool[-base_idx[1L], ])
  } else {
    pool <- rbind(data.frame(eta = base_eta, max_epochs = base_max_epochs), pool)
  }

  # Cycle through pool to fill n_searches rows
  idx <- rep_len(seq_len(nrow(pool)), n_searches)
  pool[idx, , drop = FALSE]
}

# ============================================================================
# Internal: orchestrate K independent serial Hyperband searches in parallel
# ============================================================================

.parallel.hyperband.search <- function(data, search_space, type, sfs.dims,
                                        max_epochs, eta,
                                        n_searches, cores, gpus,
                                        gpu.threshold, greedy, top_k, seed, verbose) {

  n_concurrent <- min(as.integer(cores), as.integer(n_searches))

  # Build per-search (eta, max_epochs) configs
  if (length(max_epochs) == 1L && length(eta) == 1L) {
    cfgs <- .generate.search.configs(n_searches, max_epochs, eta)
  } else {
    cfgs <- data.frame(
      eta        = rep_len(as.integer(eta), n_searches),
      max_epochs = rep_len(as.integer(max_epochs), n_searches)
    )
  }
  search_max_epochs <- cfgs$max_epochs   # integer vector, length n_searches
  search_eta        <- cfgs$eta          # integer vector, length n_searches

  if (verbose) {
    cat(sprintf("PipeMaster:: Launching %d searches (%d concurrent)\n",
                n_searches, n_concurrent))
    for (k in seq_len(n_searches))
      cat(sprintf("  Search %d: eta=%d, max_epochs=%d\n",
                  k, search_eta[k], search_max_epochs[k]))
  }

  # ==========================================================================
  # Sequential mode (cores = 1): run searches in main process
  # ==========================================================================
  if (n_concurrent <= 1L) {
    all_search_results <- list()
    # Ranked list of top-K: each entry = list(val_loss, hp, model_dir)
    topk_entries <- list()

    for (k in seq_len(n_searches)) {
      search_seed <- seed + (k - 1L) * 10000L
      if (verbose)
        cat(sprintf("\n=== Search %d/%d (seed=%d) ===\n", k, n_searches, search_seed))

      hb <- .hyperband(search_space = search_space, data = data,
                       type = type, sfs.dims = sfs.dims,
                       max_epochs = search_max_epochs[k], eta = search_eta[k],
                       seed = search_seed, verbose = verbose)

      # Retrain best config from this search
      tensorflow::tf$random$set_seed(as.integer(search_seed))
      model <- .build.nn(hb$best_hp, data, type, sfs.dims)

      weights_loaded <- tryCatch({
        keras::load_model_weights_tf(model, file.path(hb$best_weights_path, "ckpt"))
        TRUE
      }, error = function(e) FALSE)

      start_epoch <- if (weights_loaded) as.integer(hb$best_epochs) else 0L
      final_vl <- hb$best_val_loss

      if (start_epoch < search_max_epochs[k]) {
        if (verbose)
          cat(sprintf("PipeMaster:: Retraining search %d best for %d epochs (warm-start from %d)...\n",
                      k, search_max_epochs[k], start_epoch))

        retrain_history <- model |> keras::fit(
          x = data$X_train, y = data$Y_train,
          validation_data = list(data$X_val, data$Y_val),
          epochs        = as.integer(search_max_epochs[k]),
          initial_epoch = as.integer(start_epoch),
          batch_size    = as.integer(hb$best_hp$batch_size),
          callbacks     = list(
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
          final_vl <- min(final_vl, min(retrain_vl[is.finite(retrain_vl)]))
      }

      if (verbose)
        cat(sprintf("PipeMaster:: Search %d final val_loss: %.6f\n", k, final_vl))

      search_res <- hb$all_results
      search_res$search <- k
      all_search_results[[k]] <- search_res

      # Save model to temp dir and insert into ranked top-K list
      model_dir <- tempfile(sprintf("topk_model_%d_", k))
      keras::save_model_tf(model, model_dir)

      topk_entries[[length(topk_entries) + 1L]] <- list(
        val_loss  = final_vl,
        hp        = hb$best_hp,
        model_dir = model_dir
      )

      # Sort by val_loss and prune to top_k
      if (length(topk_entries) > top_k) {
        vl_vec <- vapply(topk_entries, `[[`, numeric(1), "val_loss")
        ord <- order(vl_vec)
        # Remove worst entry (last after sorting)
        worst_idx <- ord[length(ord)]
        unlink(topk_entries[[worst_idx]]$model_dir, recursive = TRUE)
        topk_entries <- topk_entries[ord[seq_len(top_k)]]
      }

      # Clean up
      unlink(hb$best_weights_path, recursive = TRUE)
      rm(model); gc()
      tryCatch(keras::k_clear_session(), error = function(e) NULL)
    }

    # Sort final entries by val_loss
    vl_vec <- vapply(topk_entries, `[[`, numeric(1), "val_loss")
    topk_entries <- topk_entries[order(vl_vec)]

    # Load all top-K models from disk
    models <- lapply(topk_entries, function(entry) {
      tryCatch(
        keras::load_model_tf(entry$model_dir),
        error = function(e) {
          warning("Could not reload model (val_loss=",
                  sprintf("%.6f", entry$val_loss), "): ",
                  conditionMessage(e), call. = FALSE)
          NULL
        }
      )
    })
    models_val_loss <- vapply(topk_entries, `[[`, numeric(1), "val_loss")

    # Clean up temp dirs
    for (entry in topk_entries) unlink(entry$model_dir, recursive = TRUE)

    # Remove NULLs from failed loads
    keep <- !vapply(models, is.null, logical(1))
    models <- models[keep]
    models_val_loss <- models_val_loss[keep]

    if (length(models) == 0L)
      stop("All models failed to reload from disk")

    combined_results <- do.call(rbind, all_search_results)
    best_hp <- topk_entries[[1L]]$hp

    # Compute R² and MPE for each model
    models_metrics <- lapply(models, .compute.model.metrics, data = data, type = type)

    if (verbose) {
      cat(sprintf("\nPipeMaster:: Best across %d searches: val_loss=%.6f\n",
                  n_searches, models_val_loss[1L]))
      cat(sprintf("PipeMaster:: Kept top %d models (val_loss: %s)\n",
                  length(models),
                  paste(sprintf("%.6f", models_val_loss), collapse = ", ")))
      cat(sprintf("PipeMaster:: Best config: %s\n", .hp.to.string(best_hp, type)))
      .print.model.metrics(models_metrics, models_val_loss, verbose)
    }

    return(list(
      best_hp         = best_hp,
      best_val_loss   = models_val_loss[1L],
      all_results     = combined_results,
      best_model      = models[[1L]],
      models          = models,
      models_val_loss = models_val_loss,
      models_metrics  = models_metrics
    ))
  }

  # ==========================================================================
  # Parallel mode: launch Rscript workers
  # ==========================================================================
  work_dir <- tempfile("hb_search_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)

  # Save shared data for workers
  X_train    <- data$X_train
  X_val      <- data$X_val
  Y_train    <- data$Y_train
  Y_val      <- data$Y_val
  n_features <- if (type %in% c("sumstat", "emulator")) ncol(data$X_train) else data$n_features
  n_bins     <- data$n_bins

  search_seeds        <- seed + (seq_len(n_searches) - 1L) * 10000L
  saved_lib_paths     <- .libPaths()
  threads_per_worker  <- .compute.threads.per.worker(n_concurrent, greedy)
  search_space_saved  <- search_space

  # Detect if PipeMaster was loaded via devtools::load_all()
  pkg_source_dir <- ""
  if ("devtools_shims" %in% search()) {
    pkg_path <- find.package("PipeMaster", quiet = TRUE)
    if (file.exists(file.path(pkg_path, "DESCRIPTION"))) pkg_source_dir <- pkg_path
  }

  save(X_train, X_val, Y_train, Y_val,
       type, sfs.dims, n_features, n_bins,
       search_seeds, saved_lib_paths, threads_per_worker,
       search_max_epochs, search_eta, search_space_saved,
       pkg_source_dir,
       file = file.path(work_dir, "shared_search.RData"))

  # Write worker script
  .write.search.worker.script(file.path(work_dir, "_search_worker.R"))

  # Build task list with per-task GPU assignment
  n_gpu_searches <- if (gpus > 0L) min(n_searches, gpu.threshold * gpus) else 0L

  tasks <- lapply(seq_len(n_searches), function(k) {
    if (gpus > 0L && k <= n_gpu_searches) {
      gpu_id <- (k - 1L) %% gpus
      task_gpu_env <- sprintf("CUDA_VISIBLE_DEVICES=%d TF_FORCE_GPU_ALLOW_GROWTH=true", gpu_id)
    } else {
      task_gpu_env <- "CUDA_VISIBLE_DEVICES=-1"
    }
    list(
      script  = "_search_worker.R",
      id      = k,
      result  = sprintf("results/search_%04d/done.txt", k),
      prefix  = "search",
      gpu_env = task_gpu_env
    )
  })

  pool_result <- .launch.rscript.pool(
    tasks, n_concurrent, work_dir,
    gpus = 0L,
    verbose = verbose, max_retries = 1L)

  # Collect results with val_loss for ranking
  all_search_results <- list()
  search_entries <- list()  # list of (val_loss, hp, search_dir)

  for (k in seq_len(n_searches)) {
    search_dir <- file.path(results_dir, sprintf("search_%04d", k))
    rds_file   <- file.path(search_dir, "result.rds")

    if (file.exists(rds_file)) {
      res <- tryCatch(readRDS(rds_file), error = function(e) NULL)
      if (!is.null(res)) {
        search_res <- res$all_results
        search_res$search <- k
        all_search_results[[length(all_search_results) + 1L]] <- search_res

        if (verbose)
          cat(sprintf("  Search %d: val_loss=%.6f\n", k, res$best_val_loss))

        search_entries[[length(search_entries) + 1L]] <- list(
          val_loss   = res$best_val_loss,
          hp         = res$best_hp,
          search_dir = search_dir
        )
      }
    } else {
      if (verbose) cat(sprintf("  Search %d: FAILED (no result)\n", k))
    }
  }

  if (length(search_entries) == 0L)
    stop("All parallel searches failed. Check worker logs in ", work_dir)

  # Sort by val_loss and take top-K
  vl_vec <- vapply(search_entries, `[[`, numeric(1), "val_loss")
  ord <- order(vl_vec)
  n_keep <- min(top_k, length(search_entries))
  search_entries <- search_entries[ord[seq_len(n_keep)]]

  # Load top-K models before cleaning work_dir (suppressWarnings: TF/reticulate
  # emits benign "NAs introduced by coercion" during model deserialization)
  models <- lapply(search_entries, function(entry) {
    tryCatch(
      suppressWarnings(keras::load_model_tf(file.path(entry$search_dir, "best_model"))),
      error = function(e) {
        warning("Could not load model (val_loss=",
                sprintf("%.6f", entry$val_loss), "): ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  })
  models_val_loss <- vapply(search_entries, `[[`, numeric(1), "val_loss")

  # Remove NULLs from failed loads
  keep <- !vapply(models, is.null, logical(1))
  models <- models[keep]
  models_val_loss <- models_val_loss[keep]

  if (length(models) == 0L)
    stop("All top-K models failed to load. Check worker logs in ", work_dir)

  best_hp <- search_entries[[1L]]$hp

  combined_results <- if (length(all_search_results) > 0)
    do.call(rbind, all_search_results) else data.frame()

  # Compute R² and MPE for each model
  models_metrics <- lapply(models, .compute.model.metrics, data = data, type = type)

  if (verbose) {
    cat(sprintf("\nPipeMaster:: Best across %d searches: val_loss=%.6f\n",
                n_searches, models_val_loss[1L]))
    cat(sprintf("PipeMaster:: Kept top %d models (val_loss: %s)\n",
                length(models),
                paste(sprintf("%.6f", models_val_loss), collapse = ", ")))
    cat(sprintf("PipeMaster:: Best config: %s\n", .hp.to.string(best_hp, type)))
    .print.model.metrics(models_metrics, models_val_loss, verbose)
  }

  Sys.sleep(1)
  unlink(work_dir, recursive = TRUE)

  list(
    best_hp         = best_hp,
    best_val_loss   = models_val_loss[1L],
    all_results     = combined_results,
    best_model      = models[[1L]],
    models          = models,
    models_val_loss = models_val_loss,
    models_metrics  = models_metrics
  )
}

# ============================================================================
# Internal: write standalone search worker Rscript
# ============================================================================

.write.search.worker.script <- function(filepath) {
  writeLines(c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    '# Load shared data',
    'load("shared_search.RData")',
    '',
    '# Threading env',
    'n_threads <- as.character(threads_per_worker)',
    'Sys.setenv(TF_NUM_INTRAOP_THREADS = n_threads,',
    '           TF_NUM_INTEROP_THREADS = "1",',
    '           OMP_NUM_THREADS = n_threads)',
    '',
    '# Restore library paths and load PipeMaster',
    '.libPaths(saved_lib_paths)',
    'if (nzchar(pkg_source_dir)) {',
    '  suppressPackageStartupMessages(devtools::load_all(pkg_source_dir, quiet = TRUE))',
    '} else {',
    '  suppressPackageStartupMessages(library(PipeMaster))',
    '}',
    '',
    'suppressPackageStartupMessages({',
    '  library(keras)',
    '  library(tensorflow)',
    '})',
    '',
    '# Enable GPU memory growth',
    'tryCatch({',
    '  tf_gpus <- tf$config$list_physical_devices("GPU")',
    '  for (gpu in tf_gpus)',
    '    tf$config$experimental$set_memory_growth(gpu, TRUE)',
    '}, error = function(e) NULL)',
    '',
    '# Create task-specific output directory',
    'out_dir <- file.path("results", sprintf("search_%04d", task_id))',
    'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)',
    '',
    '# Build data list',
    'hb_data <- list(',
    '  X_train = X_train, X_val = X_val,',
    '  Y_train = Y_train, Y_val = Y_val,',
    '  n_features = n_features, n_bins = n_bins',
    ')',
    '',
    '# Per-search config',
    'worker_seed       <- search_seeds[task_id]',
    'worker_max_epochs <- search_max_epochs[task_id]',
    'worker_eta        <- search_eta[task_id]',
    'cat(sprintf("Search %d starting (seed=%d, eta=%d, max_epochs=%d)\\n",',
    '            task_id, worker_seed, worker_eta, worker_max_epochs))',
    '',
    '# Check parent alive before starting expensive work',
    'if (file.exists(".PM_parent.pid") && !PipeMaster:::.pm.parent.alive(".PM_parent.pid")) {',
    '  cat("Parent died, worker exiting.\\n"); q("no") }',
    '',
    'hb <- PipeMaster:::.hyperband(',
    '  search_space = search_space_saved,',
    '  data = hb_data, type = type, sfs.dims = sfs.dims,',
    '  max_epochs = worker_max_epochs, eta = worker_eta,',
    '  seed = worker_seed, verbose = TRUE)',
    '',
    '# Check parent alive before retraining',
    'if (file.exists(".PM_parent.pid") && !PipeMaster:::.pm.parent.alive(".PM_parent.pid")) {',
    '  cat("Parent died, worker exiting.\\n"); q("no") }',
    '',
    '# Retrain best config to worker_max_epochs',
    'tf$random$set_seed(as.integer(worker_seed))',
    'model <- PipeMaster:::.build.nn(hb$best_hp, hb_data, type, sfs.dims)',
    '',
    'weights_loaded <- tryCatch({',
    '  load_model_weights_tf(model, file.path(hb$best_weights_path, "ckpt"))',
    '  TRUE',
    '}, error = function(e) FALSE)',
    '',
    'start_epoch <- if (weights_loaded) as.integer(hb$best_epochs) else 0L',
    'final_vl <- hb$best_val_loss',
    '',
    'if (start_epoch < worker_max_epochs) {',
    '  cat(sprintf("Search %d: retraining from epoch %d to %d\\n",',
    '              task_id, start_epoch, worker_max_epochs))',
    '  retrain_h <- model |> fit(',
    '    x = hb_data$X_train, y = hb_data$Y_train,',
    '    validation_data = list(hb_data$X_val, hb_data$Y_val),',
    '    epochs        = worker_max_epochs,',
    '    initial_epoch = as.integer(start_epoch),',
    '    batch_size    = as.integer(hb$best_hp$batch_size),',
    '    callbacks     = list(',
    '      callback_early_stopping(monitor = "val_loss", patience = 30L,',
    '                              restore_best_weights = TRUE),',
    '      callback_reduce_lr_on_plateau(monitor = "val_loss", patience = 15L,',
    '                                    factor = 0.5, min_lr = 1e-6, verbose = 0L)',
    '    ),',
    '    verbose = 0L',
    '  )',
    '  rvl <- retrain_h$metrics$val_loss',
    '  if (is.null(rvl)) rvl <- retrain_h$history$val_loss',
    '  rvl <- unlist(rvl)',
    '  if (length(rvl) > 0 && any(is.finite(rvl)))',
    '    final_vl <- min(final_vl, min(rvl[is.finite(rvl)]))',
    '}',
    '',
    '# Save model and results',
    'save_model_tf(model, file.path(out_dir, "best_model"))',
    'saveRDS(list(',
    '  best_hp       = hb$best_hp,',
    '  best_val_loss = final_vl,',
    '  all_results   = hb$all_results',
    '), file.path(out_dir, "result.rds"))',
    '',
    '# Clean up Hyperband temp weights',
    'unlink(hb$best_weights_path, recursive = TRUE)',
    '',
    'cat(sprintf("Search %d done (val_loss=%.6f)\\n", task_id, final_vl))',
    'writeLines("done", file.path(out_dir, "done.txt"))',
    'k_clear_session()'
  ), filepath)
}

# ============================================================================
# save / load tune.nn results (keras model needs special serialization)
# ============================================================================

#' Save tune.nn Result to Disk
#'
#' Saves the output of \code{tune.nn()} so it can be loaded on another machine.
#' The keras model is serialized separately using \code{keras::save_model_tf()}.
#' When multiple models are present (from \code{top_k > 1}), each model is saved
#' to a separate subdirectory under \code{models/}.
#'
#' @param tune.result list — output from \code{tune.nn()}.
#' @param path character — directory path where files will be saved (created if needed).
#'
#' @export
save.tune.result <- function(tune.result, path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  nn_backend <- .detect.backend(tune.result)

  if (nn_backend == "torch") {
    # Torch: save models as .pt files (need metadata for reconstruction)
    data <- tune.result$data
    tp <- tune.result$type
    hp <- tune.result$best_hp
    n_feat <- if (tp %in% c("sumstat", "emulator")) ncol(data$X_train) else data$n_features
    n_targ <- ncol(data$Y_train)
    n_bins <- data$n_bins  # NULL for sumstat
    sfs_dims <- tune.result$sfs.dims

    model_file <- file.path(path, "best_model.pt")
    .torch.save.model(tune.result$best_model, model_file, tp, hp,
                      n_features = n_feat, n_bins = n_bins,
                      sfs_dims = sfs_dims, n_targets = n_targ)

    models <- tune.result$models
    if (!is.null(models) && length(models) > 0L) {
      models_dir <- file.path(path, "models")
      dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
      for (i in seq_along(models)) {
        # Use per-model HP when available (top-K ensemble members may have
        # different architectures via Hyperband search). Falls back to
        # best_hp for older in-memory results without models_hp.
        model_hp <- if (!is.null(tune.result$models_hp) &&
                        i <= length(tune.result$models_hp))
          tune.result$models_hp[[i]] else hp
        .torch.save.model(models[[i]],
                          file.path(models_dir, sprintf("model_%03d.pt", i)),
                          tp, model_hp,
                          n_features = n_feat, n_bins = n_bins,
                          sfs_dims = sfs_dims, n_targets = n_targ)
      }
    }
  } else {
    # Keras: save models as TF SavedModel
    if (!requireNamespace("keras", quietly = TRUE))
      stop("save.tune.result() requires the 'keras' package for keras models.")

    model_dir <- file.path(path, "best_model")
    keras::save_model_tf(tune.result$best_model, model_dir)

    models <- tune.result$models
    if (!is.null(models) && length(models) > 0L) {
      models_dir <- file.path(path, "models")
      dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
      for (i in seq_along(models))
        keras::save_model_tf(models[[i]],
                             file.path(models_dir, sprintf("model_%03d", i)))
    }
  }

  # Save everything else as RDS (strip model objects)
  result_no_model <- tune.result
  result_no_model$best_model <- NULL
  result_no_model$models <- NULL
  saveRDS(result_no_model, file.path(path, "tune_result.rds"))

  cat(sprintf("PipeMaster:: Saved tune.nn result to %s (%d models, backend=%s)\n",
              path, length(tune.result$models), nn_backend))
}

#' Load tune.nn Result from Disk
#'
#' Loads a tune.nn() result previously saved with \code{save.tune.result()}.
#' If the saved result contains multiple models (from \code{top_k > 1}),
#' all models are loaded from the \code{models/} subdirectory. For results
#' saved before \code{top_k} support, falls back to wrapping the single
#' best model in a length-1 list.
#'
#' @param path character — directory path where files were saved.
#'
#' @return A list identical to the output of \code{tune.nn()}.
#'
#' @export
load.tune.result <- function(path) {
  rds_file <- file.path(path, "tune_result.rds")
  if (!file.exists(rds_file))
    stop("tune_result.rds not found in: ", path)

  result <- readRDS(rds_file)
  nn_backend <- if (!is.null(result$backend)) result$backend else "keras"

  if (nn_backend == "torch") {
    if (!requireNamespace("torch", quietly = TRUE))
      stop("load.tune.result() requires the 'torch' package for torch models.")

    best_pt <- file.path(path, "best_model.pt")
    if (!file.exists(best_pt))
      stop("best_model.pt not found in: ", path)
    result$best_model <- .torch.load.model(best_pt)

    models_dir <- file.path(path, "models")
    if (dir.exists(models_dir)) {
      model_files <- sort(list.files(models_dir, pattern = "\\.pt$",
                                     full.names = TRUE))
      if (length(model_files) > 0L) {
        result$models <- lapply(model_files, function(f) {
          tryCatch(.torch.load.model(f), error = function(e) {
            warning("Could not load model from ", f, ": ",
                    conditionMessage(e), call. = FALSE)
            NULL
          })
        })
        keep <- !vapply(result$models, is.null, logical(1))
        result$models <- result$models[keep]
        if (!is.null(result$models_val_loss))
          result$models_val_loss <- result$models_val_loss[keep]
      } else {
        result$models <- list(result$best_model)
      }
    } else {
      result$models <- list(result$best_model)
    }
  } else {
    if (!requireNamespace("keras", quietly = TRUE))
      stop("load.tune.result() requires the 'keras' package for keras models.")

    model_dir <- file.path(path, "best_model")
    if (!dir.exists(model_dir))
      stop("best_model/ directory not found in: ", path)
    result$best_model <- keras::load_model_tf(model_dir)

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
      } else {
        result$models <- list(result$best_model)
      }
    } else {
      result$models <- list(result$best_model)
      if (is.null(result$models_val_loss))
        result$models_val_loss <- result$best_val_loss
    }
  }

  cat(sprintf("PipeMaster:: Loaded tune.nn result from %s (%d models, backend=%s)\n",
              path, length(result$models), nn_backend))
  result
}

# ============================================================================
# nn.predict — Posterior estimation via neural network
# ============================================================================

#' Posterior Estimation via Neural Network Uncertainty Quantification
#'
#' Estimates posterior distributions for observed data using a trained neural
#' network from \code{tune.nn()}, with multiple methods for uncertainty
#' quantification.
#'
#' When the tune result contains multiple models (from \code{top_k > 1}), the
#' point estimate is computed as a proximity-weighted ensemble average, where
#' each model is weighted by how well it predicts validation samples near the
#' observed data point.
#'
#' @param tune.result list — output from \code{tune.nn()}.
#' @param observed named numeric vector or 1-row data.frame of observed summary
#'   statistics (or SFS bins).
#' @param reftable data.frame — original reference table (required for
#'   \code{"bootstrap"} and \code{"ABC_NN_regression"} methods).
#' @param param.cols character vector — parameter column names (required when
#'   reftable is provided).
#' @param type character — architecture type: \code{"sumstat"}, \code{"sfs1d"},
#'   or \code{"sfs2d"}. If NULL, uses the type stored in tune.result.
#' @param sfs.dims integer vector — for 2D CNN only. If NULL, uses tune.result.
#' @param method character — one or more of:
#'   \describe{
#'     \item{\code{"point"}}{Fast point estimate. No retraining or resampling.
#'       \code{reftable} and \code{param.cols} not required.}
#'     \item{\code{"bootstrap"}}{Locally-weighted residual bootstrap. Captures
#'       prediction uncertainty due to simulator stochasticity and emulator
#'       imprecision. Computes residuals (true - predicted) on validation data,
#'       weights them by proximity of their summary statistics to observed data
#'       via an Epanechnikov kernel, and samples \code{n_boot} residuals to add
#'       to the point estimate. Fast (single forward pass, no retraining).
#'       Note: this is NOT a Bayesian posterior — it is a frequentist prediction
#'       interval centered on the point estimate.}
#'     \item{\code{"ABC_NN_regression"}}{ABC with neural network regression
#'       adjustment (Beaumont et al. 2002). Produces a proper Bayesian posterior
#'       that incorporates both prior uncertainty and simulator stochasticity.
#'       Accepts reftable points within \code{tolerance} of the observed data
#'       (by Euclidean distance in standardized stat space), then adjusts
#'       accepted parameter values using the NN as a local regression correction:
#'       \eqn{\theta_{adj} = \theta_i - (g(S_i) - g(S_{obs}))}, where \eqn{g}
#'       is the trained neural network. Fast (single forward pass, no retraining).}
#'   }
#' @param n_boot integer — number of residual bootstrap samples (default 1000).
#' @param tolerance numeric — ABC acceptance tolerance as a quantile of pairwise
#'   distances (default 0.01, i.e. closest 1\% of reftable). Only used for
#'   \code{"ABC_NN_regression"}.
#' @param pls logical — project summary statistics into PLS space before
#'   computing distances for \code{"bootstrap"} and \code{"ABC_NN_regression"}
#'   (default FALSE). Recommended when the number of statistics exceeds ~30.
#'   The NN prediction still uses the full feature space; only the distance
#'   computation (acceptance/kernel weighting) uses the PLS projection.
#' @param n.pls integer — number of PLS components (default 20). Should be
#'   at least equal to the number of free parameters.
#' @param seed integer — random seed (default 42).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return An object of class \code{"posterior"} (a list) with:
#' \describe{
#'   \item{point_estimate}{named numeric vector — inverse-transformed point
#'     prediction from best model}
#'   \item{bootstrap}{matrix of bootstrap samples (n_boot x n_params), or NULL.
#'     Locally-weighted residual bootstrap — prediction intervals, not a
#'     Bayesian posterior.}
#'   \item{abc}{matrix of regression-adjusted ABC posterior samples
#'     (n_accepted x n_params), or NULL. Proper Bayesian posterior.}
#'   \item{prior}{matrix of prior samples from reftable, or NULL}
#'   \item{param_names}{character vector of parameter column names}
#' }
#' Use \code{summary()} for posterior quantiles, \code{plot()} for density
#' curves, \code{density()} for density objects.
#'
#' @references
#' Beaumont, M. A., Zhang, W., & Balding, D. J. (2002). Approximate Bayesian
#' computation in population genetics. Genetics, 162(4), 2025-2035.
#'
#' @export
nn.predict <- function(tune.result, observed, reftable = NULL, param.cols = NULL,
                       type = NULL, sfs.dims = NULL,
                       method = c("point", "bootstrap", "ABC_NN_regression"),
                       n_boot = 1000L, tolerance = 0.01,
                       pls = FALSE, n.pls = 20L,
                       seed = 42, verbose = TRUE) {

  # --- Detect backend ---
  nn_backend <- .detect.backend(tune.result)

  # --- Dependency check ---
  if (nn_backend == "keras") {
    if (!requireNamespace("keras", quietly = TRUE) ||
        !requireNamespace("tensorflow", quietly = TRUE))
      stop("nn.predict() with keras backend requires the 'keras' and 'tensorflow' R packages.")
  } else {
    if (!requireNamespace("torch", quietly = TRUE))
      stop("nn.predict() with torch backend requires the 'torch' R package.")
  }

  method <- match.arg(method, c("point", "bootstrap", "ABC_NN_regression"),
                      several.ok = TRUE)
  point_only <- length(method) == 1L && method == "point"

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

  needs_reftable <- any(c("bootstrap", "ABC_NN_regression") %in% method)
  if (needs_reftable && is.null(reftable))
    stop("reftable is required for bootstrap and ABC_NN_regression methods")
  if (!point_only && !is.null(reftable) && is.null(param.cols))
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
  method_labels <- c(point = "Point", bootstrap = "Residual Bootstrap",
                     ABC_NN_regression = "ABC-NN Regression")
  method_str <- paste(method_labels[method], collapse = " + ")

  if (verbose)
    cat(sprintf("PipeMaster:: nn.predict \u2014 %s\n", method_str))

  # --- Prepare observed data ---
  X_obs <- .prep.observed(observed, data, type, sfs.dims)

  if (verbose) {
    n_obs_feat <- length(observed)
    cat(sprintf("PipeMaster:: Observed: %d %s\n", n_obs_feat,
                if (type == "sumstat") "summary statistics" else "SFS bins"))
  }

  # --- Point estimate ---
  models <- tune.result$models
  models_val_loss <- tune.result$models_val_loss
  use_ensemble <- !is.null(models) && length(models) > 1L

  if (use_ensemble) {
    point_est <- .ensemble.predict(models, models_val_loss, X_obs, data, type, verbose)
    names(point_est) <- param_names
    if (verbose) {
      est_str <- paste(sprintf("%s=%.0f", param_names, point_est), collapse = " ")
      cat(sprintf("PipeMaster:: Ensemble point estimate: %s\n", est_str))
    }
  } else {
    Y_z <- .predict.nn(best_model, X_obs)
    point_est <- as.numeric(.inv.transform(Y_z, data$target_mu, data$target_sd))
    names(point_est) <- param_names
    if (verbose) {
      est_str <- paste(sprintf("%s=%.0f", param_names, point_est), collapse = " ")
      cat(sprintf("PipeMaster:: Point estimate: %s\n", est_str))
    }
  }

  # --- Early return for point-only ---
  if (point_only) {
    result <- list(
      point_estimate = point_est,
      bootstrap      = NULL,
      abc            = NULL,
      prior          = NULL,
      param_names    = param_names
    )
    class(result) <- "posterior"
    return(result)
  }

  # --- Uncertainty quantification ---
  boot_samples <- NULL
  abc_samples  <- NULL
  do_bootstrap <- "bootstrap" %in% method
  do_abc       <- "ABC_NN_regression" %in% method

  if (do_bootstrap) {
    boot_samples <- .run.residual.bootstrap(
      best_model, data, reftable, param.cols, observed,
      type, sfs.dims, n_boot, point_est, seed, verbose,
      exclude.cols = exclude.cols, nn_backend = nn_backend,
      pls = pls, n.pls = n.pls
    )
    colnames(boot_samples) <- param_names
  }

  if (do_abc) {
    abc_samples <- .run.abc.nn.regression(
      best_model, data, reftable, param.cols, observed,
      type, sfs.dims, tolerance, point_est, verbose,
      exclude.cols = exclude.cols, nn_backend = nn_backend,
      pls = pls, n.pls = n.pls
    )
    colnames(abc_samples) <- param_names
  }

  # Clip to prior range — use reftable parameter columns as bounds
  if (!is.null(reftable) && !is.null(param_names)) {
    for (j in seq_along(param_names)) {
      p <- param_names[j]
      if (p %in% colnames(reftable)) {
        lo <- min(reftable[[p]])
        hi <- max(reftable[[p]])
        if (!is.null(boot_samples)) boot_samples[, j] <- pmax(lo, pmin(hi, boot_samples[, j]))
        if (!is.null(abc_samples))  abc_samples[, j]  <- pmax(lo, pmin(hi, abc_samples[, j]))
      }
    }
  } else {
    if (!is.null(boot_samples)) boot_samples[boot_samples < 0] <- 0
    if (!is.null(abc_samples))  abc_samples[abc_samples < 0]   <- 0
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
    bootstrap      = boot_samples,
    abc            = abc_samples,
    prior          = prior_samples,
    param_names    = param_names
  )
  class(result) <- "posterior"
  result
}

#' @export
summary.posterior <- function(object, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), ...) {
  param_names <- object$param_names
  n_params <- length(param_names)
  q_labels <- paste0(formatC(probs * 100, format = "fg"), "%")

  .summarize <- function(mat) {
    tbl <- matrix(NA, nrow = n_params, ncol = length(probs) + 3)
    colnames(tbl) <- c("Mean", "Median", "Mode", q_labels)
    rownames(tbl) <- param_names
    for (j in seq_len(n_params)) {
      vals <- mat[, j]
      vals <- vals[is.finite(vals)]
      tbl[j, "Mean"]   <- mean(vals)
      tbl[j, "Median"] <- median(vals)
      d <- density(vals)
      tbl[j, "Mode"]   <- d$x[which.max(d$y)]
      tbl[j, -(1:3)]   <- quantile(vals, probs = probs)
    }
    tbl
  }

  out <- list(point_estimate = object$point_estimate)
  if (!is.null(object$bootstrap))
    out$bootstrap <- .summarize(object$bootstrap)
  if (!is.null(object$abc))
    out$abc <- .summarize(object$abc)
  class(out) <- "summary.posterior"
  out
}

#' @export
print.summary.posterior <- function(x, digits = 2, ...) {
  cat("Point estimate:\n")
  print(round(x$point_estimate, digits))
  if (!is.null(x$bootstrap)) {
    cat("\nResidual bootstrap (prediction intervals):\n")
    print(round(x$bootstrap, digits))
  }
  if (!is.null(x$abc)) {
    cat("\nABC-NN regression posterior:\n")
    print(round(x$abc, digits))
  }
  invisible(x)
}

#' @export
print.posterior <- function(x, ...) {
  methods <- c()
  if (!is.null(x$bootstrap))
    methods <- c(methods, sprintf("bootstrap (%d samples)", nrow(x$bootstrap)))
  if (!is.null(x$abc))
    methods <- c(methods, sprintf("ABC-NN regression (%d samples)", nrow(x$abc)))
  if (length(methods) == 0) methods <- "point only"
  cat(sprintf("posterior object — %s\n", paste(methods, collapse = " + ")))
  cat("Point estimate:\n")
  print(round(x$point_estimate, 2))
  cat("\nUse summary() for posterior quantiles, density() for density objects.\n")
  invisible(x)
}

#' @export
density.posterior <- function(x, method = NULL, ...) {
  param_names <- x$param_names
  sample_methods <- c("bootstrap", "abc")

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
plot.posterior <- function(x, method = NULL, col = "red", lwd = 2,
                             show_point_est = TRUE, show_prior = FALSE,
                             true_values = NULL, bw.adjust = 1, ...) {
  param_names <- x$param_names
  n_params <- length(param_names)

  all_methods <- c("prior", "bootstrap", "abc")

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
  method_cols <- c(prior = "grey50", bootstrap = "blue", abc = "steelblue")
  if (length(col) == 1 && length(post_methods) == 1)
    method_cols[post_methods] <- col

  method_labels <- c(prior = "Prior", bootstrap = "Residual bootstrap",
                     abc = "ABC-NN regression")

  sample_methods <- intersect(methods, all_methods)

  ncol <- min(n_params, 5)
  nrow <- ceiling(n_params / ncol)
  par(mfrow = c(nrow, ncol))
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
        d <- density(vals, from = prior_range[1], to = prior_range[2], adjust = bw.adjust)
      } else {
        d <- density(vals, adjust = bw.adjust)
      }
      densities[[m]]   <- d
      density_lty[[m]] <- if (m == "prior") 3 else 1
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

    if (!is.null(true_values))
      abline(v = true_values[j], col = "darkred", lwd = lwd, lty = 4)

    # Legend
    legend_labels <- c()
    legend_cols   <- c()
    legend_lty    <- c()
    for (m in sample_methods) {
      if (!is.null(densities[[m]])) {
        legend_labels <- c(legend_labels, method_labels[m])
        legend_cols   <- c(legend_cols, method_cols[m])
        legend_lty    <- c(legend_lty, if (m == "prior") 3 else 1)
      }
    }
    if (show_point_est) {
      legend_labels <- c(legend_labels, "Point est.")
      legend_cols   <- c(legend_cols, "black")
      legend_lty    <- c(legend_lty, 2)
    }
    if (!is.null(true_values)) {
      legend_labels <- c(legend_labels, "True")
      legend_cols   <- c(legend_cols, "darkred")
      legend_lty    <- c(legend_lty, 4)
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
# Internal: batch version of .prep.observed() for an N x p sim matrix
#
# Applies the same type-specific augmentation + z-scoring as .prep.observed(),
# vectorized over rows. Returns a matrix (sumstat) or 3D/4D array (sfs1d/sfs2d)
# shaped as the trained model expects. Non-finite values are set to 0.
# ============================================================================

.prep.sims <- function(S_sim, data, type, sfs.dims) {
  feat_mu <- data$feat_mu
  feat_sd <- data$feat_sd

  if (type == "sumstat") {
    sim_aug <- cbind(S_sim, log1p(abs(S_sim)))
    if (ncol(sim_aug) != length(feat_mu))
      stop(sprintf("Dimension mismatch: sims have %d features (augmented) but model expects %d.",
                   ncol(sim_aug), length(feat_mu)))
    X <- sweep(sweep(sim_aug, 2, feat_mu, "-"), 2, feat_sd, "/")
    X[!is.finite(X)] <- 0
    return(X)
  }

  if (type == "sfs1d") {
    sim_log <- log1p(S_sim)
    X <- sweep(sweep(sim_log, 2, feat_mu, "-"), 2, feat_sd, "/")
    X[!is.finite(X)] <- 0
    return(array(X, dim = c(nrow(S_sim), ncol(S_sim), 1L)))
  }

  # sfs2d
  sim_log <- log1p(S_sim)
  X <- sweep(sweep(sim_log, 2, feat_mu, "-"), 2, feat_sd, "/")
  X[!is.finite(X)] <- 0
  array(X, dim = c(nrow(S_sim), sfs.dims[1], sfs.dims[2], 1L))
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
# Internal: compute R² and MPE on validation set in original scale
# ============================================================================

.compute.model.metrics <- function(model, data, type) {
  pred_z <- predict(model, data$X_val, verbose = 0L)
  inv_fn <- if (type == "emulator") .inv.transform.emulator else .inv.transform
  true_orig <- as.matrix(inv_fn(data$Y_val, data$target_mu, data$target_sd))
  pred_orig <- as.matrix(inv_fn(pred_z, data$target_mu, data$target_sd))

  param_names <- names(data$target_mu)
  n_params <- length(param_names)
  r2  <- numeric(n_params)
  mpe <- numeric(n_params)

  for (j in seq_len(n_params)) {
    tru <- true_orig[, j]
    prd <- pred_orig[, j]
    ss_res <- sum((tru - prd)^2)
    ss_tot <- sum((tru - mean(tru))^2)
    r2[j]  <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
    mpe[j] <- mean(abs(tru - prd) / pmax(abs(tru), 1e-12)) * 100
  }
  names(r2)  <- param_names
  names(mpe) <- param_names

  list(r2 = r2, mpe = mpe)
}

.print.model.metrics <- function(models_metrics, models_val_loss, verbose) {
  if (!verbose || length(models_metrics) == 0L) return(invisible(NULL))

  param_names <- names(models_metrics[[1L]]$r2)
  n_params <- length(param_names)

  # Build header
  hdr <- sprintf("  %-5s %10s", "Model", "val_loss")
  for (p in param_names) hdr <- paste0(hdr, sprintf("  R2(%-6s)", substr(p, 1, 6)))
  for (p in param_names) hdr <- paste0(hdr, sprintf("  MPE(%-5s)", substr(p, 1, 5)))

  cat("PipeMaster:: Model metrics (validation set, original scale):\n")
  cat(hdr, "\n")

  for (i in seq_along(models_metrics)) {
    row <- sprintf("  %-5d %10.6f", i, models_val_loss[i])
    for (j in seq_len(n_params))
      row <- paste0(row, sprintf("  %10.3f", models_metrics[[i]]$r2[j]))
    for (j in seq_len(n_params))
      row <- paste0(row, sprintf("  %9.1f%%", models_metrics[[i]]$mpe[j]))
    cat(row, "\n")
  }
}

# ============================================================================
# Internal: flatten features for Dense-only architecture
# ============================================================================

.flatten.features <- function(data, type) {
  # sumstat: already 2D matrix, pass through
  # sfs1d:  (n, bins, 1) -> (n, bins)
  # sfs2d:  (n, d1, d2, 1) -> (n, d1*d2)

  flatten_one <- function(X) {
    if (type == "sumstat") {
      return(X)
    } else if (type == "sfs1d") {
      dims <- dim(X)
      return(matrix(X, nrow = dims[1], ncol = dims[2]))
    } else {
      dims <- dim(X)
      return(matrix(X, nrow = dims[1], ncol = prod(dims[2:3])))
    }
  }

  X_train_flat <- flatten_one(data$X_train)
  X_val_flat   <- flatten_one(data$X_val)

  list(
    X_train    = X_train_flat,
    X_val      = X_val_flat,
    n_feat_flat = ncol(X_train_flat)
  )
}

# ============================================================================
# Internal: proximity-weighted ensemble prediction from top-K models
# ============================================================================

.ensemble.predict <- function(models, models_val_loss, X_obs, data, type, verbose) {
  n_models <- length(models)

  # 1. Flatten X_val and X_obs to 2D for distance computation
  flat <- .flatten.features(data, type)
  X_val_flat <- flat$X_val  # (n_val, d)

  # Flatten X_obs (single observation) to match X_val_flat columns
  if (type == "sumstat") {
    x_obs_flat <- as.numeric(X_obs)
  } else if (type == "sfs1d") {
    x_obs_flat <- as.numeric(X_obs[1, , 1])
  } else {
    # sfs2d: (1, d1, d2, 1) -> vector of length d1*d2
    x_obs_flat <- as.numeric(X_obs[1, , , 1])
  }

  # 2. Euclidean distances from X_obs to each validation sample
  diffs <- sweep(X_val_flat, 2, x_obs_flat, `-`)
  dists <- sqrt(rowSums(diffs^2))

  # 3. Gaussian kernel weights (bandwidth = median distance)
  bandwidth <- median(dists)
  if (bandwidth < .Machine$double.eps) bandwidth <- 1.0
  w <- exp(-dists^2 / (2 * bandwidth^2))
  w <- w / sum(w)  # normalize

  # 4. Per-model proximity-weighted loss on validation set
  Y_val <- data$Y_val  # (n_val, n_params), Z-scored
  prox_losses <- numeric(n_models)

  for (m in seq_len(n_models)) {
    pred_z <- .predict.nn(models[[m]], data$X_val)  # (n_val, n_params)
    per_sample_mse <- rowSums((pred_z - Y_val)^2)  # (n_val,)
    prox_losses[m] <- sum(w * per_sample_mse)
  }

  # 5. Softmax ensemble weights (temperature = median proximity loss)
  temperature <- median(prox_losses)
  if (temperature < .Machine$double.eps) temperature <- 1.0
  log_weights <- -prox_losses / temperature
  log_weights <- log_weights - max(log_weights)  # shift for numerical stability
  ensemble_w <- exp(log_weights)
  ensemble_w <- ensemble_w / sum(ensemble_w)

  if (verbose) {
    w_str <- paste(sprintf("%.3f", ensemble_w), collapse = ", ")
    vl_str <- paste(sprintf("%.6f", models_val_loss), collapse = ", ")
    pl_str <- paste(sprintf("%.6f", prox_losses), collapse = ", ")
    cat(sprintf("PipeMaster:: Ensemble %d models\n", n_models))
    cat(sprintf("PipeMaster::   val_loss:    [%s]\n", vl_str))
    cat(sprintf("PipeMaster::   prox_loss:   [%s]\n", pl_str))
    cat(sprintf("PipeMaster::   weights:     [%s]\n", w_str))
  }

  # 6. Weighted average of model predictions in Z-space
  pred_z_list <- lapply(seq_len(n_models), function(m) {
    as.numeric(.predict.nn(models[[m]], X_obs))
  })

  ensemble_z <- rep(0, length(pred_z_list[[1L]]))
  for (m in seq_len(n_models)) {
    ensemble_z <- ensemble_z + ensemble_w[m] * pred_z_list[[m]]
  }

  # Inverse transform from Z-space to original scale
  point_est <- as.numeric(.inv.transform(
    matrix(ensemble_z, nrow = 1), data$target_mu, data$target_sd
  ))
  point_est
}

# ============================================================================
# Internal: locally-weighted residual bootstrap
#
# Captures prediction uncertainty due to simulator stochasticity and emulator
# imprecision. Uses the already-trained model (no retraining):
# 1. Predict theta_hat_i = g(S_i) for all validation points
# 2. Compute residuals r_i = theta_i - theta_hat_i (original scale)
# 3. Compute distances d_i = ||S_i - S_obs|| in standardized stat space
# 4. Weight residuals via Epanechnikov kernel: K(u) = 3/4(1-u^2) for |u|<=1
#    where u = d_i / bandwidth, bandwidth = quantile(distances, 0.1)
# 5. Sample n_boot residuals with replacement, weighted by kernel
# 6. Return point_est + sampled residuals
# ============================================================================

.run.residual.bootstrap <- function(best_model, data, reftable, param.cols,
                                    observed, type, sfs.dims, n_boot,
                                    point_est, seed, verbose,
                                    exclude.cols = NULL, nn_backend = "keras",
                                    pls = FALSE, n.pls = 20L) {

  nuisance    <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params    <- length(target_cols)

  if (verbose)
    cat(sprintf("\nPipeMaster:: Residual bootstrap (%d samples%s)...\n",
                n_boot, if (pls) sprintf(", PLS %d", n.pls) else ""))

  # Use validation set from tune.result$data
  X_val <- data$X_val
  Y_val <- data$Y_val  # z-scored log-targets

  # Predict on validation set
  if (nn_backend == "torch") {
    Y_pred_z <- .torch.predict(best_model, X_val)
  } else {
    Y_pred_z <- predict(best_model, X_val, verbose = 0L)
  }

  # Inverse-transform both to original scale
  Y_val_orig  <- .inv.transform(Y_val, data$target_mu, data$target_sd)
  Y_pred_orig <- .inv.transform(Y_pred_z, data$target_mu, data$target_sd)

  # Residuals in original parameter scale: true - predicted
  residuals <- Y_val_orig - Y_pred_orig  # n_val x n_params

  # Compute distances in standardized feature space
  # X_val is already z-scored, so compute obs in same space
  X_obs_z <- .prep.observed(observed, data, type, sfs.dims)

  # Flatten for distance computation
  if (type == "sumstat") {
    X_val_flat <- X_val
    X_obs_flat <- matrix(X_obs_z, nrow = 1)
  } else {
    X_val_flat <- matrix(X_val, nrow = nrow(X_val))
    X_obs_flat <- matrix(X_obs_z, nrow = 1)
  }

  # PLS projection for distance computation (NN prediction uses full space)
  if (pls) {
    # Fit PLS on validation features → parameter residuals
    pls_fit <- pls.fit(X_val_flat, Y_val_orig - Y_pred_orig, n.comp = n.pls)
    X_val_flat <- pls.project(pls_fit, X_val_flat)
    X_obs_flat <- pls.project(pls_fit, X_obs_flat)
    if (verbose)
      cat(sprintf("PipeMaster:: PLS: %d stats → %d components\n",
                  ncol(X_val), ncol(X_val_flat)))
  }

  # Euclidean distances from each validation point to observed
  diffs <- sweep(X_val_flat, 2, X_obs_flat[1, ])
  distances <- sqrt(rowSums(diffs^2))

  # Epanechnikov kernel with adaptive bandwidth (10th percentile of distances)
  bandwidth <- quantile(distances, 0.10)
  if (bandwidth < .Machine$double.eps) bandwidth <- quantile(distances, 0.25)

  u <- distances / bandwidth
  weights <- ifelse(abs(u) <= 1, 0.75 * (1 - u^2), 0)

  # Fallback: if too few non-zero weights, use Gaussian kernel
  n_nonzero <- sum(weights > 0)
  if (n_nonzero < 50) {
    if (verbose)
      cat(sprintf("PipeMaster:: Epanechnikov kernel gave %d non-zero weights, falling back to Gaussian\n",
                  n_nonzero))
    sigma_d <- median(distances)
    weights <- exp(-0.5 * (distances / sigma_d)^2)
  }

  # Normalize weights
  weights <- weights / sum(weights)

  if (verbose) {
    eff_n <- 1 / sum(weights^2)
    cat(sprintf("PipeMaster:: Effective sample size: %.0f (of %d validation points)\n",
                eff_n, length(weights)))
  }

  # Sample residuals with replacement, weighted by proximity
  set.seed(seed)
  idx <- sample(length(weights), size = n_boot, replace = TRUE, prob = weights)

  # Bootstrap samples: point estimate + sampled residual
  boot_matrix <- matrix(NA, nrow = n_boot, ncol = n_params)
  for (j in seq_len(n_params))
    boot_matrix[, j] <- point_est[j] + residuals[idx, j]

  if (verbose)
    cat(sprintf("PipeMaster:: Residual bootstrap done — %d samples\n", n_boot))

  boot_matrix
}

# ============================================================================
# Internal: ABC with NN regression adjustment (Beaumont et al. 2002)
#
# Produces a proper Bayesian posterior incorporating both prior uncertainty
# and simulator stochasticity:
# 1. Compute distances ||S_i - S_obs|| in standardized stat space for all
#    reftable rows
# 2. Accept the closest `tolerance` fraction (e.g. 1% = 1000 from 100K)
# 3. Predict theta_hat_i = g(S_i) for accepted rows and g(S_obs) for observed
# 4. Regression adjustment: theta_adj_i = theta_i - (g(S_i) - g(S_obs))
#    This corrects each accepted theta toward what the NN predicts the
#    difference should be, sharpening the ABC posterior.
# ============================================================================

.run.abc.nn.regression <- function(best_model, data, reftable, param.cols,
                                   observed, type, sfs.dims, tolerance,
                                   point_est, verbose,
                                   exclude.cols = NULL, nn_backend = "keras",
                                   pls = FALSE, n.pls = 20L) {

  nuisance    <- c("mean.rate", "sd.rate")
  target_cols <- setdiff(param.cols, nuisance)
  n_params    <- length(target_cols)
  n_total     <- nrow(reftable)

  if (verbose)
    cat(sprintf("\nPipeMaster:: ABC-NN regression (tolerance=%.3f, n=%d%s)...\n",
                tolerance, n_total, if (pls) sprintf(", PLS %d", n.pls) else ""))

  # Prepare all reftable features in the same space as the model
  all_rows   <- seq_len(n_total)
  full_split <- .prep.reftable.split(reftable, param.cols, all_rows, type, sfs.dims,
                                      exclude.cols = exclude.cols)
  features_all <- full_split$features   # raw features
  targets_all  <- full_split$targets    # raw target params (raw scale)

  # Z-score features using the model's training normalization
  feat_mu <- data$feat_mu
  feat_sd <- data$feat_sd

  X_all_z <- t((t(features_all) - feat_mu) / feat_sd)

  # Prepare observed in same z-scored space (flattened)
  X_obs_z <- .prep.observed(observed, data, type, sfs.dims)

  # Flatten for distance computation
  if (type == "sumstat") {
    X_flat <- X_all_z
    obs_flat <- matrix(X_obs_z, nrow = 1)
  } else {
    X_flat <- matrix(X_all_z, nrow = n_total)
    obs_flat <- matrix(X_obs_z, nrow = 1)
  }

  # PLS projection for distance computation (NN prediction uses full space)
  if (pls) {
    pls_fit <- pls.fit(X_flat, targets_all, n.comp = n.pls)
    X_flat_pls <- pls.project(pls_fit, X_flat)
    obs_flat_pls <- pls.project(pls_fit, obs_flat)
    if (verbose)
      cat(sprintf("PipeMaster:: PLS: %d stats → %d components\n",
                  ncol(X_flat), ncol(X_flat_pls)))
  } else {
    X_flat_pls <- X_flat
    obs_flat_pls <- obs_flat
  }

  # Euclidean distances (in PLS space if enabled, full space otherwise)
  diffs <- sweep(X_flat_pls, 2, obs_flat_pls[1, ])
  distances <- sqrt(rowSums(diffs^2))

  # Accept closest tolerance fraction
  n_accept <- max(10L, floor(tolerance * n_total))
  threshold <- sort(distances, partial = n_accept)[n_accept]
  accepted <- which(distances <= threshold)
  if (length(accepted) > n_accept) accepted <- accepted[order(distances[accepted])][1:n_accept]

  if (verbose)
    cat(sprintf("PipeMaster:: Accepted %d/%d (threshold distance=%.3f)\n",
                length(accepted), n_total, threshold))

  # Get true parameter values for accepted rows (raw scale from reftable)
  theta_accepted <- targets_all[accepted, , drop = FALSE]

  # Prepare accepted features for NN prediction — reshape for CNN if needed
  X_acc_z <- X_all_z[accepted, , drop = FALSE]
  if (type == "sfs1d") {
    dim(X_acc_z) <- c(nrow(X_acc_z), ncol(X_acc_z), 1L)
  } else if (type == "sfs2d") {
    X_acc_z <- array(X_acc_z, dim = c(nrow(X_acc_z), sfs.dims[1], sfs.dims[2], 1L))
  }

  # Predict theta_hat for accepted rows and for observed
  if (nn_backend == "torch") {
    Y_acc_pred_z <- .torch.predict(best_model, X_acc_z)
    Y_obs_pred_z <- .torch.predict(best_model, X_obs_z)
  } else {
    Y_acc_pred_z <- predict(best_model, X_acc_z, verbose = 0L)
    Y_obs_pred_z <- predict(best_model, X_obs_z, verbose = 0L)
  }

  # Inverse-transform predictions to original scale
  theta_hat_acc <- .inv.transform(Y_acc_pred_z, data$target_mu, data$target_sd)
  theta_hat_obs <- as.numeric(.inv.transform(Y_obs_pred_z, data$target_mu, data$target_sd))

  # targets_all from .prep.reftable.split() is already in original scale
  theta_true_acc <- theta_accepted

  # Regression adjustment: theta_adj = theta_true - (theta_hat_acc - theta_hat_obs)
  abc_matrix <- matrix(NA, nrow = length(accepted), ncol = n_params)
  for (j in seq_len(n_params))
    abc_matrix[, j] <- theta_true_acc[, j] - (theta_hat_acc[, j] - theta_hat_obs[j])

  if (verbose)
    cat(sprintf("PipeMaster:: ABC-NN regression done — %d posterior samples\n",
                nrow(abc_matrix)))

  abc_matrix
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
# Internal: prepare reftable split (features + targets from reftable rows)
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
# Internal: compute how many TF intra-op threads each worker should use
# ============================================================================

.compute.threads.per.worker <- function(cores, greedy = TRUE) {
  if (!greedy) return(1L)
  n_logical <- parallel::detectCores()
  n_physical <- tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(e) n_logical
  )
  if (is.na(n_physical)) n_physical <- n_logical
  if (is.na(n_physical)) return(1L)
  # detectCores(logical=FALSE) is unreliable on some systems (reports logical
  # count). Cap at half of logical cores as a conservative estimate.
  n_physical <- min(n_physical, ceiling(n_logical / 2))
  max(1L, floor(n_physical / cores))
}

# ============================================================================
# Internal: continuous worker pool with priority scheduling
# ============================================================================

.launch.rscript.pool <- function(tasks, cores, work_dir,
                                  gpus, verbose,
                                  max_retries = 2L) {
  old_wd <- getwd()
  setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)

  # Register parent PID so workers can detect if parent dies
  pid_file <- file.path(work_dir, ".PM_parent.pid")
  writeLines(as.character(Sys.getpid()), pid_file)

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

  # On exit: remove PID file + kill all workers by their PID files
  on.exit({
    tryCatch(file.remove(pid_file), error = function(e) NULL)
    pid_files <- list.files(sentinels_dir, pattern = "\\.pid$", full.names = TRUE)
    for (pf in pid_files) {
      wpid <- tryCatch(as.integer(readLines(pf, warn = FALSE)[1]), error = function(e) NA)
      if (!is.na(wpid)) {
        if (.Platform$OS.type == "unix") {
          system(sprintf("kill -9 -%d 2>/dev/null; kill -9 %d 2>/dev/null", wpid, wpid),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
        } else {
          system(sprintf("taskkill /F /T /PID %d 2>NUL", wpid),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
        }
      }
    }
  }, add = TRUE)

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

    # GPU environment: per-task override or pool-level round-robin
    if (!is.null(task$gpu_env)) {
      gpu_env <- task$gpu_env
    } else if (gpus > 0) {
      gpu_id <- gpu_counter %% gpus
      gpu_counter <<- gpu_counter + 1L
      gpu_env <- sprintf("CUDA_VISIBLE_DEVICES=%d TF_FORCE_GPU_ALLOW_GROWTH=true",
                         gpu_id)
    } else {
      gpu_env <- "CUDA_VISIBLE_DEVICES=-1"
    }

    # Shell wrapper: run Rscript, capture exit code in .done file, record PID
    # setsid creates a new process group so we can kill Rscript + children
    # echo $? stderr redirected to /dev/null to avoid "Directory nonexistent" on cleanup
    cmd <- sprintf(
      "{ setsid env %s Rscript %s %d > %s 2>&1; echo $? > %s 2>/dev/null; } & echo $! > %s",
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

  # Helper: extract one-line progress status from a worker's log file
  .worker.progress <- function(log_rel_path) {
    log_path <- file.path(work_dir, log_rel_path)
    if (!file.exists(log_path)) return("(starting...)")
    lines <- readLines(log_path, warn = FALSE)
    if (length(lines) == 0L) return("(starting...)")

    # Extract total brackets (s_max) from header line
    n_brackets <- NA_integer_
    smax_idx <- grep("s_max=", lines, fixed = TRUE)
    if (length(smax_idx) > 0L) {
      m <- regmatches(lines[smax_idx[1L]],
                      regexpr("s_max=(\\d+)", lines[smax_idx[1L]]))
      if (length(m) == 1L)
        n_brackets <- as.integer(sub("s_max=", "", m)) + 1L
    }

    # Extract latest best val_loss from Round lines (scan from bottom)
    best_loss <- NA_real_
    for (j in seq.int(length(lines), 1L)) {
      if (grepl("best val_loss=", lines[j], fixed = TRUE)) {
        m <- regmatches(lines[j],
                        regexpr("best val_loss=([0-9.]+)", lines[j]))
        if (length(m) == 1L)
          best_loss <- as.numeric(sub("best val_loss=", "", m))
        break
      }
    }

    # Scan from bottom for status (highest priority first)
    for (j in seq.int(length(lines), 1L)) {
      ln <- lines[j]

      # 1. Search done
      if (grepl("^Search .* done", ln)) {
        m <- regmatches(ln, regexpr("val_loss=([0-9.]+)", ln))
        vl <- if (length(m) == 1L) sub("val_loss=", "", m) else "?"
        return(sprintf("done (val_loss=%s)", vl))
      }

      # 2. Retraining phase
      if (grepl("retraining from epoch", ln, fixed = TRUE))
        return("retraining best model...")

      # 3. Bracket done with ETA
      if (grepl("Bracket.*done.*elapsed", ln)) {
        m_bracket <- regmatches(ln, regexpr("Bracket (\\d+)", ln))
        bnum <- if (length(m_bracket) == 1L)
          as.integer(sub("Bracket ", "", m_bracket)) else NA_integer_

        m_eta <- regmatches(ln, regexpr("~(\\d+) min remaining", ln))
        eta_str <- if (length(m_eta) == 1L) {
          paste0("~", sub(" min remaining", "m left", m_eta))
        } else NULL

        bracket_str <- if (!is.na(bnum) && !is.na(n_brackets)) {
          # brackets count down from s_max to 0, so done count = s_max - bnum + 1
          # but the last bracket done is bnum, meaning brackets s_max..bnum are done
          done_count <- n_brackets - bnum
          sprintf("bracket %d/%d", done_count, n_brackets)
        } else if (!is.na(bnum)) {
          sprintf("after bracket %d", bnum)
        } else "between brackets"

        parts <- bracket_str
        if (!is.na(best_loss)) parts <- paste0(parts, sprintf(", best=%.4f", best_loss))
        if (!is.null(eta_str)) parts <- paste0(parts, ", ", eta_str)
        return(parts)
      }

      # 4. Currently inside a bracket (not done yet)
      if (grepl("Bracket.*\\|.*configs", ln)) {
        m_bracket <- regmatches(ln, regexpr("Bracket (\\d+)", ln))
        bnum <- if (length(m_bracket) == 1L)
          as.integer(sub("Bracket ", "", m_bracket)) else NA_integer_

        bracket_str <- if (!is.na(bnum) && !is.na(n_brackets)) {
          in_count <- n_brackets - bnum - 1L
          sprintf("in bracket %d/%d", in_count + 1L, n_brackets)
        } else if (!is.na(bnum)) {
          sprintf("in bracket %d", bnum)
        } else "running bracket"

        if (!is.na(best_loss))
          bracket_str <- paste0(bracket_str, sprintf(", best=%.4f", best_loss))
        return(bracket_str)
      }
    }
    "(starting...)"
  }

  # Fill initial slots (conformal tasks come first = priority)
  while (length(active) < cores && length(pending) > 0) {
    tidx <- pending[1]
    pending <- pending[-1]
    active[[length(active) + 1]] <- launch_one(tidx)
  }

  if (verbose) {
    cat(sprintf("  [pool] Launched %d / %d tasks (%d pending)\n",
                length(active), n_tasks, length(pending)))
    flush.console()
  }

  # Poll loop — check every second, 4-way detection
  last_progress <- Sys.time()
  while (length(active) > 0) {
    Sys.sleep(1)

    still_active <- list()
    for (a in active) {
      task <- tasks[[a$task_idx]]
      elapsed_sec <- as.numeric(difftime(Sys.time(), a$start_time, units = "secs"))

      if (file.exists(task$result)) {
        # --- SUCCESS: result CSV exists ---
        completed <- c(completed, a$task_idx)
        if (verbose) {
          elapsed_min <- elapsed_sec / 60
          cat(sprintf("  [pool] %s %d done (%d/%d, %.1f min)\n",
                      task$prefix, task$id,
                      length(completed), n_tasks, elapsed_min))
          flush.console()
        }

      } else if (file.exists(a$done_file)) {
        # --- CRASH: process exited but no result CSV ---
        exit_code <- suppressWarnings(tryCatch({
          raw <- readLines(a$done_file, n = 1L, warn = FALSE)
          if (length(raw) == 0L || nchar(trimws(raw)) == 0L) NA_integer_
          else as.integer(trimws(raw))
        }, error = function(e) NA_integer_))
        tail_lines <- .log_tail(a$log_file)

        if (retry_count[a$task_idx] < max_retries) {
          retry_count[a$task_idx] <- retry_count[a$task_idx] + 1L
          if (verbose) {
            cat(sprintf("  [pool] %s %d CRASHED (exit %s, %.0fs) — retry %d/%d\n",
                        task$prefix, task$id,
                        if (is.na(exit_code)) "?" else as.character(exit_code),
                        elapsed_sec,
                        retry_count[a$task_idx], max_retries))
            flush.console()
          }
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

      } else {
        # --- RUNNING ---
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

    # Periodic progress report from worker logs
    if (verbose && length(active) > 0 &&
        as.numeric(difftime(Sys.time(), last_progress, units = "secs")) >= 60) {
      last_progress <- Sys.time()
      cat(sprintf("  [pool] Progress (%s):\n",
                  format(Sys.time(), "%H:%M:%S")))
      for (a in active) {
        task <- tasks[[a$task_idx]]
        elapsed_min <- as.numeric(difftime(Sys.time(), a$start_time,
                                           units = "mins"))
        status <- .worker.progress(a$log_file)
        cat(sprintf("    %s %d: %.0fm \u2014 %s\n",
                    task$prefix, task$id, elapsed_min, status))
      }
      flush.console()
    }
  }

  list(completed = completed, failed = failed)
}
# ============================================================================
# Torch backend for tune.nn
#
# Dispatched from tune.nn() when backend = "torch". Uses .torch.hyperband()
# and .torch.train.model() from torch_training.R, and .build.nn.torch()
# from torch_modules.R.
# ============================================================================

#' @keywords internal
.tune.nn.torch <- function(reftable, param.cols, type, sfs.dims,
                           max_epochs, eta, search_space, exclude.cols,
                           val.frac, n_searches = 1L, cores = 1L,
                           gpus = 0L, gpu.threshold = 4L, greedy = TRUE,
                           top_k, seed, verbose) {

  if (!requireNamespace("torch", quietly = TRUE))
    stop("tune.nn(backend='torch') requires the 'torch' R package.\n",
         "Install with: install.packages('torch')\n",
         "Then run: torch::install_torch()")

  type <- match.arg(type, c("sumstat", "sfs1d", "sfs2d", "emulator"))

  # Note: torch backend currently supports sumstat/emulator (ResNet) only
  if (type %in% c("sfs1d", "sfs2d"))
    stop("Torch backend does not yet support CNN architectures (sfs1d/sfs2d). ",
         "Use backend='keras' for SFS-based models.")

  top_k <- as.integer(top_k)
  if (top_k < 1L) stop("top_k must be >= 1")
  n_searches <- as.integer(n_searches)
  n_concurrent <- min(as.integer(cores), n_searches)

  # --- Prepare data (reuse existing prep) ---
  if (verbose) cat(sprintf("PipeMaster:: tune.nn \u2014 Hyperband (ResNet)\n"))
  data <- .prep.data(reftable, param.cols, type, sfs.dims, exclude.cols,
                     val.frac, seed)
  reftable <- NULL  # release reference, data has everything we need

  n_feat <- ncol(data$X_train)
  n_targ <- ncol(data$Y_train)

  if (verbose) cat(sprintf("PipeMaster:: %d features, %d targets | %d train, %d val\n",
                           n_feat, n_targ, nrow(data$X_train), nrow(data$X_val)))

  # --- Search space ---
  ss <- if (!is.null(search_space)) search_space else .emulator.search.space()

  # ==========================================================================
  # Sequential mode (single search, no workers)
  # ==========================================================================
  if (n_searches <= 1L || n_concurrent <= 1L) {

    device <- "cpu"
    if (torch::cuda_is_available()) {
      device <- "cuda"
      if (verbose) cat("PipeMaster:: Using CUDA GPU for torch training\n")
    }

    torch::torch_manual_seed(as.integer(seed))
    hb <- .torch.hyperband(ss, data, type, sfs.dims = sfs.dims,
                           max_epochs = as.integer(max_epochs),
                           eta = as.integer(eta),
                           seed = as.integer(seed),
                           verbose = verbose,
                           device = device)

    if (verbose)
      cat(sprintf("\nPipeMaster:: Best config: %s\n", .hp.to.string(hb$best_hp, type)))

    torch::torch_manual_seed(as.integer(seed))
    best_model <- .build.nn.torch(hb$best_hp, data, type)
    best_model$to(device = torch::torch_device(device))

    if (!is.null(hb$best_state))
      best_model$load_state_dict(hb$best_state)

    start_ep <- if (!is.null(hb$best_state)) hb$best_epochs else 0L
    retrain <- .torch.train.model(
      best_model, data$X_train, data$Y_train, data$X_val, data$Y_val,
      hp = hb$best_hp, type = type,
      epochs = as.integer(max_epochs), initial_epoch = start_ep,
      patience = 30L, lr_patience = 15L,
      device = device
    )

    if (verbose)
      cat(sprintf("PipeMaster:: Retrained: val_loss=%.6f, epochs=%d\n",
                  retrain$val_loss, retrain$epochs_trained))

    metrics <- .torch.compute.model.metrics(best_model, data, type, device = device)

    if (verbose) {
      cat(sprintf("  R\u00b2:  %s\n", paste(sprintf("%s=%.4f", names(metrics$r2), metrics$r2), collapse = "  ")))
      cat(sprintf("  MPE: %s\n", paste(sprintf("%s=%.1f%%", names(metrics$mpe), metrics$mpe), collapse = "  ")))
    }

    return(list(
      best_hp         = hb$best_hp,
      best_val_loss   = retrain$val_loss,
      all_results     = hb$all_results,
      best_model      = best_model,
      models          = list(best_model),
      models_hp       = list(hb$best_hp),
      models_val_loss = retrain$val_loss,
      models_metrics  = list(metrics),
      data            = data,
      type            = type,
      sfs.dims        = sfs.dims,
      exclude.cols    = exclude.cols,
      backend         = "torch"
    ))
  }

  # ==========================================================================
  # Parallel mode: save data as torch tensors, launch Rscript workers
  # ==========================================================================
  work_dir <- tempfile("hb_torch_search_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)

  # Save data for workers. Two strategies:
  # (a) bigmemory: file-backed matrices shared via mmap — ~0 RSS per worker
  # (b) fallback:  .rds files — each worker loads a full copy
  use_bigmemory <- requireNamespace("bigmemory", quietly = TRUE)

  if (use_bigmemory) {
    .save_bigmatrix <- function(mat, name) {
      bm <- bigmemory::as.big.matrix(mat, type = "double",
              backingfile = paste0(name, ".bin"),
              descriptorfile = paste0(name, ".desc"),
              backingpath = work_dir)
      bigmemory::flush(bm)
      bm
    }
    .save_bigmatrix(data$X_train, "X_train")
    .save_bigmatrix(data$Y_train, "Y_train")
    .save_bigmatrix(data$X_val,   "X_val")
    .save_bigmatrix(data$Y_val,   "Y_val")

    if (verbose) {
      bin_mb <- sum(file.size(file.path(work_dir,
        c("X_train.bin", "Y_train.bin", "X_val.bin", "Y_val.bin")))) / 1e6
      cat(sprintf("PipeMaster:: Saved shared data via bigmemory (%.1f MB, mmap)\n", bin_mb))
    }
  } else {
    saveRDS(data$X_train, file.path(work_dir, "X_train.rds"), compress = FALSE)
    saveRDS(data$Y_train, file.path(work_dir, "Y_train.rds"), compress = FALSE)
    saveRDS(data$X_val,   file.path(work_dir, "X_val.rds"),   compress = FALSE)
    saveRDS(data$Y_val,   file.path(work_dir, "Y_val.rds"),   compress = FALSE)

    if (verbose) {
      rds_mb <- sum(file.size(file.path(work_dir,
        c("X_train.rds", "Y_train.rds", "X_val.rds", "Y_val.rds")))) / 1e6
      cat(sprintf("PipeMaster:: Saved data files (%.1f MB, per-worker copies)\n", rds_mb))
    }
  }

  # Save metadata only (small .RData — no data matrices)
  n_features <- n_feat
  n_targets  <- n_targ
  n_bins     <- data$n_bins

  cfgs <- if (length(max_epochs) == 1L && length(eta) == 1L) {
    .generate.search.configs(n_searches, max_epochs, eta)
  } else {
    data.frame(eta = rep_len(as.integer(eta), n_searches),
               max_epochs = rep_len(as.integer(max_epochs), n_searches))
  }
  search_max_epochs  <- cfgs$max_epochs
  search_eta         <- cfgs$eta
  search_seeds       <- seed + (seq_len(n_searches) - 1L) * 10000L
  saved_lib_paths    <- .libPaths()
  threads_per_worker <- .compute.threads.per.worker(n_concurrent, greedy)
  search_space_saved <- ss

  pkg_source_dir <- ""
  if ("devtools_shims" %in% search()) {
    pkg_path <- find.package("PipeMaster", quiet = TRUE)
    if (file.exists(file.path(pkg_path, "DESCRIPTION"))) pkg_source_dir <- pkg_path
  }

  save(type, sfs.dims, n_features, n_targets, n_bins,
       use_bigmemory,
       search_seeds, saved_lib_paths, threads_per_worker,
       search_max_epochs, search_eta, search_space_saved,
       pkg_source_dir,
       file = file.path(work_dir, "shared_search_meta.RData"))

  # Free parent memory: save data for later metrics, drop heavy objects
  data_rds <- file.path(work_dir, "prep_data.rds")
  saveRDS(data, data_rds, compress = FALSE)
  rm(data, reftable); gc()

  # Write torch worker script
  .torch.write.search.worker.script(file.path(work_dir, "_search_worker.R"))

  if (verbose) {
    cat(sprintf("PipeMaster:: Launching %d searches (%d concurrent)\n",
                n_searches, n_concurrent))
    for (k in seq_len(n_searches))
      cat(sprintf("  Search %d: eta=%d, max_epochs=%d\n",
                  k, search_eta[k], search_max_epochs[k]))
  }

  # Build task list
  n_gpu_searches <- if (gpus > 0L) min(n_searches, gpu.threshold * gpus) else 0L
  tasks <- lapply(seq_len(n_searches), function(k) {
    if (gpus > 0L && k <= n_gpu_searches) {
      gpu_id <- (k - 1L) %% gpus
      task_gpu_env <- sprintf("CUDA_VISIBLE_DEVICES=%d", gpu_id)
    } else {
      task_gpu_env <- "CUDA_VISIBLE_DEVICES=-1"
    }
    list(
      script  = "_search_worker.R",
      id      = k,
      result  = sprintf("results/search_%04d/done.txt", k),
      prefix  = "search",
      gpu_env = task_gpu_env
    )
  })

  pool_result <- .launch.rscript.pool(
    tasks, n_concurrent, work_dir,
    gpus = 0L,
    verbose = verbose, max_retries = 1L)

  # Collect results and rank by val_loss
  all_search_results <- list()
  search_entries <- list()

  for (k in seq_len(n_searches)) {
    search_dir <- file.path(results_dir, sprintf("search_%04d", k))
    rds_file   <- file.path(search_dir, "result.rds")

    if (file.exists(rds_file)) {
      res <- tryCatch(readRDS(rds_file), error = function(e) NULL)
      if (!is.null(res)) {
        search_res <- res$all_results
        search_res$search <- k
        all_search_results[[length(all_search_results) + 1L]] <- search_res

        if (verbose)
          cat(sprintf("  Search %d: val_loss=%.6f\n", k, res$best_val_loss))

        search_entries[[length(search_entries) + 1L]] <- list(
          val_loss   = res$best_val_loss,
          hp         = res$best_hp,
          search_dir = search_dir
        )
      }
    } else {
      if (verbose) cat(sprintf("  Search %d: FAILED (no result)\n", k))
    }
  }

  if (length(search_entries) == 0L)
    stop("All parallel searches failed. Check worker logs in ", work_dir)

  # Sort by val_loss and take top-K
  vl_vec <- vapply(search_entries, `[[`, numeric(1), "val_loss")
  ord <- order(vl_vec)
  n_keep <- min(top_k, length(search_entries))
  search_entries <- search_entries[ord[seq_len(n_keep)]]

  # Load top-K torch models
  models <- lapply(search_entries, function(entry) {
    tryCatch(
      .torch.load.model(file.path(entry$search_dir, "best_model.pt")),
      error = function(e) {
        warning("Could not load torch model (val_loss=",
                sprintf("%.6f", entry$val_loss), "): ",
                conditionMessage(e), call. = FALSE)
        NULL
      }
    )
  })
  models_val_loss <- vapply(search_entries, `[[`, numeric(1), "val_loss")

  keep <- !vapply(models, is.null, logical(1))
  models <- models[keep]
  models_val_loss <- models_val_loss[keep]
  # Per-model HP: each parallel search produces its own architecture via
  # Hyperband, so top-K models can have different unit counts / block counts.
  # .torch.save.tune.result() needs these to reconstruct the correct shapes
  # on reload — otherwise state_dict load fails with tensor-size mismatch.
  models_hp <- lapply(search_entries[keep], `[[`, "hp")

  if (length(models) == 0L)
    stop("All top-K models failed to load. Check worker logs in ", work_dir)

  best_hp <- search_entries[[1L]]$hp

  combined_results <- if (length(all_search_results) > 0)
    do.call(rbind, all_search_results) else data.frame()

  # Reload data for metrics computation (was freed before worker launch)
  data <- readRDS(data_rds)

  # Compute R² and MPE for each model
  models_metrics <- lapply(models, .torch.compute.model.metrics,
                           data = data, type = type, device = "cpu")

  if (verbose) {
    cat(sprintf("\nPipeMaster:: Best across %d searches: val_loss=%.6f\n",
                n_searches, models_val_loss[1L]))
    cat(sprintf("PipeMaster:: Kept top %d models (val_loss: %s)\n",
                length(models),
                paste(sprintf("%.6f", models_val_loss), collapse = ", ")))
    cat(sprintf("PipeMaster:: Best config: %s\n", .hp.to.string(best_hp, type)))
    .print.model.metrics(models_metrics, models_val_loss, verbose)
  }

  Sys.sleep(1)
  unlink(work_dir, recursive = TRUE)

  list(
    best_hp         = best_hp,
    best_val_loss   = models_val_loss[1L],
    all_results     = combined_results,
    best_model      = models[[1L]],
    models          = models,
    models_hp       = models_hp,
    models_val_loss = models_val_loss,
    models_metrics  = models_metrics,
    data            = data,
    type            = type,
    sfs.dims        = sfs.dims,
    exclude.cols    = exclude.cols,
    backend         = "torch"
  )
}

# ============================================================================
# Unified model predict dispatcher
#
# Detects keras vs torch models and dispatches to the right predict function.
# Used by nn.predict(), OOD.posttrain(), .ensemble.predict(), etc.
# ============================================================================

#' @keywords internal
.predict.nn <- function(model, X) {
  if (inherits(model, "nn_module")) {
    .torch.predict(model, X)
  } else {
    predict(model, X, verbose = 0L)
  }
}

# ============================================================================
# Detect backend from tune result or model object
# ============================================================================

#' @keywords internal
.detect.backend <- function(tune.result) {
  # Explicit backend field (new results)
  if (!is.null(tune.result$backend)) return(tune.result$backend)
  # Detect from model type
  if (inherits(tune.result$best_model, "nn_module")) return("torch")
  "keras"
}


# ============================================================================
# Sensitivity / Inverse Jacobian for tune.nn (regression NN: stats -> params)
#
# Mirror of .torch.compute.jacobian() in torch_emulator.R, which is for the
# forward emulator (params -> stats). Here the chain rule runs the other way:
#
#   x_raw ── augment ──> x_aug = [x_raw, log1p|x_raw|]
#         ── z-score ──> x_z   = (x_aug - feat_mu) / feat_sd
#         ── NN ────────> y_z
#         ── un-z + exp> y_orig = exp(y_z * target_sd + target_mu)
#
# d y_orig_j / d x_raw_i =
#   y_orig_j * target_sd_j *
#   [ J_z[j, i_raw] / feat_sd_raw[i]
#   + J_z[j, i_log] * sign(x_raw_i)/(1 + |x_raw_i|) / feat_sd_log[i] ]
# ============================================================================

#' @keywords internal
.torch.compute.jacobian.inverse <- function(model, data, observed_raw,
                                            device = "cpu") {
  dev <- torch::torch_device(device)
  feat_mu <- data$feat_mu;    feat_sd <- data$feat_sd
  target_mu <- data$target_mu; target_sd <- data$target_sd

  n_stats <- length(observed_raw)
  obs_aug <- c(observed_raw, log1p(abs(observed_raw)))
  obs_z   <- (obs_aug - feat_mu) / feat_sd

  x <- torch::torch_tensor(matrix(obs_z, nrow = 1L),
                           dtype = torch::torch_float(),
                           device = dev, requires_grad = TRUE)
  model$eval()
  pred_z <- model(x)
  n_targets <- pred_z$size(2)
  n_aug     <- length(obs_z)

  # Per-target backward via autograd_grad (same pattern as forward jacobian)
  J_z <- matrix(0, n_targets, n_aug)
  for (j in seq_len(n_targets)) {
    grad_out <- torch::torch_zeros_like(pred_z)
    grad_out[1, j] <- 1.0
    g <- torch::autograd_grad(outputs = pred_z, inputs = x,
                              grad_outputs = grad_out,
                              retain_graph = TRUE, create_graph = FALSE)[[1]]
    J_z[j, ] <- as.numeric(g$cpu())
  }

  pred_z_np <- as.numeric(as.matrix(pred_z$detach()$cpu()))
  pred_log  <- pred_z_np * target_sd + target_mu
  pred_orig <- exp(pred_log)

  # Decompose into raw and log-augmented partial derivatives
  J_raw <- J_z[, seq_len(n_stats), drop = FALSE]
  J_log <- J_z[, n_stats + seq_len(n_stats), drop = FALSE]

  feat_sd_raw <- feat_sd[seq_len(n_stats)]
  feat_sd_log <- feat_sd[n_stats + seq_len(n_stats)]

  # Chain rule for log1p(|x|): derivative is sign(x)/(1+|x|); 0 at x=0
  log_chain <- sign(observed_raw) / (1 + abs(observed_raw))
  log_chain[!is.finite(log_chain)] <- 0

  # dY_z / dx_raw_i
  J_z_x <- sweep(J_raw, 2, 1 / feat_sd_raw, "*") +
           sweep(J_log, 2, log_chain / feat_sd_log, "*")

  # Output chain: dY_orig / dY_z = pred_orig * target_sd
  out_scale <- pred_orig * target_sd
  J_orig <- sweep(J_z_x, 1, out_scale, "*")

  # Elasticity (dimensionless): E[j,i] = (x_raw_i / Y_j) * J_orig[j,i]
  obs_safe <- ifelse(abs(observed_raw) < 1e-12, 1e-12, observed_raw)
  ratio    <- outer(1 / pred_orig, obs_safe)
  E_orig   <- J_orig * ratio

  rownames(J_orig) <- rownames(E_orig) <- names(target_mu)
  colnames(J_orig) <- colnames(E_orig) <- data$stat_cols

  list(jacobian = J_orig, elasticity = E_orig, predictions = pred_orig)
}


#' Sensitivity Analysis for `tune.nn()` (Inverse Jacobian)
#'
#' Computes \eqn{J[j, i] = d\theta_j / dS_i} of the trained NN at one or more
#' observed datasets — i.e., how each parameter prediction responds to each
#' summary statistic. Mirror of \code{emulator.sensitivity()}, which computes
#' the opposite direction \eqn{dS_k / d\theta_j} for the forward emulator.
#'
#' Two outputs:
#' \itemize{
#'   \item \strong{Jacobian}: in original parameter and statistic units.
#'         Useful when stats and params share interpretable scales.
#'   \item \strong{Elasticity}: dimensionless,
#'         \eqn{E[j, i] = (S_i / \theta_j) \cdot J[j, i]} — comparable across
#'         stats with very different magnitudes.
#' }
#'
#' Ensemble averaging follows the conventions of \code{nn.predict()}:
#' weights are inverse validation loss (default) or uniform.
#'
#' @param tune.result list — output from \code{tune.nn()}.
#' @param observed numeric vector or 1-row data.frame — observed summary stats.
#' @param ensemble_weights character — \code{"inv_val_loss"} (default) weights
#'   each ensemble member by the inverse of its validation loss; \code{"uniform"}
#'   gives equal weight.
#' @param aggregate character — \code{"mean"} (default), \code{"median"}, or
#'   \code{"none"} (return per-model array).
#' @param device character — \code{"cpu"} or \code{"cuda"} (default chooses cuda
#'   when available).
#' @param verbose logical — print progress (default TRUE).
#'
#' @return A list of class \code{"tune_nn_sensitivity"} with:
#' \describe{
#'   \item{jacobian}{matrix (n_params x n_stats) when aggregated, or array
#'     (n_params x n_stats x K) when \code{aggregate = "none"}}
#'   \item{elasticity}{same shape as \code{jacobian}, dimensionless}
#'   \item{predictions}{ensemble parameter predictions at obs}
#'   \item{stat_cols, param_cols, observed, ensemble_weights, aggregate}{}
#' }
#'
#' @seealso \code{\link{emulator.sensitivity}} for the forward direction.
#' @export
tune.nn.sensitivity <- function(tune.result, observed,
                                ensemble_weights = c("inv_val_loss", "uniform"),
                                aggregate = c("mean", "median", "none"),
                                device = NULL, verbose = TRUE) {
  ensemble_weights <- match.arg(ensemble_weights)
  aggregate <- match.arg(aggregate)

  backend <- .detect.backend(tune.result)
  if (backend != "torch")
    stop("tune.nn.sensitivity currently supports torch backend only.")

  data <- tune.result$data
  stat_cols <- data$stat_cols
  if (is.data.frame(observed) || is.matrix(observed)) {
    obs_raw <- as.numeric(observed[1, stat_cols])
  } else if (!is.null(names(observed))) {
    obs_raw <- as.numeric(observed[stat_cols])
  } else {
    obs_raw <- as.numeric(observed)
    if (length(obs_raw) != length(stat_cols))
      stop(sprintf("observed has %d entries but tune.nn used %d stats",
                   length(obs_raw), length(stat_cols)))
  }
  if (any(!is.finite(obs_raw)))
    stop("observed contains non-finite values for the trained stat columns.")

  models <- tune.result$models %||% list(tune.result$best_model)
  K <- length(models)
  if (K == 0L) stop("tune.result has no models")

  val_loss <- tune.result$models_val_loss
  if (ensemble_weights == "inv_val_loss" && !is.null(val_loss) &&
      all(is.finite(val_loss)) && length(val_loss) == K) {
    w <- (1 / val_loss) / sum(1 / val_loss)
  } else {
    w <- rep(1 / K, K)
  }

  # Detect the actual device the loaded models live on; fall back to cpu.
  if (is.null(device)) {
    device <- tryCatch({
      d <- models[[1]]$parameters[[1]]$device
      if (grepl("cuda", as.character(d), fixed = TRUE)) "cuda" else "cpu"
    }, error = function(e) "cpu")
  }
  # Move models to the requested device if needed
  for (k in seq_len(K)) models[[k]]$to(device = torch::torch_device(device))

  if (verbose)
    cat(sprintf("PipeMaster:: tune.nn.sensitivity — %d ensemble models, weights=%s, aggregate=%s, device=%s\n",
                K, ensemble_weights, aggregate, device))

  J_list <- vector("list", K); E_list <- vector("list", K)
  P_list <- vector("list", K)
  for (k in seq_len(K)) {
    res <- .torch.compute.jacobian.inverse(models[[k]], data, obs_raw,
                                            device = device)
    J_list[[k]] <- res$jacobian
    E_list[[k]] <- res$elasticity
    P_list[[k]] <- res$predictions
  }

  if (aggregate == "none") {
    J_out <- simplify2array(J_list)
    E_out <- simplify2array(E_list)
    P_out <- do.call(rbind, P_list)
  } else {
    P_mat <- do.call(rbind, P_list)
    if (aggregate == "median") {
      J_out <- apply(simplify2array(J_list), c(1, 2), median)
      E_out <- apply(simplify2array(E_list), c(1, 2), median)
      P_out <- apply(P_mat, 2, median)
    } else {
      J_out <- Reduce("+", Map("*", J_list, w))
      E_out <- Reduce("+", Map("*", E_list, w))
      P_out <- as.numeric(crossprod(w, P_mat))
    }
    rownames(J_out) <- rownames(E_out) <- names(data$target_mu)
    colnames(J_out) <- colnames(E_out) <- stat_cols
    names(P_out)    <- names(data$target_mu)
  }

  result <- list(jacobian = J_out, elasticity = E_out,
                 predictions = P_out,
                 stat_cols = stat_cols,
                 param_cols = names(data$target_mu),
                 observed = obs_raw,
                 ensemble_weights = w,
                 aggregate = aggregate)
  class(result) <- "tune_nn_sensitivity"

  if (verbose && aggregate != "none") {
    cat(sprintf("PipeMaster:: |dtheta/dS| range: [%.4g, %.4g]\n",
                min(abs(result$jacobian)), max(abs(result$jacobian))))
    cat(sprintf("PipeMaster:: |elasticity| range: [%.4g, %.4g]\n",
                min(abs(result$elasticity)), max(abs(result$elasticity))))
  }
  result
}

# Local "or-else" helper; tune.nn.R doesn't define this elsewhere
`%||%` <- function(a, b) if (is.null(a)) b else a


#' @export
print.tune_nn_sensitivity <- function(x, top = 10L, ...) {
  K <- length(x$param_cols)
  cat(sprintf("tune.nn.sensitivity — %d params, %d stats, aggregate=%s\n",
              K, length(x$stat_cols), x$aggregate))
  if (x$aggregate == "none") {
    cat("Per-model Jacobians stored in $jacobian (3D array). ",
        "Use aggregate='mean' for ensemble summary.\n", sep = "")
    return(invisible(x))
  }
  J <- x$jacobian; E <- x$elasticity
  for (j in seq_len(K)) {
    ord <- order(-abs(J[j, ]))[seq_len(min(top, ncol(J)))]
    cat(sprintf("\n%s — top %d stats:\n", x$param_cols[j], length(ord)))
    for (i in ord)
      cat(sprintf("  %-30s J=%+11.4g  E=%+8.4f\n",
                  x$stat_cols[i], J[j, i], E[j, i]))
  }
  invisible(x)
}


#' @export
summary.tune_nn_sensitivity <- function(object, top = 5L, ...) {
  if (object$aggregate == "none") {
    cat("Use aggregate='mean'/'median' for summary.\n"); return(invisible(NULL))
  }
  J <- object$jacobian; E <- object$elasticity
  out_list <- list()
  for (j in seq_along(object$param_cols)) {
    ord <- order(-abs(J[j, ]))[seq_len(min(top, ncol(J)))]
    out_list[[object$param_cols[j]]] <- data.frame(
      stat       = object$stat_cols[ord],
      jacobian   = J[j, ord],
      elasticity = E[j, ord],
      stringsAsFactors = FALSE
    )
  }
  print(out_list)
  invisible(out_list)
}


#' @export
plot.tune_nn_sensitivity <- function(x, top = 8L,
                                     use = c("jacobian", "elasticity"), ...) {
  use <- match.arg(use)
  if (x$aggregate == "none")
    stop("aggregate must be mean/median for plotting.")
  M <- if (use == "elasticity") x$elasticity else x$jacobian
  K <- length(x$param_cols)
  old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
  par(mfrow = c(K, 1), mar = c(3, 9, 2, 1))
  for (j in seq_len(K)) {
    ord  <- order(-abs(M[j, ]))[seq_len(min(top, ncol(M)))]
    vals <- M[j, ord]
    cols <- ifelse(vals > 0, "#1f77b4", "#d62728")
    barplot(rev(vals), names.arg = rev(x$stat_cols[ord]),
            horiz = TRUE, las = 1, col = rev(cols),
            main = sprintf("%s - top %d stats by |%s|",
                           x$param_cols[j], length(ord), use),
            xlab = sprintf("%s (signed)", use), cex.names = 0.8)
    abline(v = 0, col = "grey50")
  }
  invisible(NULL)
}
