# ============================================================================
# Torch Neural Network Classifier for PipeMaster Model Selection
#
# Multi-class classification on summary statistics. Reuses the ResNet backbone
# from torch_modules.R; replaces the regression head with K-class logits and
# the Huber loss with cross-entropy.
#
# Public API:
#   tune.nn.classify(reftable, model_col, ...)   — train classifier
#   nn.predict.classify(classifier, observed, ...) — predict class probs on obs
#
# Internal helpers (.torch.* convention):
#   .prep.classify.sumstat()   — feature/label/split prep for sumstat backbone
#   .torch.classify.train()    — training loop with CE loss + accuracy
#   .torch.classify.hyperband()— Hyperband HP search using CE loss
#   .torch.classify.metrics()  — accuracy + per-class precision/recall on val
#   .torch.classify.predict()  — forward pass → softmax probabilities
#
# S3 class "model_selection" with plot/summary/print methods.
# ============================================================================

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

#' @keywords internal
.softmax <- function(logits) {
  if (is.matrix(logits)) {
    m <- apply(logits, 1, max)
    e <- exp(logits - m)
    e / rowSums(e)
  } else {
    e <- exp(logits - max(logits))
    e / sum(e)
  }
}

#' @keywords internal
.classify.feature.cols <- function(reftable, model_col, exclude.cols = NULL) {
  param_pattern <- "^(Ne[0-9]*\\.|Ne\\.anc|t\\.Ne|alpha|join[0-9]+|mig[0-9]+|mean\\.rate|sd\\.rate|mean\\.rec\\.rate|sd\\.rec\\.rate|rec\\.rate)"
  param_cols <- grep(param_pattern, colnames(reftable), value = TRUE)
  setdiff(colnames(reftable), c(model_col, param_cols, exclude.cols))
}

# Stratified sample of indices (one per class, val.frac of each class in val)
#' @keywords internal
.stratified.split <- function(Y_int, val.frac, seed) {
  set.seed(seed)
  n <- length(Y_int)
  classes <- sort(unique(Y_int))
  tr <- integer(0); va <- integer(0)
  for (k in classes) {
    idx_k <- which(Y_int == k)
    idx_k <- sample(idx_k)
    n_v <- max(1L, floor(val.frac * length(idx_k)))
    va <- c(va, idx_k[seq_len(n_v)])
    tr <- c(tr, idx_k[(n_v + 1L):length(idx_k)])
  }
  list(tr = sample(tr), va = sample(va))
}

# ----------------------------------------------------------------------------
# Data prep (sumstat backbone — z-scored stats + log1p augmentation)
# ----------------------------------------------------------------------------

#' @keywords internal
.prep.classify.sumstat <- function(reftable, model_col, exclude.cols,
                                   val.frac, seed) {
  if (!model_col %in% colnames(reftable))
    stop(sprintf(".prep.classify.sumstat: model_col '%s' not in reftable",
                 model_col))

  raw_labels  <- as.character(reftable[[model_col]])
  class_names <- sort(unique(raw_labels))
  K <- length(class_names)
  if (K < 2L) stop(".prep.classify.sumstat: need >= 2 classes; got ", K)

  Y_int <- match(raw_labels, class_names)  # 1..K

  feat_cols <- .classify.feature.cols(reftable, model_col, exclude.cols)
  if (length(feat_cols) == 0L)
    stop(".prep.classify.sumstat: no feature columns left after exclusions")

  features_raw <- as.matrix(reftable[, feat_cols, drop = FALSE])
  features <- cbind(features_raw, log1p(abs(features_raw)))

  # Drop rows with non-finite features
  bad <- apply(features, 1, function(x) any(!is.finite(x)))
  if (any(bad)) {
    features <- features[!bad, , drop = FALSE]
    Y_int    <- Y_int[!bad]
  }

  split <- .stratified.split(Y_int, val.frac, seed)
  tr <- split$tr; va <- split$va

  feat_mu <- colMeans(features[tr, , drop = FALSE])
  feat_sd <- apply(features[tr, , drop = FALSE], 2, sd)
  feat_sd[feat_sd == 0] <- 1
  zscore <- function(X) t((t(X) - feat_mu) / feat_sd)

  X_train <- zscore(features[tr, , drop = FALSE])
  X_val   <- zscore(features[va, , drop = FALSE])
  Y_train <- Y_int[tr]
  Y_val   <- Y_int[va]

  list(
    X_train = X_train, X_val = X_val,
    Y_train = Y_train, Y_val = Y_val,
    n_features  = ncol(X_train),
    n_classes   = K,
    class_names = class_names,
    feat_mu = feat_mu, feat_sd = feat_sd,
    stat_cols = feat_cols
  )
}

# ----------------------------------------------------------------------------
# Training loop (cross-entropy loss; tracks accuracy on validation)
# ----------------------------------------------------------------------------

#' @keywords internal
.torch.classify.train <- function(model, X_train, Y_train, X_val, Y_val,
                                  hp, epochs = 500L, initial_epoch = 0L,
                                  patience = 30L, lr_patience = 15L,
                                  lr_factor = 0.5, min_lr = 1e-6,
                                  device = "cpu") {
  dev <- torch::torch_device(device)

  X_tr <- torch::torch_tensor(as.matrix(X_train), dtype = torch::torch_float(), device = dev)
  Y_tr <- torch::torch_tensor(as.integer(Y_train), dtype = torch::torch_long(), device = dev)
  X_va <- torch::torch_tensor(as.matrix(X_val),   dtype = torch::torch_float(), device = dev)
  Y_va <- torch::torch_tensor(as.integer(Y_val),  dtype = torch::torch_long(), device = dev)

  n_train    <- X_tr$size(1)
  batch_size <- as.integer(hp$batch_size)
  l2_reg     <- if (!is.null(hp$l2_reg)) hp$l2_reg else 0

  optimizer <- torch::optim_adam(model$parameters, lr = hp$learning_rate,
                                 weight_decay = 2 * l2_reg)
  scheduler <- torch::lr_reduce_on_plateau(optimizer, mode = "min",
                                           factor = lr_factor,
                                           patience = as.integer(lr_patience),
                                           min_lr = min_lr)

  best_val_loss <- Inf
  best_state    <- NULL
  best_epoch    <- initial_epoch
  wait <- 0L

  for (epoch in seq(initial_epoch + 1L, epochs)) {
    model$train()
    perm <- sample.int(n_train)
    n_batches <- ceiling(n_train / batch_size)

    for (b in seq_len(n_batches)) {
      s <- (b - 1L) * batch_size + 1L
      e <- min(b * batch_size, n_train)
      idx <- perm[s:e]
      X_b <- X_tr[idx]
      Y_b <- Y_tr[idx]

      optimizer$zero_grad()
      logits <- model(X_b)
      loss <- torch::nnf_cross_entropy(logits, Y_b)
      loss$backward()
      optimizer$step()
    }

    model$eval()
    val_loss <- torch::with_no_grad({
      logits_v <- model(X_va)
      torch::nnf_cross_entropy(logits_v, Y_va)$item()
    })

    scheduler$step(val_loss)

    if (val_loss < best_val_loss) {
      best_val_loss <- val_loss
      best_state <- lapply(model$state_dict(), function(t) t$clone())
      best_epoch <- epoch
      wait <- 0L
    } else {
      wait <- wait + 1L
      if (wait >= patience) break
    }
  }

  if (!is.null(best_state)) model$load_state_dict(best_state)
  list(val_loss = best_val_loss, epochs_trained = best_epoch)
}

# ----------------------------------------------------------------------------
# Hyperband (cross-entropy loss)
# ----------------------------------------------------------------------------

#' @keywords internal
.torch.classify.hyperband <- function(search_space, data,
                                      max_epochs, eta, seed, verbose,
                                      device = "cpu") {
  s_max <- min(floor(log(max_epochs) / log(eta)), 3L)

  if (verbose)
    cat(sprintf("PipeMaster:: max_epochs=%d, eta=%d, s_max=%d, %d brackets\n\n",
                max_epochs, eta, s_max, s_max + 1))

  all_results <- data.frame(hp_string = character(), val_loss = numeric(),
                            bracket = integer(), round = integer(),
                            stringsAsFactors = FALSE)
  global_best_loss <- Inf
  global_best_hp   <- NULL
  global_best_epochs <- 0L
  best_state <- NULL

  for (s in s_max:0) {
    n <- ceiling((s_max + 1) / (s + 1)) * as.integer(eta^s)
    r <- max_epochs / eta^s

    if (verbose) cat(sprintf("  Bracket %d | %d configs × %d epochs\n",
                             s, n, round(r)))

    set.seed(seed + s)
    configs <- lapply(seq_len(n), function(i) .sample.config(search_space))

    prev_epochs   <- rep(0L, n)
    config_states <- vector("list", n)

    for (i in 0:s) {
      r_i <- round(r * eta^i)
      n_i <- max(1, floor(n / eta^i))
      n_keep <- max(1, ceiling(n_i / eta))

      val_losses <- rep(Inf, length(configs))

      for (j in seq_along(configs)) {
        tryCatch({
          torch::torch_manual_seed(as.integer(seed + s + j))

          hp <- configs[[j]]
          model <- .PipeMasterResNet(data$n_features, data$n_classes, hp)
          model$to(device = torch::torch_device(device))

          if (prev_epochs[j] > 0 && !is.null(config_states[[j]]))
            model$load_state_dict(config_states[[j]])

          result <- .torch.classify.train(
            model, data$X_train, data$Y_train, data$X_val, data$Y_val,
            hp = hp, epochs = as.integer(r_i),
            initial_epoch = as.integer(prev_epochs[j]),
            patience = 10L, lr_patience = 5L,
            device = device
          )
          val_losses[j] <- result$val_loss
          config_states[[j]] <- lapply(model$state_dict(),
                                       function(t) t$clone()$cpu())

          if (val_losses[j] < global_best_loss) {
            global_best_loss   <- val_losses[j]
            global_best_hp     <- hp
            global_best_epochs <- as.integer(r_i)
            best_state <- lapply(model$state_dict(), function(t) t$clone()$cpu())
            if (verbose) cat(sprintf("    ★ new best: val_CE=%.4f (bracket %d, round %d)\n",
                                     val_losses[j], s, i))
          }

          rm(model); gc()
          if (device != "cpu") torch::cuda_empty_cache()
        }, error = function(e) {
          if (verbose) cat(sprintf("    [warn] config %d error: %s\n",
                                   j, conditionMessage(e)))
          val_losses[j] <<- Inf
        })
      }

      prev_epochs[seq_along(configs)] <- r_i

      for (j in seq_along(configs)) {
        if (is.finite(val_losses[j])) {
          all_results <- rbind(all_results, data.frame(
            hp_string = .hp.to.string(configs[[j]], "sumstat"),
            val_loss  = val_losses[j],
            bracket = s, round = i,
            stringsAsFactors = FALSE
          ))
        }
      }

      if (i < s) {
        ranking <- order(val_losses)
        keep <- ranking[seq_len(min(n_keep, length(ranking)))]
        if (verbose) cat(sprintf("    Round %d: best val_CE=%.4f | pruning to %d\n",
                                 i, min(val_losses[is.finite(val_losses)]),
                                 length(keep)))
        configs       <- configs[keep]
        prev_epochs   <- prev_epochs[keep]
        config_states <- config_states[keep]
      } else {
        if (verbose) cat(sprintf("    Round %d: best val_CE=%.4f\n",
                                 i, min(val_losses[is.finite(val_losses)])))
      }
    }
  }

  list(best_hp = global_best_hp,
       best_val_loss = global_best_loss,
       best_epochs = global_best_epochs,
       best_state = best_state,
       all_results = all_results)
}

# ----------------------------------------------------------------------------
# Validation metrics — accuracy + per-class precision/recall + confusion matrix
# ----------------------------------------------------------------------------

#' @keywords internal
.torch.classify.metrics <- function(model, data, device = "cpu") {
  logits <- .torch.predict(model, data$X_val, device = device)
  pred   <- max.col(logits, ties.method = "first")
  truth  <- data$Y_val
  K <- data$n_classes

  acc <- mean(pred == truth)
  conf <- table(factor(truth, levels = seq_len(K)),
                factor(pred,  levels = seq_len(K)))
  rownames(conf) <- paste0("true:", data$class_names)
  colnames(conf) <- paste0("pred:", data$class_names)

  precision <- numeric(K); recall <- numeric(K); f1 <- numeric(K)
  for (k in seq_len(K)) {
    tp <- conf[k, k]
    fp <- sum(conf[, k]) - tp
    fn <- sum(conf[k, ]) - tp
    precision[k] <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    recall[k]    <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    f1[k]        <- if (!is.na(precision[k]) && !is.na(recall[k]) &&
                        (precision[k] + recall[k]) > 0)
                      2 * precision[k] * recall[k] / (precision[k] + recall[k])
                    else NA_real_
  }
  names(precision) <- names(recall) <- names(f1) <- data$class_names

  list(accuracy = acc, confusion = conf,
       precision = precision, recall = recall, f1 = f1)
}

# ----------------------------------------------------------------------------
# Worker script for parallel classifier Hyperband (one Hyperband search per
# Rscript worker, distributed across GPUs via CUDA_VISIBLE_DEVICES). Mirrors
# .torch.write.search.worker.script() in torch_training.R for the regression
# path. Each worker reads the shared bigmemory/rds data, runs an independent
# .torch.classify.hyperband(), and writes its best state_dict + result.
# ----------------------------------------------------------------------------

#' @keywords internal
.torch.write.classify.search.worker.script <- function(filepath) {
  writeLines(c(
    '#!/usr/bin/env Rscript',
    'args <- commandArgs(trailingOnly = TRUE)',
    'task_id <- as.integer(args[1])',
    '',
    '# Load metadata (small — no data matrices)',
    'load("shared_search_meta.RData")',
    '',
    '# Threading env',
    'n_threads <- as.character(threads_per_worker)',
    'Sys.setenv(OMP_NUM_THREADS = n_threads,',
    '           MKL_NUM_THREADS = n_threads,',
    '           TORCH_NUM_THREADS = n_threads)',
    '',
    '# Restore library paths and load PipeMaster',
    '.libPaths(saved_lib_paths)',
    'if (nzchar(pkg_source_dir)) {',
    '  suppressPackageStartupMessages(devtools::load_all(pkg_source_dir, quiet = TRUE))',
    '} else {',
    '  suppressPackageStartupMessages(library(PipeMaster))',
    '}',
    '',
    'suppressPackageStartupMessages(library(torch))',
    'torch::torch_set_num_threads(as.integer(n_threads))',
    'torch::torch_set_num_interop_threads(1L)',
    '',
    '# GPU setup — CUDA_VISIBLE_DEVICES pre-set by parent in the worker env',
    'device <- "cpu"',
    'if (Sys.getenv("CUDA_VISIBLE_DEVICES", "-1") != "-1") {',
    '  if (cuda_is_available()) device <- "cuda"',
    '}',
    '',
    '# Load data — bigmemory (shared mmap) or rds (per-worker copy)',
    'if (use_bigmemory) {',
    '  suppressPackageStartupMessages(library(bigmemory))',
    '  X_train <- attach.big.matrix("X_train.desc")[,]',
    '  X_val   <- attach.big.matrix("X_val.desc")[,]',
    '  Y_train <- as.integer(attach.big.matrix("Y_train.desc")[, 1L])',
    '  Y_val   <- as.integer(attach.big.matrix("Y_val.desc")[, 1L])',
    '} else {',
    '  X_train <- readRDS("X_train.rds")',
    '  X_val   <- readRDS("X_val.rds")',
    '  Y_train <- as.integer(readRDS("Y_train.rds"))',
    '  Y_val   <- as.integer(readRDS("Y_val.rds"))',
    '}',
    '',
    '# Reconstruct hb_data shape expected by .torch.classify.hyperband',
    'hb_data <- list(',
    '  X_train = X_train, Y_train = Y_train,',
    '  X_val   = X_val,   Y_val   = Y_val,',
    '  feat_mu = feat_mu, feat_sd = feat_sd,',
    '  stat_cols = stat_cols,',
    '  class_names = class_names,',
    '  n_features = n_features, n_classes = n_classes',
    ')',
    '',
    'out_dir <- file.path("results", sprintf("search_%04d", task_id))',
    'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)',
    '',
    'worker_seed       <- search_seeds[task_id]',
    'worker_max_epochs <- search_max_epochs[task_id]',
    'worker_eta        <- search_eta[task_id]',
    'cat(sprintf("Search %d starting (seed=%d, eta=%d, max_epochs=%d, device=%s)\\n",',
    '            task_id, worker_seed, worker_eta, worker_max_epochs, device))',
    '',
    '# Parent-alive watchdog',
    'if (file.exists(".PM_parent.pid") && !PipeMaster:::.pm.parent.alive(".PM_parent.pid")) {',
    '  cat("Parent died, worker exiting.\\n"); q("no") }',
    '',
    'hb <- PipeMaster:::.torch.classify.hyperband(',
    '  search_space = search_space_saved,',
    '  data = hb_data,',
    '  max_epochs = worker_max_epochs, eta = worker_eta,',
    '  seed = worker_seed, verbose = TRUE, device = device)',
    '',
    'if (file.exists(".PM_parent.pid") && !PipeMaster:::.pm.parent.alive(".PM_parent.pid")) {',
    '  cat("Parent died, worker exiting.\\n"); q("no") }',
    '',
    '# Convert state dict (list of torch tensors) to plain numeric arrays.',
    '# Torch tensors are external pointers and do NOT survive saveRDS round-trip',
    '# (you get back stale pointers that segfault on use). Plain arrays do.',
    'best_state_arr <- lapply(hb$best_state, function(t) as.array(t$cpu()))',
    '',
    'saveRDS(list(',
    '  best_hp       = hb$best_hp,',
    '  best_val_loss = hb$best_val_loss,',
    '  all_results   = hb$all_results,',
    '  best_state    = best_state_arr',
    '), file.path(out_dir, "result.rds"))',
    '',
    'cat(sprintf("Search %d done (val_CE=%.6f)\\n", task_id, hb$best_val_loss))',
    'writeLines("done", file.path(out_dir, "done.txt"))'
  ), filepath)
}

# ----------------------------------------------------------------------------
# Public API: tune.nn.classify
# ----------------------------------------------------------------------------

#' Train a NN classifier for demographic-model selection
#'
#' Multi-class classification on summary statistics. Stacks reftables labelled
#' by `model_col` and learns to discriminate between models from sim stats.
#'
#' @param reftable Combined reftable: rows = sims, columns = stats + a label
#'   column (`model_col`). Demographic params are auto-detected and excluded
#'   from features.
#' @param model_col Column name with class labels (character or integer).
#' @param max_epochs Maximum training epochs per Hyperband config (default 500).
#' @param eta Hyperband halving factor (default 3).
#' @param search_space Named list of hyperparameter sampling fns. Default
#'   uses the same ResNet search space as the regression `tune.nn()`.
#' @param exclude.cols Additional column names to exclude from features.
#' @param val.frac Fraction of sims for validation (stratified per class).
#' @param n_searches Integer — number of independent Hyperband searches to run
#'   (default 1L). Mirrors `tune.nn(n_searches=)`: each search runs a full
#'   independent Hyperband with its own seed; the best across all searches is
#'   returned. With `cores > 1`, searches run in parallel via Rscript workers.
#' @param cores Integer — maximum number of concurrent search workers
#'   (default 1L). Ignored when `n_searches = 1`. Per-worker RAM scales with
#'   the data size; bigmemory (if available) shares the data via mmap.
#' @param gpus Integer — number of GPUs to distribute searches across
#'   (default 0L, CPU only). With `gpus > 0`, the first
#'   `min(n_searches, gpu.threshold*gpus)` workers are pinned to GPUs
#'   round-robin via `CUDA_VISIBLE_DEVICES`. With `n_searches = 1`, the single
#'   search uses CUDA directly if available.
#' @param gpu.threshold Integer — maximum concurrent searches per GPU
#'   (default 4L). Total GPU searches = `min(n_searches, gpu.threshold*gpus)`;
#'   excess searches run on CPU. Ignored when `gpus = 0`.
#' @param greedy Logical — thread allocation policy across concurrent workers
#'   (default TRUE). Passed through to `.compute.threads.per.worker()`.
#' @param backend Backend selector. Currently torch only (keras path TODO).
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A `model_selection_classifier` object: list with
#'   `model` (best torch nn_module), `data` (z-score params + class names),
#'   `class_names`, `metrics` (accuracy/confusion/F1), `best_hp`, `best_val_loss`,
#'   `all_results` (Hyperband trace), `feature.cols`, `n_features`, `n_classes`.
#'
#' @export
tune.nn.classify <- function(reftable, model_col,
                             max_epochs = 500L, eta = 3L,
                             search_space = NULL,
                             exclude.cols = NULL,
                             val.frac = 0.1,
                             n_searches = 1L, cores = 1L, gpus = 0L,
                             gpu.threshold = 4L, greedy = TRUE,
                             backend = c("torch", "keras"),
                             seed = 42, verbose = TRUE) {

  backend <- match.arg(backend)
  if (backend == "keras")
    stop("tune.nn.classify keras backend not yet implemented; use backend='torch'.")

  if (!requireNamespace("torch", quietly = TRUE) ||
      !torch::torch_is_installed())
    stop("torch backend requested but torch is not installed.\n",
         "Install with: install.packages('torch'); torch::install_torch()")

  n_searches   <- as.integer(n_searches)
  cores        <- as.integer(cores)
  gpus         <- as.integer(gpus)
  n_concurrent <- min(cores, n_searches)

  ss <- if (is.null(search_space)) .default.search.space("sumstat") else search_space

  if (verbose)
    cat(sprintf("PipeMaster:: tune.nn.classify — Hyperband ResNet (cross-entropy)\n"))

  data <- .prep.classify.sumstat(reftable, model_col, exclude.cols, val.frac, seed)
  reftable <- NULL

  if (verbose) {
    counts <- table(data$Y_train)
    cat(sprintf("PipeMaster:: %d features, %d classes, %d train, %d val\n",
                data$n_features, data$n_classes,
                length(data$Y_train), length(data$Y_val)))
    cat(sprintf("PipeMaster:: classes: %s\n",
                paste(sprintf("%s(n=%d)", data$class_names, counts), collapse = ", ")))
  }

  # ==========================================================================
  # Sequential mode (single search, no workers) — preserves old behavior
  # ==========================================================================
  if (n_searches <= 1L || n_concurrent <= 1L) {

    device <- if (torch::cuda_is_available()) "cuda" else "cpu"
    if (verbose) cat(sprintf("PipeMaster:: device: %s\n", device))

    hb <- .torch.classify.hyperband(
      search_space = ss, data = data,
      max_epochs = max_epochs, eta = eta,
      seed = seed, verbose = verbose,
      device = device
    )

    best_model <- .PipeMasterResNet(data$n_features, data$n_classes, hb$best_hp)
    best_model$to(device = torch::torch_device(device))
    best_model$load_state_dict(hb$best_state)

    metrics <- .torch.classify.metrics(best_model, data, device = device)

    if (verbose) {
      cat(sprintf("\nPipeMaster:: Best val_CE = %.4f, val accuracy = %.3f\n",
                  hb$best_val_loss, metrics$accuracy))
      cat(sprintf("PipeMaster:: Confusion matrix:\n"))
      print(metrics$confusion)
    }

    result <- list(
      model         = best_model,
      data          = data,
      class_names   = data$class_names,
      n_features    = data$n_features,
      n_classes     = data$n_classes,
      feature.cols  = data$stat_cols,
      feat_mu       = data$feat_mu,
      feat_sd       = data$feat_sd,
      best_hp       = hb$best_hp,
      best_val_loss = hb$best_val_loss,
      all_results   = hb$all_results,
      metrics       = metrics,
      backend       = "torch",
      device        = device
    )
    class(result) <- c("model_selection_classifier", "tune_classify")
    return(result)
  }

  # ==========================================================================
  # Parallel mode: save data as bigmemory/rds, launch Rscript workers
  # (mirrors .tune.nn.torch parallel path; one Hyperband search per worker,
  # round-robin across GPUs via CUDA_VISIBLE_DEVICES)
  # ==========================================================================
  work_dir <- tempfile("hb_classify_search_")
  dir.create(work_dir, recursive = TRUE)
  results_dir <- file.path(work_dir, "results")
  dir.create(results_dir)

  use_bigmemory <- requireNamespace("bigmemory", quietly = TRUE)

  if (use_bigmemory) {
    .save_bigmatrix <- function(mat, name) {
      bm <- bigmemory::as.big.matrix(as.matrix(mat), type = "double",
              backingfile = paste0(name, ".bin"),
              descriptorfile = paste0(name, ".desc"),
              backingpath = work_dir)
      bigmemory::flush(bm)
      bm
    }
    .save_bigmatrix(data$X_train, "X_train")
    .save_bigmatrix(data$X_val,   "X_val")
    # Y_train / Y_val are integer vectors — store as 1-col double matrices
    .save_bigmatrix(matrix(as.numeric(data$Y_train), ncol = 1L), "Y_train")
    .save_bigmatrix(matrix(as.numeric(data$Y_val),   ncol = 1L), "Y_val")

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

  # Save metadata only (small)
  n_features  <- data$n_features
  n_classes   <- data$n_classes
  class_names <- data$class_names
  feat_mu     <- data$feat_mu
  feat_sd     <- data$feat_sd
  stat_cols   <- data$stat_cols

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

  save(n_features, n_classes, class_names, feat_mu, feat_sd, stat_cols,
       use_bigmemory,
       search_seeds, saved_lib_paths, threads_per_worker,
       search_max_epochs, search_eta, search_space_saved,
       pkg_source_dir,
       file = file.path(work_dir, "shared_search_meta.RData"))

  # Write classifier worker script
  .torch.write.classify.search.worker.script(file.path(work_dir, "_search_worker.R"))

  if (verbose) {
    cat(sprintf("PipeMaster:: Launching %d searches (%d concurrent, %d GPUs)\n",
                n_searches, n_concurrent, gpus))
    for (k in seq_len(n_searches))
      cat(sprintf("  Search %d: eta=%d, max_epochs=%d\n",
                  k, search_eta[k], search_max_epochs[k]))
  }

  # Build task list with GPU assignment
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

  # Collect results
  all_search_results <- list()
  search_entries     <- list()

  for (k in seq_len(n_searches)) {
    search_dir <- file.path(results_dir, sprintf("search_%04d", k))
    rds_file   <- file.path(search_dir, "result.rds")

    if (file.exists(rds_file)) {
      res <- tryCatch(readRDS(rds_file), error = function(e) NULL)
      if (!is.null(res)) {
        search_res <- res$all_results
        if (!is.null(search_res) && nrow(search_res) > 0) {
          search_res$search <- k
          all_search_results[[length(all_search_results) + 1L]] <- search_res
        }

        if (verbose)
          cat(sprintf("  Search %d: val_CE=%.6f\n", k, res$best_val_loss))

        search_entries[[length(search_entries) + 1L]] <- list(
          val_loss   = res$best_val_loss,
          hp         = res$best_hp,
          best_state = res$best_state
        )
      }
    } else {
      if (verbose) cat(sprintf("  Search %d: FAILED (no result)\n", k))
    }
  }

  if (length(search_entries) == 0L)
    stop("All parallel classifier searches failed. Check worker logs in ", work_dir)

  # Pick the single best across searches
  vl_vec <- vapply(search_entries, `[[`, numeric(1), "val_loss")
  best_idx <- which.min(vl_vec)
  best_entry <- search_entries[[best_idx]]

  # Reconstruct best model on target device. Convert the worker-saved plain
  # numeric arrays back to torch tensors before load_state_dict.
  device <- if (gpus > 0L && torch::cuda_is_available()) "cuda" else "cpu"
  best_model <- .PipeMasterResNet(n_features, n_classes, best_entry$hp)
  best_model$to(device = torch::torch_device(device))
  best_state_tensors <- lapply(best_entry$best_state,
                                function(arr) torch::torch_tensor(arr))
  best_model$load_state_dict(best_state_tensors)

  metrics <- .torch.classify.metrics(best_model, data, device = device)

  combined_results <- if (length(all_search_results) > 0)
    do.call(rbind, all_search_results) else data.frame()

  if (verbose) {
    cat(sprintf("\nPipeMaster:: Best across %d searches: val_CE=%.6f (search %d)\n",
                length(search_entries), best_entry$val_loss, best_idx))
    cat(sprintf("PipeMaster:: Val accuracy = %.3f\n", metrics$accuracy))
    cat(sprintf("PipeMaster:: Confusion matrix:\n"))
    print(metrics$confusion)
  }

  Sys.sleep(1)
  unlink(work_dir, recursive = TRUE)

  result <- list(
    model         = best_model,
    data          = data,
    class_names   = class_names,
    n_features    = n_features,
    n_classes     = n_classes,
    feature.cols  = stat_cols,
    feat_mu       = feat_mu,
    feat_sd       = feat_sd,
    best_hp       = best_entry$hp,
    best_val_loss = best_entry$val_loss,
    all_results   = combined_results,
    metrics       = metrics,
    backend       = "torch",
    device        = device
  )
  class(result) <- c("model_selection_classifier", "tune_classify")
  result
}

# ----------------------------------------------------------------------------
# Public API: nn.predict.classify
# ----------------------------------------------------------------------------

#' Predict the best demographic model for an observed dataset
#'
#' Forward pass through the trained classifier with one of three output
#' rules:
#' \itemize{
#'   \item \code{"softmax"} (default): the classifier's native softmax over
#'         logits. Measures class-direction alignment in latent space; can
#'         disagree with class-manifold distance when obs is OOD.
#'   \item \code{"mahalanobis"}: Mahalanobis distance from obs's penultimate
#'         latent to each class centroid using the pooled within-class
#'         covariance (Lee et al. 2018, NeurIPS). Probability via
#'         \code{softmax(-d/2)}. Density-aware; better for OOD obs.
#'   \item \code{"latent_nn"}: per-class nearest-neighbor distance in the
#'         penultimate latent space; probability via \code{softmax(-d/2)}.
#'         Pure-geometric; predicts the class with the closest sim.
#' }
#'
#' @param classifier Object returned by \code{tune.nn.classify()}.
#' @param observed Numeric vector or 1-row data.frame/matrix of observed
#'   summary statistics (column names matching the training feature set
#'   *before* log1p augmentation).
#' @param method One of \code{"softmax"}, \code{"mahalanobis"}, or
#'   \code{"latent_nn"}. Default \code{"softmax"}.
#' @param reftable Reference table containing classifier features and a
#'   label column. Required for \code{"mahalanobis"} and \code{"latent_nn"};
#'   ignored for \code{"softmax"}.
#' @param model_col Column name with class labels in \code{reftable}
#'   (default \code{"model_id"}).
#' @param subsample Per-class subsample size for distance-based methods
#'   (default \code{NULL} = use all sims).
#' @param mc_samples Number of MC-dropout forward passes for CIs (default 0
#'   = single deterministic prediction). Only applies to \code{"softmax"}.
#' @param ci_probs Two-element numeric, lower/upper percentiles for CI.
#' @param seed Random seed for subsampling.
#' @param verbose Print progress.
#'
#' @return A \code{model_selection} object: list with
#'   \code{method}, \code{prob_mean}, \code{prob_ci}, \code{prob_samples},
#'   \code{distances} (NULL for softmax), \code{best_model},
#'   \code{best_prob}, \code{class_names}, \code{mc_samples}.
#'
#' @export
nn.predict.classify <- function(classifier, observed,
                                 method = c("softmax", "mahalanobis", "latent_nn"),
                                 reftable = NULL,
                                 model_col = "model_id",
                                 subsample = NULL,
                                 mc_samples = 0L,
                                 ci_probs = c(0.025, 0.975),
                                 seed = 42,
                                 verbose = TRUE) {
  if (!inherits(classifier, "model_selection_classifier"))
    stop("classifier must be an object from tune.nn.classify()")
  method <- match.arg(method)

  device      <- classifier$device %||% "cpu"
  feat_cols   <- classifier$feature.cols
  feat_mu     <- classifier$feat_mu
  feat_sd     <- classifier$feat_sd
  K           <- classifier$n_classes
  class_names <- classifier$class_names

  # ---- Coerce observed to named numeric vector ----
  obs <- if (is.data.frame(observed) || is.matrix(observed)) {
    if (nrow(observed) > 1L)
      stop("observed must be a single row (one observed dataset)")
    setNames(as.numeric(observed[1, ]), colnames(observed))
  } else if (!is.null(names(observed))) {
    as.numeric(observed)
  } else {
    stop("observed must have named entries (statistic names)")
  }
  if (is.null(names(obs))) names(obs) <- colnames(observed)

  missing_cols <- setdiff(feat_cols, names(obs))
  if (length(missing_cols) > 0L)
    stop("observed is missing feature columns: ",
         paste(missing_cols, collapse = ", "))
  obs_raw <- obs[feat_cols]
  obs_aug <- c(obs_raw, log1p(abs(obs_raw)))
  obs_z_mat <- matrix((as.numeric(obs_aug) - feat_mu) / feat_sd, nrow = 1L)

  # ---- Method dispatch ----
  if (method == "softmax") {
    return(.predict.classify.softmax(classifier, obs_z_mat, K, class_names,
                                      device, mc_samples, ci_probs, verbose))
  }

  # mahalanobis / latent_nn
  if (mc_samples > 0L)
    warning(sprintf("mc_samples ignored for method='%s'", method))
  if (is.null(reftable))
    stop(sprintf("method='%s' requires the 'reftable' argument", method))
  if (!model_col %in% colnames(reftable))
    stop(sprintf("model_col '%s' not in reftable", model_col))

  raw_labels <- as.character(reftable[[model_col]])
  miss_lbl <- setdiff(unique(raw_labels), class_names)
  if (length(miss_lbl) > 0L)
    stop("reftable contains labels not seen by classifier: ",
         paste(miss_lbl, collapse = ", "))
  model_id_int <- match(raw_labels, class_names)

  set.seed(seed)
  if (!is.null(subsample)) {
    counts <- table(model_id_int)
    if (any(counts > subsample)) {
      keep_idx <- unlist(lapply(seq_len(K), function(k) {
        r <- which(model_id_int == k)
        if (length(r) <= subsample) r else sample(r, subsample)
      }))
      reftable <- reftable[keep_idx, , drop = FALSE]
      model_id_int <- model_id_int[keep_idx]
    }
  }

  S_raw <- as.matrix(reftable[, feat_cols, drop = FALSE])
  bad <- apply(S_raw, 1, function(x) any(!is.finite(x)))
  if (any(bad)) {
    if (verbose)
      cat(sprintf("PipeMaster:: dropping %d sims with non-finite raw stats\n",
                  sum(bad)))
    S_raw <- S_raw[!bad, , drop = FALSE]
    model_id_int <- model_id_int[!bad]
  }
  S_aug <- cbind(S_raw, log1p(abs(S_raw)))
  S_z   <- sweep(sweep(S_aug, 2, feat_mu, "-"), 2, feat_sd, "/")
  bad_z <- apply(S_z, 1, function(x) any(!is.finite(x)))
  if (any(bad_z)) {
    if (verbose)
      cat(sprintf("PipeMaster:: dropping %d sims with non-finite z-scored features\n",
                  sum(bad_z)))
    S_z <- S_z[!bad_z, , drop = FALSE]
    model_id_int <- model_id_int[!bad_z]
  }

  if (verbose)
    cat(sprintf("PipeMaster:: extracting latent (%d sims)\n", nrow(S_z)))
  dev <- torch::torch_device(device)
  Z_sim <- .torch.penultimate(classifier$model, S_z, device = dev)
  Z_obs <- .torch.penultimate(classifier$model, obs_z_mat, device = dev)

  # Filter latent rows that came back non-finite (rare numerical edge case
  # in the trained net for extreme prior-tail inputs).
  bad_lat <- apply(Z_sim, 1, function(x) any(!is.finite(x)))
  if (any(bad_lat)) {
    if (verbose)
      cat(sprintf("PipeMaster:: dropping %d sims with non-finite latent\n",
                  sum(bad_lat)))
    Z_sim <- Z_sim[!bad_lat, , drop = FALSE]
    model_id_int <- model_id_int[!bad_lat]
  }
  if (any(!is.finite(Z_obs)))
    stop("Observed produces non-finite latent activations. ",
         "Inspect obs (extreme/missing values?) and the classifier's eval-mode forward pass.")
  if (nrow(Z_sim) < K * 2L)
    stop(sprintf("After filtering, only %d valid latent rows for %d classes; ",
                 nrow(Z_sim), K),
         "cannot estimate per-class density.")

  # Z-score latent globally so cross-class distances share metric
  mu_z <- colMeans(Z_sim)
  sd_z <- apply(Z_sim, 2, sd); sd_z[sd_z == 0 | !is.finite(sd_z)] <- 1
  Z_n     <- sweep(sweep(Z_sim, 2, mu_z, "-"), 2, sd_z, "/")
  z_obs_n <- as.numeric((Z_obs[1, ] - mu_z) / sd_z)

  if (method == "latent_nn") {
    if (verbose) cat("PipeMaster:: per-class NN distance (latent_nn)\n")
    distances <- vapply(seq_len(K), function(k) {
      rows <- which(model_id_int == k)
      if (length(rows) < 1L) return(NA_real_)
      Z_k <- Z_n[rows, , drop = FALSE]
      .ood.knn.obs(Z_k, matrix(z_obs_n, nrow = 1L))
    }, numeric(1))
  } else {
    if (verbose) cat("PipeMaster:: Mahalanobis with pooled within-class covariance\n")
    class_means <- t(vapply(seq_len(K), function(k) {
      rows <- which(model_id_int == k)
      if (length(rows) < 1L) rep(NA_real_, ncol(Z_n))
      else colMeans(Z_n[rows, , drop = FALSE])
    }, numeric(ncol(Z_n))))
    Z_centered <- Z_n
    for (k in seq_len(K)) {
      rows <- which(model_id_int == k)
      if (length(rows) > 0L)
        Z_centered[rows, ] <- sweep(Z_n[rows, , drop = FALSE], 2,
                                     class_means[k, ], "-")
    }
    Sigma_pooled <- crossprod(Z_centered) / max(1L, nrow(Z_centered) - K)
    ridge <- max(1e-6, 1e-3 * mean(diag(Sigma_pooled)))
    Sigma_inv <- tryCatch(
      solve(Sigma_pooled + diag(ridge, ncol(Z_n))),
      error = function(e) {
        warning("Singular pooled covariance; using ginv")
        if (!requireNamespace("MASS", quietly = TRUE))
          stop("Mahalanobis fallback needs MASS::ginv; install MASS")
        MASS::ginv(Sigma_pooled + diag(ridge, ncol(Z_n)))
      }
    )
    distances <- vapply(seq_len(K), function(k) {
      diff <- z_obs_n - class_means[k, ]
      as.numeric(diff %*% Sigma_inv %*% diff)
    }, numeric(1))
  }
  names(distances) <- class_names

  # softmax(-d/2) for relative probability (Gaussian-likelihood interpretation)
  log_p <- -distances / 2
  log_p <- log_p - max(log_p, na.rm = TRUE)
  prob_mean <- exp(log_p) / sum(exp(log_p), na.rm = TRUE)
  names(prob_mean) <- class_names

  best_idx <- which.min(distances)
  result <- list(
    method        = method,
    prob_mean     = prob_mean,
    prob_ci       = matrix(NA_real_, nrow = 2L, ncol = K,
                            dimnames = list(sprintf("%g%%", 100 * ci_probs),
                                            class_names)),
    prob_samples  = NULL,
    distances     = distances,
    best_model    = class_names[best_idx],
    best_prob     = unname(prob_mean[best_idx]),
    class_names   = class_names,
    mc_samples    = 0L
  )
  class(result) <- "model_selection"
  result
}


# ----------------------------------------------------------------------------
# Internal: native softmax prediction (factored from nn.predict.classify so the
# distance-based methods can call the new top-level easily).
# ----------------------------------------------------------------------------

#' @keywords internal
.predict.classify.softmax <- function(classifier, obs_z, K, class_names,
                                       device, mc_samples, ci_probs, verbose) {
  # MC-dropout: set the whole model to eval (so BatchNorm uses running stats
  # and isn't corrupted by batch-of-1 inputs), then selectively enable only
  # the dropout layers. This is the standard MC-dropout pattern -- without it,
  # train()-mode BN updates running stats from single-sample batches and the
  # running_var collapses toward 0, breaking later eval-mode forward passes.
  predict_once <- function(use_dropout) {
    classifier$model$eval()  # BN stays in eval, dropout disabled by default
    if (use_dropout) {
      for (nm in names(classifier$model$modules)) {
        m <- classifier$model$modules[[nm]]
        if (inherits(m, "nn_dropout") || inherits(m, "nn_dropout_") ||
            inherits(m, "MCDropout"))
          m$train()
      }
    }
    on.exit(classifier$model$eval(), add = TRUE)
    Xt <- torch::torch_tensor(obs_z, dtype = torch::torch_float(),
                              device = torch::torch_device(device))
    torch::with_no_grad({
      logits <- as.matrix(classifier$model(Xt)$cpu())
    })
    as.numeric(.softmax(logits))
  }

  if (mc_samples >= 1L) {
    if (verbose)
      cat(sprintf("PipeMaster:: MC-dropout sampling (%d passes)\n", mc_samples))
    samples <- matrix(NA_real_, nrow = mc_samples, ncol = K)
    for (i in seq_len(mc_samples)) samples[i, ] <- predict_once(TRUE)
    prob_mean <- colMeans(samples)
    prob_ci <- apply(samples, 2, quantile, probs = ci_probs)
    rownames(prob_ci) <- sprintf("%g%%", 100 * ci_probs)
    colnames(prob_ci) <- class_names
    names(prob_mean) <- class_names
  } else {
    samples   <- NULL
    prob_mean <- setNames(predict_once(FALSE), class_names)
    prob_ci   <- matrix(NA_real_, nrow = 2L, ncol = K,
                        dimnames = list(sprintf("%g%%", 100 * ci_probs),
                                        class_names))
  }

  best_idx <- which.max(prob_mean)
  result <- list(
    method       = "softmax",
    prob_mean    = prob_mean,
    prob_ci      = prob_ci,
    prob_samples = samples,
    distances    = NULL,
    best_model   = class_names[best_idx],
    best_prob    = unname(prob_mean[best_idx]),
    class_names  = class_names,
    mc_samples   = as.integer(mc_samples)
  )
  class(result) <- "model_selection"
  result
}

# Internal "or-else" helper (used in this file only)
`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# Save / load
# ----------------------------------------------------------------------------

#' Save a classifier produced by `tune.nn.classify()`
#' @param classifier the object returned by `tune.nn.classify()`
#' @param path file path (`.pt` extension recommended)
#' @export
save.classifier <- function(classifier, path) {
  if (!inherits(classifier, "model_selection_classifier"))
    stop("classifier must come from tune.nn.classify()")

  meta <- list(
    state_dict   = classifier$model$state_dict(),
    hp           = classifier$best_hp,
    n_features   = classifier$n_features,
    n_classes    = classifier$n_classes,
    class_names  = classifier$class_names,
    feature.cols = classifier$feature.cols,
    feat_mu      = classifier$feat_mu,
    feat_sd      = classifier$feat_sd,
    best_val_loss = classifier$best_val_loss,
    metrics      = classifier$metrics,
    all_results  = classifier$all_results,
    backend      = classifier$backend
  )
  torch::torch_save(meta, path)
  invisible(path)
}

#' Load a classifier produced by `save.classifier()`
#' @param path path written by `save.classifier()`
#' @param device "cpu" or "cuda"; default chooses cuda when available
#' @export
load.classifier <- function(path, device = NULL) {
  meta <- torch::torch_load(path)
  if (is.null(device))
    device <- if (torch::cuda_is_available()) "cuda" else "cpu"

  model <- .PipeMasterResNet(meta$n_features, meta$n_classes, meta$hp)
  model$to(device = torch::torch_device(device))
  model$load_state_dict(meta$state_dict)
  model$eval()

  result <- list(
    model        = model,
    data         = NULL,
    class_names  = meta$class_names,
    n_features   = meta$n_features,
    n_classes    = meta$n_classes,
    feature.cols = meta$feature.cols,
    feat_mu      = meta$feat_mu,
    feat_sd      = meta$feat_sd,
    best_hp      = meta$hp,
    best_val_loss = meta$best_val_loss,
    metrics      = meta$metrics,
    all_results  = meta$all_results,
    backend      = meta$backend %||% "torch",
    device       = device
  )
  class(result) <- c("model_selection_classifier", "tune_classify")
  result
}

# ----------------------------------------------------------------------------
# S3 methods for model_selection
# ----------------------------------------------------------------------------

#' @export
print.model_selection <- function(x, ...) {
  meth <- if (is.null(x$method)) "softmax" else x$method
  cat(sprintf("Model selection [%s]: best = %s (P = %.3f)\n",
              meth, x$best_model, x$best_prob))
  cat("All class probabilities:\n")
  print(round(x$prob_mean, 4))
  if (!is.null(x$distances)) {
    cat("\nPer-class distances (smaller = closer):\n")
    print(round(x$distances, 3))
  }
  if (x$mc_samples > 0L) {
    cat(sprintf("\n%s-%s CIs from %d MC-dropout passes:\n",
                rownames(x$prob_ci)[1],
                rownames(x$prob_ci)[2],
                x$mc_samples))
    print(round(x$prob_ci, 4))
  }
  invisible(x)
}

#' @export
summary.model_selection <- function(object, ...) {
  K <- length(object$class_names)
  meth <- if (is.null(object$method)) "softmax" else object$method
  out <- data.frame(
    model = object$class_names,
    prob  = round(unname(object$prob_mean), 4),
    stringsAsFactors = FALSE
  )
  if (!is.null(object$distances))
    out$distance <- round(unname(object$distances), 3)
  if (object$mc_samples > 0L) {
    out$ci_lo <- round(object$prob_ci[1, ], 4)
    out$ci_hi <- round(object$prob_ci[2, ], 4)
  }
  # Sort by closeness: distance ascending if available, else prob descending
  if (!is.null(object$distances))
    out <- out[order(out$distance), ]
  else
    out <- out[order(-out$prob), ]
  rownames(out) <- NULL
  cat(sprintf("Best model [%s]: %s (P = %.3f)\n\n",
              meth, object$best_model, object$best_prob))
  print(out)
  invisible(out)
}

#' @export
plot.model_selection <- function(x, col = NULL, las = 1, ...) {
  K <- length(x$class_names)
  meth <- if (is.null(x$method)) "softmax" else x$method
  if (is.null(col))
    col <- ifelse(seq_len(K) == which.max(x$prob_mean),
                  "#d62728", "#9ecae1")

  if (!is.null(x$distances)) {
    # Two-panel: probabilities (top) + distances (bottom)
    old_par <- par(no.readonly = TRUE); on.exit(par(old_par))
    par(mfrow = c(2, 1), mar = c(4, 4, 3, 1))

    ord <- order(-x$prob_mean)
    barplot(x$prob_mean[ord], names.arg = x$class_names[ord],
            col = col[ord], ylim = c(0, max(x$prob_mean) * 1.1),
            ylab = "Relative probability",
            main = sprintf("Model selection [%s]: best = %s (P = %.3f)",
                           meth, x$best_model, x$best_prob),
            las = las, ...)
    abline(h = 1 / K, lty = 2, col = "grey50")

    ord_d <- order(x$distances)
    bar_col <- ifelse(seq_along(x$distances) == which.min(x$distances),
                      "#d62728", "#9ecae1")
    barplot(x$distances[ord_d], names.arg = x$class_names[ord_d],
            col = bar_col[ord_d],
            ylab = sprintf("Distance to class (%s)", meth),
            main = "Per-class distance (smaller = closer)",
            las = las)
    return(invisible(NULL))
  }

  ord <- order(-x$prob_mean)
  bp <- barplot(x$prob_mean[ord], names.arg = x$class_names[ord],
                col = col[ord], ylim = c(0, 1.05),
                ylab = "Posterior probability",
                main = sprintf("Model selection [%s]: best = %s (P = %.3f)",
                               meth, x$best_model, x$best_prob),
                las = las, ...)
  if (x$mc_samples > 0L) {
    arrows(bp[, 1], x$prob_ci[1, ord],
           bp[, 1], x$prob_ci[2, ord],
           code = 3, angle = 90, length = 0.05, lwd = 1.5)
  }
  abline(h = 1 / K, lty = 2, col = "grey50")
  text(par("usr")[1] + 0.15, 1 / K + 0.025, "uniform", col = "grey40",
       cex = 0.7, adj = 0)
  invisible(bp)
}
