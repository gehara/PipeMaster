# ============================================================================
# Torch Training Infrastructure for PipeMaster
#
# Custom training loop, Hyperband, save/load — all in pure R torch.
#
# Functions:
#   .torch.train.model()  — custom training loop (early stopping, LR scheduler)
#   .torch.hyperband()    — Hyperband HP search
#   .torch.compute.model.metrics() — R² and MPE on validation set
# ============================================================================

# ============================================================================
# Custom training loop:
#   - Mini-batch SGD via torch::dataloader()
#   - Early stopping with restore_best_weights
#   - ReduceLROnPlateau
#   - Warm-start from initial_epoch (for Hyperband successive halving)
#   - Huber or MSE loss
#
# Returns: list(val_loss = best validation loss, epochs_trained = integer)
# ============================================================================

#' @keywords internal
.torch.train.model <- function(model, X_train, Y_train, X_val, Y_val,
                               hp, type,
                               epochs = 500L, initial_epoch = 0L,
                               patience = 30L, lr_patience = 15L,
                               lr_factor = 0.5, min_lr = 1e-6,
                               verbose = 0L, device = "cpu") {

  dev <- torch::torch_device(device)

  # Data loading strategy:
  # - big.matrix: materialize from mmap to CPU tensor once (avoids per-batch
  #   disk I/O), transfer only batches to GPU. Keeps GPU memory low so many
  #   workers can share GPUs.
  # - tensor/matrix: move to GPU directly (small datasets, single-worker)
  is_bigmem <- inherits(X_train, "big.matrix")

  .to_tensor <- function(x, target_dev) {
    if (inherits(x, "torch_tensor")) return(x$to(dtype = torch::torch_float(), device = target_dev))
    if (inherits(x, "big.matrix")) x <- x[,]
    torch::torch_tensor(x, dtype = torch::torch_float(), device = target_dev)
  }

  if (is_bigmem) {
    # CPU tensors for training data, GPU for validation (small)
    X_tr_t <- .to_tensor(X_train, "cpu")
    Y_tr_t <- .to_tensor(Y_train, "cpu")
  } else {
    X_tr_t <- .to_tensor(X_train, dev)
    Y_tr_t <- .to_tensor(Y_train, dev)
  }
  X_va_t <- .to_tensor(X_val, dev)
  Y_va_t <- .to_tensor(Y_val, dev)
  n_train <- X_tr_t$size(1)

  batch_size <- as.integer(hp$batch_size)

  # Loss function
  loss_type <- if (!is.null(hp$loss)) hp$loss else "huber"
  huber_delta <- if (!is.null(hp$huber_delta)) hp$huber_delta else 1.0

  loss_fn <- if (identical(loss_type, "mse")) {
    function(pred, target) torch::nnf_mse_loss(pred, target)
  } else {
    function(pred, target) .torch.huber.loss(pred, target, delta = huber_delta)
  }

  # Optimizer with L2 regularization (weight_decay).
  # Conventional L2 regularizer_l2(l) adds l*sum(w^2) to loss → gradient adds 2*l*w.
  # Torch's weight_decay adds wd*w to gradient directly, so pass 2*l2_reg.
  l2_reg <- if (!is.null(hp$l2_reg)) hp$l2_reg else 0
  optimizer <- torch::optim_adam(model$parameters, lr = hp$learning_rate,
                                weight_decay = 2 * l2_reg)

  # LR scheduler
  scheduler <- torch::lr_reduce_on_plateau(optimizer, mode = "min",
                                           factor = lr_factor,
                                           patience = as.integer(lr_patience),
                                           min_lr = min_lr)

  # Early stopping state
  best_val_loss <- Inf
  best_state_dict <- NULL
  wait <- 0L
  best_epoch <- initial_epoch

  # Epoch loop
  for (epoch in seq(initial_epoch + 1L, epochs)) {
    model$train()

    # Shuffle indices
    perm <- sample.int(n_train)
    n_batches <- ceiling(n_train / batch_size)
    epoch_loss <- 0.0

    for (b in seq_len(n_batches)) {
      start_idx <- (b - 1L) * batch_size + 1L
      end_idx   <- min(b * batch_size, n_train)
      batch_idx <- perm[start_idx:end_idx]

      if (is_bigmem) {
        X_batch <- X_tr_t[batch_idx]$to(device = dev)
        Y_batch <- Y_tr_t[batch_idx]$to(device = dev)
      } else {
        X_batch <- X_tr_t[batch_idx]
        Y_batch <- Y_tr_t[batch_idx]
      }

      optimizer$zero_grad()
      pred <- model(X_batch)
      loss <- loss_fn(pred, Y_batch)
      loss$backward()
      optimizer$step()

      epoch_loss <- epoch_loss + loss$item() * (end_idx - start_idx + 1L)
    }
    epoch_loss <- epoch_loss / n_train

    # Validation loss
    model$eval()
    val_loss <- torch::with_no_grad({
      pred_val <- model(X_va_t)
      loss_fn(pred_val, Y_va_t)$item()
    })

    # LR scheduler step
    scheduler$step(val_loss)

    # Early stopping check
    if (val_loss < best_val_loss) {
      best_val_loss <- val_loss
      best_state_dict <- lapply(model$state_dict(), function(t) t$clone())
      wait <- 0L
      best_epoch <- epoch
    } else {
      wait <- wait + 1L
      if (wait >= patience) break
    }
  }

  # Restore best weights
  if (!is.null(best_state_dict))
    model$load_state_dict(best_state_dict)

  list(val_loss = best_val_loss, epochs_trained = best_epoch)
}

# ============================================================================
# Hyperband algorithm
#
# Inputs: search_space, data, type, sfs.dims, max_epochs, eta, seed.
# ============================================================================

#' @keywords internal
.torch.hyperband <- function(search_space, data, type, sfs.dims,
                             max_epochs, eta, seed, verbose,
                             device = "cpu") {

  s_max <- min(floor(log(max_epochs) / log(eta)), 3L)

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
  best_state         <- NULL  # torch state_dict (in-memory, no temp files)

  # Precompute total budget for ETA
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

    # Track state_dicts for warm-starting between rounds
    prev_epochs <- rep(0L, n)
    config_states <- vector("list", n)  # list of state_dicts

    for (i in 0:s) {
      r_i <- round(r * eta^i)
      n_i <- max(1, floor(n / eta^i))
      n_keep <- max(1, ceiling(n_i / eta))

      val_losses <- rep(Inf, length(configs))
      cfg_times  <- rep(NA_real_, length(configs))

      for (j in seq_along(configs)) {
        tryCatch({
          cfg_t0 <- proc.time()[3]
          torch::torch_manual_seed(as.integer(seed + s + j))

          model <- .build.nn.torch(configs[[j]], data, type, sfs.dims,
                                   device = device)

          # Load state from previous round (if any)
          if (prev_epochs[j] > 0 && !is.null(config_states[[j]]))
            model$load_state_dict(config_states[[j]])

          result <- .torch.train.model(
            model, data$X_train, data$Y_train, data$X_val, data$Y_val,
            hp = configs[[j]], type = type,
            epochs = as.integer(r_i),
            initial_epoch = as.integer(prev_epochs[j]),
            patience = 10L, lr_patience = 5L,
            device = device
          )

          val_losses[j] <- result$val_loss
          cfg_times[j]  <- proc.time()[3] - cfg_t0

          # Save state for next round
          config_states[[j]] <- lapply(model$state_dict(), function(t) t$clone()$cpu())

          # Update global best
          if (val_losses[j] < global_best_loss) {
            global_best_loss   <- val_losses[j]
            global_best_hp     <- configs[[j]]
            global_best_epochs <- as.integer(r_i)
            best_state <- lapply(model$state_dict(), function(t) t$clone()$cpu())

            if (verbose) cat(sprintf("    \u2605 new best: val_loss=%.4f (bracket %d, round %d)\n",
                                     val_losses[j], s, i))
          }

          rm(model); gc()
          if (device != "cpu") torch::cuda_empty_cache()

        }, error = function(e) {
          if (verbose) cat(sprintf("    [warn] config %d error: %s\n",
                                   j, conditionMessage(e)))
          val_losses[j] <<- Inf
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

        configs       <- configs[keep]
        prev_epochs   <- prev_epochs[keep]
        config_states <- config_states[keep]
      } else {
        if (verbose) cat("\n")
      }
    }

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
    best_state        = best_state,
    all_results       = all_results
  )
}

# ============================================================================
# Compute R² and MPE on validation set (torch backend)
#
# Drop-in replacement for .compute.model.metrics()
# ============================================================================

#' @keywords internal
.torch.compute.model.metrics <- function(model, data, type, device = "cpu") {
  pred_z <- .torch.predict(model, data$X_val, device = device)
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

# ============================================================================
# Ensemble predict (torch backend)
#
# Drop-in replacement for .ensemble.predict()
# ============================================================================

#' @keywords internal
.torch.ensemble.predict <- function(models, models_val_loss, X_obs, data,
                                    type, verbose, device = "cpu") {
  n_models <- length(models)

  # 1. Flatten X_val and X_obs for distance computation
  flat <- .flatten.features(data, type)
  X_val_flat <- flat$X_val

  if (type == "sumstat") {
    x_obs_flat <- as.numeric(X_obs)
  } else if (type == "sfs1d") {
    x_obs_flat <- as.numeric(X_obs[1, , 1])
  } else {
    x_obs_flat <- as.numeric(X_obs[1, , , 1])
  }

  # 2. Euclidean distances
  diffs <- sweep(X_val_flat, 2, x_obs_flat, `-`)
  dists <- sqrt(rowSums(diffs^2))

  # 3. Gaussian kernel weights
  bandwidth <- median(dists)
  if (bandwidth < .Machine$double.eps) bandwidth <- 1.0
  w <- exp(-dists^2 / (2 * bandwidth^2))
  w <- w / sum(w)

  # 4. Per-model proximity-weighted loss
  Y_val <- data$Y_val
  prox_losses <- numeric(n_models)

  for (m in seq_len(n_models)) {
    pred_z <- .torch.predict(models[[m]], data$X_val, device = device)
    per_sample_mse <- rowSums((pred_z - Y_val)^2)
    prox_losses[m] <- sum(w * per_sample_mse)
  }

  # 5. Softmax ensemble weights
  temperature <- median(prox_losses)
  if (temperature < .Machine$double.eps) temperature <- 1.0
  log_weights <- -prox_losses / temperature
  log_weights <- log_weights - max(log_weights)
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

  # 6. Weighted average in Z-space
  pred_z_list <- lapply(seq_len(n_models), function(m) {
    as.numeric(.torch.predict(models[[m]], X_obs, device = device))
  })

  ensemble_z <- rep(0, length(pred_z_list[[1L]]))
  for (m in seq_len(n_models))
    ensemble_z <- ensemble_z + ensemble_w[m] * pred_z_list[[m]]

  as.numeric(.inv.transform(
    matrix(ensemble_z, nrow = 1), data$target_mu, data$target_sd
  ))
}

# ============================================================================
# Write standalone search worker script (torch backend)
#
# Each worker runs a complete Hyperband search in a separate R process.
# ============================================================================

#' @keywords internal
.torch.write.search.worker.script <- function(filepath) {
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
    '# GPU setup',
    'device <- "cpu"',
    'if (Sys.getenv("CUDA_VISIBLE_DEVICES", "-1") != "-1") {',
    '  if (cuda_is_available()) device <- "cuda"',
    '}',
    '',
    '# Load data — bigmemory (shared mmap) or rds (per-worker copy)',
    'if (use_bigmemory) {',
    '  suppressPackageStartupMessages(library(bigmemory))',
    '  X_train <- attach.big.matrix("X_train.desc")',
    '  Y_train <- attach.big.matrix("Y_train.desc")',
    '  X_val   <- torch_tensor(attach.big.matrix("X_val.desc")[,], dtype = torch_float())',
    '  Y_val   <- torch_tensor(attach.big.matrix("Y_val.desc")[,], dtype = torch_float())',
    '} else {',
    '  .rds_to_tensor <- function(path) {',
    '    m <- readRDS(path)',
    '    t <- torch_tensor(m, dtype = torch_float())',
    '    rm(m); gc()',
    '    t',
    '  }',
    '  X_train <- .rds_to_tensor("X_train.rds")',
    '  Y_train <- .rds_to_tensor("Y_train.rds")',
    '  X_val   <- .rds_to_tensor("X_val.rds")',
    '  Y_val   <- .rds_to_tensor("Y_val.rds")',
    '}',
    '',
    '# Create task-specific output directory',
    'out_dir <- file.path("results", sprintf("search_%04d", task_id))',
    'dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)',
    '',
    '# Build data list',
    'hb_data <- list(',
    '  X_train = X_train, X_val = X_val,',
    '  Y_train = Y_train, Y_val = Y_val,',
    '  n_features = n_features, n_targets = n_targets, n_bins = n_bins',
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
    'hb <- PipeMaster:::.torch.hyperband(',
    '  search_space = search_space_saved,',
    '  data = hb_data, type = type, sfs.dims = sfs.dims,',
    '  max_epochs = worker_max_epochs, eta = worker_eta,',
    '  seed = worker_seed, verbose = TRUE, device = device)',
    '',
    '# Check parent alive before retraining',
    'if (file.exists(".PM_parent.pid") && !PipeMaster:::.pm.parent.alive(".PM_parent.pid")) {',
    '  cat("Parent died, worker exiting.\\n"); q("no") }',
    '',
    '# Retrain best config to worker_max_epochs',
    'torch_manual_seed(as.integer(worker_seed))',
    'model <- PipeMaster:::.build.nn.torch(hb$best_hp, hb_data, type, sfs.dims,',
    '                                       device = device)',
    '',
    'if (!is.null(hb$best_state))',
    '  model$load_state_dict(hb$best_state)',
    '',
    'start_epoch <- if (!is.null(hb$best_state)) as.integer(hb$best_epochs) else 0L',
    'final_vl <- hb$best_val_loss',
    '',
    'if (start_epoch < worker_max_epochs) {',
    '  cat(sprintf("Search %d: retraining from epoch %d to %d\\n",',
    '              task_id, start_epoch, worker_max_epochs))',
    '  retrain <- PipeMaster:::.torch.train.model(',
    '    model, hb_data$X_train, hb_data$Y_train,',
    '    hb_data$X_val, hb_data$Y_val,',
    '    hp = hb$best_hp, type = type,',
    '    epochs = worker_max_epochs,',
    '    initial_epoch = start_epoch,',
    '    patience = 30L, lr_patience = 15L,',
    '    device = device)',
    '  final_vl <- min(final_vl, retrain$val_loss)',
    '}',
    '',
    '# Save model and results',
    'PipeMaster:::.torch.save.model(model, file.path(out_dir, "best_model.pt"),',
    '  type = type, hp = hb$best_hp,',
    '  n_features = n_features,',
    '  n_bins = hb_data$n_bins,',
    '  sfs_dims = sfs.dims,',
    '  n_targets = n_targets)',
    '',
    'saveRDS(list(',
    '  best_hp       = hb$best_hp,',
    '  best_val_loss = final_vl,',
    '  all_results   = hb$all_results',
    '), file.path(out_dir, "result.rds"))',
    '',
    'cat(sprintf("Search %d done (val_loss=%.6f)\\n", task_id, final_vl))',
    'writeLines("done", file.path(out_dir, "done.txt"))'
  ), filepath)
}

# ============================================================================
# Write standalone model-builder script for parallel workers.
# Used by conformal/bootstrap workers.
# ============================================================================

#' @keywords internal
.torch.write.builder.script <- function(filepath, type) {
  lines <- c(
    '# Auto-generated torch model builder for nn.predict parallel workers',
    'build_nn_torch <- function(hp, data, type, sfs.dims) {',
    '  PipeMaster:::.build.nn.torch(hp, data, type, sfs.dims)',
    '}',
    '',
    'train_model_torch <- function(model, data, hp, type, max_epochs, device = "cpu") {',
    '  PipeMaster:::.torch.train.model(',
    '    model, data$X_train, data$Y_train, data$X_val, data$Y_val,',
    '    hp = hp, type = type, epochs = as.integer(max_epochs),',
    '    patience = 30L, lr_patience = 15L, device = device)',
    '}',
    '',
    'predict_torch <- function(model, X, device = "cpu") {',
    '  PipeMaster:::.torch.predict(model, X, device = device)',
    '}'
  )
  writeLines(lines, filepath)
}

# ============================================================================
# save / load tune.nn results (torch backend)
# ============================================================================

#' @keywords internal
.torch.save.tune.result <- function(tune.result, path) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  type <- tune.result$type
  hp   <- tune.result$best_hp
  data <- tune.result$data
  n_features <- if (type %in% c("sumstat", "emulator")) ncol(data$X_train) else NULL
  n_bins     <- data$n_bins
  sfs_dims   <- tune.result$sfs.dims
  n_targets  <- ncol(data$Y_train)

  # Save best_model
  .torch.save.model(tune.result$best_model,
                    file.path(path, "best_model.pt"),
                    type = type, hp = hp,
                    n_features = n_features, n_bins = n_bins,
                    sfs_dims = sfs_dims, n_targets = n_targets)

  # Save all top-K models
  models <- tune.result$models
  if (!is.null(models) && length(models) > 0L) {
    models_dir <- file.path(path, "models")
    dir.create(models_dir, showWarnings = FALSE, recursive = TRUE)
    for (i in seq_along(models)) {
      model_hp <- if (!is.null(tune.result$models_hp) && i <= length(tune.result$models_hp))
        tune.result$models_hp[[i]] else hp
      .torch.save.model(models[[i]],
                        file.path(models_dir, sprintf("model_%03d.pt", i)),
                        type = type, hp = model_hp,
                        n_features = n_features, n_bins = n_bins,
                        sfs_dims = sfs_dims, n_targets = n_targets)
    }
  }

  # Save everything else as RDS (strip torch models)
  result_no_model <- tune.result
  result_no_model$best_model <- NULL
  result_no_model$models <- NULL
  saveRDS(result_no_model, file.path(path, "tune_result.rds"))

  cat(sprintf("PipeMaster:: Saved tune.nn result to %s (%d models, torch format)\n",
              path, length(models)))
}

#' @keywords internal
.torch.load.tune.result <- function(path, device = "cpu") {
  rds_file <- file.path(path, "tune_result.rds")

  if (!file.exists(rds_file))
    stop("tune_result.rds not found in: ", path)

  result <- readRDS(rds_file)

  # Load best_model
  best_pt <- file.path(path, "best_model.pt")
  if (file.exists(best_pt)) {
    result$best_model <- .torch.load.model(best_pt, device = device)
  } else {
    stop("best_model.pt not found in: ", path)
  }

  # Load top-K models
  models_dir <- file.path(path, "models")
  if (dir.exists(models_dir)) {
    model_files <- sort(list.files(models_dir, pattern = "\\.pt$", full.names = TRUE))
    if (length(model_files) > 0L) {
      result$models <- lapply(model_files, function(f) {
        tryCatch(.torch.load.model(f, device = device), error = function(e) {
          warning("Could not load model from ", f, ": ",
                  conditionMessage(e), call. = FALSE)
          NULL
        })
      })
      keep <- !vapply(result$models, is.null, logical(1))
      result$models <- result$models[keep]
      if (!is.null(result$models_val_loss))
        result$models_val_loss <- result$models_val_loss[keep]
      cat(sprintf("PipeMaster:: Loaded tune.nn result from %s (%d models, torch)\n",
                  path, length(result$models)))
    } else {
      result$models <- list(result$best_model)
      cat(sprintf("PipeMaster:: Loaded tune.nn result from %s (torch)\n", path))
    }
  } else {
    result$models <- list(result$best_model)
    if (is.null(result$models_val_loss))
      result$models_val_loss <- result$best_val_loss
    cat(sprintf("PipeMaster:: Loaded tune.nn result from %s (torch)\n", path))
  }

  result
}
