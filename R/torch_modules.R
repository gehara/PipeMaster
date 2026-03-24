# ============================================================================
# Torch Neural Network Modules for PipeMaster
#
# Pure R torch implementations of the keras architectures used in tune.nn()
# and train.emulator(). Eliminates the Python/TensorFlow/keras dependency.
#
# Modules:
#   .MCDropout         — always-on dropout (for MC dropout inference)
#   .ResBlock          — single residual block (dense → BN → dense → BN + skip)
#   .PipeMasterResNet  — full ResNet for sumstat / emulator
#   .ConvBlock1D       — conv1d → BN → ReLU (+ optional residual)
#   .PipeMasterCNN1D   — full 1D CNN for single-pop SFS
#   .ConvBlock2D       — conv2d → BN → ReLU (+ optional residual)
#   .PipeMasterCNN2D   — full 2D CNN for joint SFS
#
# Builders:
#   .build.nn.torch()  — unified dispatcher (replaces .build.nn())
#   .torch.predict()   — forward pass wrapper (replaces predict(keras_model, ...))
#
# ============================================================================

# ============================================================================
# MC Dropout: always-on dropout (ignores eval/train mode)
# ============================================================================

.MCDropout <- torch::nn_module(
  "MCDropout",
  initialize = function(p = 0.5) {
    self$p <- p
  },
  forward = function(x) {
    torch::nnf_dropout(x, p = self$p, training = TRUE)
  }
)

# ============================================================================
# Residual Block: Dense → BN → Dense → BN + skip → ReLU
# ============================================================================

.ResBlock <- torch::nn_module(
  "ResBlock",
  initialize = function(units) {
    self$dense1 <- torch::nn_linear(units, units)
    self$bn1    <- torch::nn_batch_norm1d(units)
    self$dense2 <- torch::nn_linear(units, units)
    self$bn2    <- torch::nn_batch_norm1d(units)
  },
  forward = function(x) {
    skip <- x
    x <- self$dense1(x)
    x <- torch::nnf_relu(x)
    x <- self$bn1(x)
    x <- self$dense2(x)
    x <- self$bn2(x)
    x <- x + skip
    torch::nnf_relu(x)
  }
)

# ============================================================================
# ResNet for summary statistics / emulator
#
# Architecture:
#   input → dense(units_1) → BN → [ResBlock(units_1)] × n_resblocks_1
#         → dense(units_2) → BN → [dropout] → [ResBlock(units_2)] × n_resblocks_2
#         → dense(units_3) → BN → [dropout] → dense(n_targets)
# ============================================================================

.PipeMasterResNet <- torch::nn_module(
  "PipeMasterResNet",

  initialize = function(n_features, n_targets, hp, mc_dropout = FALSE) {
    # Store config for serialization
    self$hp <- hp
    self$n_features <- n_features
    self$n_targets <- n_targets
    self$mc_dropout_mode <- mc_dropout

    use_drop <- isTRUE(hp$use_dropout) && hp$dropout > 0

    # First group: projection → BN → residual blocks
    self$proj1  <- torch::nn_linear(n_features, as.integer(hp$units_1))
    self$bn1    <- torch::nn_batch_norm1d(as.integer(hp$units_1))
    self$res1   <- torch::nn_module_list()
    for (i in seq_len(hp$n_resblocks_1))
      self$res1$append(.ResBlock(as.integer(hp$units_1)))

    # Middle transition
    self$proj2  <- torch::nn_linear(as.integer(hp$units_1), as.integer(hp$units_2))
    self$bn2    <- torch::nn_batch_norm1d(as.integer(hp$units_2))
    self$drop1  <- if (use_drop) {
      if (mc_dropout) .MCDropout(hp$dropout) else torch::nn_dropout(p = hp$dropout)
    } else NULL

    # Second residual group
    self$res2 <- torch::nn_module_list()
    if (hp$n_resblocks_2 > 0)
      for (i in seq_len(hp$n_resblocks_2))
        self$res2$append(.ResBlock(as.integer(hp$units_2)))

    # Final dense
    self$proj3 <- torch::nn_linear(as.integer(hp$units_2), as.integer(hp$units_3))
    self$bn3   <- torch::nn_batch_norm1d(as.integer(hp$units_3))
    self$drop2 <- if (use_drop) {
      if (mc_dropout) .MCDropout(hp$dropout) else torch::nn_dropout(p = hp$dropout)
    } else NULL

    # Output
    self$output_layer <- torch::nn_linear(as.integer(hp$units_3), n_targets)
  },

  forward = function(x) {
    # First group
    x <- torch::nnf_relu(self$bn1(self$proj1(x)))
    for (i in seq_along(self$res1))
      x <- self$res1[[i]](x)

    # Middle transition
    x <- torch::nnf_relu(self$bn2(self$proj2(x)))
    if (!is.null(self$drop1))
      x <- self$drop1(x)

    # Second residual group
    for (i in seq_along(self$res2))
      x <- self$res2[[i]](x)

    # Final
    x <- torch::nnf_relu(self$bn3(self$proj3(x)))
    if (!is.null(self$drop2))
      x <- self$drop2(x)

    self$output_layer(x)
  }
)

# ============================================================================
# 1D CNN for single-population SFS
#
# Architecture:
#   input [N, n_bins, 1] → ConvBlocks → GlobalAvgPool1D → Dense head → output
# ============================================================================

.ConvBlock1D <- torch::nn_module(
  "ConvBlock1D",

  initialize = function(in_channels, out_channels, kernel_size,
                        use_residual = FALSE) {
    self$use_residual <- use_residual
    pad <- (kernel_size - 1L) %/% 2L  # "same" padding

    self$conv1 <- torch::nn_conv1d(in_channels, out_channels, kernel_size,
                                   padding = pad)
    self$bn1   <- torch::nn_batch_norm1d(out_channels)

    if (use_residual) {
      self$conv2 <- torch::nn_conv1d(out_channels, out_channels, kernel_size,
                                     padding = pad)
      self$bn2   <- torch::nn_batch_norm1d(out_channels)
      self$conv3 <- torch::nn_conv1d(out_channels, out_channels, kernel_size,
                                     padding = pad)
      self$bn3   <- torch::nn_batch_norm1d(out_channels)
    }
  },

  forward = function(x) {
    x <- torch::nnf_relu(self$bn1(self$conv1(x)))

    if (self$use_residual) {
      skip <- x
      x <- torch::nnf_relu(self$bn2(self$conv2(x)))
      x <- self$bn3(self$conv3(x))
      x <- torch::nnf_relu(x + skip)
    }
    x
  }
)

.PipeMasterCNN1D <- torch::nn_module(
  "PipeMasterCNN1D",

  initialize = function(n_bins, n_targets, hp, mc_dropout = FALSE) {
    self$hp <- hp
    self$n_bins <- n_bins
    self$n_targets <- n_targets
    self$mc_dropout_mode <- mc_dropout

    # Build conv blocks
    self$conv_blocks <- torch::nn_module_list()
    in_ch <- 1L
    for (b in seq_len(hp$n_blocks)) {
      filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))
      ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))
      use_res <- isTRUE(hp$use_residual) && b > 1 && b < hp$n_blocks

      self$conv_blocks$append(.ConvBlock1D(in_ch, filters, ks,
                                           use_residual = use_res))
      in_ch <- filters
    }
    # last_filters = in_ch after all conv blocks

    # Dense head
    self$dense_layers <- torch::nn_module_list()
    self$dense_bns    <- torch::nn_module_list()
    self$dense_drops  <- torch::nn_module_list()
    prev_units <- in_ch  # after global avg pool
    for (d in seq_len(hp$n_dense)) {
      units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
      self$dense_layers$append(torch::nn_linear(prev_units, units))
      self$dense_bns$append(torch::nn_batch_norm1d(units))
      if (hp$dropout > 0) {
        self$dense_drops$append(
          if (mc_dropout) .MCDropout(hp$dropout) else torch::nn_dropout(p = hp$dropout)
        )
      } else {
        self$dense_drops$append(NULL)
      }
      prev_units <- units
    }

    self$output_layer <- torch::nn_linear(prev_units, n_targets)
  },

  forward = function(x) {
    # x shape: [N, n_bins, 1] from R — transpose to [N, 1, n_bins] for conv1d
    if (x$dim() == 3L && x$size(3) == 1L)
      x <- x$permute(c(1L, 3L, 2L))

    for (i in seq_along(self$conv_blocks))
      x <- self$conv_blocks[[i]](x)

    # Global average pooling: [N, C, L] → [N, C]
    x <- x$mean(dim = 3L)

    # Dense head
    for (d in seq_along(self$dense_layers)) {
      x <- torch::nnf_relu(self$dense_bns[[d]](self$dense_layers[[d]](x)))
      if (!is.null(self$dense_drops[[d]]))
        x <- self$dense_drops[[d]](x)
    }

    self$output_layer(x)
  }
)

# ============================================================================
# 2D CNN for joint SFS
#
# Architecture:
#   input [N, dim1, dim2, 1] → ConvBlocks (+ MaxPool) → GlobalAvgPool2D
#     → Dense head → output
# ============================================================================

.ConvBlock2D <- torch::nn_module(
  "ConvBlock2D",

  initialize = function(in_channels, out_channels, kernel_size,
                        use_residual = FALSE) {
    self$use_residual <- use_residual
    pad <- (kernel_size - 1L) %/% 2L

    self$conv1 <- torch::nn_conv2d(in_channels, out_channels, kernel_size,
                                   padding = pad)
    self$bn1   <- torch::nn_batch_norm2d(out_channels)

    if (use_residual) {
      self$conv2 <- torch::nn_conv2d(out_channels, out_channels, kernel_size,
                                     padding = pad)
      self$bn2   <- torch::nn_batch_norm2d(out_channels)
      self$conv3 <- torch::nn_conv2d(out_channels, out_channels, kernel_size,
                                     padding = pad)
      self$bn3   <- torch::nn_batch_norm2d(out_channels)
    }
  },

  forward = function(x) {
    x <- torch::nnf_relu(self$bn1(self$conv1(x)))

    if (self$use_residual) {
      skip <- x
      x <- torch::nnf_relu(self$bn2(self$conv2(x)))
      x <- self$bn3(self$conv3(x))
      x <- torch::nnf_relu(x + skip)
    }
    x
  }
)

.PipeMasterCNN2D <- torch::nn_module(
  "PipeMasterCNN2D",

  initialize = function(sfs_dims, n_targets, hp, mc_dropout = FALSE) {
    self$hp <- hp
    self$sfs_dims <- sfs_dims
    self$n_targets <- n_targets
    self$mc_dropout_mode <- mc_dropout
    n_blocks <- hp$n_blocks

    # Build conv blocks
    self$conv_blocks <- torch::nn_module_list()
    # Track which blocks get max pooling (first and middle)
    self$pool_after <- c(1L, n_blocks %/% 2L + 1L)

    in_ch <- 1L
    for (b in seq_len(n_blocks)) {
      filters <- as.integer(hp$base_filters * (2 ^ min(b - 1, 2)))
      ks <- as.integer(max(3L, hp$kernel_start - (b - 1) * 2))
      use_res <- isTRUE(hp$use_residual) && b > 1 && b < n_blocks

      self$conv_blocks$append(.ConvBlock2D(in_ch, filters, ks,
                                           use_residual = use_res))
      in_ch <- filters
    }

    # Dense head
    self$dense_layers <- torch::nn_module_list()
    self$dense_bns    <- torch::nn_module_list()
    self$dense_drops  <- torch::nn_module_list()
    prev_units <- in_ch
    for (d in seq_len(hp$n_dense)) {
      units <- as.integer(hp$dense_units / (2 ^ (d - 1)))
      self$dense_layers$append(torch::nn_linear(prev_units, units))
      self$dense_bns$append(torch::nn_batch_norm1d(units))
      if (hp$dropout > 0) {
        self$dense_drops$append(
          if (mc_dropout) .MCDropout(hp$dropout) else torch::nn_dropout(p = hp$dropout)
        )
      } else {
        self$dense_drops$append(NULL)
      }
      prev_units <- units
    }

    self$output_layer <- torch::nn_linear(prev_units, n_targets)
  },

  forward = function(x) {
    # x shape from R: [N, dim1, dim2, 1] — permute to [N, 1, dim1, dim2]
    if (x$dim() == 4L && x$size(4) == 1L)
      x <- x$permute(c(1L, 4L, 2L, 3L))

    for (b in seq_along(self$conv_blocks)) {
      x <- self$conv_blocks[[b]](x)
      if (b %in% self$pool_after)
        x <- torch::nnf_max_pool2d(x, kernel_size = 2L)
    }

    # Global average pooling: [N, C, H, W] → [N, C]
    x <- x$mean(dim = c(3L, 4L))

    # Dense head
    for (d in seq_along(self$dense_layers)) {
      x <- torch::nnf_relu(self$dense_bns[[d]](self$dense_layers[[d]](x)))
      if (!is.null(self$dense_drops[[d]]))
        x <- self$dense_drops[[d]](x)
    }

    self$output_layer(x)
  }
)

# ============================================================================
# Unified builder: replaces .build.nn() for torch backend
# ============================================================================

#' @keywords internal
.build.nn.torch <- function(hp, data, type, sfs.dims = NULL, device = "cpu") {

  .ncol2 <- function(x) if (inherits(x, "torch_tensor")) x$size(2) else ncol(x)
  n_targets <- .ncol2(data$Y_train)

  model <- switch(type,
    sumstat = , emulator = {
      n_features <- .ncol2(data$X_train)
      .PipeMasterResNet(n_features, n_targets, hp)
    },
    sfs1d = {
      .PipeMasterCNN1D(data$n_bins, n_targets, hp)
    },
    sfs2d = {
      .PipeMasterCNN2D(sfs.dims, n_targets, hp)
    },
    stop(sprintf("Unknown type: %s", type))
  )

  model$to(device = torch::torch_device(device))
}

# ============================================================================
# Torch-native predict: replaces predict(keras_model, X, verbose = 0L)
#
# Returns an R matrix (same as keras predict).
# ============================================================================

#' @keywords internal
.torch.predict <- function(model, X, device = NULL) {
  was_training <- model$training
  model$eval()
  on.exit(if (was_training) model$train(), add = TRUE)

  # Auto-detect model device if not specified
  if (is.null(device))
    device <- model$parameters[[1]]$device

  if (!inherits(X, "torch_tensor"))
    X <- torch::torch_tensor(X, dtype = torch::torch_float(), device = device)

  torch::with_no_grad({
    pred <- model(X)
    as.matrix(pred$cpu())
  })
}

# ============================================================================
# Huber loss function (matches keras loss_huber with configurable delta)
#
# R torch's nnf_smooth_l1_loss has no delta parameter, so we implement it
# directly:  L = 0.5*(x/d)^2  if |x| < d,  |x| - 0.5*d  otherwise
# where x = input - target, d = delta
# ============================================================================

#' @keywords internal
.torch.huber.loss <- function(input, target, delta = 1.0) {
  diff <- input - target
  abs_diff <- diff$abs()
  quadratic <- torch::torch_clamp(abs_diff, max = delta)
  linear <- abs_diff - quadratic
  loss <- 0.5 * quadratic$pow(2) / delta + delta * linear
  loss$mean()
}

# ============================================================================
# Save / load a torch model with its HP config and architecture metadata
#
# Saves: list(state_dict, hp, type, n_features/n_bins/sfs_dims, n_targets,
#             mc_dropout)
# ============================================================================

#' @keywords internal
.torch.save.model <- function(model, path, type, hp,
                              n_features = NULL, n_bins = NULL,
                              sfs_dims = NULL, n_targets = NULL) {
  meta <- list(
    state_dict = model$state_dict(),
    hp         = hp,
    type       = type,
    n_features = n_features,
    n_bins     = n_bins,
    sfs_dims   = sfs_dims,
    n_targets  = n_targets,
    mc_dropout = isTRUE(model$mc_dropout_mode)
  )
  torch::torch_save(meta, path)
}

#' @keywords internal
.torch.load.model <- function(path, device = "cpu") {
  meta <- torch::torch_load(path)

  hp      <- meta$hp
  type    <- meta$type
  mc_drop <- isTRUE(meta$mc_dropout)

  model <- switch(type,
    sumstat = , emulator = {
      .PipeMasterResNet(meta$n_features, meta$n_targets, hp,
                        mc_dropout = mc_drop)
    },
    sfs1d = {
      .PipeMasterCNN1D(meta$n_bins, meta$n_targets, hp,
                       mc_dropout = mc_drop)
    },
    sfs2d = {
      .PipeMasterCNN2D(meta$sfs_dims, meta$n_targets, hp,
                       mc_dropout = mc_drop)
    },
    stop(sprintf("Unknown model type: %s", type))
  )

  model$load_state_dict(meta$state_dict)
  model$to(device = torch::torch_device(device))
  model
}
