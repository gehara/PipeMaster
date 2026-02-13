#!/usr/bin/env Rscript
#
# Neural Network Regression for Demographic Parameter Estimation
# Trains Keras NNs on coalescent simulation data (from sim.sumstat())
# and predicts parameters from observed summary statistics.
# Compares NN point estimates vs ABC rejection vs true values.
#
# Models: Vaquita2Epoch (1-pop bottleneck), OutOfAfrica_3G09 (2-pop with migration)
#
# Architecture: Residual blocks with batch normalization, Huber loss, ensemble of 5

suppressPackageStartupMessages({
  library(keras)
  library(tensorflow)
})

set.seed(42)
N_ENSEMBLE <- 5

# ============================================================================
# Load data
# ============================================================================

cat("Loading data...\n")
load("data_to_test/sim_sumstat_results.RData")     # sim_vaq, sim_ooa
load("data_to_test/observed_msABC_sumstats.RData")  # observed_msABC_Vaquita2Epoch, observed_msABC_OutOfAfrica_3G09
load("data_to_test/test_models.RData")              # model objects with true_params

# True params in PipeMaster scale (from model attributes)
extract_true_params <- function(model_obj) {
  tp <- attr(model_obj, "true_params")
  is_num <- sapply(tp, is.numeric)
  unlist(tp[is_num])
}
true_vaq <- extract_true_params(Vaquita2Epoch)
true_ooa <- extract_true_params(OutOfAfrica_3G09)

# ============================================================================
# ABC rejection (standardized Euclidean distance, 0.1% tolerance)
# ============================================================================

abc_reject <- function(sim_df, obs_vals, stat_cols, tol_frac) {
  sim_matrix <- as.matrix(sim_df[, stat_cols, drop = FALSE])
  obs_vector <- obs_vals[stat_cols]
  sim_matrix[!is.finite(sim_matrix)] <- 0
  obs_vector[!is.finite(obs_vector)] <- 0
  sim_means <- colMeans(sim_matrix)
  sim_sd <- apply(sim_matrix, 2, sd, na.rm = TRUE)
  sim_sd[sim_sd == 0] <- 1
  sim_scaled <- sweep(sim_matrix, 2, sim_means) / rep(sim_sd, each = nrow(sim_matrix))
  obs_scaled <- (obs_vector - sim_means) / sim_sd
  distances <- sqrt(rowSums(sweep(sim_scaled, 2, obs_scaled)^2))
  tol <- quantile(distances, tol_frac)
  which(distances <= tol)
}

run_abc <- function(sim_df, obs_mat, param_cols, stat_cols, model_name) {
  cat(sprintf("\n  ABC rejection (0.1%% tolerance = %d accepted)...\n",
              floor(nrow(sim_df) * 0.001)))

  obs_vals <- as.numeric(obs_mat[1, ])
  names(obs_vals) <- colnames(obs_mat)

  in_range <- character(0)
  for (col in stat_cols) {
    ov <- obs_vals[col]
    if (is.na(ov) || !is.finite(ov)) next
    sr <- range(sim_df[, col], na.rm = TRUE)
    if (ov >= sr[1] && ov <= sr[2]) in_range <- c(in_range, col)
  }
  cat(sprintf("  Using %d / %d in-range stats\n", length(in_range), length(stat_cols)))

  accepted <- abc_reject(sim_df, obs_vals, in_range, 0.001)
  cat(sprintf("  Accepted %d simulations\n", length(accepted)))

  abc_est <- sapply(param_cols, function(p) median(sim_df[accepted, p]))
  abc_est
}

# ============================================================================
# Helper: prepare data with feature augmentation
# ============================================================================

prepare_data <- function(sim_df, obs_mat, param_names, model_name) {

  cat(sprintf("\n=== %s ===\n", model_name))

  nuisance <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(sim_df), c(param_names, nuisance))
  target_cols <- setdiff(param_names, nuisance)

  targets <- as.matrix(sim_df[, target_cols])
  features_raw <- as.matrix(sim_df[, stat_cols])

  # Feature augmentation: add log1p(abs(features)) columns
  features_log <- log1p(abs(features_raw))
  features <- cbind(features_raw, features_log)

  cat(sprintf("  Simulations: %d | Raw features: %d | Augmented features: %d | Targets: %d\n",
              nrow(features), ncol(features_raw), ncol(features), ncol(targets)))
  cat(sprintf("  Target params: %s\n", paste(target_cols, collapse = ", ")))

  # Remove rows with NA/Inf in features
  bad_rows <- apply(features, 1, function(x) any(is.na(x) | is.infinite(x)))
  n_bad <- sum(bad_rows)
  if (n_bad > 0) {
    cat(sprintf("  Removing %d rows with NA/Inf in features\n", n_bad))
    features <- features[!bad_rows, ]
    targets  <- targets[!bad_rows, ]
  }

  # 80/10/10 split
  n <- nrow(features)
  idx <- sample(n)
  n_train <- floor(0.8 * n)
  n_val   <- floor(0.1 * n)

  train_idx <- idx[1:n_train]
  val_idx   <- idx[(n_train + 1):(n_train + n_val)]
  test_idx  <- idx[(n_train + n_val + 1):n]

  X_train <- features[train_idx, ]
  X_val   <- features[val_idx, ]
  X_test  <- features[test_idx, ]

  Y_train_raw <- targets[train_idx, ]
  Y_val_raw   <- targets[val_idx, ]
  Y_test_raw  <- targets[test_idx, ]

  # Standardize features (z-score from training set)
  feat_mean <- colMeans(X_train)
  feat_sd   <- apply(X_train, 2, sd)
  feat_sd[feat_sd == 0] <- 1

  scale_features <- function(X) t((t(X) - feat_mean) / feat_sd)

  X_train <- scale_features(X_train)
  X_val   <- scale_features(X_val)
  X_test  <- scale_features(X_test)

  # Prepare observed data: augment + scale
  obs_vec <- as.numeric(obs_mat[1, ])
  names(obs_vec) <- colnames(obs_mat)
  obs_raw <- obs_vec[stat_cols]
  obs_log <- log1p(abs(obs_raw))
  obs_aug <- c(obs_raw, obs_log)
  X_obs <- matrix((obs_aug - feat_mean) / feat_sd, nrow = 1)
  X_obs[is.na(X_obs) | is.nan(X_obs)] <- 0

  # Log-transform targets, then standardize
  Y_train_log <- log(Y_train_raw)
  Y_val_log   <- log(Y_val_raw)
  Y_test_log  <- log(Y_test_raw)

  target_mean <- colMeans(Y_train_log)
  target_sd   <- apply(Y_train_log, 2, sd)
  target_sd[target_sd == 0] <- 1

  scale_targets <- function(Y) t((t(Y) - target_mean) / target_sd)

  Y_train <- scale_targets(Y_train_log)
  Y_val   <- scale_targets(Y_val_log)
  Y_test  <- scale_targets(Y_test_log)

  cat(sprintf("  Train: %d | Val: %d | Test: %d\n",
              nrow(X_train), nrow(X_val), nrow(X_test)))

  list(
    X_train = X_train, X_val = X_val, X_test = X_test, X_obs = X_obs,
    Y_train = Y_train, Y_val = Y_val, Y_test = Y_test,
    Y_test_raw = Y_test_raw,
    target_cols = target_cols, stat_cols = stat_cols,
    target_mean = target_mean, target_sd = target_sd,
    feat_mean = feat_mean, feat_sd = feat_sd
  )
}

# ============================================================================
# Helper: build residual model (functional API)
# ============================================================================

build_model <- function(n_features, n_targets) {

  l2_reg <- regularizer_l2(1e-4)

  # Residual block: Dense(relu)->BN->Dense->BN + skip -> relu
  res_block <- function(x, units) {
    skip <- x
    x <- x %>%
      layer_dense(units = units, activation = "relu",
                  kernel_regularizer = l2_reg) %>%
      layer_batch_normalization() %>%
      layer_dense(units = units, activation = "linear",
                  kernel_regularizer = l2_reg) %>%
      layer_batch_normalization()
    x <- layer_add(list(x, skip))
    x <- layer_activation(x, activation = "relu")
    x
  }

  input <- layer_input(shape = n_features)

  # Projection to 256
  x <- input %>%
    layer_dense(units = 256, activation = "relu",
                kernel_regularizer = l2_reg) %>%
    layer_batch_normalization()

  # Two residual blocks at 256
  x <- res_block(x, 256)
  x <- res_block(x, 256)

  # Bottleneck to 128
  x <- x %>%
    layer_dense(units = 128, activation = "relu",
                kernel_regularizer = l2_reg) %>%
    layer_batch_normalization()

  # One residual block at 128
  x <- res_block(x, 128)

  # Final layers
  x <- x %>%
    layer_dense(units = 64, activation = "relu",
                kernel_regularizer = l2_reg) %>%
    layer_batch_normalization()

  output <- x %>%
    layer_dense(units = n_targets, activation = "linear")

  model <- keras_model(input, output)

  model %>% compile(
    loss = loss_huber(delta = 1.0),
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = list("mae")
  )

  model
}

# ============================================================================
# Helper: train ensemble of models
# ============================================================================

train_ensemble <- function(data, model_name, n_models = N_ENSEMBLE) {

  n_features <- ncol(data$X_train)
  n_targets  <- ncol(data$Y_train)

  cat(sprintf("  Architecture: Input(%d) -> 256[res x2] -> 128[res x1] -> 64 -> %d\n",
              n_features, n_targets))
  cat(sprintf("  Loss: Huber | Ensemble: %d models\n", n_models))

  models <- list()

  for (i in seq_len(n_models)) {
    cat(sprintf("\n  --- Training model %d/%d (%s) ---\n", i, n_models, model_name))

    # Set different seed for each model (affects weight init)
    tensorflow::tf$random$set_seed(as.integer(42 + i * 7))

    model <- build_model(n_features, n_targets)

    callbacks <- list(
      callback_early_stopping(
        monitor = "val_loss",
        patience = 30,
        restore_best_weights = TRUE
      ),
      callback_reduce_lr_on_plateau(
        monitor = "val_loss",
        patience = 15,
        factor = 0.5,
        min_lr = 1e-6,
        verbose = 1
      )
    )

    history <- model %>% fit(
      x = data$X_train,
      y = data$Y_train,
      validation_data = list(data$X_val, data$Y_val),
      epochs = 1000,
      batch_size = 512,
      callbacks = callbacks,
      verbose = 0
    )

    h <- history$metrics
    n_epochs <- length(h$loss)
    cat(sprintf("  Stopped at epoch %d | Train loss: %.6f | Val loss: %.6f\n",
                n_epochs, h$loss[n_epochs], h$val_loss[n_epochs]))

    models[[i]] <- model
  }

  models
}

# ============================================================================
# Helper: predict with ensemble (average)
# ============================================================================

predict_ensemble <- function(models, X, data) {
  # Returns: list(mean = matrix, sd = matrix) in original scale
  preds_scaled <- lapply(models, function(m) predict(m, X, verbose = 0))
  preds_array <- array(unlist(preds_scaled),
                       dim = c(nrow(X), ncol(data$Y_train), length(models)))

  # Inverse transform each: scaled -> log -> original
  inverse_transform <- function(Y_scaled) {
    Y_log <- t(t(Y_scaled) * data$target_sd + data$target_mean)
    exp(Y_log)
  }

  preds_orig <- lapply(seq_along(models), function(i) {
    inverse_transform(preds_scaled[[i]])
  })
  preds_orig_array <- array(unlist(preds_orig),
                            dim = c(nrow(X), ncol(data$Y_train), length(models)))

  # Mean and SD across models
  ens_mean <- apply(preds_orig_array, c(1, 2), mean)
  ens_sd   <- apply(preds_orig_array, c(1, 2), sd)

  colnames(ens_mean) <- data$target_cols
  colnames(ens_sd)   <- data$target_cols

  list(mean = ens_mean, sd = ens_sd)
}

# ============================================================================
# Helper: evaluate ensemble + ABC and compare
# ============================================================================

evaluate_and_compare <- function(models, data, sim_df, obs_mat, true_params, model_name) {

  # --- Ensemble predictions on test set ---
  ens_test <- predict_ensemble(models, data$X_test, data)

  cat(sprintf("\n  --- %s: Test Set Evaluation (ensemble of %d) ---\n",
              model_name, length(models)))
  cat(sprintf("  %-15s %8s %12s\n", "Parameter", "R\u00b2", "RMSE"))
  cat(sprintf("  %s\n", paste(rep("-", 37), collapse = "")))

  for (j in seq_along(data$target_cols)) {
    y_true <- data$Y_test_raw[, j]
    y_pred <- ens_test$mean[, j]
    ss_res <- sum((y_true - y_pred)^2)
    ss_tot <- sum((y_true - mean(y_true))^2)
    r2 <- 1 - ss_res / ss_tot
    rmse <- sqrt(mean((y_true - y_pred)^2))
    cat(sprintf("  %-15s %8.4f %12.2f\n", data$target_cols[j], r2, rmse))
  }

  # --- Ensemble prediction on observed ---
  ens_obs <- predict_ensemble(models, data$X_obs, data)
  nn_est <- as.numeric(ens_obs$mean)
  nn_sd  <- as.numeric(ens_obs$sd)
  names(nn_est) <- data$target_cols
  names(nn_sd)  <- data$target_cols

  cat(sprintf("\n  Ensemble prediction uncertainty (SD across %d models):\n", length(models)))
  for (p in data$target_cols) {
    cat(sprintf("    %-15s  mean: %12.4f  sd: %12.4f  cv: %.1f%%\n",
                p, nn_est[p], nn_sd[p],
                ifelse(nn_est[p] != 0, abs(nn_sd[p] / nn_est[p]) * 100, 0)))
  }

  # --- ABC rejection ---
  abc_est <- run_abc(sim_df, obs_mat, data$target_cols, data$stat_cols, model_name)

  # --- Comparison table ---
  cat(sprintf("\n  --- %s: True vs ABC vs NN (ensemble) ---\n", model_name))
  cat(sprintf("  %-15s %12s %12s %12s %10s %10s\n",
              "Parameter", "True", "ABC", "NN", "ABC err%", "NN err%"))
  cat(sprintf("  %s\n", paste(rep("-", 73), collapse = "")))

  for (p in data$target_cols) {
    tv <- if (p %in% names(true_params)) true_params[p] else NA
    av <- abc_est[p]
    nv <- nn_est[p]

    abc_err <- if (!is.na(tv) && tv != 0) abs(av - tv) / abs(tv) * 100 else NA
    nn_err  <- if (!is.na(tv) && tv != 0) abs(nv - tv) / abs(tv) * 100 else NA

    if (!is.na(tv) && abs(tv) < 1) {
      cat(sprintf("  %-15s %12.4f %12.4f %12.4f %9.1f%% %9.1f%%\n",
                  p, tv, av, nv, abc_err, nn_err))
    } else {
      cat(sprintf("  %-15s %12.1f %12.1f %12.1f %9.1f%% %9.1f%%\n",
                  p, tv, av, nv, abc_err, nn_err))
    }
  }

  list(nn = nn_est, nn_sd = nn_sd, abc = abc_est)
}

# ============================================================================
# Run for Vaquita2Epoch
# ============================================================================

vaq_param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1", "mean.rate", "sd.rate")

data_vaq <- prepare_data(
  sim_df = sim_vaq,
  obs_mat = observed_msABC_Vaquita2Epoch,
  param_names = vaq_param_cols,
  model_name = "Vaquita2Epoch"
)

models_vaq <- train_ensemble(data_vaq, "Vaquita2Epoch")
res_vaq <- evaluate_and_compare(
  models_vaq, data_vaq, sim_vaq,
  observed_msABC_Vaquita2Epoch, true_vaq, "Vaquita2Epoch"
)

# ============================================================================
# Run for OutOfAfrica_3G09
# ============================================================================

ooa_param_cols <- c("Ne0.pop1", "Ne0.pop2", "Ne1.pop2", "Ne1.pop1", "join1",
                    "t.Ne1.pop2", "t.Ne1.pop1", "mig0.1_2", "mig0.2_1",
                    "mean.rate", "sd.rate")

data_ooa <- prepare_data(
  sim_df = sim_ooa,
  obs_mat = observed_msABC_OutOfAfrica_3G09,
  param_names = ooa_param_cols,
  model_name = "OutOfAfrica_3G09"
)

models_ooa <- train_ensemble(data_ooa, "OutOfAfrica_3G09")
res_ooa <- evaluate_and_compare(
  models_ooa, data_ooa, sim_ooa,
  observed_msABC_OutOfAfrica_3G09, true_ooa, "OutOfAfrica_3G09"
)

# ============================================================================
# Final summary
# ============================================================================

cat("\n\n========================================================\n")
cat("  FINAL SUMMARY: ABC vs NN ensemble (which is closer to true?)\n")
cat("========================================================\n")

summarize_winner <- function(res, true_params, target_cols, model_name) {
  cat(sprintf("\n%s:\n", model_name))
  abc_wins <- 0; nn_wins <- 0
  for (p in target_cols) {
    tv <- if (p %in% names(true_params)) true_params[p] else next
    abc_err <- abs(res$abc[p] - tv) / abs(tv) * 100
    nn_err  <- abs(res$nn[p] - tv) / abs(tv) * 100
    winner <- if (abc_err < nn_err) "ABC" else "NN"
    if (abc_err < nn_err) abc_wins <- abc_wins + 1 else nn_wins <- nn_wins + 1
    cat(sprintf("  %-15s  ABC: %6.1f%%  NN: %6.1f%%  -> %s\n", p, abc_err, nn_err, winner))
  }
  cat(sprintf("  Score: ABC %d - NN %d\n", abc_wins, nn_wins))
}

summarize_winner(res_vaq, true_vaq, data_vaq$target_cols, "Vaquita2Epoch")
summarize_winner(res_ooa, true_ooa, data_ooa$target_cols, "OutOfAfrica_3G09")

cat("\nDone.\n")
