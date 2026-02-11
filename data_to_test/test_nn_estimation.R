#!/usr/bin/env Rscript
#
# Neural Network Regression for Demographic Parameter Estimation
# Trains Keras NNs on coalescent simulation data (from sim.sumstat())
# and predicts parameters from observed summary statistics.
# Compares NN point estimates vs ABC rejection vs true values.
#
# Models: Vaquita2Epoch (1-pop bottleneck), OutOfAfrica_3G09 (2-pop with migration)

suppressPackageStartupMessages({
  library(keras)
  library(tensorflow)
})

set.seed(42)

# ============================================================================
# Load data
# ============================================================================

cat("Loading data...\n")
load("data_to_test/sim_sumstat_results.RData")     # sim_vaq, sim_ooa
load("data_to_test/observed_msABC_sumstats.RData")  # observed_msABC_Vaquita2Epoch, observed_msABC_OutOfAfrica_3G09
load("data_to_test/test_models.RData")              # model objects with true_params

# True params in PipeMaster scale (from model attributes)
# Filter to numeric entries only (exclude 'model' string)
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
  # Replace non-finite with 0
  sim_matrix[!is.finite(sim_matrix)] <- 0
  obs_vector[!is.finite(obs_vector)] <- 0
  # Standardize
  sim_means <- colMeans(sim_matrix)
  sim_sd <- apply(sim_matrix, 2, sd, na.rm = TRUE)
  sim_sd[sim_sd == 0] <- 1
  sim_scaled <- sweep(sim_matrix, 2, sim_means) / rep(sim_sd, each = nrow(sim_matrix))
  obs_scaled <- (obs_vector - sim_means) / sim_sd
  # Euclidean distance
  distances <- sqrt(rowSums(sweep(sim_scaled, 2, obs_scaled)^2))
  tol <- quantile(distances, tol_frac)
  which(distances <= tol)
}

run_abc <- function(sim_df, obs_mat, param_cols, stat_cols, model_name) {
  cat(sprintf("\n  ABC rejection (0.1%% tolerance = %d accepted)...\n",
              floor(nrow(sim_df) * 0.001)))

  obs_vals <- as.numeric(obs_mat[1, ])
  names(obs_vals) <- colnames(obs_mat)

  # Find in-range stats (exclude stats where observed is outside sim range)
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

  # Posterior medians for target params
  abc_est <- sapply(param_cols, function(p) median(sim_df[accepted, p]))
  abc_est
}

# ============================================================================
# Helper: prepare data for a model
# ============================================================================

prepare_data <- function(sim_df, obs_mat, param_names, model_name) {

  cat(sprintf("\n=== %s ===\n", model_name))

  # Separate targets and features
  nuisance <- c("mean.rate", "sd.rate")
  stat_cols <- setdiff(colnames(sim_df), c(param_names, nuisance))
  target_cols <- setdiff(param_names, nuisance)

  targets <- as.matrix(sim_df[, target_cols])
  features <- as.matrix(sim_df[, stat_cols])

  cat(sprintf("  Simulations: %d | Features: %d | Targets: %d\n",
              nrow(features), ncol(features), ncol(targets)))
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

  # Prepare observed data with same scaling
  obs_vec <- as.numeric(obs_mat[1, ])
  names(obs_vec) <- colnames(obs_mat)
  obs_vec <- obs_vec[stat_cols]
  X_obs <- matrix((obs_vec - feat_mean) / feat_sd, nrow = 1)
  # Impute NaN/NA in observed data with 0 (= training mean after z-scoring)
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
# Helper: build and train model
# ============================================================================

build_and_train <- function(data, model_name) {

  n_features <- ncol(data$X_train)
  n_targets  <- ncol(data$Y_train)

  cat(sprintf("  Building model: %d -> 128 -> 64 -> 32 -> %d\n", n_features, n_targets))

  model <- keras_model_sequential() %>%
    layer_dense(units = 128, activation = "relu", input_shape = n_features) %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 64, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = 32, activation = "relu") %>%
    layer_batch_normalization() %>%
    layer_dropout(rate = 0.1) %>%
    layer_dense(units = n_targets, activation = "linear")

  model %>% compile(
    loss = "mse",
    optimizer = optimizer_adam(learning_rate = 0.001),
    metrics = list("mae")
  )

  callbacks <- list(
    callback_early_stopping(
      monitor = "val_loss",
      patience = 20,
      restore_best_weights = TRUE
    ),
    callback_reduce_lr_on_plateau(
      monitor = "val_loss",
      patience = 10,
      factor = 0.5,
      verbose = 1
    )
  )

  cat("  Training...\n")
  history <- model %>% fit(
    x = data$X_train,
    y = data$Y_train,
    validation_data = list(data$X_val, data$Y_val),
    epochs = 500,
    batch_size = 256,
    callbacks = callbacks,
    verbose = 0
  )

  h <- history$metrics
  n_epochs <- length(h$loss)
  cat(sprintf("  Stopped at epoch %d | Train loss: %.6f | Val loss: %.6f\n",
              n_epochs, h$loss[n_epochs], h$val_loss[n_epochs]))

  list(model = model, history = history)
}

# ============================================================================
# Helper: evaluate NN + ABC and compare
# ============================================================================

evaluate_and_compare <- function(model, data, sim_df, obs_mat, true_params, model_name) {

  # --- NN predictions ---
  inverse_transform <- function(Y_scaled) {
    Y_log <- t(t(Y_scaled) * data$target_sd + data$target_mean)
    exp(Y_log)
  }

  # Test set R² and RMSE
  Y_test_pred_scaled <- predict(model, data$X_test, verbose = 0)
  Y_test_pred <- inverse_transform(Y_test_pred_scaled)

  cat(sprintf("\n  --- %s: Test Set Evaluation ---\n", model_name))
  cat(sprintf("  %-15s %8s %12s\n", "Parameter", "R²", "RMSE"))
  cat(sprintf("  %s\n", paste(rep("-", 37), collapse = "")))

  for (j in seq_along(data$target_cols)) {
    y_true <- data$Y_test_raw[, j]
    y_pred <- Y_test_pred[, j]
    ss_res <- sum((y_true - y_pred)^2)
    ss_tot <- sum((y_true - mean(y_true))^2)
    r2 <- 1 - ss_res / ss_tot
    rmse <- sqrt(mean((y_true - y_pred)^2))
    cat(sprintf("  %-15s %8.4f %12.2f\n", data$target_cols[j], r2, rmse))
  }

  # NN on observed
  Y_obs_pred_scaled <- predict(model, data$X_obs, verbose = 0)
  Y_obs_pred <- inverse_transform(Y_obs_pred_scaled)
  nn_est <- as.numeric(Y_obs_pred)
  names(nn_est) <- data$target_cols

  # --- ABC rejection ---
  abc_est <- run_abc(sim_df, obs_mat, data$target_cols, data$stat_cols, model_name)

  # --- Comparison table ---
  cat(sprintf("\n  --- %s: True vs ABC vs NN ---\n", model_name))
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

  list(nn = nn_est, abc = abc_est)
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

fit_vaq <- build_and_train(data_vaq, "Vaquita2Epoch")
res_vaq <- evaluate_and_compare(
  fit_vaq$model, data_vaq, sim_vaq,
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

fit_ooa <- build_and_train(data_ooa, "OutOfAfrica_3G09")
res_ooa <- evaluate_and_compare(
  fit_ooa$model, data_ooa, sim_ooa,
  observed_msABC_OutOfAfrica_3G09, true_ooa, "OutOfAfrica_3G09"
)

# ============================================================================
# Final summary
# ============================================================================

cat("\n\n========================================================\n")
cat("  FINAL SUMMARY: ABC vs NN (which is closer to true?)\n")
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
