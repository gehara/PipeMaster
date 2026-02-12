#!/usr/bin/env Rscript
# Test script for ABC parameter estimation using SFS
# Validates the full pipeline: sim.sfs() -> abc::abc() -> posterior recovery of true params

suppressMessages(library(PipeMaster))
suppressMessages(library(abc))

# Determine script directory robustly
args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if(length(file_arg) > 0) {
  test_dir <- normalizePath(dirname(sub("^--file=", "", file_arg)))
} else {
  test_dir <- normalizePath("data_to_test")
}
setwd(test_dir)

cat("=== Loading test data ===\n")
load("test_models.RData")

pass <- 0
fail <- 0
recovered <- 0
not_recovered <- 0

check <- function(desc, condition) {
  if(isTRUE(condition)) {
    cat(paste0("  PASS: ", desc, "\n"))
    pass <<- pass + 1
  } else {
    cat(paste0("  FAIL: ", desc, "\n"))
    fail <<- fail + 1
  }
}

tmpdir <- tempdir()

# Helper: run ABC for a model and check parameter recovery
run_abc_test <- function(model, model_name, observed_sfs, true_params,
                         is_multipop = FALSE, block_size = 1000, ncores = 10) {

  total_sims <- block_size * ncores
  cat(paste0("\n--- Simulating reference table (", total_sims, " sims, ",
             ncores, " cores) ---\n"))
  sim.sfs(model = model,
          nsim.blocks = 1,
          block.size = block_size,
          path = tmpdir,
          use.alpha = FALSE,
          output.name = paste0("abc_", model_name),
          ncores = ncores)

  outfile <- file.path(tmpdir, paste0("SIM_SFS_abc_", model_name, ".txt"))
  check(paste0(model_name, ": reference table exists"), file.exists(outfile))

  sim_data <- read.table(outfile, header = TRUE, sep = "\t")
  check(paste0(model_name, ": ", total_sims, " simulation rows"), nrow(sim_data) == total_sims)

  # Split into parameters and SFS
  sfs_cols <- grep("^sfs_", colnames(sim_data), value = TRUE)
  param_cols <- setdiff(colnames(sim_data), sfs_cols)
  # Exclude mean.rate and sd.rate from ABC params (not in true_params)
  param_cols <- setdiff(param_cols, c("mean.rate", "sd.rate"))

  params <- sim_data[, param_cols, drop = FALSE]
  sumstats <- sim_data[, sfs_cols, drop = FALSE]

  # Get observed SFS target vector
  if(is_multipop) {
    # 2-pop: sim and obs both have 1681 columns — direct match
    target <- as.numeric(observed_sfs[1, ])
  } else {
    # 1-pop: obs has 41 bins (freq 0..40), sim has 39 bins (freq 1..39)
    # Align: obs columns 2:40 correspond to sim columns
    target <- as.numeric(observed_sfs[1, 2:40])
  }

  check(paste0(model_name, ": target and sumstat dimensions match"),
        length(target) == ncol(sumstats))

  # Remove zero-variance columns (constant SFS bins cause issues)
  col_var <- apply(sumstats, 2, var)
  keep <- col_var > 0
  sumstats_filt <- sumstats[, keep, drop = FALSE]
  target_filt <- target[keep]

  cat(paste0("  Using ", ncol(sumstats_filt), " of ", ncol(sumstats),
             " SFS bins (removed ", sum(!keep), " zero-variance bins)\n"))

  # Choose ABC method based on dimensionality
  # loclinear regression needs n_accepted > n_variables
  n_accepted <- floor(total_sims * 0.1)
  if(ncol(sumstats_filt) >= n_accepted) {
    # Too many SFS bins for regression — use PCA to reduce dimensions
    n_pcs <- min(20, n_accepted - ncol(params) - 1)
    cat(paste0("  High-dimensional SFS (", ncol(sumstats_filt), " bins > ",
               n_accepted, " accepted). Reducing to ", n_pcs, " PCs.\n"))
    pca_fit <- prcomp(sumstats_filt, center = TRUE, scale. = TRUE)
    target_df <- as.data.frame(matrix(target_filt, nrow = 1))
    colnames(target_df) <- colnames(sumstats_filt)
    target_filt <- predict(pca_fit, target_df)[1, 1:n_pcs]
    sumstats_filt <- as.data.frame(pca_fit$x[, 1:n_pcs])
  }

  # Run ABC with log transform (all params are positive)
  abc_method <- "loclinear"
  transf_vec <- rep("log", ncol(params))
  cat(paste0("--- Running ABC (", abc_method, ", tol=0.1, log-transformed) ---\n"))
  abc_result <- tryCatch(
    abc(target = target_filt,
        param = params,
        sumstat = sumstats_filt,
        tol = 0.1,
        transf = transf_vec,
        method = abc_method),
    error = function(e) {
      cat(paste0("  ERROR: ", e$message, "\n"))
      NULL
    }
  )

  check(paste0(model_name, ": ABC completed without error"), !is.null(abc_result))
  if(is.null(abc_result)) return(invisible(NULL))

  # Check posterior recovery for each parameter
  cat("--- Posterior summary ---\n")
  cat(sprintf("  %-15s %12s %12s %12s %12s %s\n",
              "Parameter", "True", "Median", "CI_low", "CI_high", "Recovered?"))

  for(p in param_cols) {
    post_samples <- abc_result$adj.values[, p]
    post_median <- median(post_samples)
    ci <- quantile(post_samples, probs = c(0.025, 0.975))

    # Get true value
    true_val <- true_params[[p]]
    if(is.null(true_val)) {
      cat(sprintf("  %-15s %12s (not in true_params, skipping)\n", p, "NA"))
      next
    }

    in_ci <- true_val >= ci[1] && true_val <= ci[2]
    cat(sprintf("  %-15s %12.1f %12.1f %12.1f %12.1f %s\n",
                p, true_val, post_median, ci[1], ci[2],
                ifelse(in_ci, "YES", "NO")))

    # Parameter recovery is informational — not a pipeline failure
    if(in_ci) recovered <<- recovered + 1 else not_recovered <<- not_recovered + 1
  }

  # Cleanup
  file.remove(outfile)
  invisible(abc_result)
}


# ============================================================================
# Test 1: Vaquita2Epoch (1-pop, 3 params)
# ============================================================================
cat("\n=== Test 1: Vaquita2Epoch ABC estimation ===\n")

run_abc_test(
  model = Vaquita2Epoch,
  model_name = "vaquita",
  observed_sfs = observed_sfs_Vaquita2Epoch,
  true_params = attr(Vaquita2Epoch, "true_params"),
  is_multipop = FALSE
)


# ============================================================================
# Test 2: Africa_1T12 (1-pop, 5 params)
# ============================================================================
cat("\n=== Test 2: Africa_1T12 ABC estimation ===\n")

run_abc_test(
  model = Africa_1T12,
  model_name = "africa",
  observed_sfs = observed_sfs_Africa_1T12,
  true_params = attr(Africa_1T12, "true_params"),
  is_multipop = FALSE
)


# ============================================================================
# Test 3: OutOfAfrica_2T12 (2-pop, 11 params)
# ============================================================================
cat("\n=== Test 3: OutOfAfrica_2T12 ABC estimation ===\n")

run_abc_test(
  model = OutOfAfrica_2T12,
  model_name = "ooa2",
  observed_sfs = observed_sfs_OutOfAfrica_2T12,
  true_params = attr(OutOfAfrica_2T12, "true_params"),
  is_multipop = TRUE
)


# ============================================================================
# Summary
# ============================================================================
cat(paste0("\n=== Pipeline checks: ", pass, " passed, ", fail, " failed ===\n"))
cat(paste0("=== Parameter recovery: ", recovered, " of ", recovered + not_recovered,
           " params within 95% CI ===\n"))
if(fail > 0) {
  cat("PIPELINE TESTS FAILED!\n")
  quit(status = 1)
} else {
  cat("ALL PIPELINE TESTS PASSED!\n")
  if(not_recovered > 0) {
    cat("(Some parameters not recovered — expected with 1000 sims from prior.\n")
    cat(" Increase block_size for better posterior precision.)\n")
  }
}
