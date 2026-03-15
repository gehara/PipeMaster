#' Simulate summary statistics using scrm (SMC') coalescent engine
#'
#' Uses the vendored scrm C++ engine for fast coalescent simulation with
#' recombination, combined with PipeMaster's C summary statistic functions.
#' All computation stays in C/C++ memory — no R matrix overhead. This is
#' ~5x faster than msABC for WGS-scale loci (e.g. 100kb with recombination).
#'
#' @param model A PipeMaster model object (from main.menu or build functions).
#' @param nsim.blocks Number of simulation blocks per worker.
#' @param block.size Number of simulations per block.
#' @param mu.rates Mutation rate per base per generation (single numeric value).
#' @param rec.rates Recombination rate per base per generation (single numeric value).
#' @param use.alpha Logical or vector. If TRUE or c(TRUE, pop1, pop2, ...),
#'   uses exponential growth for specified populations.
#' @param skip.zns Logical. If TRUE (default), skip ZnS computation
#'   (O(segsites^2), very slow for large loci).
#' @param ncores Number of parallel worker processes.
#' @param path Output directory (default: current directory).
#' @param output.name Base name for output file (will be prefixed with SIMS_).
#' @param variable_samples Logical. If TRUE, allows variable sample sizes across
#'   loci but skips SFS computation. Note: variable sample sizes in scrm require
#'   separate C calls per sample-size group, which is slower. For RADseq data with
#'   variable samples, consider \code{sim.sumstats()} instead.
#'   Default is FALSE.
#' @param append.sims Logical. If TRUE, append to existing file.
#' @param verbose Logical. Print progress messages.
#' @return Invisibly returns the output file path. Writes results to
#'   SIMS_<output.name>.txt in the specified path.
#'
#' @details
#' Total simulations = nsim.blocks * block.size * ncores.
#'
#' For WGS-scale simulations (many long loci with recombination), this function
#' is significantly faster than \code{sim.sumstats()} because:
#' \itemize{
#'   \item scrm uses the SMC' approximation, which is O(n*L) instead of O(n*L*rho)
#'   \item Summary statistics are computed in C directly from scrm's in-memory data
#'   \item No R matrix allocation or text parsing overhead
#' }
#'
#' Requires fixed mutation and recombination rates (same for all loci).
#' All loci must have the same length.
#'
#' Multi-core parallelism uses independent Rscript worker processes
#' (same approach as \code{sim.sumstats()}).
#'
#' @export
sim.scrm.sumstats <- function(model, nsim.blocks, block.size,
                           mu.rates, rec.rates,
                           use.alpha = FALSE, skip.zns = TRUE,
                           ncores = 1, path = ".", output.name = "scrm",
                           variable_samples = FALSE,
                           append.sims = FALSE, verbose = TRUE) {

  if (!is.numeric(mu.rates) || length(mu.rates) != 1)
    stop("mu.rates must be a single numeric value")
  if (!is.numeric(rec.rates) || length(rec.rates) != 1)
    stop("rec.rates must be a single numeric value")

  npop <- as.integer(model$I[1, 3])
  pop_cols <- 4:(3 + npop)
  config <- as.integer(model$I[1, pop_cols])
  nsam <- sum(config)
  nloci <- nrow(model$loci)
  locus_lengths <- as.numeric(model$loci[, 2])
  uniform_len <- length(unique(locus_lengths)) == 1

  # Check per-population sample sizes across loci
  pop_sizes_mat <- matrix(as.numeric(model$I[, pop_cols]), ncol = npop)
  uniform_samples <- nrow(unique(pop_sizes_mat)) == 1

  if(!uniform_samples && !variable_samples) {
    stop("Per-population sample sizes vary across loci. SFS cannot be computed.\n",
         "  Either use optimize.sfs.model() to downsample to uniform sizes,\n",
         "  or set variable_samples = TRUE to skip SFS (summary statistics only).")
  }
  if(!uniform_samples && variable_samples) {
    stop("Variable sample sizes are not yet supported in the scrm engine.\n",
         "  The scrm C backend requires uniform sample sizes across loci.\n",
         "  Use sim.sumstats() with variable_samples = TRUE instead,\n",
         "  or use optimize.sfs.model() to downsample to uniform sizes.")
  }

  # Reference Ne for coalescent scaling (same as msABC.commander)
  Ne0 <- 100000
  ms_scalar <- 4 * Ne0

  # Build per-group scrm commands (one group per unique locus length)
  # For uniform lengths, this is just one group
  unique_lens <- sort(unique(locus_lengths))
  n_groups <- length(unique_lens)
  group_nloci <- as.integer(table(factor(locus_lengths, levels = unique_lens)))

  base_cmds <- character(n_groups)
  for (g in seq_along(unique_lens)) {
    gl <- unique_lens[g]
    theta_g <- ms_scalar * mu.rates * gl
    rho_g <- ms_scalar * rec.rates * gl
    cmd <- sprintf("%d %d -t %g -r %g %d",
                   nsam, group_nloci[g], theta_g, rho_g, as.integer(gl))
    if (npop > 1) {
      pop_str <- paste(c("-I", npop, config), collapse = " ")
      cmd <- paste(cmd, pop_str)
    }
    cmd <- paste(cmd, "-l", group_nloci[g])
    base_cmds[g] <- cmd
  }

  # Output setup
  abs_path <- normalizePath(path, mustWork = TRUE)
  outfile <- file.path(abs_path, paste0("SIMS_", output.name, ".txt"))

  total_sims <- nsim.blocks * block.size * ncores

  if (verbose) {
    if (uniform_len) {
      cat(sprintf("PipeMaster:: scrm engine: %d sims (%d blocks x %d x %d cores), %d loci x %d bp, %d pops\n",
                  total_sims, nsim.blocks, block.size, ncores,
                  nloci, as.integer(unique_lens[1]), npop))
    } else {
      cat(sprintf("PipeMaster:: scrm engine: %d sims (%d blocks x %d x %d cores), %d loci (%d length groups, %d-%d bp), %d pops\n",
                  total_sims, nsim.blocks, block.size, ncores,
                  nloci, n_groups, as.integer(min(unique_lens)), as.integer(max(unique_lens)), npop))
    }
  }

  # Write header (build column names deterministically, no sim needed)
  if (!append.sims || !file.exists(outfile)) {
    col_names <- .scrm.col.names(model, npop, config, use.alpha)
    writeLines(paste(col_names, collapse = "\t"), outfile)
  }

  ############### Multi-core path: independent worker processes
  if (ncores > 1) {

    # Clean up stale worker dirs and done files
    for (w in 1:ncores) {
      unlink(file.path(abs_path, paste0(".scrm_worker_", w)), recursive = TRUE)
      f <- file.path(abs_path, paste0(".scrm_worker_", w, ".done"))
      if (file.exists(f)) file.remove(f)
    }

    # Save parameters for workers
    save(model, nsim.blocks, block.size, use.alpha, skip.zns,
         mu.rates, rec.rates, output.name, variable_samples,
         file = file.path(abs_path, ".PM_scrm_worker_params.RData"))

    # Worker script
    worker_script <- paste(
      'args <- commandArgs(TRUE)',
      'worker_id <- as.integer(args[1])',
      'suppressMessages(library(PipeMaster))',
      sprintf('base_path <- "%s"', abs_path),
      'load(file.path(base_path, ".PM_scrm_worker_params.RData"))',
      'worker_dir <- file.path(base_path, paste0(".scrm_worker_", worker_id))',
      'dir.create(worker_dir, showWarnings = FALSE)',
      'sim.scrm.sumstats(model = model, nsim.blocks = nsim.blocks,',
      '               block.size = block.size, mu.rates = mu.rates,',
      '               rec.rates = rec.rates, use.alpha = use.alpha,',
      '               skip.zns = skip.zns, output.name = output.name,',
      '               path = worker_dir, variable_samples = variable_samples,',
      '               ncores = 1, append.sims = TRUE)',
      'write("done", file.path(base_path, paste0(".scrm_worker_", worker_id, ".done")))',
      'quit(save = "no")',
      sep = "\n")
    writeLines(worker_script, file.path(abs_path, ".PM_scrm_worker.R"))

    start_time <- Sys.time()
    for (w in 1:ncores) {
      system(paste("Rscript", file.path(abs_path, ".PM_scrm_worker.R"), w,
                   ">", file.path(abs_path, paste0(".scrm_worker_", w, ".log")), "2>&1"),
             wait = FALSE)
    }
    cat(sprintf("PipeMaster:: Launched %d worker processes\n", ncores))

    total_expected <- nsim.blocks * block.size * ncores
    prev_total_sims <- -1
    prev_done_count <- -1
    while (TRUE) {
      Sys.sleep(5)
      done_count <- sum(file.exists(file.path(abs_path,
                         paste0(".scrm_worker_", 1:ncores, ".done"))))

      # Count sims across all workers
      total_sims_done <- 0
      for (w in 1:ncores) {
        wf <- file.path(abs_path, paste0(".scrm_worker_", w),
                        paste0("SIMS_", output.name, ".txt"))
        if (file.exists(wf)) {
          n <- as.integer(system(paste("wc -l <", shQuote(wf)), intern = TRUE))
          if (!is.na(n) && n > 1) total_sims_done <- total_sims_done + (n - 1)  # subtract header
        }
      }

      if (total_sims_done != prev_total_sims || done_count != prev_done_count) {
        elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
        if (elapsed_h > 0.001 && total_sims_done > 0) {
          rate <- round(total_sims_done / elapsed_h)
          remaining <- round(max(0, (total_expected - total_sims_done) / rate), 2)
        } else {
          rate <- "..."
          remaining <- "..."
        }
        cat(sprintf("PipeMaster:: %d/%d sims (~%s sims/h) | ~%s h remaining | %d/%d workers done\n",
                    total_sims_done, total_expected, rate, remaining,
                    done_count, ncores))
        prev_total_sims <- total_sims_done
        prev_done_count <- done_count
      }

      if (done_count >= ncores) break
    }

    # Compile results from workers
    cat("PipeMaster:: Compiling results from workers\n")
    for (w in 1:ncores) {
      wf <- file.path(abs_path, paste0(".scrm_worker_", w),
                      paste0("SIMS_", output.name, ".txt"))
      if (file.exists(wf)) {
        worker_data <- readLines(wf)
        # Skip header line from each worker
        if (length(worker_data) > 1) {
          cat(paste(worker_data[-1], collapse = "\n"), "\n",
              file = outfile, append = TRUE, sep = "")
        }
      }
      unlink(file.path(abs_path, paste0(".scrm_worker_", w)), recursive = TRUE)
      f <- file.path(abs_path, paste0(".scrm_worker_", w, ".done"))
      if (file.exists(f)) file.remove(f)
    }
    file.remove(file.path(abs_path, ".PM_scrm_worker_params.RData"))
    file.remove(file.path(abs_path, ".PM_scrm_worker.R"))
    for (w in 1:ncores) {
      f <- file.path(abs_path, paste0(".scrm_worker_", w, ".log"))
      if (file.exists(f)) file.remove(f)
    }

    end_time <- Sys.time()
    elapsed_h <- as.numeric(difftime(end_time, start_time, units = "hours"))
    cat(sprintf("PipeMaster:: Done! %d simulations in %.3f hours (~%d sims/h)\n",
                total_expected, elapsed_h, round(total_expected / elapsed_h)))

  } else {
    ############### Single-core path
    start_time <- Sys.time()
    total_done <- 0

    for (j in 1:nsim.blocks) {
      block_results <- vector("list", block.size)
      for (i in 1:block.size) {
        block_results[[i]] <- .scrm.run.one(model, base_cmds, config, npop,
                                            use.alpha, skip.zns, mu.rates)
      }

      block_mat <- do.call(rbind, block_results)
      write.table(block_mat, file = outfile, append = TRUE, quote = FALSE,
                  row.names = FALSE, col.names = FALSE, sep = "\t")

      total_done <- total_done + block.size
      if (verbose) {
        elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
        nsim_total <- nsim.blocks * block.size
        rate <- if (elapsed_h > 0.001) round(total_done / elapsed_h) else "..."
        remaining <- if (elapsed_h > 0.001) {
          round((nsim_total - total_done) / (total_done / elapsed_h), 3)
        } else "..."
        cat(sprintf("PipeMaster:: %d/%d sims (~%s sims/h) | ~%s h remaining\n",
                    total_done, nsim_total, rate, remaining))
      }
    }

    elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
    if (verbose)
      cat(sprintf("PipeMaster:: Done! %d simulations in %.3f hours (~%d sims/h)\n",
                  total_done, elapsed_h, round(total_done / elapsed_h)))
  }

  invisible(outfile)
}


# Internal: run one scrm simulation and return named numeric vector
# base_cmds: character vector of scrm commands (one per length group)
.scrm.run.one <- function(model, base_cmds, config, npop,
                          use.alpha, skip.zns, mu.rates) {
  nloci <- nrow(model$loci)

  # Sample parameters from priors
  cmd_result <- msABC.commander(model, use.alpha = use.alpha, arg = 1)
  ms_command <- cmd_result[[1]]
  params <- cmd_result[[2]]

  # Extract demographic flags from msABC command
  demog_part <- ms_command

  # Remove the leading "nsam nloci" part
  demog_part <- sub("^\\s*\\d+\\s+\\d+\\s*", "", demog_part)

  # Remove the --frag-begin ... --frag-end part
  demog_part <- sub("\\s*--frag-begin.*$", "", demog_part)

  # Remove the -I npop n1 n2 ... part (already in base_cmds)
  demog_part <- sub("^\\s*-I\\s+\\d+(\\s+\\d+)+\\s*", "", demog_part)

  # Build full scrm commands (one per length group, same demog flags)
  scrm_cmds <- paste(base_cmds, demog_part)

  if (length(scrm_cmds) == 1L) {
    # Uniform length: use original single-command call
    result <- .Call("scrm_stats_call", scrm_cmds, config, as.integer(npop),
                    as.logical(skip.zns), PACKAGE = "PipeMaster")
  } else {
    # Variable lengths: use multi-command call with shared accumulators
    result <- .Call("scrm_stats_multi_call", scrm_cmds, config, as.integer(npop),
                    as.logical(skip.zns), as.integer(nloci),
                    PACKAGE = "PipeMaster")
  }

  # Build output row: params + stats + sfs
  par_vec <- as.numeric(params[2, ])
  names(par_vec) <- params[1, ]

  # Add mu rate info
  par_vec <- c(par_vec, mean.rate = mu.rates, sd.rate = 0)

  # Fold and name the SFS entries
  # C code already folds per-site to minor allele; trim to correct length
  sfs_vec <- result$sfs
  nsam <- sum(config)
  if (npop == 1) {
    # 1-pop: C output is length nsam-1 but only first floor(nsam/2) bins populated
    sfs_len <- floor(nsam / 2)
    sfs_vec <- sfs_vec[1:sfs_len]
    names(sfs_vec) <- paste0("sfs_fold_", seq(0, sfs_len - 1))
  } else {
    # Multi-pop joint SFS: already folded to minor allele in C
    # Use expand.grid naming convention (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(config, function(n) 0:n))
    names(sfs_vec) <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_vec, result$stats, sfs_vec)
}


# Internal: build column names deterministically (no simulation needed)
.scrm.col.names <- function(model, npop, config, use.alpha) {
  # Parameter names from model
  size_pars <- rbind(model$flags$n, model$flags$en$size)
  mig_pars <- rbind(model$flags$m, model$flags$em$size)
  time_pars <- rbind(model$flags$ej, model$flags$en$time, model$flags$em$time)

  par_names <- c(size_pars[, 1], time_pars[, 1])
  if (!is.null(mig_pars)) par_names <- c(par_names, mig_pars[, 1])
  par_names <- c(par_names, "mean.rate", "sd.rate")

  # Stat names: by-stat-type layout matching scrm_stats.cpp make_stat_names()
  nsam <- sum(config)
  main_stats <- c("segs", "pi", "thetaW", "tajd", "ZnS")
  main_names <- unlist(lapply(main_stats, function(s) {
    c(paste0(s, "_", 1:npop), s)
  }))
  pair_block <- character(0)
  if (npop > 1) {
    pair_block <- "Fst"
    for (i in 1:(npop - 1)) {
      for (j in (i + 1):npop) {
        pair_block <- c(pair_block,
                        paste0(c("shared", "private", "fixed", "Fst"), "_", i, "_", j))
      }
    }
  }
  hap_names <- c(unlist(lapply(1:npop, function(p) paste0(c("nhap", "Hd"), "_", p))),
                 "nhap", "Hd")
  base_stat_names <- c(main_names, pair_block, hap_names)
  prefixes <- c("s_mean_", "s_var_", "s_skew_", "s_kurt_")
  stat_names <- unlist(lapply(prefixes, function(p) paste0(p, base_stat_names)))

  # SFS names (folded, matching sim.sumstats() convention)
  if (npop == 1) {
    sfs_len <- floor(nsam / 2)
    sfs_names <- paste0("sfs_fold_", seq(0, sfs_len - 1))
  } else {
    # Joint SFS: expand.grid naming (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(config, function(n) 0:n))
    sfs_names <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_names, stat_names, sfs_names)
}
