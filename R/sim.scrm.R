#' Simulate summary statistics using scrm (SMC') coalescent engine
#'
#' Uses the vendored scrm C++ engine for fast coalescent simulation with
#' recombination, combined with PipeMaster's C summary statistic functions.
#' All computation stays in C/C++ memory — no R matrix overhead. This is
#' ~5x faster than msABC for WGS-scale loci (e.g. 100kb with recombination).
#'
#' @param model A PipeMaster model object (from main.menu or build functions).
#' @param nsims Total number of simulations to run. The actual number may be
#'   slightly higher due to rounding to batch.size * ncores.
#' @param batch.size Number of simulations per batch (controls R overhead).
#'   Default is 32.
#' @param mu.rates Mutation rate per base per generation. Either:
#'   (a) A single numeric value: applied uniformly to all loci, all sims (back-compat default).
#'   (b) A list with distribution spec for per-locus per-sim sampling. The
#'       distribution must be \code{"lognormal"} (the only currently supported
#'       shape); \code{median} is fixed across sims; the dispersion parameter
#'       \code{sigma_log} can be either FIXED or SAMPLED per sim:
#'       \itemize{
#'         \item Fixed:  \code{list(distribution = "lognormal", median = 5.83e-9, sigma_log = 0.3)}
#'         \item Sampled: \code{list(distribution = "lognormal", median = 5.83e-9, sigma_log_range = c(0.05, 0.5))}
#'           — per sim, \code{sigma_log ~ Uniform(lo, hi)}; the realized value
#'             is recorded in the reftable as \code{sigma_log_mu} so it can be
#'             treated as a learnable nuisance parameter by \code{tune.nn()}.
#'       }
#'       Per-locus mu heterogeneity is essentially free (mu only affects
#'       post-hoc SegSites mutation placement, not the ARG). Output reftable's
#'       \code{mean.rate}/\code{sd.rate} record the realized per-sim mean/sd
#'       of the sampled per-locus rates.
#' @param rec.rates Recombination rate per base per generation. Same shape as
#'   \code{mu.rates}: scalar OR \code{list(distribution="lognormal", median,
#'   sigma_log)} for fixed dispersion OR \code{list(..., sigma_log_range =
#'   c(lo, hi))} for sampled dispersion. Realized \code{sigma_log_rec} is
#'   recorded in the reftable.
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
#' The actual number of simulations is \code{ceiling(nsims / (ncores * batch.size)) * batch.size * ncores},
#' which may slightly exceed \code{nsims}.
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
sim.scrm.sumstats <- function(model, nsims, batch.size = 32,
                           mu.rates, rec.rates,
                           skip.zns = TRUE,
                           ncores = 1, path = ".", output.name = "scrm",
                           variable_samples = FALSE,
                           append.sims = FALSE, verbose = TRUE,
                           .parent.pid.file = NULL) {

  if (is.null(model$use.alpha))
    stop("model$use.alpha is missing. Set it on the model (e.g. model$use.alpha <- FALSE).")

  nsim.blocks <- ceiling(nsims / (ncores * batch.size))

  # Internal helper: validate a rates spec (scalar OR distribution list with
  # either fixed `sigma_log` or sampled `sigma_log_range = c(lo, hi)` prior).
  .validate_rate_spec <- function(spec, name) {
    if (!is.list(spec) || is.null(spec$distribution)) {
      if (!is.numeric(spec) || length(spec) != 1)
        stop(sprintf("%s must be a single numeric value OR a list with ",
                     "$distribution + $median + (sigma_log XOR sigma_log_range)", name))
      return(invisible(NULL))
    }
    if (!(spec$distribution %in% c("lognormal")))
      stop(sprintf("%s$distribution: only 'lognormal' is currently supported", name))
    if (is.null(spec$median) || !is.numeric(spec$median) ||
        length(spec$median) != 1 || spec$median <= 0)
      stop(sprintf("%s$median must be a positive numeric scalar", name))
    # Use exact name match — R's $ does partial-prefix matching, which would
    # treat sigma_log as matching sigma_log_range.
    has_fixed <- "sigma_log"       %in% names(spec)
    has_range <- "sigma_log_range" %in% names(spec)
    if (!has_fixed && !has_range)
      stop(sprintf("%s list must specify sigma_log (fixed) OR sigma_log_range (sampled)", name))
    if (has_fixed && has_range)
      stop(sprintf("%s: provide sigma_log OR sigma_log_range, not both", name))
    if (has_fixed) {
      v <- spec[["sigma_log"]]
      if (!is.numeric(v) || length(v) != 1 || v < 0)
        stop(sprintf("%s$sigma_log must be a non-negative numeric scalar", name))
    } else {
      r <- spec[["sigma_log_range"]]
      if (!is.numeric(r) || length(r) != 2 || r[1] < 0 || r[2] < r[1])
        stop(sprintf("%s$sigma_log_range must be c(lo, hi) with 0 <= lo <= hi", name))
    }
    invisible(NULL)
  }
  .validate_rate_spec(mu.rates,  "mu.rates")
  .validate_rate_spec(rec.rates, "rec.rates")

  mu_is_distribution  <- is.list(mu.rates)  && !is.null(mu.rates$distribution)
  rec_is_distribution <- is.list(rec.rates) && !is.null(rec.rates$distribution)

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

  # Per-locus rec rate path needs a "placeholder" rho in the args
  # because scrm's parser requires -r to be valid syntactically.
  # The C side overrides per-locus via Model::setRecombinationRate
  # before each tree build, so the placeholder value is unused.
  rec_scalar_for_args <- if (rec_is_distribution) rec.rates$median
                         else                     rec.rates
  # Same idea for mu: -t needs a scalar placeholder; per-locus override
  # happens on the C side via Model::setMutationRate before each locus.
  mu_scalar_for_args  <- if (mu_is_distribution)  mu.rates$median
                         else                     mu.rates

  # Build per-group scrm commands (one group per unique locus length)
  # For uniform lengths, this is just one group
  unique_lens <- sort(unique(locus_lengths))
  n_groups <- length(unique_lens)
  group_nloci <- as.integer(table(factor(locus_lengths, levels = unique_lens)))

  base_cmds <- character(n_groups)
  for (g in seq_along(unique_lens)) {
    gl <- unique_lens[g]
    theta_g <- ms_scalar * mu_scalar_for_args * gl
    rho_g <- ms_scalar * rec_scalar_for_args * gl
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

  total_sims <- nsim.blocks * batch.size * ncores

  if (verbose) {
    if (total_sims != nsims)
      cat(sprintf("PipeMaster:: Requested %d sims, running %d (rounded to batch.size=%d x ncores=%d)\n",
                  nsims, total_sims, batch.size, ncores))
    if (uniform_len) {
      cat(sprintf("PipeMaster:: scrm engine: %d sims, %d loci x %d bp, %d pops\n",
                  total_sims, nloci, as.integer(unique_lens[1]), npop))
    } else {
      cat(sprintf("PipeMaster:: scrm engine: %d sims, %d loci (%d length groups, %d-%d bp), %d pops\n",
                  total_sims, nloci, n_groups, as.integer(min(unique_lens)), as.integer(max(unique_lens)), npop))
    }
  }

  # Write header (build column names deterministically, no sim needed)
  if (!append.sims || !file.exists(outfile)) {
    col_names <- .scrm.col.names(model, npop, config)
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

    # Register parent PID + on.exit cleanup (kills workers on exit/error/interrupt)
    pid_file <- file.path(abs_path, ".PM_parent.pid")
    worker_pids_env <- new.env(parent = emptyenv())
    worker_pids_env$pids <- integer(0)
    .pm.register.parent(pid_file, worker_pids_env)

    worker_nsims <- nsim.blocks * batch.size
    save(model, worker_nsims, batch.size, skip.zns,
         mu.rates, rec.rates, output.name, variable_samples,
         file = file.path(abs_path, ".PM_scrm_worker_params.RData"))

    # Worker script
    worker_script <- paste(
      'args <- commandArgs(TRUE)',
      'worker_id <- as.integer(args[1])',
      'suppressMessages(library(PipeMaster))',
      sprintf('base_path <- "%s"', abs_path),
      'pid_file <- file.path(base_path, ".PM_parent.pid")',
      'load(file.path(base_path, ".PM_scrm_worker_params.RData"))',
      'worker_dir <- file.path(base_path, paste0(".scrm_worker_", worker_id))',
      'dir.create(worker_dir, showWarnings = FALSE)',
      'sim.scrm.sumstats(model = model, nsims = worker_nsims,',
      '               batch.size = batch.size, mu.rates = mu.rates,',
      '               rec.rates = rec.rates, skip.zns = skip.zns,',
      '               output.name = output.name,',
      '               path = worker_dir, variable_samples = variable_samples,',
      '               ncores = 1, append.sims = TRUE,',
      '               .parent.pid.file = pid_file)',
      'write("done", file.path(base_path, paste0(".scrm_worker_", worker_id, ".done")))',
      'quit(save = "no")',
      sep = "\n")
    writeLines(worker_script, file.path(abs_path, ".PM_scrm_worker.R"))

    start_time <- Sys.time()
    for (w in 1:ncores) {
      pid <- system(paste("Rscript", file.path(abs_path, ".PM_scrm_worker.R"), w,
                   ">", file.path(abs_path, paste0(".scrm_worker_", w, ".log")), "2>&1 & echo $!"),
             intern = TRUE)
      wpid <- suppressWarnings(as.integer(pid[length(pid)]))
      if (!is.na(wpid)) worker_pids_env$pids <- c(worker_pids_env$pids, wpid)
    }
    cat(sprintf("PipeMaster:: Launched %d worker processes\n", ncores))

    total_expected <- nsim.blocks * batch.size * ncores
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
      # Check if parent process is still alive (workers only)
      if (!is.null(.parent.pid.file) && !.pm.parent.alive(.parent.pid.file)) {
        cat("PipeMaster:: Parent process died, worker exiting.\n")
        return(invisible(NULL))
      }
      block_results <- vector("list", batch.size)
      for (i in 1:batch.size) {
        block_results[[i]] <- .scrm.run.one(model, base_cmds, config, npop,
                                            skip.zns, mu.rates, rec.rates)
      }

      block_mat <- do.call(rbind, block_results)
      write.table(block_mat, file = outfile, append = TRUE, quote = FALSE,
                  row.names = FALSE, col.names = FALSE, sep = "\t")

      total_done <- total_done + batch.size
      if (verbose) {
        elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
        nsim_total <- nsim.blocks * batch.size
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
# rec.rates: scalar (back-compat) or list(distribution, median, sigma_log)
.scrm.run.one <- function(model, base_cmds, config, npop,
                          skip.zns, mu.rates, rec.rates = NULL) {
  nloci <- nrow(model$loci)

  # Sample parameters from priors
  cmd_result <- msABC.commander(model, arg = 1)
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

  # Per-locus rate sampling (A.2 heterogeneity).
  #
  # The per-locus rates are i.i.d. lognormal draws. The C side reads the
  # vector in length-group order, but since the draws are exchangeable
  # there is no per-group structure to preserve — a single rlnorm(nloci,
  # ...) call produces a statistically equivalent vector.
  #
  # sigma_log may be FIXED (spec$sigma_log) or SAMPLED per sim from a
  # uniform prior (spec$sigma_log_range = c(lo, hi)). The realized
  # sigma_log is recorded in the reftable as sigma_log_mu / sigma_log_rec
  # so tune.nn() can treat it as a regression target if desired.
  .draw_sigma_log <- function(spec) {
    # Use exact name lookup — avoid R's $ prefix matching.
    if ("sigma_log_range" %in% names(spec)) {
      r <- spec[["sigma_log_range"]]
      runif(1, r[1], r[2])
    } else {
      spec[["sigma_log"]]
    }
  }

  rec_per_locus    <- NULL
  sigma_log_rec_used <- NA_real_
  mean_rec_rate <- if (is.list(rec.rates)) rec.rates$median else
                   if (is.numeric(rec.rates)) rec.rates else NA_real_
  sd_rec_rate <- 0
  if (is.list(rec.rates) && !is.null(rec.rates$distribution)) {
    sigma_log_rec_used <- .draw_sigma_log(rec.rates)
    rec_per_locus <- rlnorm(nloci, meanlog = log(rec.rates$median),
                                   sdlog   = sigma_log_rec_used)
    mean_rec_rate <- mean(rec_per_locus)
    sd_rec_rate   <- stats::sd(rec_per_locus)
  }

  # Same pattern for mu. Even cheaper because mu does NOT affect the ARG
  # (only post-hoc SegSites mutation placement) — per-locus override is
  # essentially free.
  mu_per_locus     <- NULL
  sigma_log_mu_used <- NA_real_
  mean_mu_rate <- if (is.list(mu.rates)) mu.rates$median else
                  if (is.numeric(mu.rates)) mu.rates else NA_real_
  sd_mu_rate <- 0
  if (is.list(mu.rates) && !is.null(mu.rates$distribution)) {
    sigma_log_mu_used <- .draw_sigma_log(mu.rates)
    mu_per_locus <- rlnorm(nloci, meanlog = log(mu.rates$median),
                                  sdlog   = sigma_log_mu_used)
    mean_mu_rate <- mean(mu_per_locus)
    sd_mu_rate   <- stats::sd(mu_per_locus)
  }

  # scrm has a hardcoded default_pop_size_ = 10000 (src/scrm/model.h:76).
  # When scrm parses -t theta and -r rho (theta/rho built by PipeMaster with
  # Ne0 = 100000), it divides by 4*default_pop_size_ -- so internally scrm
  # stores mu and rec scaled by Ne0_PM / Ne0_scrm = 100000/10000 = 10. All
  # quantities (Ne, times, mu, rec) get the same 10x rescaling, leaving
  # theta and t/Ne preserved; the simulation is consistent.
  # Per-locus overrides via setMutationRate / setRecombinationRate bypass
  # that scaling, so we apply it here before passing to C. Skipping this
  # would make per-locus rates effectively 10x lower than nominal.
  scrm_ne_scaling <- 100000 / 10000   # PipeMaster Ne0 / scrm default_pop_size_
  rec_per_locus_c <- if (!is.null(rec_per_locus)) rec_per_locus * scrm_ne_scaling else NULL
  mu_per_locus_c  <- if (!is.null(mu_per_locus))  mu_per_locus  * scrm_ne_scaling else NULL

  if (length(scrm_cmds) == 1L) {
    # Uniform length: use original single-command call
    result <- .Call("scrm_stats_call", scrm_cmds, config, as.integer(npop),
                    as.logical(skip.zns), rec_per_locus_c, mu_per_locus_c,
                    PACKAGE = "PipeMaster")
  } else {
    # Variable lengths: use multi-command call with shared accumulators
    result <- .Call("scrm_stats_multi_call", scrm_cmds, config, as.integer(npop),
                    as.logical(skip.zns), as.integer(nloci),
                    rec_per_locus_c, mu_per_locus_c,
                    PACKAGE = "PipeMaster")
  }

  # Build output row: params + stats + sfs
  par_vec <- as.numeric(params[2, ])
  names(par_vec) <- params[1, ]

  # Add mu and rec rate info (mean.rate / sd.rate reflect the realized
  # per-sim mean and sd of the per-locus rates). sigma_log_mu / sigma_log_rec
  # record the drawn distribution-shape parameter when sigma_log_range is
  # used (NA when fixed sigma_log or scalar). When sigma_log_range is set,
  # these become learnable nuisance parameters (e.g. tune.nn() target).
  par_vec <- c(par_vec,
               mean.rate     = mean_mu_rate,
               sd.rate       = sd_mu_rate,
               sigma_log_mu  = sigma_log_mu_used,
               mean.rec.rate = mean_rec_rate,
               sd.rec.rate   = sd_rec_rate,
               sigma_log_rec = sigma_log_rec_used)

  # Fold and name the SFS entries
  # C code already folds per-site to minor allele; trim to correct length
  sfs_vec <- result$sfs
  nsam <- sum(config)
  if (npop == 1) {
    # 1-pop: C output is length nsam-1 but only first floor(nsam/2) bins populated
    sfs_len <- floor(nsam / 2)
    sfs_vec <- sfs_vec[1:sfs_len]
    names(sfs_vec) <- paste0("sfs_", seq(0, sfs_len - 1))
  } else {
    # Multi-pop joint SFS: already folded to minor allele in C
    # Use expand.grid naming convention (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(config, function(n) 0:n))
    names(sfs_vec) <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_vec, result$stats, sfs_vec)
}


# Internal: build column names deterministically (no simulation needed)
.scrm.col.names <- function(model, npop, config) {
  # Parameter names from model
  size_pars <- rbind(model$flags$n, model$flags$en$size)
  mig_pars <- rbind(model$flags$m, model$flags$em$size)
  time_pars <- rbind(model$flags$ej, model$flags$en$time, model$flags$em$time)

  par_names <- c(size_pars[, 1], time_pars[, 1])
  if (!is.null(mig_pars)) par_names <- c(par_names, mig_pars[, 1])
  par_names <- c(par_names,
                 "mean.rate", "sd.rate", "sigma_log_mu",
                 "mean.rec.rate", "sd.rec.rate", "sigma_log_rec")

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
    sfs_names <- paste0("sfs_", seq(0, sfs_len - 1))
  } else {
    # Joint SFS: expand.grid naming (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(config, function(n) 0:n))
    sfs_names <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_names, stat_names, sfs_names)
}
