#' Simulate summary statistics using scrm (SMC') coalescent engine
#'
#' Uses the vendored scrm C++ engine for fast coalescent simulation with
#' recombination, combined with PipeMaster's C summary statistic functions.
#' All computation stays in C/C++ memory -- no R matrix overhead. This is
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
#'       shape). The location is either \code{median} (fixed across sims) or
#'       \code{median_range = c(lo, hi)}, drawn per sim from a LOG-uniform prior
#'       and recorded in the reftable as \code{median_mu} / \code{median_rec};
#'       give exactly one. Note \code{meanlog = log(median)}, so the MEAN
#'       per-locus rate is \code{median * exp(sigma_log^2/2)} -- with a fixed
#'       median the mean rate is unestimable and \code{sigma_log} is the only
#'       parameter that can raise it. The dispersion parameter
#'       \code{sigma_log} can be either FIXED or SAMPLED per sim:
#'       \itemize{
#'         \item Fixed:  \code{list(distribution = "lognormal", median = 5.83e-9, sigma_log = 0.3)}
#'         \item Sampled: \code{list(distribution = "lognormal", median = 5.83e-9, sigma_log_range = c(0.05, 0.5))}
#'           -- per sim, \code{sigma_log ~ Uniform(lo, hi)}; the realized value
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
#' @param skip.sfs Logical. If TRUE, skip the site frequency spectrum entirely
#'   and return a zero-length \code{sfs}. The joint SFS has
#'   \code{prod(config + 1)} bins, which is allocated, zeroed and accumulated
#'   PER LOCUS regardless of how many are non-zero, so it dominates runtime for
#'   designs with many populations (10 populations x 10 haplotypes is 2.6e10
#'   bins). Default FALSE, so existing behaviour is unchanged. Summary
#'   statistics are bit-identical either way.
#' @param stat.config Integer vector of per-population sample sizes used for
#'   SUMMARY STATISTIC calculation, or NULL (default). Must sum to the total
#'   number of haplotypes. NULL means the statistics follow the simulated deme
#'   structure in \code{model$I} -- the previous behaviour, unchanged.
#'
#'   The demography is always simulated on all \code{npop} demes of
#'   \code{model$I}; this argument controls only how the resulting haplotypes
#'   are partitioned when statistics are computed, so that models with different
#'   population structures emit identical, directly comparable statistic
#'   vectors:
#'   \itemize{
#'     \item \code{-I 1 50} with \code{stat.config = c(28,20,2)} -- splits 1 into 3
#'     \item \code{-I 3 28 20 2} with \code{stat.config = c(28,20,2)} -- identity
#'     \item \code{-I 7 12 4 6 6 12 8 2} with \code{stat.config = c(28,20,2)} -- pools 7 into 3
#'   }
#'   All three then produce the same columns and can be compared by
#'   \code{tune.nn.classify()} or \code{OOD.pretrain.classify()}. The first is
#'   the way to build a genuinely panmictic null: simulate one population and
#'   summarise it as three, rather than faking panmixia with joins at t = 0.
#'
#'   Haplotypes are assigned in order: the first \code{stat.config[1]}
#'   haplotypes to statistics population 1, and so on. scrm emits haplotypes in
#'   \code{-I} deme order, so pooling demes requires them to be ADJACENT in
#'   \code{model$I}; a \code{stat.config} whose boundaries cut through a deme
#'   is allowed (it is what the panmictic null needs) but reported, since it is
#'   otherwise usually a mistake.
#'
#'   The observed data needs no counterpart argument: pass
#'   \code{observed.sumstats()} a \code{pop.assign} with the pooled
#'   populations and a model whose \code{I} carries the pooled sample sizes.
#'
#'   Affects the joint SFS too -- its dimensions follow \code{stat.config}, so
#'   this is also the way to keep \code{prod(config + 1)} tractable for
#'   many-deme designs.
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
#' @param stall.seconds Numeric or NULL. Per-worker stall timeout in seconds
#'   for the head-node watchdog. If a worker's log file mtime has not advanced
#'   for more than \code{stall.seconds}, the head kills the worker via
#'   \code{SIGKILL} and treats its already-completed blocks as the final output
#'   (recovered from the worker subdir). \code{NULL} (default) disables the
#'   watchdog. See \code{Details} for sensible values per data regime.
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
#'   \item scrm walks the sequence updating a local tree, rather than building
#'     the full ARG up front. The SMC' approximation window is set with
#'     \code{-l 500r} (500 recombination events, scrm's own default). Our loci
#'     carry rho = 4*Ne*r*L of order 40-400, below that window, so the result
#'     is effectively exact; the window only bites for much longer loci.
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
#' Under wide priors, a small fraction of parameter draws produce pathological
#' ARGs that can take orders of magnitude longer than typical sims to complete.
#' One stuck sim blocks an entire worker (each worker is single-threaded over
#' its assigned blocks), and the head waits for all workers to finish before
#' compiling the final reftable. The \code{stall.seconds} watchdog runs in the
#' head: if a worker's per-block progress (tracked via its log mtime) stops
#' advancing for more than \code{stall.seconds}, the worker is killed and its
#' already-completed blocks are salvaged from its subdir during the compile
#' step. The reftable will contain fewer rows than requested when this fires.
#'
#' Sensible \code{stall.seconds} values depend on the per-block wall time
#' (= \code{batch.size} x typical-sim-time). For reference, batch.size = 32 with:
#' \itemize{
#'   \item RAD / short-locus (typical sim ~ seconds): try 300 (5 min)
#'   \item pseudo-WGS 100 kb informative (typical block ~ minutes): try 1800
#'   \item real-WGS informative (typical block ~ 90 min): try 18000 (5 h)
#'   \item real-WGS diffuse, OoA-scale (typical block ~ 3.3 h): try 36000 (10 h)
#' }
#' Tighten when the typical block time is well below the threshold; loosen if
#' you see false positives on the slow tail.
#'
#' @export
sim.scrm.sumstats <- function(model, nsims, batch.size = 32,
                           mu.rates, rec.rates,
                           skip.zns = TRUE,
                           skip.sfs = FALSE,
                           stat.config = NULL,
                           ncores = 1, path = ".", output.name = "scrm",
                           variable_samples = FALSE,
                           append.sims = FALSE, verbose = TRUE,
                           stall.seconds = NULL,
                           .parent.pid.file = NULL) {

  if (is.null(model$use.alpha))
    stop("model$use.alpha is missing. Set it on the model (e.g. model$use.alpha <- FALSE).")

  nsim.blocks <- ceiling(nsims / (ncores * batch.size))

  # Internal helper: validate a rates spec (scalar OR distribution list with
  # either fixed `sigma_log` or sampled `sigma_log_range = c(lo, hi)` prior).
  .validate_rate_spec <- function(spec, name) {
    if (!is.list(spec) || is.null(spec$distribution)) {
      if (!is.numeric(spec) || length(spec) != 1)
        stop(sprintf(paste("%s must be a single numeric value OR a list with",
                           "$distribution + (median XOR median_range) +",
                           "(sigma_log XOR sigma_log_range)"), name))
      return(invisible(NULL))
    }
    if (!(spec$distribution %in% c("lognormal")))
      stop(sprintf("%s$distribution: only 'lognormal' is currently supported", name))
    # PM-MEDIAN-20260903: the median may be FIXED (spec$median) or SAMPLED per
    # sim from a log-uniform prior (spec$median_range = c(lo, hi)). Exactly one.
    # Rationale: with meanlog pinned to log(median), the MEAN per-locus rate is
    # median * exp(sigma_log^2/2), so a fixed median leaves the mean rate
    # unestimable and sigma_log is the only lever that can raise it -- which is
    # why sigma_log_rec pinned at its prior ceiling on PonAbe real WGS.
    # PM-MEANSPEC-20260905: the location may be given as a MEDIAN or as a MEAN,
    # fixed or sampled -- exactly one of median / median_range / mean / mean_range.
    # With meanlog = log(median) the realised MEAN is median*exp(sigma^2/2), so a
    # median spec lets sigma silently rescale the rate (measured 1.69x at
    # sigma ~ U(0.05,1.60)) and hence every Ne, since Ne scales as 1/mu. A mean
    # spec sets meanlog = log(mean) - sigma^2/2 instead, so sigma becomes pure
    # dispersion and the mean rate is whatever was asked for.
    locs <- c("median","median_range","mean","mean_range")
    have <- locs[locs %in% names(spec)]
    if (length(have) != 1)
      stop(sprintf("%s must specify exactly one of median / median_range / mean / mean_range (got %d)",
                   name, length(have)))
    v <- spec[[have]]
    if (grepl("_range$", have)) {
      if (!is.numeric(v) || length(v) != 2 || v[1] <= 0 || v[2] < v[1])
        stop(sprintf("%s$%s must be c(lo, hi) with 0 < lo <= hi", name, have))
    } else {
      if (!is.numeric(v) || length(v) != 1 || v <= 0)
        stop(sprintf("%s$%s must be a positive numeric scalar", name, have))
    }
    # Use exact name match -- R's $ does partial-prefix matching, which would
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

  if (!is.null(stall.seconds) && (!is.numeric(stall.seconds) ||
                                  length(stall.seconds) != 1L ||
                                  !is.finite(stall.seconds) ||
                                  stall.seconds <= 0))
    stop("stall.seconds must be a single positive number, or NULL to disable.")

  mu_is_distribution  <- is.list(mu.rates)  && !is.null(mu.rates$distribution)
  rec_is_distribution <- is.list(rec.rates) && !is.null(rec.rates$distribution)

  npop <- as.integer(model$I[1, 3])
  pop_cols <- 4:(3 + npop)
  config <- as.integer(model$I[1, pop_cols])
  nsam <- sum(config)

  # Statistics population structure. The demography is simulated on all `npop`
  # demes; `stat_config` controls only how the emitted haplotypes are
  # partitioned when statistics are computed, so that models with different
  # population structures produce identical, comparable statistic vectors.
  # It can pool demes (7 -> 3), split a single population (1 -> 3), or be the
  # identity. See the stat.config documentation.
  if (is.null(stat.config)) {
    stat_config <- config
    stat_npop   <- npop
  } else {
    stat_config <- as.integer(stat.config)
    if (anyNA(stat_config) || length(stat_config) < 1L || any(stat_config < 1L))
      stop("stat.config must be a vector of positive integers.")
    if (sum(stat_config) != nsam)
      stop(sprintf(paste0("stat.config must sum to the total number of haplotypes.\n",
                          "  stat.config sums to %d, model$I has %d haplotypes (%s)."),
                   sum(stat_config), nsam, paste(config, collapse = " + ")))
    stat_npop <- length(stat_config)
    # Boundaries that cut through a simulated deme are legitimate (that is what
    # a panmictic null needs) but are usually an arithmetic slip, so report them.
    if (npop > 1L && stat_npop > 1L) {
      deme_edges <- cumsum(config)[-npop]
      stat_edges <- cumsum(stat_config)[-stat_npop]
      if (!all(stat_edges %in% deme_edges))
        cat(sprintf(paste0("PipeMaster:: NOTE: stat.config boundaries (%s) do not align with\n",
                           "  deme boundaries (%s); at least one deme is split across statistics\n",
                           "  populations. Intended for panmictic nulls -- check if not.\n"),
                    paste(stat_edges, collapse = ", "), paste(deme_edges, collapse = ", ")))
    }
  }
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
  # PM-MEDIAN-20260903: with median_range the spec has no $median, so take the
  # geometric centre of the range. This is a PLACEHOLDER only -- the C side
  # overrides it per locus -- but it must be a valid scalar: a NULL here makes
  # rho_g numeric(0) and silently mangles every per-group scrm command.
  .placeholder_rate <- function(spec) {
    for (k in c("median_range","mean_range"))
      if (k %in% names(spec)) return(exp(mean(log(spec[[k]]))))
    for (k in c("median","mean"))
      if (k %in% names(spec)) return(spec[[k]])
    stop("no location field in rate spec")
  }
  rec_scalar_for_args <- if (rec_is_distribution) .placeholder_rate(rec.rates)
                         else                     rec.rates
  # Same idea for mu: -t needs a scalar placeholder; per-locus override
  # happens on the C side via Model::setMutationRate before each locus.
  mu_scalar_for_args  <- if (mu_is_distribution)  .placeholder_rate(mu.rates)
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
    # PM-FIX-20260819: scrm's -l sets the SMC' approximation window. The plain
    # form is a length in BASE PAIRS (param.cc:314 -> set_window_length_seq);
    # a trailing "r" makes it a count of RECOMBINATION EVENTS
    # (param.cc:313 -> set_window_length_rec).
    #
    # This previously passed `-l group_nloci[g]` -- the LOCUS COUNT -- as a
    # base-pair window, so the accuracy/speed dial was set by an unrelated
    # quantity:
    #   * pseudo-WGS (one group of 10,000 loci)  -> 10,000 bp window on 100 kb
    #   * real WGS   (thousands of length groups) -> 1-127 bp windows, varying
    #     BETWEEN GROUPS WITHIN ONE MODEL, on loci of 1 kb - 18 Mb
    # Measured effect: first moments are unaffected (s_mean_* shifts < 0.01%),
    # but s_var_* moves up to 5.8% (z = 33) between 10,000 and 50,000 bp
    # windows, so the bp window was not in the saturated regime even for the
    # pseudo-WGS cells.
    #
    # scrm's own default is set_window_length_rec(500) (model.cc:40), i.e. 500
    # recombination events -- scale-free, so it adapts to any Ne, rec rate and
    # locus length instead of depending on how many loci happen to be
    # simulated. Our loci carry rho = 4*Ne*r*L of order 40-400 events, well
    # under 500, so the whole locus is retained and the ARG is effectively
    # exact. Passed explicitly rather than relying on the vendored default.
    cmd <- paste(cmd, "-l", "500r")
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
    if (!identical(stat_config, config))
      cat(sprintf("PipeMaster::   %d demes simulated (%s), statistics on %d populations (%s)\n",
                  npop, paste(config, collapse = " "),
                  stat_npop, paste(stat_config, collapse = " ")))
  }

  # Write header (build column names deterministically, no sim needed)
  if (!append.sims || !file.exists(outfile)) {
    col_names <- .scrm.col.names(model, stat_npop, stat_config, skip.sfs)
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
    save(model, worker_nsims, batch.size, skip.zns, skip.sfs, stat.config,
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
      '               skip.sfs = skip.sfs, stat.config = stat.config,',
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
    killed_workers <- integer(0)
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

      # Watchdog: kill workers whose log file mtime has not advanced for
      # more than stall.seconds. The worker's already-completed blocks are
      # left on disk for the compile step to pick up.
      if (!is.null(stall.seconds)) {
        for (w in setdiff(seq_len(ncores), killed_workers)) {
          done_f <- file.path(abs_path, sprintf(".scrm_worker_%d.done", w))
          if (file.exists(done_f)) next
          log_f  <- file.path(abs_path, sprintf(".scrm_worker_%d.log", w))
          if (!file.exists(log_f)) next
          age <- as.numeric(difftime(Sys.time(), file.mtime(log_f),
                                     units = "secs"))
          if (age > stall.seconds) {
            pid <- if (length(worker_pids_env$pids) >= w)
                     worker_pids_env$pids[w] else NA_integer_
            cat(sprintf("PipeMaster:: worker %d stalled %.0f sec (>%.0f), killing PID %s\n",
                        w, age, stall.seconds,
                        if (is.na(pid)) "?" else as.character(pid)))
            if (!is.na(pid))
              try(tools::pskill(pid, signal = tools::SIGKILL), silent = TRUE)
            writeLines("killed", done_f)
            killed_workers <- c(killed_workers, w)
          }
        }
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
    if (length(killed_workers) > 0) {
      cat(sprintf("PipeMaster:: %d worker(s) killed by stall watchdog (>%g sec): %s\n",
                  length(killed_workers), stall.seconds,
                  paste(killed_workers, collapse = ", ")))
      cat("PipeMaster:: their completed blocks were salvaged into the final reftable.\n")
    }

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
                                            skip.zns, mu.rates, rec.rates, skip.sfs,
                                            stat_config, stat_npop)
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
                          skip.zns, mu.rates, rec.rates = NULL,
                          skip.sfs = FALSE,
                          stat_config = NULL, stat_npop = NULL) {
  # config/npop drive the simulation (the -I flag); stat_config/stat_npop drive
  # the statistic calculation. They differ whenever stat.config was supplied.
  if (is.null(stat_config)) stat_config <- config
  if (is.null(stat_npop))   stat_npop   <- length(stat_config)
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
  # there is no per-group structure to preserve -- a single rlnorm(nloci,
  # ...) call produces a statistically equivalent vector.
  #
  # sigma_log may be FIXED (spec$sigma_log) or SAMPLED per sim from a
  # uniform prior (spec$sigma_log_range = c(lo, hi)). The realized
  # sigma_log is recorded in the reftable as sigma_log_mu / sigma_log_rec
  # so tune.nn() can treat it as a regression target if desired.
  # PM-MEDIAN-20260903: log-uniform draw of the median rate when median_range is
  # given (a rate is a scale parameter, so the prior is uniform on log scale).
  # Returns the fixed median unchanged otherwise, so the value recorded in the
  # reftable is always the median actually used.
  .draw_median <- function(spec, sigma_log) {
    # returns the MEDIAN of the lognormal actually used, i.e. exp(meanlog)
    if ("median_range" %in% names(spec)) {
      r <- spec[["median_range"]]; exp(runif(1, log(r[1]), log(r[2])))
    } else if ("median" %in% names(spec)) {
      spec[["median"]]
    } else {
      # mean spec: meanlog = log(mean) - sigma^2/2  =>  median = mean*exp(-sigma^2/2)
      m <- if ("mean_range" %in% names(spec)) {
             r <- spec[["mean_range"]]; exp(runif(1, log(r[1]), log(r[2])))
           } else spec[["mean"]]
      m * exp(-sigma_log^2 / 2)
    }
  }

  .draw_sigma_log <- function(spec) {
    # Use exact name lookup -- avoid R's $ prefix matching.
    if ("sigma_log_range" %in% names(spec)) {
      r <- spec[["sigma_log_range"]]
      runif(1, r[1], r[2])
    } else {
      spec[["sigma_log"]]
    }
  }

  rec_per_locus    <- NULL
  sigma_log_rec_used <- NA_real_
  median_rec_used <- if (is.list(rec.rates)) NA_real_ else
                     if (is.numeric(rec.rates)) rec.rates else NA_real_
  mean_rec_rate <- if (is.list(rec.rates)) NA_real_ else
                   if (is.numeric(rec.rates)) rec.rates else NA_real_
  sd_rec_rate <- 0
  if (is.list(rec.rates) && !is.null(rec.rates$distribution)) {
    sigma_log_rec_used <- .draw_sigma_log(rec.rates)
    median_rec_used    <- .draw_median(rec.rates, sigma_log_rec_used)
    rec_per_locus <- rlnorm(nloci, meanlog = log(median_rec_used),
                                   sdlog   = sigma_log_rec_used)
    mean_rec_rate <- mean(rec_per_locus)
    sd_rec_rate   <- stats::sd(rec_per_locus)
  }

  # Same pattern for mu. Even cheaper because mu does NOT affect the ARG
  # (only post-hoc SegSites mutation placement) -- per-locus override is
  # essentially free.
  mu_per_locus     <- NULL
  sigma_log_mu_used <- NA_real_
  median_mu_used <- if (is.list(mu.rates)) NA_real_ else
                    if (is.numeric(mu.rates)) mu.rates else NA_real_
  mean_mu_rate <- if (is.list(mu.rates)) NA_real_ else
                  if (is.numeric(mu.rates)) mu.rates else NA_real_
  sd_mu_rate <- 0
  if (is.list(mu.rates) && !is.null(mu.rates$distribution)) {
    sigma_log_mu_used <- .draw_sigma_log(mu.rates)
    median_mu_used    <- .draw_median(mu.rates, sigma_log_mu_used)
    mu_per_locus <- rlnorm(nloci, meanlog = log(median_mu_used),
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
    result <- .Call("scrm_stats_call", scrm_cmds, stat_config, as.integer(stat_npop),
                    as.logical(skip.zns), rec_per_locus_c, mu_per_locus_c,
                    as.logical(skip.sfs),
                    PACKAGE = "PipeMaster")
  } else {
    # Variable lengths: use multi-command call with shared accumulators
    result <- .Call("scrm_stats_multi_call", scrm_cmds, stat_config, as.integer(stat_npop),
                    as.logical(skip.zns), as.integer(nloci),
                    rec_per_locus_c, mu_per_locus_c,
                    as.logical(skip.sfs),
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
               median_mu     = median_mu_used,
               mean.rec.rate = mean_rec_rate,
               sd.rec.rate   = sd_rec_rate,
               sigma_log_rec = sigma_log_rec_used,
               median_rec    = median_rec_used)

  # Fold and name the SFS entries
  # C code already folds per-site to minor allele; trim to correct length
  sfs_vec <- result$sfs
  nsam <- sum(stat_config)
  if (length(sfs_vec) == 0L) {
    # skip.sfs = TRUE: C returned a zero-length sfs; emit no SFS columns.
    # The header builder below must agree, so it takes skip.sfs too.
    sfs_vec <- numeric(0)
  } else if (stat_npop == 1) {
    # 1-pop: C output is length nsam-1 but only first floor(nsam/2) bins populated
    sfs_len <- floor(nsam / 2)
    sfs_vec <- sfs_vec[1:sfs_len]
    names(sfs_vec) <- paste0("sfs_", seq(0, sfs_len - 1))
  } else {
    # Multi-pop joint SFS: already folded to minor allele in C
    # Use expand.grid naming convention (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(stat_config, function(n) 0:n))
    names(sfs_vec) <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_vec, result$stats, sfs_vec)
}


# Internal: build column names deterministically (no simulation needed)
.scrm.col.names <- function(model, npop, config, skip.sfs = FALSE) {
  # Parameter names from model
  size_pars <- rbind(model$flags$n, model$flags$en$size)
  mig_pars <- rbind(model$flags$m, model$flags$em$size)
  time_pars <- rbind(model$flags$ej, model$flags$en$time, model$flags$em$time)

  par_names <- c(size_pars[, 1], time_pars[, 1])
  if (!is.null(mig_pars)) par_names <- c(par_names, mig_pars[, 1])
  par_names <- c(par_names,
                 "mean.rate", "sd.rate", "sigma_log_mu", "median_mu",
                 "mean.rec.rate", "sd.rec.rate", "sigma_log_rec", "median_rec")

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
  if (isTRUE(skip.sfs)) {
    sfs_names <- character(0)
  } else if (npop == 1) {
    sfs_len <- floor(nsam / 2)
    sfs_names <- paste0("sfs_", seq(0, sfs_len - 1))
  } else {
    # Joint SFS: expand.grid naming (sfs_0_0, sfs_1_0, ...)
    idx_grid <- expand.grid(lapply(config, function(n) 0:n))
    sfs_names <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  c(par_names, stat_names, sfs_names)
}
