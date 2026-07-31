# Fast VCF-to-ms converter (C backend)
# Streams through VCF line by line in C, detecting locus boundaries on CHROM change.
# Each locus becomes one ms block with segsites/positions/haplotype strings.
# Expects diploid GT field (e.g. "0|1" or "0/1"), outputs phased haplotypes.
# Monomorphic loci (no variants in VCF) are emitted as "segsites: 0" blocks.
# CHROM names must be "locus_N" (1-based numbering) for gap detection.
vcf_to_ms_file <- function(filepath, pop.assign, contig_length, n_loci,
                           output_file, verbose = FALSE) {
  pops <- data.frame(pop.assign)
  pops <- pops[order(pops[,2]), ]

  if(verbose) cat(sprintf("PipeMaster::   VCF->ms (C): %d samples, %d expected loci\n",
                          nrow(pops), n_loci))
  t0 <- proc.time()

  loci_written <- .Call("vcf_to_ms_call",
                        path.expand(filepath),
                        as.character(pops[,1]),
                        as.integer(contig_length),
                        as.integer(n_loci),
                        output_file,
                        PACKAGE = "PipeMaster")

  if(verbose) {
    ms_size <- file.info(output_file)$size
    cat(sprintf("PipeMaster::   %d loci written (%.1f sec, ms file: %.1f MB)\n",
                loci_written, (proc.time() - t0)[3], ms_size / 1e6))
  }
  invisible(loci_written)
}

#' Subsample one SNP per locus from a multi-locus VCF or PHYLIP file
#' @description Reads a multi-locus VCF or PHYLIP file and writes a VCF
#'              containing exactly one randomly sampled segregating site per
#'              locus.  Loci with no segregating sites are skipped.  Output is
#'              always VCF regardless of input format.
#'              Uses reservoir sampling in C for speed.
#' @param path.to.vcf Path to input VCF file. Each CHROM value defines a locus.
#'        Exactly one of \code{path.to.vcf} or \code{path.to.phylip} must be provided.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file.
#' @param output Path to output VCF file.
#' @param seed Optional integer seed for reproducibility (passed to \code{set.seed}).
#' @param verbose Logical. If TRUE prints summary. Default is TRUE.
#' @return Invisibly returns the number of loci written (those with >= 1 seg site).
#' @export
one.snp.per.locus <- function(path.to.vcf = NULL, path.to.phylip = NULL,
                               output, seed = NULL, verbose = TRUE) {
  n_inputs <- sum(!is.null(path.to.vcf), !is.null(path.to.phylip))
  if (n_inputs != 1)
    stop("Exactly one of path.to.vcf or path.to.phylip must be provided")

  if (!is.null(seed)) set.seed(seed)

  input_path <- if (!is.null(path.to.vcf)) path.to.vcf else path.to.phylip
  input_type <- if (!is.null(path.to.vcf)) "VCF" else "PHYLIP"

  if (verbose) cat(sprintf("PipeMaster:: one.snp.per.locus -- input: %s (%s)\n",
                            input_path, input_type))
  t0 <- proc.time()

  if (!is.null(path.to.vcf)) {
    n_loci <- .Call("vcf_one_snp_call",
                    path.expand(path.to.vcf),
                    path.expand(output),
                    PACKAGE = "PipeMaster")
  } else {
    n_loci <- .Call("phylip_one_snp_call",
                    path.expand(path.to.phylip),
                    path.expand(output),
                    PACKAGE = "PipeMaster")
  }

  if (verbose) {
    sz <- file.info(path.expand(output))$size
    cat(sprintf("PipeMaster::   %d loci written (%.1f sec, %.1f MB)\n",
                n_loci, (proc.time() - t0)[3], sz / 1e6))
  }
  invisible(n_loci)
}

# Fast PHYLIP-to-ms converter (internal helper)
# Converts a multi-locus sequential PHYLIP file directly to a concatenated
# ms-format file, avoiding character matrices and per-nucleotide strsplit.
phylip_to_ms_file <- function(filepath, pop.assign, output_file, verbose = FALSE) {
  lines <- readLines(filepath)
  pops <- data.frame(pop.assign)
  pops <- pops[order(pops[,2]), ]

  # Pre-scan to count loci (lines parseable as "ntax nchar")
  nloci <- 0L
  ii <- 1L
  while(ii <= length(lines)) {
    if(nchar(trimws(lines[ii])) == 0) { ii <- ii + 1L; next }
    hdr <- as.integer(strsplit(trimws(lines[ii]), "\\s+")[[1]])
    nloci <- nloci + 1L
    ii <- ii + 1L + hdr[1]  # skip ntax sequence lines
  }

  if(verbose) {
    pb <- txtProgressBar(min = 0, max = nloci, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  con <- file(output_file, "w")
  on.exit(close(con), add = TRUE)

  i <- 1
  locus <- 0L
  while(i <= length(lines)) {
    # Skip blank lines
    if(nchar(trimws(lines[i])) == 0) { i <- i + 1; next }

    # Parse header: "ntax nchar"
    header <- as.integer(strsplit(trimws(lines[i]), "\\s+")[[1]])
    ntax <- header[1]; nchar_seq <- header[2]
    i <- i + 1

    # Read sequences as strings (NOT character matrices)
    # Detect tab-delimited (long names) vs fixed-width (10-char names)
    names_vec <- character(ntax)
    seqs <- character(ntax)
    for(j in 1:ntax) {
      if (grepl("\t", lines[i])) {
        parts <- strsplit(lines[i], "\t")[[1]]
        names_vec[j] <- trimws(parts[1])
        seqs[j] <- gsub("\\s", "", parts[2])
      } else {
        names_vec[j] <- trimws(substr(lines[i], 1, 10))
        seqs[j] <- gsub("\\s", "", substr(lines[i], 11, nchar(lines[i])))
      }
      i <- i + 1
    }

    # Reorder by pop.assign
    ord <- match(pops[,1], names_vec)
    ord <- ord[!is.na(ord)]
    seqs <- seqs[ord]
    n <- length(seqs)

    # Convert to raw matrix for vectorized SNP detection
    raw_mat <- do.call(rbind, lapply(seqs, function(s) charToRaw(tolower(s))))

    # Find variable sites (exclude gaps '-'=0x2d, N 'n'=0x6e)
    gap_mask <- colSums(raw_mat == as.raw(0x2d)) > 0
    n_mask <- colSums(raw_mat == as.raw(0x6e)) > 0
    ref <- raw_mat[1, ]
    var_mask <- colSums(raw_mat != rep(ref, each = n)) > 0
    snp_cols <- which(var_mask & !gap_mask & !n_mask)
    ss <- length(snp_cols)

    if(ss > 0) {
      # Encode: major allele -> 0, others -> 1
      snp_raw <- raw_mat[, snp_cols, drop = FALSE]
      for(col in 1:ss) {
        col_bytes <- snp_raw[, col]
        freq <- tabulate(as.integer(col_bytes), nbins = 122)
        major <- as.raw(which.max(freq))
        encoded <- rep(as.raw(0x31), length(col_bytes))
        encoded[col_bytes == major] <- as.raw(0x30)
        snp_raw[, col] <- encoded
      }
      hap_strings <- character(n)
      for(k in 1:n) hap_strings[k] <- rawToChar(snp_raw[k, ])
      pos <- snp_cols / nchar_seq

      writeLines(paste("segsites:", ss), con)
      writeLines(paste("positions:", paste(pos, collapse = " ")), con)
      writeLines(hap_strings, con)
    } else {
      writeLines("segsites: 0", con)
      writeLines("positions:", con)
    }
    writeLines("//", con)

    locus <- locus + 1L
    if(verbose) setTxtProgressBar(pb, locus)
  }
}

#' Observed summary statistics and SFS from a PHYLIP or VCF file
#' @description Calculates observed summary statistics and the Site Frequency Spectrum
#'              from a multi-locus sequential PHYLIP file or a VCF file. The output
#'              columns match \code{sim.sumstats()} (summary statistics followed by
#'              SFS bins), so the result can be used directly for ABC or emulator inference.
#' @param model A model object built by the main.menu.gui() function.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file.
#'        Exactly one of path.to.phylip or path.to.vcf must be provided.
#' @param path.to.vcf Path to a VCF file where each CHROM represents a locus.
#'        The VCF should contain phased or unphased diploid genotypes (GT field).
#'        Sample names in pop.assign must match VCF header sample columns.
#' @param pop.assign A two-column data frame: sample names in column 1,
#'        population number in column 2.
#' @param one.snp Logical. If TRUE, one segregating site is randomly sampled
#'        per locus for the SFS. Default is FALSE.
#' @param monomorphic Logical. If TRUE a monomorphic bin is prepended to the SFS.
#'        Default is FALSE.
#' @param verbose Logical. If TRUE prints progress messages. Default is TRUE.
#' @return A one-row matrix with columns: summary statistics (4 moments) then SFS bins,
#'         matching the non-parameter columns of \code{sim.sumstats()}.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
observed.sumstats <- function(model, path.to.phylip = NULL, path.to.vcf = NULL,
                          pop.assign,
                          one.snp = FALSE,
                          monomorphic = FALSE, variable_samples = FALSE,
                          verbose = TRUE) {
  WD <- getwd()

  # Input validation
  n_inputs <- sum(!is.null(path.to.phylip), !is.null(path.to.vcf))
  if(n_inputs == 0)
    stop("One of path.to.phylip or path.to.vcf must be provided")
  if(n_inputs > 1)
    stop("Only one of path.to.phylip or path.to.vcf should be provided")

  if(class(model) != "Model")
    stop("First argument should be an object of class Model generated by the main.menu.gui() function.")
  if(model$I[1,1] == "genomic")
    stop("Your model is incomplete. You need to get the data structure first (get.data.structure function)")
  if(is.null(nrow(model$loci)))
    stop("Your model is incomplete. Go through the gene menu first and then get.data.structure().")

  if (is.character(pop.assign) && length(pop.assign) == 1L && file.exists(pop.assign))
    pop.assign <- read.table(pop.assign, header = TRUE)
  pop.assign <- data.frame(pop.assign)
  if(ncol(pop.assign) < 2) stop("pop.assign must have at least 2 columns")
  if(length(which(pop.assign[,2] %in% 1:10 == FALSE)) > 0)
    stop("Population names should be numbers")

  # Population structure and SFS dimensions
  npop <- as.numeric(model$I[1, 3])
  pop_cols <- 4:(3 + npop)
  pop_sizes_mat <- matrix(as.numeric(model$I[, pop_cols]), ncol = npop)

  # Validate full population representation per locus
  zero_rep_mask <- apply(pop_sizes_mat, 1, function(x) any(x == 0))
  n_no_rep <- sum(zero_rep_mask)
  if (n_no_rep > 0) {
    bad_per_pop <- vapply(seq_len(ncol(pop_sizes_mat)),
                          function(p) sum(pop_sizes_mat[, p] == 0), integer(1))
    affected <- which(bad_per_pop > 0)
    msg <- sprintf("%d of %d loci have no representation in at least one population:\n",
                   n_no_rep, nrow(pop_sizes_mat))
    for (p in affected) {
      msg <- paste0(msg, sprintf("    pop %d: %d loci with 0 samples\n", p, bad_per_pop[p]))
    }
    msg <- paste0(msg,
                  "  Remove these loci, drop the affected population(s), or add more samples\n",
                  "  before running observed.sumstats().")
    stop(msg, call. = FALSE)
  }

  uniform_samples <- nrow(unique(pop_sizes_mat)) == 1

  if (!uniform_samples && !variable_samples)
    stop("Sample sizes vary across loci. SFS cannot be computed.\n",
         "  Use optimize.sfs.model() to downsample to uniform sizes, or\n",
         "  set variable_samples = TRUE to compute sumstats without SFS.")

  if (!uniform_samples && variable_samples)
    cat("PipeMaster:: WARNING: Variable sample sizes -- computing sumstats only (no SFS).\n",
        "  Use optimize.sfs.model() to downsample to uniform sizes for SFS.\n")

  pop_sizes_vec <- as.integer(pop_sizes_mat[1, ])
  nsam <- sum(pop_sizes_vec)

  input_path <- if(!is.null(path.to.vcf)) path.to.vcf else path.to.phylip
  input_type <- if(!is.null(path.to.vcf)) "VCF" else "PHYLIP"
  contig_length <- as.integer(model$loci[1, 2])

  t_total <- proc.time()
  if(verbose) {
    cat(sprintf("PipeMaster:: observed.sumstats -- input: %s (%s)\n", input_path, input_type))
    cat(sprintf("PipeMaster::   samples: %d, populations: %d\n",
                nrow(pop.assign), length(unique(pop.assign[,2]))))
  }

  tmpdir <- tempfile("obs_allstats_")
  dir.create(tmpdir)

  # Step 1: Convert input -> concatenated ms file
  n_steps <- if(uniform_samples) 4L else 3L
  if(verbose) cat(sprintf("PipeMaster:: [1/%d] Converting %s to ms format\n", n_steps, input_type))
  t0 <- proc.time()
  ms_file <- file.path(tmpdir, "combined.ms")
  if(!is.null(path.to.vcf)) {
    vcf_to_ms_file(path.to.vcf, pop.assign, contig_length, nrow(model$I), ms_file, verbose = verbose)
  } else {
    phylip_to_ms_file(path.to.phylip, pop.assign, ms_file, verbose = verbose)
  }
  if(verbose) {
    ms_size <- file.info(ms_file)$size
    cat(sprintf("PipeMaster::   done (%.1f sec, ms file: %.1f MB)\n",
                (proc.time() - t0)[3], ms_size / 1e6))
  }

  # Step 2: Create locfile and build command with --obs
  if(verbose) cat(sprintf("PipeMaster:: [2/%d] Setting up msABC command...", n_steps))
  locfile <- PipeMaster:::get.locfile(model)
  nloci_rows <- nrow(locfile)
  write.table(locfile, file.path(tmpdir, ".1locfile.txt"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")

  saved_wd <- getwd()
  setwd(tmpdir)
  com <- PipeMaster:::msABC.commander(model, arg = 1)
  command <- paste(com[[1]], "--obs", ms_file)
  if(verbose) cat(" done\n")

  # =========================================================================
  # Variable samples path: sumstats only (no SFS)
  # =========================================================================
  if (!uniform_samples) {

    if(verbose) cat(sprintf("PipeMaster:: [3/%d] Computing summary statistics via msABC...", n_steps))
    t0 <- proc.time()
    result <- run.msABC(command)
    setwd(saved_wd)
    if(verbose) cat(sprintf(" done (%.1f sec)\n", (proc.time() - t0)[3]))

    # Parse output
    x <- strsplit(result, "\t")
    frag_nam <- x[[1]]
    values <- as.numeric(x[[2]])

    # Filter thomson, FayWuH/fwh columns (ZnS retained)
    cols_remove <- grep("thomson|FayWuH|fwh", frag_nam)
    if(length(cols_remove) > 0) {
      frag_nam <- frag_nam[-cols_remove]
      values <- values[-cols_remove]
    }

    observed <- t(data.frame(values))
    colnames(observed) <- frag_nam
    rownames(observed) <- NULL

    unlink(tmpdir, recursive = TRUE)
    elapsed <- (proc.time() - t_total)[3]
    if(verbose) cat(sprintf("PipeMaster:: Finished -- %d summary statistics in %.1f sec (no SFS)\n",
                            ncol(observed), elapsed))
    setwd(WD)
    return(observed)
  }

  # =========================================================================
  # Uniform samples path: sumstats + SFS via compute_sumstats (same C as scrm)
  # =========================================================================

  # Step 3: Parse ms file into haplotype matrices
  if(verbose) cat("PipeMaster:: [3/5] Parsing ms file into haplotype matrices...")
  t0 <- proc.time()
  setwd(saved_wd)

  ms_lines <- readLines(ms_file)
  n_ms <- length(ms_lines)
  haplist <- vector("list", nrow(model$I))
  loc <- 0L
  i <- 1L

  while (i <= n_ms) {
    line <- ms_lines[i]
    if (!grepl("^segsites:", line)) { i <- i + 1L; next }

    loc <- loc + 1L
    ss <- as.integer(sub("segsites: ", "", line))
    i <- i + 2L  # skip positions line

    if (ss == 0) {
      # Leave a NULL placeholder -- compute_sumstats treats NULL entries as
      # monomorphic loci and uses length(haplist) as the per-locus denominator.
      # haplist[[loc]] <- NULL would DELETE the element and shrink the list,
      # deflating that denominator and inflating every moment statistic.
      haplist[loc] <- list(NULL)
      next
    }

    hap_strings <- ms_lines[i:(i + nsam - 1L)]
    i <- i + nsam

    # Convert to integer matrix (nsam x segsites)
    mat <- matrix(0L, nrow = nsam, ncol = ss)
    for (h in seq_len(nsam))
      mat[h, ] <- as.integer(strsplit(hap_strings[h], "", fixed = TRUE)[[1]])

    # one.snp: subsample one random column for SFS (stats use full matrix)
    if (one.snp && ss > 1L) {
      col_idx <- sample.int(ss, 1L)
      haplist[[loc]] <- mat[, col_idx, drop = FALSE]
    } else {
      haplist[[loc]] <- mat
    }
  }
  if(verbose) cat(sprintf(" done (%d loci, %.1f sec)\n", loc, (proc.time() - t0)[3]))

  # Step 4: Compute stats + SFS via compute_sumstats (same C code as scrm)
  if(verbose) cat("PipeMaster:: [4/5] Computing summary statistics and SFS...")
  t0 <- proc.time()

  cs_result <- .Call("compute_sumstats_batch_call",
                     haplist,
                     as.integer(pop_sizes_vec),
                     as.integer(npop),
                     FALSE,  # skip_zns = FALSE for observed data
                     PACKAGE = "PipeMaster")
  if(verbose) cat(sprintf(" done (%.1f sec)\n", (proc.time() - t0)[3]))

  # Step 5: Build output
  if(verbose) cat("PipeMaster:: [5/5] Building output...")

  # --- Summary statistics (already named: s_mean_*, s_var_*, s_skew_*, s_kurt_*) ---
  stat_values <- as.numeric(cs_result$stats)
  stat_names  <- names(cs_result$stats)

  # Filter thomson, FayWuH/fwh (same as sim.sumstats; ZnS retained)
  cols_remove <- grep("thomson|FayWuH|fwh", stat_names)
  if(length(cols_remove) > 0) {
    stat_values <- stat_values[-cols_remove]
    stat_names  <- stat_names[-cols_remove]
  }

  # --- SFS (folded, same convention as scrm: sfs[folded-1]) ---
  sfs_vec <- as.numeric(cs_result$sfs)

  if(npop == 1) {
    sfs_len <- floor(nsam / 2)
    sfs_vec <- sfs_vec[1:sfs_len]
    sfs.names <- paste0("sfs_", seq(0, sfs_len - 1))
  } else {
    idx_grid <- expand.grid(lapply(pop_sizes_vec, function(n) 0:n))
    sfs.names <- apply(idx_grid, 1, function(x) paste0("sfs_", paste(x, collapse = "_")))
  }

  if(monomorphic) {
    total_sites <- sum(as.numeric(model$loci[, 2]))
    mono <- total_sites - sum(sfs_vec)
    sfs_vec <- c(mono, sfs_vec)
    sfs.names <- c("sfs_mono", sfs.names)
  }

  # Combine: sumstats | SFS (same column order as sim.sumstats, minus parameters)
  observed <- c(stat_values, sfs_vec)
  names_all <- c(stat_names, sfs.names)
  result <- t(data.frame(observed))
  colnames(result) <- names_all
  rownames(result) <- NULL

  unlink(tmpdir, recursive = TRUE)

  elapsed <- (proc.time() - t_total)[3]
  if(verbose) cat(sprintf(" done\nPipeMaster:: Finished -- %d sumstats + %d SFS bins in %.1f sec\n",
                          length(stat_names), length(sfs.names), elapsed))

  setwd(WD)
  return(result)
}

