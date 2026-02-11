#' Observed Site Frequency Spectrum from FASTA Data
#' @description Computes the observed Site Frequency Spectrum (SFS) from empirical FASTA alignments.
#'              For single-population data, returns the 1D SFS. For multi-population data, returns
#'              the joint SFS. Output format matches sim.sfs() for direct comparison.
#' @param model A model object built by the main.menu function. Any model with the same number of
#'              populations as your data will work (used for consistency checks).
#' @param path.to.fasta Path to the folder containing FASTA alignments (.fas files).
#'                      Invariable sites must be included. Alignments must contain phased data.
#' @param pop.assign A two-column data frame with sample names in the first column and population
#'                   numbers in the second column. Numbers should match population ordering in the model.
#' @param folded Logical. If TRUE the SFS will be folded (1-pop only). Default is FALSE.
#' @param ncores Number of cores for parallel execution. When ncores > 1, loci are processed
#'               in parallel using mclapply. Default is 1.
#' @return A one-row data frame with SFS bins as columns. Column names follow sim.sfs() convention:
#'         sfs_0, sfs_1, ... for 1-pop; sfs_0_0, sfs_1_0, ... for multi-pop (expand.grid order).
#'         For 2-pop models, the result has an "sfs_matrix" attribute for use with plot.2D.sfs().
#' @note Polarization: fasta.snp.2ms() assigns the major allele as ancestral (0).
#'       This differs from true ancestral polarization, so the unfolded SFS may not match
#'       simulation-based SFS that use true ancestral alleles. Consider using folded=TRUE
#'       when comparing with such SFS, or when no outgroup is available.
#' @examples
#' \dontrun{
#' # Load a model and population assignment
#' pop_assign <- read.table("pop_assign.txt", header = FALSE)
#'
#' # Compute observed 1D SFS
#' obs <- obs.sfs(model = my_model,
#'                path.to.fasta = "path/to/fastas",
#'                pop.assign = pop_assign)
#'
#' # Compute folded SFS
#' obs_folded <- obs.sfs(model = my_model,
#'                       path.to.fasta = "path/to/fastas",
#'                       pop.assign = pop_assign,
#'                       folded = TRUE)
#' }
#' @author Marcelo Gehara
#' @export
obs.sfs <- function(model, path.to.fasta, pop.assign, folded = FALSE, ncores = 1) {

  WD <- getwd()
  on.exit(setwd(WD))

  # Validate inputs
  pop.assign <- data.frame(pop.assign)
  if(ncol(pop.assign) < 2) stop("pop.assign must have at least 2 columns")
  if(length(which(pop.assign[,2] %in% c(1:10) == FALSE)) > 0)
    stop("Population numbers should be integers")

  npop <- length(unique(pop.assign[,2]))
  pop_sizes <- as.integer(table(pop.assign[,2]))
  nsam <- sum(pop_sizes)

  # List FASTA files
  setwd(path.to.fasta)
  fasta.files <- list.files(pattern = "\\.fa")
  if(length(fasta.files) == 0) stop("No FASTA files found in ", path.to.fasta)
  n_loci <- length(fasta.files)

  cat(paste0("PipeMaster:: Computing observed SFS from ", n_loci,
             " loci with ", ncores, " cores\n"))

  # Process one locus: convert FASTA -> ms format -> parse -> compute SFS
  process.one.locus <- function(i) {
    ms.output <- PipeMaster:::fasta.snp.2ms(path.to.fasta, fasta.files[i],
                                              write.file = FALSE, pop.assign)
    ms_lines <- ms.output[[1]]

    # Parse segregating sites count
    ss <- as.integer(strsplit(ms_lines[3], ": ")[[1]][2])

    if(is.na(ss) || ss == 0) {
      if(npop == 1) return(numeric(nsam + 1))
      else return(array(0L, dim = pop_sizes + 1L))
    }

    # Parse haplotype strings into binary matrix (samples x SNPs)
    hap_strings <- ms_lines[5:length(ms_lines)]
    hap_strings <- hap_strings[nchar(hap_strings) > 0]
    n_haps <- length(hap_strings)
    seg_sites <- matrix(0L, nrow = n_haps, ncol = ss)
    for(h in seq_along(hap_strings)) {
      seg_sites[h, ] <- as.integer(strsplit(hap_strings[h], "")[[1]])
    }

    # Compute per-locus SFS
    if(npop == 1) {
      sfs <- tabulate(colSums(seg_sites) + 1L, nbins = nsam + 1L)
    } else {
      sfs <- compute_joint_sfs(seg_sites, pop_sizes)
    }
    sfs
  }

  # Process all loci (with optional parallelization)
  if(ncores > 1) {
    all_sfs <- list()
    batch_size <- ncores * 2
    for(batch_start in seq(1, n_loci, by = batch_size)) {
      batch_end <- min(batch_start + batch_size - 1, n_loci)
      batch_idx <- batch_start:batch_end
      batch_results <- parallel::mclapply(batch_idx, process.one.locus,
                                           mc.cores = ncores)
      all_sfs <- c(all_sfs, batch_results)
      cat(paste0("PipeMaster:: ", batch_end, " of ", n_loci, " loci processed\n"))
    }
  } else {
    all_sfs <- lapply(seq_len(n_loci), function(i) {
      result <- process.one.locus(i)
      if(i %% 50 == 0)
        cat(paste0("PipeMaster:: ", i, " of ", n_loci, " loci processed\n"))
      result
    })
  }

  # Sum SFS across all loci
  if(npop == 1) {
    total_sfs <- numeric(nsam + 1)
  } else {
    total_sfs <- array(0L, dim = pop_sizes + 1L)
  }
  for(s in all_sfs) {
    total_sfs <- total_sfs + s
  }

  # Format output as data frame (matching sim.sfs() convention)
  if(npop == 1) {
    if(folded) {
      total_sfs <- fold_sfs(total_sfs)
      sfs.names <- paste0("sfs_fold_", seq(0, length(total_sfs) - 1))
    } else {
      sfs.names <- paste0("sfs_", seq(0, length(total_sfs) - 1))
    }
    result <- data.frame(matrix(total_sfs, nrow = 1))
    colnames(result) <- sfs.names
  } else {
    # Multi-pop: flatten using expand.grid order (pop1 varies fastest)
    idx_grid <- expand.grid(lapply(pop_sizes, function(n) 0:n))
    sfs.names <- apply(idx_grid, 1, function(x)
      paste0("sfs_", paste(x, collapse = "_")))
    sfs_vec <- as.vector(total_sfs)
    result <- data.frame(matrix(sfs_vec, nrow = 1))
    colnames(result) <- sfs.names

    # Attach matrix attribute for 2-pop models (for plot.2D.sfs)
    if(npop == 2) {
      attr(result, "sfs_matrix") <- matrix(as.integer(total_sfs),
                                            nrow = pop_sizes[1] + 1,
                                            ncol = pop_sizes[2] + 1)
    }
  }

  cat("PipeMaster:: Done!\n")
  result
}
