#' Observed Site Frequency Spectrum from FASTA or PHYLIP Data
#' @description Computes the observed Site Frequency Spectrum (SFS) from empirical alignments.
#'              For single-population data, returns the 1D SFS. For multi-population data, returns
#'              the joint SFS. Output format matches sim.sfs() for direct comparison.
#'              SFS computation is performed natively in C for speed.
#' @param model A model object built by the main.menu function. Any model with the same number of
#'              populations as your data will work (used for consistency checks).
#' @param path.to.fasta Path to the folder containing FASTA alignments (.fas files).
#'                      Invariable sites must be included. Alignments must contain phased data.
#'                      Either path.to.fasta or path.to.phylip must be provided.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file. Each locus block starts
#'                       with a "ntax nchar" header followed by ntax sequence lines. Alternative
#'                       to path.to.fasta for datasets with many loci.
#' @param pop.assign A two-column data frame with sample names in the first column and population
#'                   numbers in the second column. Numbers should match population ordering in the model.
#' @param one.snp Logical. If TRUE, one segregating site is randomly sampled per locus before
#'                computing the SFS. This reduces variance inflation caused by linkage among
#'                sites within the same locus. Recommended for short loci where recombination
#'                is negligible. Default is FALSE.
#' @param folded Logical. If TRUE the SFS will be folded (1-pop only). Default is FALSE.
#' @return A one-row data frame with SFS bins as columns. Column names follow sim.sfs() convention:
#'         sfs_0, sfs_1, ... for 1-pop; sfs_0_0, sfs_1_0, ... for multi-pop (expand.grid order).
#'         For 2-pop models, the result has an "sfs_matrix" attribute for use with plot.2D.sfs().
#' @note Polarization: the major allele (most frequent) is assigned as ancestral (0).
#'       This differs from true ancestral polarization, so the unfolded SFS may not match
#'       simulation-based SFS that use true ancestral alleles. Consider using folded=TRUE
#'       when comparing with such SFS, or when no outgroup is available.
#' @examples
#' \dontrun{
#' # Load a model and population assignment
#' pop_assign <- read.table("pop_assign.txt", header = FALSE)
#'
#' # Compute observed 1D SFS from FASTA files
#' obs <- obs.sfs(model = my_model,
#'                path.to.fasta = "path/to/fastas",
#'                pop.assign = pop_assign)
#'
#' # Compute observed 1D SFS from PHYLIP file
#' obs <- obs.sfs(model = my_model,
#'                path.to.phylip = "path/to/data.phy",
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
obs.sfs <- function(model, path.to.fasta = NULL, path.to.phylip = NULL,
                    pop.assign, one.snp = FALSE, folded = FALSE) {

  # Validate inputs
  if(is.null(path.to.fasta) && is.null(path.to.phylip))
    stop("Either path.to.fasta or path.to.phylip must be provided")

  pop.assign <- data.frame(pop.assign)
  if(ncol(pop.assign) < 2) stop("pop.assign must have at least 2 columns")
  if(length(which(pop.assign[,2] %in% c(1:10) == FALSE)) > 0)
    stop("Population numbers should be integers")

  npop <- length(unique(pop.assign[,2]))
  pop_sizes <- as.integer(table(pop.assign[,2]))
  nsam <- sum(pop_sizes)

  # Read loci into list of character matrices
  if(!is.null(path.to.phylip)) {
    cat("PipeMaster:: Reading multi-locus PHYLIP file...\n")
    loci <- read.phylip.loci(path.to.phylip)
  } else {
    fasta.files <- list.files(path.to.fasta, pattern = "\\.fa")
    if(length(fasta.files) == 0) stop("No FASTA files found in ", path.to.fasta)
    cat(paste0("PipeMaster:: Reading ", length(fasta.files), " FASTA files...\n"))
    loci <- lapply(fasta.files, function(f) {
      dna <- ape::read.dna(file.path(path.to.fasta, f), format = "fasta")
      as.character(dna)
    })
  }

  n_loci <- length(loci)
  cat(paste0("PipeMaster:: Computing observed SFS from ", n_loci, " loci (native C)\n"))

  # Compute sample index mapping: pop-ordered position -> matrix row (0-based)
  pop.assign <- pop.assign[order(pop.assign[,2]), ]
  sample_names <- rownames(loci[[1]])
  sample_indices <- match(pop.assign[,1], sample_names) - 1L

  if(any(is.na(sample_indices)))
    stop("Some sample names in pop.assign not found in alignment data")

  # Native C SFS computation (all loci in one call)
  total_sfs <- .Call("obs_sfs_call", loci, sample_indices, pop_sizes, one.snp,
                     PACKAGE = "PipeMaster")

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
    result <- data.frame(matrix(total_sfs, nrow = 1))
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
