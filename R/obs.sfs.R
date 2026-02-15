#' Observed Site Frequency Spectrum from FASTA, PHYLIP, or VCF Data
#' @description Computes the observed Site Frequency Spectrum (SFS) from empirical alignments
#'              or VCF files. For single-population data, returns the 1D SFS. For multi-population
#'              data, returns the joint SFS. Output format matches sim.sfs() for direct comparison.
#'              SFS computation is performed natively in C for speed.
#' @param model A model object built by the main.menu.gui() function. Any model with the same number of
#'              populations as your data will work (used for consistency checks). When using VCF input,
#'              the model should first be configured with get.data.structure() so that model$loci
#'              contains per-chromosome bp values.
#' @param path.to.fasta Path to the folder containing FASTA alignments (.fas files).
#'                      Invariable sites must be included. Alignments must contain phased data.
#'                      Exactly one of path.to.fasta, path.to.phylip, or path.to.vcf must be provided.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file. Each locus block starts
#'                       with a "ntax nchar" header followed by ntax sequence lines. Alternative
#'                       to path.to.fasta for datasets with many loci.
#' @param path.to.vcf Path to a VCF file (plain text, whole-genome data). The VCF should contain
#'                    biallelic SNPs. Multi-allelic sites and sites with missing genotypes are
#'                    automatically skipped. Genotypes are treated as diploid (2 alleles per sample).
#'                    The monomorphic bin is computed automatically from model$loci bp values.
#' @param pop.assign A two-column data frame with sample names in the first column and population
#'                   numbers in the second column. Numbers should match population ordering in the model.
#'                   For VCF input, sample names must match VCF header sample columns.
#' @param monomorphic Logical. If TRUE, a monomorphic bin is included in the SFS: sfs_mono =
#'                    total_sites - sum(segregating bins), where total_sites is derived from
#'                    sum(model$loci[,2]). Automatically set to TRUE for VCF input. For
#'                    FASTA/PHYLIP, default is FALSE (standard multi-locus SFS without monomorphic bin).
#' @param one.snp Logical. If TRUE, one segregating site is randomly sampled per locus before
#'                computing the SFS. Not applicable for VCF input. Default is FALSE.
#' @param folded Logical. If TRUE the SFS will be folded (1-pop only). Default is FALSE.
#' @return A one-row data frame with SFS bins as columns. Column names follow sim.sfs() convention:
#'         sfs_0, sfs_1, ... for 1-pop; sfs_0_0, sfs_1_0, ... for multi-pop (expand.grid order).
#'         When monomorphic bin is included, the first column is sfs_mono.
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
#' # Compute observed SFS from VCF file (monomorphic bin added automatically)
#' obs <- obs.sfs(model = my_model,
#'                path.to.vcf = "path/to/data.vcf",
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
                    path.to.vcf = NULL,
                    pop.assign, monomorphic = FALSE,
                    one.snp = FALSE, folded = FALSE) {

  # Validate inputs
  n_inputs <- sum(!is.null(path.to.fasta), !is.null(path.to.phylip), !is.null(path.to.vcf))
  if(n_inputs == 0)
    stop("One of path.to.fasta, path.to.phylip, or path.to.vcf must be provided")
  if(n_inputs > 1)
    stop("Only one of path.to.fasta, path.to.phylip, or path.to.vcf should be provided")

  # VCF always includes monomorphic bin
  if(!is.null(path.to.vcf)) monomorphic <- TRUE

  pop.assign <- data.frame(pop.assign)
  if(ncol(pop.assign) < 2) stop("pop.assign must have at least 2 columns")
  if(length(which(pop.assign[,2] %in% c(1:10) == FALSE)) > 0)
    stop("Population numbers should be integers")

  npop <- length(unique(pop.assign[,2]))
  pop_sizes <- as.integer(table(pop.assign[,2]))
  nsam <- sum(pop_sizes)

  # ---- VCF path ----
  if(!is.null(path.to.vcf)) {
    if(!file.exists(path.to.vcf))
      stop("VCF file not found: ", path.to.vcf)

    cat("PipeMaster:: Reading VCF header to map samples...\n")

    # Read VCF header line to get sample column names
    con <- file(path.to.vcf, "r")
    vcf_header <- NULL
    while(TRUE) {
      line <- readLines(con, n = 1)
      if(length(line) == 0) break
      if(startsWith(line, "##")) next
      if(startsWith(line, "#CHROM")) {
        vcf_header <- strsplit(line, "\t")[[1]]
        break
      }
    }
    close(con)

    if(is.null(vcf_header))
      stop("Could not find #CHROM header line in VCF file")

    # VCF sample columns start at position 10 (1-based), i.e., index 10+
    vcf_samples <- vcf_header[10:length(vcf_header)]

    # Order pop.assign by population number
    pop.assign <- pop.assign[order(pop.assign[,2]), ]

    # Map pop.assign sample names to VCF column indices (0-based from first sample column)
    sample_col_indices <- match(pop.assign[,1], vcf_samples) - 1L
    if(any(is.na(sample_col_indices)))
      stop("Some sample names in pop.assign not found in VCF header: ",
           paste(pop.assign[is.na(sample_col_indices), 1], collapse = ", "))

    # Population map: which pop (1-based) each target sample belongs to
    sample_pop_map <- as.integer(pop.assign[,2])

    # Haploid sample sizes (2 * diploid count per pop)
    haploid_pop_sizes <- as.integer(2L * pop_sizes)

    cat(paste0("PipeMaster:: Computing observed SFS from VCF (native C)\n"))

    total_sfs <- .Call("obs_sfs_vcf_call", path.to.vcf, sample_pop_map,
                       sample_col_indices, haploid_pop_sizes, as.integer(npop),
                       PACKAGE = "PipeMaster")

  } else {
    # ---- FASTA/PHYLIP path ----

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
  }

  # ---- Add monomorphic bin if requested ----
  # Derive total sites from model loci bp values
  if(monomorphic) {
    total_sites <- sum(as.numeric(model$loci[, 2]))
    mono_count <- total_sites - sum(total_sfs)
    cat(paste0("PipeMaster:: Monomorphic bin: ", mono_count,
               " (total sites: ", total_sites, ")\n"))
  }

  # ---- Format output ----
  # For VCF path, use haploid pop sizes for SFS dimensions
  if(!is.null(path.to.vcf)) {
    sfs_pop_sizes <- as.integer(2L * pop_sizes)
  } else {
    sfs_pop_sizes <- pop_sizes
  }

  if(npop == 1) {
    if(folded) {
      total_sfs <- fold_sfs(total_sfs)
      sfs.names <- paste0("sfs_fold_", seq(0, length(total_sfs) - 1))
    } else {
      sfs.names <- paste0("sfs_", seq(0, length(total_sfs) - 1))
    }
    if(monomorphic) {
      total_sfs <- c(mono_count, total_sfs)
      sfs.names <- c("sfs_mono", sfs.names)
    }
    result <- data.frame(matrix(total_sfs, nrow = 1))
    colnames(result) <- sfs.names
  } else {
    # Multi-pop: flatten using expand.grid order (pop1 varies fastest)
    idx_grid <- expand.grid(lapply(sfs_pop_sizes, function(n) 0:n))
    sfs.names <- apply(idx_grid, 1, function(x)
      paste0("sfs_", paste(x, collapse = "_")))
    if(monomorphic) {
      total_sfs <- c(mono_count, total_sfs)
      sfs.names <- c("sfs_mono", sfs.names)
    }
    result <- data.frame(matrix(total_sfs, nrow = 1))
    colnames(result) <- sfs.names

    # Attach matrix attribute for 2-pop models (for plot.2D.sfs)
    if(npop == 2) {
      sfs_for_matrix <- if(monomorphic) total_sfs[-1] else total_sfs
      attr(result, "sfs_matrix") <- matrix(as.integer(sfs_for_matrix),
                                            nrow = sfs_pop_sizes[1] + 1,
                                            ncol = sfs_pop_sizes[2] + 1)
    }
  }

  cat("PipeMaster:: Done!\n")
  result
}
