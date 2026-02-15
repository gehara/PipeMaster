#' Read the observed data to get the simulation parameters; base pairs and number of individuals per population.
#' @param model A model object generated with main.menu.gui function.
#' @param path.to.fasta Path to the fasta alignments directory. Exactly one of path.to.fasta,
#'                      path.to.phylip, or path.to.vcf must be provided.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file.
#' @param path.to.vcf Path to a VCF file (whole-genome data only). The CHROM column is used to
#'                    define loci (one locus per chromosome/scaffold). Requires chrom.sizes.
#' @param pop.assign A two-column data frame with population assignment of individuals. First column
#'                   must have the sample names, the second column the population in numbers (1, 2, etc.).
#' @param chrom.sizes A two-column data frame with chromosome/scaffold names in the first column and
#'                    total callable base pairs in the second column. Required when using VCF input.
#'                    Only chromosomes present in both chrom.sizes and the VCF are included.
#' @param sanger logical. If TRUE the inheritance scalar and mutation rates set up in the main.menu are kept. If FALSE the mutation rate is propagated to all loci. Default is FALSE.
#' @author Marcelo Gehara
#' @return Model object with updated gene parameters.
#' @export
get.data.structure <- function(model, path.to.fasta = NULL, path.to.phylip = NULL,
                               path.to.vcf = NULL, pop.assign, chrom.sizes = NULL,
                               sanger = FALSE) {

  # Validate inputs
  n_inputs <- sum(!is.null(path.to.fasta), !is.null(path.to.phylip), !is.null(path.to.vcf))
  if(n_inputs == 0)
    stop("One of path.to.fasta, path.to.phylip, or path.to.vcf must be provided")
  if(n_inputs > 1)
    stop("Only one of path.to.fasta, path.to.phylip, or path.to.vcf should be provided")

  if(!is.null(path.to.vcf) && is.null(chrom.sizes))
    stop("chrom.sizes is required when using VCF input")

  pop.assign <- data.frame(pop.assign)
  pops <- pop.assign
  pops <- pops[with(pops, order(pops[,2])), ]
  n_pops <- unique(pops[,2])

  # all samples per pop in different files
  pops_samples <- list()
  for(i in n_pops){
    pops_samples[[i]] <- pops[which(pops[,2] == i), ]
  }

  if(!is.null(path.to.vcf)) {
    # ---- VCF path (whole-genome data) ----
    chrom.sizes <- data.frame(chrom.sizes)
    if(ncol(chrom.sizes) < 2)
      stop("chrom.sizes must have at least 2 columns (chrom name, bp)")
    colnames(chrom.sizes)[1:2] <- c("chrom", "bp")
    chrom.sizes$bp <- as.numeric(chrom.sizes$bp)

    # Read VCF header to get sample names
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

    vcf_samples <- vcf_header[10:length(vcf_header)]

    # Discover which chromosomes are present in the VCF
    cat("PipeMaster:: Scanning VCF for chromosomes...\n")
    vcf_chroms <- unique(read.table(path.to.vcf, comment.char = "#", sep = "\t",
                                     colClasses = c("character", rep("NULL", 8)),
                                     nrows = -1)[, 1])

    # Keep only chromosomes present in both VCF and chrom.sizes
    shared_chroms <- intersect(as.character(chrom.sizes$chrom), vcf_chroms)
    if(length(shared_chroms) == 0)
      stop("No chromosomes in common between VCF and chrom.sizes")

    chrom.sizes <- chrom.sizes[chrom.sizes$chrom %in% shared_chroms, ]
    chrom.sizes <- chrom.sizes[match(shared_chroms, chrom.sizes$chrom), ]

    n_loci <- nrow(chrom.sizes)
    base_pairs <- as.integer(chrom.sizes$bp)

    # Count samples per population (diploid individuals -> haploid = 2x)
    pop_str <- integer(length(n_pops))
    for(j in seq_along(pops_samples)) {
      matched <- sum(pops_samples[[j]][,1] %in% vcf_samples)
      if(matched == 0)
        warning(paste("No samples from population", j, "found in VCF header"))
      pop_str[j] <- matched * 2L  # haploid sample size
    }

    LOCI <- cbind(paste0("chrom_", shared_chroms),
                  base_pairs,
                  rep(1, n_loci),
                  rep(model$loci[1, 4], n_loci),
                  rep(model$loci[1, 5], n_loci),
                  rep(model$loci[1, 6], n_loci))
    colnames(LOCI) <- NULL

    I <- cbind(paste0("chrom_", shared_chroms),
               rep("-I", n_loci),
               rep(model$I[1, 3], n_loci),
               matrix(rep(pop_str, each = n_loci), nrow = n_loci))

    model$loci <- LOCI
    model$I <- I

    cat(paste0("PipeMaster:: VCF data structure: ", n_loci, " chromosomes, ",
               sum(base_pairs), " total sites, ",
               paste(pop_str, collapse = "/"), " haploid samples per pop\n"))

  } else if(!is.null(path.to.phylip)) {
    # ---- PHYLIP path ----
    cat("PipeMaster:: Reading multi-locus PHYLIP file...\n")
    loci <- read.phylip.loci(path.to.phylip)
    n_loci <- length(loci)

    base_pairs <- NULL
    pop_str <- NULL
    for(i in 1:n_loci) {
      mat <- loci[[i]]
      bp <- ncol(mat)
      sample_names <- rownames(mat)

      npop_counts <- list()
      for(j in seq_along(pops_samples)) {
        npop_counts[[j]] <- length(na.omit(match(sample_names,
                                                  as.character(pops_samples[[j]][,1]))))
      }

      if(sum(unlist(npop_counts)) != nrow(mat)) {
        stop(paste("one or more samples in locus", i,
                   "have no assignment in your pop.assign file"))
      }

      pop_str <- rbind(pop_str, unlist(npop_counts))
      base_pairs <- c(base_pairs, bp)
      cat(paste(i, "loci"), "\n")
    }

    LOCI <- cbind(paste("locus", 1:length(base_pairs), sep = ""),
                  base_pairs,
                  rep(1, length(base_pairs)),
                  rep(model$loci[1, 4], length(base_pairs)),
                  rep(model$loci[1, 5], length(base_pairs)),
                  rep(model$loci[1, 6], length(base_pairs)))
    colnames(LOCI) <- NULL

    I <- cbind(paste("locus", 1:length(base_pairs), sep = ""),
               rep("-I", length(base_pairs)),
               rep(model$I[1, 3], length(base_pairs)),
               pop_str)

    if(sanger == TRUE) {
      model$loci[, 2] <- LOCI[, 2]
      model$I[, 4:ncol(I)] <- I[, 4:ncol(I)]
    } else {
      model$loci <- LOCI
      model$I <- I
    }

  } else {
    # ---- FASTA path (original behavior) ----
    orig.path <- getwd()
    setwd(path.to.fasta)
    on.exit(setwd(orig.path))

    # get all fasta files in WD
    fasta <- list.files(pattern = ".fa")
    fasta <- fasta[grep(".fa", fasta, fixed = TRUE)]

    #### get the population structure for every locus
    base_pairs <- NULL
    pop_str <- NULL
    for(i in 1:length(fasta)) {
      seq <- read.dna(fasta[i], format = "fasta")

      bp <- ncol(seq)

      npop <- list()
      for(j in 1:length(pops_samples)) {
        npop[[j]] <- length(na.omit(match(rownames(seq),
                                           as.character(pops_samples[[j]][,1]))))
      }

      if(sum(unlist(npop)) != nrow(seq)) {
        stop(paste("one or more samples in locus", fasta[i],
                   "have no assignment in your pop.assign file"))
      }

      pop_str <- rbind(pop_str, unlist(npop))
      base_pairs <- c(base_pairs, bp)
      cat(paste(i, "loci"), "\n")
    }

    LOCI <- cbind(paste("locus", 1:length(base_pairs), sep = ""),
                  base_pairs,
                  rep(1, length(base_pairs)),
                  rep(model$loci[1, 4], length(base_pairs)),
                  rep(model$loci[1, 5], length(base_pairs)),
                  rep(model$loci[1, 6], length(base_pairs)))
    colnames(LOCI) <- NULL

    I <- cbind(paste("locus", 1:length(base_pairs), sep = ""),
               rep("-I", length(base_pairs)),
               rep(model$I[1, 3], length(base_pairs)),
               pop_str)

    if(sanger == TRUE) {
      model$loci[, 2] <- LOCI[, 2]
      model$I[, 4:5] <- I[, 4:5]
    } else {
      model$loci <- LOCI
      model$I <- I
    }
  }

  return(model)
}
