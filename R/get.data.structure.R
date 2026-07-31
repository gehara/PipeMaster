#' Read observed data and set simulation parameters (bp and samples per population)
#' @param model A model object generated with main.menu.gui function.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file.
#'        Exactly one of path.to.phylip or path.to.vcf must be provided.
#' @param path.to.vcf Path to a VCF file. One locus per chromosome/scaffold. The locus set and
#'                    their lengths are read from the \code{##contig=<ID=...,length=...>} header
#'                    lines, so loci that were sequenced but carry no variants are retained --
#'                    they belong in the per-locus denominator of the summary statistics. Supply
#'                    \code{chrom.sizes} only for VCFs written without contig headers, or to
#'                    override them.
#' @param pop.assign A two-column data frame with population assignment of individuals. First column
#'                   must have the sample names, the second column the population in numbers (1, 2, etc.).
#' @param chrom.sizes A two-column data frame with chromosome/scaffold names in the first column and
#'                    total callable base pairs in the second column. Optional: needed only when the
#'                    VCF has no \code{##contig} header lines, or to override them. When supplied it
#'                    defines the callable set as given -- contigs are no longer intersected with
#'                    those carrying variants.
#' @param verbose Logical. If TRUE prints progress messages. Default is TRUE.
#' @author Marcelo Gehara
#' @return Model object with updated gene parameters.
#' @export
get.data.structure <- function(model, path.to.phylip = NULL,
                               path.to.vcf = NULL, pop.assign, chrom.sizes = NULL,
                               verbose = TRUE) {

  # Validate inputs
  n_inputs <- sum(!is.null(path.to.phylip), !is.null(path.to.vcf))
  if(n_inputs == 0)
    stop("One of path.to.phylip or path.to.vcf must be provided")
  if(n_inputs > 1)
    stop("Only one of path.to.phylip or path.to.vcf should be provided")

  pop.assign <- data.frame(pop.assign)
  pops <- pop.assign
  pops <- pops[order(pops[, 2]), ]
  n_pops <- unique(pops[,2])

  # all samples per pop in different files
  pops_samples <- list()
  for(i in n_pops){
    pops_samples[[i]] <- pops[which(pops[,2] == i), ]
  }

  if(!is.null(path.to.vcf)) {
    # ---- VCF path ----
    # Read the header once: sample names plus any ##contig declarations.
    # The contig header is the authoritative record of what was sequenced --
    # it lists loci with no variants, which the CHROM column cannot.
    con <- file(path.to.vcf, "r")
    vcf_header <- NULL
    contig_lines <- list()
    nblock <- 0L
    repeat {
      block <- readLines(con, n = 10000L)
      if(length(block) == 0) break
      hit <- which(startsWith(block, "#CHROM"))
      if(length(hit) > 0) {
        vcf_header <- strsplit(block[hit[1]], "\t")[[1]]
        block <- if(hit[1] > 1L) block[seq_len(hit[1] - 1L)] else character(0)
      }
      cl <- block[startsWith(block, "##contig=")]
      if(length(cl) > 0) {
        nblock <- nblock + 1L
        contig_lines[[nblock]] <- cl
      }
      if(!is.null(vcf_header)) break
    }
    close(con)

    if(is.null(vcf_header))
      stop("Could not find #CHROM header line in VCF file")

    vcf_samples <- vcf_header[10:length(vcf_header)]

    # Parse ##contig=<ID=name,length=bp> declarations
    contig_lines <- unlist(contig_lines, use.names = FALSE)
    header_chroms <- NULL
    if(length(contig_lines) > 0) {
      inner <- sub("^##contig=<(.*)>[[:space:]]*$", "\\1", contig_lines)
      has_id  <- grepl("(^|,)ID=",     inner, perl = TRUE)
      has_len <- grepl("(^|,)length=", inner, perl = TRUE)
      if(any(has_id & has_len)) {
        keep <- has_id & has_len
        header_chroms <- data.frame(
          chrom = sub("^.*?(?:^|,)ID=([^,]+).*$",     "\\1", inner[keep], perl = TRUE),
          bp    = as.numeric(sub("^.*?(?:^|,)length=([0-9]+).*$", "\\1",
                                 inner[keep], perl = TRUE)),
          stringsAsFactors = FALSE)
      }
    }

    user_sizes <- !is.null(chrom.sizes)
    if(user_sizes) {
      # Explicit chrom.sizes wins outright -- the caller is declaring the
      # callable set, so it is used as given (no intersection with the VCF).
      chrom.sizes <- data.frame(chrom.sizes)
      if(ncol(chrom.sizes) < 2)
        stop("chrom.sizes must have at least 2 columns (chrom name, bp)")
      colnames(chrom.sizes)[1:2] <- c("chrom", "bp")
      chrom.sizes$chrom <- as.character(chrom.sizes$chrom)
      chrom.sizes$bp    <- as.numeric(chrom.sizes$bp)
    } else if(!is.null(header_chroms)) {
      chrom.sizes <- header_chroms
    } else {
      stop("This VCF has no usable ##contig=<ID=...,length=...> header lines.\n",
           "  Supply chrom.sizes (chromosome/scaffold name, callable bp) instead.")
    }

    if(nrow(chrom.sizes) == 0)
      stop("No contigs found for the VCF")
    if(any(is.na(chrom.sizes$bp)))
      stop("Missing contig lengths for: ",
           paste(utils::head(chrom.sizes$chrom[is.na(chrom.sizes$bp)], 5),
                 collapse = ", "))

    shared_chroms <- chrom.sizes$chrom
    n_loci <- nrow(chrom.sizes)
    base_pairs <- as.integer(chrom.sizes$bp)

    if(verbose) {
      src <- if(user_sizes) "chrom.sizes argument" else "VCF ##contig header"
      cat(sprintf("PipeMaster:: %d loci declared by the %s\n", n_loci, src))
      cat("PipeMaster:: Scanning VCF for chromosomes carrying variants...\n")
      vcf_chroms <- unique(read.table(path.to.vcf, comment.char = "#", sep = "\t",
                                       colClasses = c("character", rep("NULL", 8)),
                                       nrows = -1)[, 1])
      n_var <- sum(shared_chroms %in% vcf_chroms)
      cat(sprintf("PipeMaster::   %d with variants, %d invariant (retained)\n",
                  n_var, n_loci - n_var))
      unknown <- setdiff(vcf_chroms, shared_chroms)
      if(length(unknown) > 0)
        warning(sprintf(paste0("%d chromosome(s) carry variants but are not in the ",
                               "declared locus set and will be IGNORED: %s%s"),
                        length(unknown), paste(utils::head(unknown, 5), collapse = ", "),
                        if(length(unknown) > 5) ", ..." else ""), call. = FALSE)
    }

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

    if(verbose) cat(paste0("PipeMaster:: VCF data structure: ", n_loci, " chromosomes, ",
                          sum(base_pairs), " total sites, ",
                          paste(pop_str, collapse = "/"), " haploid samples per pop\n"))

  } else {
    # ---- PHYLIP path ----
    if(verbose) cat("PipeMaster:: Reading multi-locus PHYLIP file...\n")
    loci <- read.phylip.loci(path.to.phylip)
    n_loci <- length(loci)

    if(verbose) pb <- txtProgressBar(min = 0, max = n_loci, style = 3)

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
      if(verbose) setTxtProgressBar(pb, i)
    }
    if(verbose) { close(pb); cat("\n") }

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

    model$loci <- LOCI
    model$I <- I
  }

  return(model)
}
