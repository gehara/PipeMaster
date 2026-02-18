# Fast PHYLIP-to-ms converter (internal helper)
# Converts a multi-locus sequential PHYLIP file directly to a concatenated
# ms-format file, avoiding character matrices and per-nucleotide strsplit.
phylip_to_ms_file <- function(filepath, pop.assign, output_file) {
  lines <- readLines(filepath)
  pops <- data.frame(pop.assign)
  pops <- pops[order(pops[,2]), ]

  con <- file(output_file, "w")
  on.exit(close(con))

  i <- 1
  while(i <= length(lines)) {
    # Skip blank lines
    if(nchar(trimws(lines[i])) == 0) { i <- i + 1; next }

    # Parse header: "ntax nchar"
    header <- as.integer(strsplit(trimws(lines[i]), "\\s+")[[1]])
    ntax <- header[1]; nchar_seq <- header[2]
    i <- i + 1

    # Read sequences as strings (NOT character matrices)
    names_vec <- character(ntax)
    seqs <- character(ntax)
    for(j in 1:ntax) {
      names_vec[j] <- trimws(substr(lines[i], 1, 10))
      seqs[j] <- gsub("\\s", "", substr(lines[i], 11, nchar(lines[i])))
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
  }
}

#' Observed summary statistics for nexgen data
#' @description This function calculates the observed summary statistics from an empirical data.
#'              This summary statistics are the same as those simulated by the sim.sumstat function.
#'              It is optimized for nexgen data.
#' @param model A model object built by the main.menu.gui() function. Any model with the same number of populations of your empirical data will work. This is just to build the sumstats names correctly.
#' @param path.to.fasta Path to the folder containing all fastas to be included in the calculation.
#'                      Invariable sites must be included in the fasta alignments.
#'                      Invariable loci must also be included. Alignments must contain phased data.
#'                      Either path.to.fasta or path.to.phylip must be provided.
#' @param path.to.phylip Path to a multi-locus sequential PHYLIP file. Each locus block starts
#'                       with a "ntax nchar" header followed by ntax sequence lines. Alternative
#'                       to path.to.fasta for datasets with many loci.
#' @param pop.assign A two-column data frame with sample names in the first column and the corresponding population membership preferably as numbers in the second column (If you have a single population the numbers won't matter). The numbers should match the population number in the model object.
#' @param moments Logical. If TRUE computes the four moments (mean, variance, skewness, kurtosis) of each summary statistic across loci. Default is TRUE.
#' @param ncores Number of cores for parallel execution. When ncores > 1, loci are processed in parallel using mclapply. Default is 1.
#' @return A table containing the observed summary stats.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
obs.sumstat.ngs<-function(model=NULL,path.to.fasta=NULL,path.to.phylip=NULL,
                          pop.assign,moments=TRUE,ncores=1){

   WD <- getwd()
  if(is.null(model)){
    print("You did not specify the model. Sumstat calculations will output a table with the specific stat of each locus.")
    x<-readline("Would you like to continue? (Yes or No)")
    if(x %in% c("Y","y","Yes","YES","yes")){}else{stop()}
    }

   if(is.null(path.to.fasta) && is.null(path.to.phylip))
     stop("Either path.to.fasta or path.to.phylip must be provided")

   if(is.null(nrow(model$loci))) stop("Your model is incomplete. Go through the gene menu first (main.menu function) and then get the data structure (get.data.structure function).")

   if(model$loci[1,1]=="genomic") stop("Your model is incomplete. You need to get the data structure first (get.data.structure function)")

   if(ncol(pop.assign) < 2) stop ("Your pop.assign file has more than 2 columns")

   pop.assign <- data.frame(pop.assign)
   if(length(which(pop.assign[,2] %in% c(1:10) == F)) > 0) stop ("Your population names should be numbers")

  # Determine input mode
  use_phylip <- !is.null(path.to.phylip)

  if(use_phylip) {
    tmpdir <- tempfile("obs_sumstat_")
    dir.create(tmpdir)

    # Step 1: Convert PHYLIP -> concatenated ms file
    cat("PipeMaster:: Converting PHYLIP to ms format...\n")
    ms_file <- file.path(tmpdir, "combined.ms")
    phylip_to_ms_file(path.to.phylip, pop.assign, ms_file)

    # Step 2: Create locfile
    locfile <- PipeMaster:::get.locfile(model)
    locfile_path <- file.path(tmpdir, ".1locfile.txt")
    write.table(locfile, locfile_path, row.names=FALSE, col.names=TRUE, quote=FALSE, sep=" ")

    # Step 3: Build fragment-mode command with --obs
    saved_wd <- getwd()
    setwd(tmpdir)
    com <- PipeMaster:::msABC.commander(model, use.alpha=FALSE, arg=1)
    command <- paste(com[[1]], "--obs", ms_file)

    # Step 4: Single msABC call -> returns header + moments
    cat("PipeMaster:: Computing summary statistics via msABC...\n")
    result <- run.msABC(command)
    setwd(saved_wd)

    # Step 5: Parse output
    x <- strsplit(result, "\t")
    frag_nam <- x[[1]]  # header names
    values <- as.numeric(x[[2]])  # moment values

    # Filter thomson, ZnS, FayWuH/fwh columns
    cols_remove <- grep("thomson|ZnS|FayWuH|fwh", frag_nam)
    if(length(cols_remove) > 0) {
      frag_nam <- frag_nam[-cols_remove]
      values <- values[-cols_remove]
    }

    observed <- t(data.frame(values))
    colnames(observed) <- frag_nam

    unlink(tmpdir, recursive=TRUE)
    setwd(WD)
    return(observed)
  }

  # FASTA path
  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fa",fasta.files,fixed=T)]
  n_loci <- length(fasta.files)

  # Per-locus: convert FASTA alignment -> ms -> run msABC --obs
  process.one.locus <- function(i) {
    ms.output <- PipeMaster:::fasta.snp.2ms(path.to.fasta, fasta.files[i], write.file=TRUE, pop.assign)
    locus.name <- strsplit(fasta.files[i], ".", fixed=TRUE)[[1]][1]
    xx <- strsplit(ms.output[[1]][1], " ")
    xx <- xx[[1]][2:length(xx[[1]])]
    xx <- paste(xx, collapse=" ")
    obs_result <- run.msABC(paste(xx, "--obs", paste0(locus.name, ".ms")))
    writeLines(obs_result, paste0(locus.name, ".out"))
    result <- read.table(text=obs_result, header=TRUE)
    snps <- grep("segs", names(result))
    message(paste(i, "  ", result[snps[length(snps)]], "SNPs"))
    result
  }

  cat(paste0("PipeMaster:: Processing ", n_loci, " loci with ", ncores, " cores\n"))

  if(ncores > 1) {
    observed <- list()
    batch_size <- ncores * 2
    for(batch_start in seq(1, n_loci, by=batch_size)) {
      batch_end <- min(batch_start + batch_size - 1, n_loci)
      batch_idx <- batch_start:batch_end
      batch_results <- parallel::mclapply(batch_idx, process.one.locus, mc.cores = ncores)
      observed <- c(observed, batch_results)
      cat(paste0("PipeMaster:: ", batch_end, " of ", n_loci, " loci processed\n"))
    }
  } else {
    observed <- lapply(seq_len(n_loci), process.one.locus)
  }

  observed <- matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

  # Get column names via msABC.commander (fragment mode)
  locfile <- PipeMaster:::get.locfile(model)
  write.table(locfile, ".1locfile.txt", row.names=F, col.names=T, quote=F, sep=" ")
  com <- PipeMaster:::msABC.commander(model, use.alpha=F, arg=1)
  x <- strsplit(run.msABC(com[[1]]), "\t")
  frag_nam <- x[[1]]

  # Derive per-locus column names from fragment-mode header
  # Fragment mode outputs 4 moments per stat: stat_mean, stat_var, stat_skew, stat_kurt
  frag_nam <- frag_nam[!grepl("^X$", frag_nam) & nchar(frag_nam) > 0]
  perlocus_names <- sub("_mean$", "", frag_nam[seq(1, length(frag_nam), by=4)])
  colnames(observed) <- perlocus_names

  # Filter thomson, ZnS, and FayWuH/fwh columns (match sim.sumstat filtering)
  # FayWuH/fwh are excluded because they require knowing ancestral vs derived
  # alleles (polarization), which fasta.snp.2ms() cannot determine without an
  # outgroup — it assigns the major allele as ancestral (0), folding the spectrum.
  cols <- grep("thomson", colnames(observed))
  cols <- c(cols, grep("ZnS", colnames(observed)))
  cols <- c(cols, grep("FayWuH", colnames(observed)))
  cols <- c(cols, grep("fwh", colnames(observed)))
  if(length(cols) != 0) observed <- observed[, -cols, drop=FALSE]

  if(!is.null(model)){
    # Compute 4 moments across loci (mean, var, skewness, kurtosis)
    obs_mean <- colMeans(observed, na.rm=TRUE)
    obs_var  <- apply(observed, 2, var, na.rm=TRUE)
    obs_skew <- apply(observed, 2, e1071::skewness, na.rm=TRUE)
    # msABC computes regular kurtosis (m4/m2^2, normal=3) while e1071::kurtosis()
    # returns excess kurtosis (normal=0), so we add 3 to match.
    obs_kurt <- apply(observed, 2, e1071::kurtosis, na.rm=TRUE) + 3

    # Interleave: for each stat -> [mean, var, skew, kurt]
    moments_mat <- rbind(obs_mean, obs_var, obs_skew, obs_kurt)
    observed <- as.vector(moments_mat)

    # Filter same columns from fragment-mode names for final column naming
    cols <- grep("thomson", frag_nam)
    cols <- c(cols, grep("ZnS", frag_nam))
    cols <- c(cols, grep("FayWuH", frag_nam))
    cols <- c(cols, grep("fwh", frag_nam))
    if(length(cols) != 0) frag_nam <- frag_nam[-cols]

    observed <- t(data.frame(observed))
    colnames(observed) <- frag_nam
  }
  setwd(WD)
  return(observed)
}
