#' Observed coexpansion hyper-summary statistics (NGS)
#'
#' @description Computes the 24 hyper-summary statistics from empirical PHYLIP
#'   files (one per species/population), matching the output of
#'   \code{\link{sim.coexp.ngs}}. For each species the per-locus summary
#'   statistics are computed via msABC in fragment mode (\code{--obs}), then the
#'   4 moments (mean, variance, kurtosis, skewness) of 6 statistics across
#'   species yield 24 hyper-stats.
#'
#' @param hmodel An hModel object (from \code{h.menu.gui()}). If provided, PHYLIP
#'   paths are taken from \code{hmodel$phylip_paths}.
#' @param phylip_paths Character vector of paths to multi-locus sequential PHYLIP
#'   files, one per species/population. Each file should contain phased data with
#'   invariable sites included.
#'
#' @return A named numeric vector of length 24 with hyper-summary statistics.
#'   The per-species \code{s_mean_*} matrix is attached as the
#'   \code{"species_stats"} attribute.
#'
#' @details At least one of \code{hmodel} or \code{phylip_paths} must be given.
#'   At least 2 species are required for hyper-stats computation.
#'
#'   The 6 per-locus statistics (computed by msABC) are: segs, pi, thetaW, tajd,
#'   nhap, Hd. For each species the mean across loci (\code{s_mean_*}) is extracted,
#'   then the 4 moments across species produce the 24 output values.
#'
#' @author Marcelo Gehara
#' @export
obs.coexp.ngs <- function(hmodel = NULL, phylip_paths = NULL) {

  # Input validation
  if (is.null(hmodel) && is.null(phylip_paths))
    stop("Either 'hmodel' or 'phylip_paths' must be provided")

  if (!is.null(hmodel)) {
    phylip_paths <- hmodel$phylip_paths
    if (is.null(phylip_paths))
      stop("hmodel does not contain 'phylip_paths'")
  }

  if (length(phylip_paths) < 2)
    stop("At least 2 species/populations are needed for coexpansion hyper-stats")

  nspecies <- length(phylip_paths)
  species_stats <- matrix(NA_real_, nrow = nspecies, ncol = 6)
  stat_names <- NULL

  cat("PipeMaster:: Computing observed coexpansion hyper-stats for", nspecies, "species\n")

  for (sp in 1:nspecies) {
    phylip_path <- phylip_paths[sp]
    if (!file.exists(phylip_path))
      stop(paste("PHYLIP file not found:", phylip_path))

    cat("PipeMaster::", sp, "of", nspecies, "-", basename(phylip_path), "\n")

    # --- Parse PHYLIP headers: get per-locus ntax/nchar and all sample names ---
    lines <- readLines(phylip_path)
    locus_info <- list()
    all_samples <- character()
    i <- 1
    locus_id <- 0

    while (i <= length(lines)) {
      if (nchar(trimws(lines[i])) == 0) { i <- i + 1; next }

      header <- as.integer(strsplit(trimws(lines[i]), "\\s+")[[1]])
      ntax <- header[1]; nchar_seq <- header[2]
      locus_id <- locus_id + 1
      i <- i + 1

      for (j in 1:ntax) {
        sname <- trimws(substr(lines[i], 1, 10))
        all_samples <- c(all_samples, sname)
        i <- i + 1
      }

      locus_info[[locus_id]] <- c(ntax = ntax, length = nchar_seq)
    }

    nloci <- length(locus_info)
    all_samples <- unique(all_samples)

    # --- Build pop_assign: all samples -> pop 1 ---
    pop_assign <- data.frame(sample = all_samples, pop = 1,
                             stringsAsFactors = FALSE)

    # --- Build locfile ---
    locfile <- data.frame(
      id     = paste0("loc", seq_len(nloci)),
      n      = vapply(locus_info, function(x) x["ntax"], integer(1)),
      pop    = 1L,
      length = vapply(locus_info, function(x) x["length"], integer(1)),
      mu     = 0,
      rec    = 0,
      stringsAsFactors = FALSE
    )

    # --- Set up temp directory ---
    tmpdir <- tempfile("obs_coexp_")
    dir.create(tmpdir)

    # --- Convert PHYLIP -> ms format ---
    ms_file <- file.path(tmpdir, "combined.ms")
    phylip_to_ms_file(phylip_path, pop_assign, ms_file)

    # --- Write locfile ---
    locfile_path <- file.path(tmpdir, ".1locfile.txt")
    write.table(locfile, locfile_path, row.names = FALSE, col.names = TRUE,
                quote = FALSE, sep = " ")

    # --- Build and run msABC command ---
    ntax_cmd <- locus_info[[1]]["ntax"]
    command <- paste(ntax_cmd, "1 -eN 0 1",
                     "--frag-begin --finp", locfile_path,
                     "--N 100000 --frag-end",
                     "--obs", ms_file)

    result <- run.msABC(command)

    # --- Parse output: header (line 1) + values (line 2) ---
    x <- strsplit(result, "\t")
    frag_nam <- x[[1]]
    values <- as.numeric(x[[2]])

    # Keep only s_mean_* columns (mean of each stat across loci)
    keep <- grep("^s_mean_", frag_nam)
    frag_nam <- frag_nam[keep]
    values <- values[keep]

    # Remove thomson, FayWuH (same filtering as sim.coexp.ngs; ZnS retained)
    cols_rm <- c(grep("thomson", frag_nam),
                 grep("FayWuH", frag_nam))
    if (length(cols_rm) > 0) {
      frag_nam <- frag_nam[-cols_rm]
      values <- values[-cols_rm]
    }

    species_stats[sp, ] <- values
    if (is.null(stat_names)) stat_names <- frag_nam

    # Clean up
    unlink(tmpdir, recursive = TRUE)
  }

  colnames(species_stats) <- stat_names
  rownames(species_stats) <- basename(phylip_paths)

  # --- Compute hyper-summary stats across species (same as coexp.msABC.batch) ---
  average <- colMeans(species_stats)
  vari    <- apply(species_stats, 2, var, na.rm = TRUE)
  kur     <- apply(species_stats, 2, e1071::kurtosis, na.rm = TRUE)
  skew    <- apply(species_stats, 2, e1071::skewness, na.rm = TRUE)

  hyper <- c(average, vari, kur, skew)

  base_names <- sub("^s_mean_", "", stat_names)
  names(hyper) <- c(
    paste0("s_mean_", base_names),
    paste0("s_var_",  base_names),
    paste0("s_kurt_", base_names),
    paste0("s_skew_", base_names)
  )

  attr(hyper, "species_stats") <- species_stats

  cat("PipeMaster:: Done. 24 hyper-summary statistics computed.\n")
  return(hyper)
}
