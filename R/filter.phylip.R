#' Filter and downsample a multi-locus PHYLIP file to match a uniform-sample model
#'
#' @description Writes a new PHYLIP file containing only the loci retained by
#'   \code{optimize.sfs.model()}, with each population downsampled to the
#'   per-population target size encoded in the model. The output PHYLIP has
#'   uniform sample sizes across all loci and is suitable for the
#'   \code{uniform_samples} path of \code{observed.sumstats()}.
#'
#' @param path.in  Path to the input PHYLIP file (multi-locus, interleaved).
#' @param path.out Path to write the filtered PHYLIP file.
#' @param model    A model object produced by \code{optimize.sfs.model()}.
#'   Must have uniform per-population sample sizes across all loci.
#' @param pop.assign A two-column data frame (or path to a tab/whitespace file
#'   with a header) mapping sample names to population numbers.
#' @param method   Downsampling method: \code{"random"} (default) picks
#'   \code{target[p]} samples uniformly at random per population per locus;
#'   \code{"first"} picks the first \code{target[p]} matching samples in the
#'   order they appear in the input PHYLIP.
#' @param seed     Optional integer seed for \code{method = "random"}.
#' @param verbose  Print a progress bar (default \code{TRUE}).
#' @return The output path, invisibly.
#' @author Marcelo Gehara
#' @export
filter.phylip <- function(path.in, path.out, model, pop.assign,
                          method = c("random", "first"), seed = NULL,
                          verbose = TRUE) {

  method <- match.arg(method)
  if (!file.exists(path.in)) stop("path.in does not exist: ", path.in)
  if (!inherits(model, "Model"))
    stop("model must be of class 'Model' (e.g., from optimize.sfs.model()).")

  npop <- as.integer(model$I[1, 3])
  pop_cols <- 4:(3 + npop)
  pop_sizes <- matrix(as.numeric(model$I[, pop_cols]), ncol = npop)
  if (nrow(unique(pop_sizes)) != 1L)
    stop("model$I has non-uniform pop sizes across loci. Run optimize.sfs.model() first.")
  target <- as.integer(pop_sizes[1, ])

  locus_names <- model$loci[, 1]
  locus_idx <- suppressWarnings(as.integer(sub("^locus", "", locus_names)))
  if (any(is.na(locus_idx)))
    stop("model$loci[,1] must contain locus names like 'locus<N>' (PHYLIP-derived models).")
  keep_set <- sort(unique(locus_idx))

  if (is.character(pop.assign) && length(pop.assign) == 1L && file.exists(pop.assign))
    pop.assign <- utils::read.table(pop.assign, header = TRUE, stringsAsFactors = FALSE)
  pop.assign <- data.frame(pop.assign, stringsAsFactors = FALSE)
  if (ncol(pop.assign) < 2L) stop("pop.assign must have at least 2 columns.")
  pop.assign[[1]] <- trimws(pop.assign[[1]])
  pop.assign[[2]] <- as.integer(trimws(as.character(pop.assign[[2]])))
  pop_samples <- lapply(seq_len(npop),
                        function(p) pop.assign[pop.assign[[2]] == p, 1])

  if (!is.null(seed)) set.seed(seed)

  lines <- readLines(path.in)
  n_lines <- length(lines)
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = length(keep_set), style = 3)
    on.exit(close(pb), add = TRUE)
  }

  out_con <- file(path.out, "w")
  on.exit(close(out_con), add = TRUE)

  i <- 1L
  cur_locus <- 0L
  n_written <- 0L
  skipped_short <- integer(0)

  while (i <= n_lines) {
    if (nchar(trimws(lines[i])) == 0L) { i <- i + 1L; next }

    hdr <- as.integer(strsplit(trimws(lines[i]), "\\s+")[[1]])
    ntax <- hdr[1]
    nchar_seq <- hdr[2]
    i <- i + 1L
    cur_locus <- cur_locus + 1L

    if (!(cur_locus %in% keep_set)) {
      i <- i + ntax
      next
    }

    names_vec <- character(ntax)
    seqs <- character(ntax)
    for (j in seq_len(ntax)) {
      line <- lines[i]
      if (grepl("\t", line, fixed = TRUE)) {
        parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
        names_vec[j] <- trimws(parts[1])
        seqs[j] <- gsub("\\s", "", parts[2])
      } else {
        names_vec[j] <- trimws(substr(line, 1, 10))
        seqs[j] <- gsub("\\s", "", substr(line, 11, nchar(line)))
      }
      i <- i + 1L
    }

    sel_names <- character(0)
    sel_seqs  <- character(0)
    short_locus <- FALSE
    for (p in seq_len(npop)) {
      matched <- which(names_vec %in% pop_samples[[p]])
      if (length(matched) < target[p]) {
        short_locus <- TRUE
        break
      }
      keep_idx <- if (method == "random")
        sample(matched, target[p]) else matched[seq_len(target[p])]
      sel_names <- c(sel_names, names_vec[keep_idx])
      sel_seqs  <- c(sel_seqs,  seqs[keep_idx])
    }
    if (short_locus) {
      skipped_short <- c(skipped_short, cur_locus)
      next
    }

    writeLines(sprintf("%d %d", length(sel_names), nchar_seq), out_con)
    for (k in seq_along(sel_names))
      writeLines(paste(sel_names[k], sel_seqs[k], sep = "\t"), out_con)

    n_written <- n_written + 1L
    if (verbose) utils::setTxtProgressBar(pb, n_written)
  }

  if (verbose) cat("\n")
  cat(sprintf("PipeMaster:: filter.phylip wrote %d loci (target per pop: %s) to %s\n",
              n_written, paste(target, collapse = "/"), path.out))
  if (length(skipped_short) > 0L)
    warning(sprintf("%d loci skipped because at least one pop had fewer samples than target (locus indices: %s%s)",
                    length(skipped_short),
                    paste(utils::head(skipped_short, 10), collapse = ", "),
                    if (length(skipped_short) > 10) ", ..." else ""))

  invisible(path.out)
}
