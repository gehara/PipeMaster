#' Reduce a multi-population joint SFS by marginalization or binning
#'
#' Converts a high-dimensional joint SFS (e.g., 3-pop with 69K bins) to a
#' lower-dimensional representation suitable for ABC or SML.
#'
#' Accepts either an R object (numeric vector or matrix) or a file path to a
#' reftable. When a file path is provided, the reftable is processed in batches
#' to avoid loading the full file into memory.
#'
#' @param sfs Numeric vector (single SFS), matrix (nsim x n_bins), or
#'   character file path to a tab-separated reftable with header.
#'   For vectors/matrices, column order must match the
#'   expand.grid(0:n1, 0:n2, ...) convention (pop1 index varies fastest),
#'   as produced by sim.sumstats() and sim.scrm.sumstats().
#' @param sample_sizes Integer vector of haploid sample sizes per population.
#' @param method Character: "marginal2d" for pairwise 2D marginals, or
#'   "binned" for coarsened SFS with adjacent bins grouped.
#' @param bin_size Integer bin width for "binned" method (default 4).
#' @param output Character path to output file or directory. Only used when
#'   \code{sfs} is a file path. If a directory, output is named
#'   \code{SIMS_reftable_reduced.txt} inside it.
#' @param batch_size Number of rows per batch when reading from file
#'   (default 1000).
#' @param overwrite Logical: if FALSE (default) and output file exists,
#'   return cached result. Only used when \code{sfs} is a file path.
#' @param verbose Logical: print progress (default TRUE).
#'
#' @return For vector input: named numeric vector. For matrix input: matrix
#'   with column names. For file input: path to the reduced reftable file
#'   (invisibly).
#'
#' @details
#' \strong{marginal2d}: For npop populations, computes C(npop,2) pairwise
#' 2D marginal SFS by summing over all other dimensions. Requires npop >= 3.
#'
#' \strong{binned}: Groups adjacent allele count bins in each dimension.
#' Works for any number of populations.
#'
#' When \code{sfs} is a file path, non-SFS columns (parameters, summary
#' statistics) are preserved. Memory usage per batch:
#' \code{batch_size * n_sfs_cols * 8} bytes.
#'
#' When \code{sfs} is an R object larger than 2 GB, the function stops with
#' a message suggesting file-based input instead.
#'
#' @examples
#' \dontrun{
#' # In-memory: single observation
#' obs_reduced <- reduce.sfs(obs_sfs_vec, sample_sizes = c(40, 40, 40),
#'                           method = "marginal2d")
#'
#' # In-memory: matrix
#' reduced <- reduce.sfs(sfs_mat, sample_sizes = c(40, 40, 40),
#'                       method = "marginal2d")
#'
#' # File-based: large reftable (69K SFS cols)
#' reduced_file <- reduce.sfs(
#'   "data/SIMS_reftable_100K.txt",
#'   sample_sizes = c(40, 40, 40),
#'   method = "marginal2d",
#'   output = "results/"
#' )
#' reftable <- read.table(reduced_file, header = TRUE)
#' }
#'
#' @export
reduce.sfs <- function(sfs, sample_sizes, method = c("marginal2d", "binned"),
                        bin_size = 4L, output = NULL, batch_size = 1000L,
                        overwrite = FALSE, verbose = TRUE) {
  method <- match.arg(method)

  # --- Dispatch: file path vs R object ---
  if (is.character(sfs) && length(sfs) == 1L) {
    return(.reduce.sfs.from.file(
      input = sfs, output = output, sample_sizes = sample_sizes,
      method = method, bin_size = bin_size, batch_size = batch_size,
      overwrite = overwrite))
  }

  # --- R object path ---
  # Size check: stop if object > 2 GB
  obj_bytes <- object.size(sfs)
  if (obj_bytes > 2e9)
    stop(sprintf(
      "SFS object is %.1f GB — too large for in-memory reduction. Pass a file path instead.",
      as.numeric(obj_bytes) / 1e9))

  npop <- length(sample_sizes)
  dims <- as.integer(sample_sizes + 1L)
  n_total <- prod(dims)

  # Handle vector input (single observation)
  single <- is.null(dim(sfs))
  if (single) sfs <- matrix(sfs, nrow = 1)

  if (ncol(sfs) != n_total)
    stop(sprintf("SFS has %d columns but sample_sizes imply %d bins (prod(%s))",
                 ncol(sfs), n_total, paste(dims, collapse = " x ")))

  nsim <- nrow(sfs)

  if (method == "marginal2d") {
    out <- .reduce.sfs.marginal2d(sfs, nsim, npop, dims, n_total, verbose)
  } else {
    out <- .reduce.sfs.binned(sfs, nsim, npop, dims, n_total, bin_size, verbose)
  }

  if (single) return(setNames(as.numeric(out), colnames(out)))
  return(out)
}


# ============================================================================
# Internal: marginal 2D reduction
# ============================================================================
.reduce.sfs.marginal2d <- function(sfs, nsim, npop, dims, n_total, verbose) {
  if (npop < 3) {
    if (verbose) cat("PipeMaster:: reduce.sfs: marginal2d requires >= 3 pops, returning unchanged\n")
    return(sfs)
  }

  grid <- do.call(expand.grid, lapply(dims, function(d) seq(0L, d - 1L)))
  pairs <- combn(npop, 2)
  results <- vector("list", ncol(pairs))

  if (verbose) {
    n_groups_total <- sum(sapply(seq_len(ncol(pairs)), function(p)
      dims[pairs[1, p]] * dims[pairs[2, p]]))
    pb <- txtProgressBar(min = 0, max = n_groups_total, style = 3)
    pb_count <- 0L
  }

  for (p in seq_len(ncol(pairs))) {
    a <- pairs[1, p]; b <- pairs[2, p]
    Da <- dims[a]; Db <- dims[b]
    n_2d <- Da * Db

    map_2d <- grid[, a] + grid[, b] * Da
    groups <- split(seq_len(n_total), map_2d)

    marginal <- matrix(0, nrow = nsim, ncol = n_2d)
    for (g in names(groups)) {
      cols <- groups[[g]]
      idx <- as.integer(g) + 1L
      if (length(cols) == 1L) {
        marginal[, idx] <- sfs[, cols]
      } else {
        marginal[, idx] <- rowSums(sfs[, cols, drop = FALSE])
      }
      if (verbose) { pb_count <- pb_count + 1L; setTxtProgressBar(pb, pb_count) }
    }

    grid_2d <- expand.grid(seq(0L, Da - 1L), seq(0L, Db - 1L))
    colnames(marginal) <- paste0("sfs2d_", a, "_", b, "_",
                                 grid_2d[, 1], "_", grid_2d[, 2])
    results[[p]] <- marginal
  }

  if (verbose) { close(pb); cat("\n") }

  out <- do.call(cbind, results)
  if (verbose)
    cat(sprintf("PipeMaster:: reduce.sfs: %dD SFS (%d bins) -> %d pairwise 2D marginals (%d bins)\n",
                npop, n_total, ncol(pairs), ncol(out)))
  out
}


# ============================================================================
# Internal: binned reduction
# ============================================================================
.reduce.sfs.binned <- function(sfs, nsim, npop, dims, n_total, bin_size, verbose) {
  bin_size <- as.integer(bin_size)
  new_dims <- as.integer(ceiling(dims / bin_size))
  n_new <- prod(new_dims)

  grid <- do.call(expand.grid, lapply(dims, function(d) seq(0L, d - 1L)))
  binned_grid <- floor(grid / bin_size)

  strides <- c(1L, cumprod(new_dims[-npop]))
  map_bin <- as.integer(rowSums(t(t(binned_grid) * strides)))

  groups <- split(seq_len(n_total), map_bin)

  if (verbose) pb <- txtProgressBar(min = 0, max = length(groups), style = 3)

  out <- matrix(0, nrow = nsim, ncol = n_new)
  gi <- 0L
  for (g in names(groups)) {
    cols <- groups[[g]]
    idx <- as.integer(g) + 1L
    if (length(cols) == 1L) {
      out[, idx] <- sfs[, cols]
    } else {
      out[, idx] <- rowSums(sfs[, cols, drop = FALSE])
    }
    if (verbose) { gi <- gi + 1L; setTxtProgressBar(pb, gi) }
  }

  if (verbose) { close(pb); cat("\n") }

  binned_unique <- do.call(expand.grid, lapply(new_dims, function(d) seq(0L, d - 1L)))
  colnames(out) <- apply(binned_unique, 1, function(x)
    paste0("sfs_bin_", paste(x, collapse = "_")))

  if (verbose)
    cat(sprintf("PipeMaster:: reduce.sfs: %dD SFS (%d bins) -> binned (%s = %d bins, bin_size=%d)\n",
                npop, n_total, paste(new_dims, collapse = " x "), n_new, bin_size))
  out
}


# ============================================================================
# Internal: file-based batched reduction
# ============================================================================
.reduce.sfs.from.file <- function(input, output, sample_sizes, method,
                                   bin_size, batch_size, overwrite) {

  input <- path.expand(input)
  if (!file.exists(input))
    stop(sprintf("Input file not found: %s", input))

  if (is.null(output))
    stop("'output' path is required when 'sfs' is a file path")

  output <- path.expand(output)
  if (dir.exists(output)) {
    output_file <- file.path(output, "SIMS_reftable_reduced.txt")
  } else {
    output_file <- output
  }

  if (!overwrite && file.exists(output_file)) {
    cat(sprintf("PipeMaster:: reduce.sfs: using cached %s\n", output_file))
    return(invisible(output_file))
  }

  # Read header
  header_line <- readLines(input, n = 1)
  all_cols <- strsplit(header_line, "\t")[[1]]
  sfs_cols <- grep("^sfs_\\d", all_cols, value = TRUE)
  nonsfs_cols <- setdiff(all_cols, sfs_cols)

  if (length(sfs_cols) == 0) {
    cat("PipeMaster:: reduce.sfs: no SFS columns found, copying file unchanged\n")
    file.copy(input, output_file, overwrite = TRUE)
    return(invisible(output_file))
  }

  n_sfs <- length(sfs_cols)
  n_expected <- prod(sample_sizes + 1L)
  if (n_sfs != n_expected)
    stop(sprintf("Found %d SFS columns but sample_sizes imply %d bins", n_sfs, n_expected))

  n_total <- as.integer(system(sprintf("wc -l < '%s'", input), intern = TRUE)) - 1L
  cat(sprintf("PipeMaster:: reduce.sfs: %s total cols (%d SFS + %d other), %s rows\n",
              format(length(all_cols), big.mark = ","), n_sfs, length(nonsfs_cols),
              format(n_total, big.mark = ",")))

  if (file.exists(output_file)) file.remove(output_file)

  sfs_idx <- as.integer(match(sfs_cols, all_cols))
  nonsfs_idx <- as.integer(match(nonsfs_cols, all_cols))

  # Process in batches using C reader for speed
  n_batches <- ceiling(n_total / batch_size)
  cat(sprintf("PipeMaster:: reduce.sfs: %d batches x %d rows, method=%s\n",
              n_batches, batch_size, method))

  offset <- 0L
  batch_num <- 0L
  t0 <- proc.time()

  while (offset < n_total) {
    n_read <- min(batch_size, n_total - offset)

    # Read SFS columns via C (skip = 1 + offset for header + previous rows)
    sfs_mat <- .Call("read_tsv_call", input, sfs_idx,
                     as.integer(1L + offset), as.integer(n_read),
                     PACKAGE = "PipeMaster")

    # Read non-SFS columns via C
    nonsfs_mat <- .Call("read_tsv_call", input, nonsfs_idx,
                        as.integer(1L + offset), as.integer(n_read),
                        PACKAGE = "PipeMaster")

    # Replace NAs with 0
    sfs_mat[is.na(sfs_mat)] <- 0
    nonsfs_mat[is.na(nonsfs_mat)] <- 0

    # Reduce SFS
    reduced <- reduce.sfs(sfs_mat, sample_sizes, method = method,
                          bin_size = bin_size, verbose = FALSE)
    if (is.null(dim(reduced))) reduced <- matrix(reduced, nrow = 1)

    # Combine and write
    out <- cbind(as.data.frame(nonsfs_mat), reduced)
    colnames(out) <- c(nonsfs_cols, colnames(reduced))

    suppressWarnings(
      write.table(out, file = output_file, append = (batch_num > 0L),
                  col.names = (batch_num == 0L), row.names = FALSE,
                  quote = FALSE, sep = "\t")
    )

    offset <- offset + n_read
    batch_num <- batch_num + 1L

    elapsed_so_far <- (proc.time() - t0)[3]
    rate <- offset / elapsed_so_far
    remaining <- (n_total - offset) / rate
    pct <- 100 * offset / n_total

    cat(sprintf("\rPipeMaster:: reduce.sfs: %s/%s rows (%.0f%%) | %.0f rows/h | ~%.2f h remaining    ",
                format(offset, big.mark = ","),
                format(n_total, big.mark = ","),
                pct, rate * 3600, remaining / 3600))
    flush(stdout())
  }

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("\nPipeMaster:: reduce.sfs: done — %s rows in %.2f h (%.0f rows/h) -> %s\n",
              format(offset, big.mark = ","), elapsed / 3600, offset / elapsed * 3600, output_file))

  return(invisible(output_file))
}
