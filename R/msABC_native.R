#' Run msABC coalescent simulator natively
#'
#' Calls the msABC simulator directly via compiled C code, avoiding
#' the overhead of spawning an external process via system2().
#'
#' @param command Character string with msABC arguments (e.g., "80 1 -I 2 40 40 -t 10 -r 5 1000").
#'   This should NOT include the program name — just the arguments.
#' @param seed Optional integer vector of length 3 for reproducible simulations.
#'   If NULL (default), the internal seed state continues from previous calls.
#'
#' @return A character vector of output lines from the msABC simulator.
#'
#' @examples
#' \dontrun{
#' # Basic single-population simulation
#' result <- run.msABC("10 1 -t 5.0")
#' cat(result, sep = "\n")
#'
#' # Reproducible simulation with seed
#' r1 <- run.msABC("10 2 -t 5.0", seed = c(123L, 456L, 789L))
#' r2 <- run.msABC("10 2 -t 5.0", seed = c(123L, 456L, 789L))
#' identical(r1, r2)  # TRUE
#'
#' # Multi-population simulation
#' run.msABC("20 1 -I 2 10 10 -t 10 -ej 1.0 1 2")
#'
#' # Fragment mode (multi-locus) — returns mean, variance, skewness
#' # and kurtosis of each summary statistic across loci
#' locfile <- system.file("extdata", "example_locfile.txt", package = "PipeMaster")
#' cmd <- paste("10 2 -t 5.0 --frag-begin --finp", locfile, "--N 10000 --frag-end")
#' result <- run.msABC(cmd, seed = c(123L, 456L, 789L))
#' tab <- read.table(text = paste(result, collapse = "\n"), header = TRUE, sep = "\t")
#' head(tab)
#' }
#'
#' @useDynLib PipeMaster, msABC_call
#' @export
run.msABC <- function(command, seed = NULL) {
  if (!is.character(command) || length(command) != 1) {
    stop("'command' must be a single character string")
  }
  if (!is.null(seed)) {
    seed <- as.integer(seed)
    if (length(seed) != 3 || any(is.na(seed))) {
      stop("'seed' must be an integer vector of length 3, or NULL")
    }
  }
  output <- .Call("msABC_call", command, seed, PACKAGE = "PipeMaster")
  strsplit(output, "\n", fixed = TRUE)[[1]]
}
