#' Readable labels for model parameters (display only)
#'
#' Builds a human-readable label for every parameter of a Model, using the
#' population names in \code{model$labels$pops} and the population indices stored
#' in each flag row. Internal parameter names are NOT changed anywhere: they stay
#' numeric (\code{Ne0.pop1}, \code{mig0.5_7}, ...), which is what simulation,
#' conditions, reference tables and trained networks rely on. The labels are only
#' for people reading tables and figures.
#'
#' Conventions (ms semantics):
#' \itemize{
#'   \item \code{Ne0.pop1} -> "Ne0 BO1"; \code{t.Ne1.pop4} -> "t.Ne1 BO4".
#'   \item \code{Ne.anc_4_5} -> "Ne.anc BO4+ST" (the join it belongs to).
#'   \item \code{join4_5} -> "join BO4->ST": \code{-ej t 4 5} moves pop 4 into pop 5.
#'   \item \code{mig0.5_7} -> "mig0 ST<-NT": \code{-m 5 7} = pop 5 RECEIVES migrants
#'         from pop 7 (forward in time, individuals move NT -> ST).
#'   \item \code{t.mig1.4_5} -> "t.mig1 BO4<-ST".
#' }
#' A population without a name in \code{model$labels$pops} is shown as "popN".
#' Parameters not found in the model flags keep their own name.
#'
#' @param model A Model object.
#' @param names Optional character vector of parameter names; if given, the result
#'   is returned for exactly these names (unknown names map to themselves).
#' @return A named character vector: names are internal parameter names, values
#'   are display labels.
#' @examples
#' \dontrun{
#' lab <- param.labels(model)
#' plot(posterior, labels = lab)
#' }
#' @export
param.labels <- function(model, names = NULL) {
  pops <- model$labels$pops
  lab <- function(i) {
    i <- as.character(i)
    v <- if (!is.null(pops) && i %in% names(pops)) as.character(pops[[i]]) else ""
    if (is.na(v) || !nzchar(v)) paste0("pop", i) else v
  }
  pair <- function(idx, op) {
    p <- strsplit(trimws(idx), "\\s+")[[1]]
    if (length(p) == 2) paste0(lab(p[1]), op, lab(p[2])) else lab(p[1])
  }
  rows <- function(m) if (is.null(m) || !is.matrix(m) || nrow(m) == 0) NULL else m
  fl <- model$flags
  n  <- rows(fl$n); m <- rows(fl$m); ej <- rows(fl$ej)
  ens <- rows(fl$en$size); ent <- rows(fl$en$time)
  ems <- rows(fl$em$size); emt <- rows(fl$em$time)

  out <- character(0)
  add <- function(nm, v) out[nm] <<- v

  # joins: "join BO4->ST" ; also a lookup pair-key -> "BO4+ST" for Ne.anc rows
  join_pair <- character(0)
  if (!is.null(ej)) for (r in seq_len(nrow(ej))) {
    stem <- sub("_?[0-9]+(_[0-9]+)?$", "", ej[r, 1]); if (!nzchar(stem)) stem <- ej[r, 1]
    add(ej[r, 1], paste(stem, pair(ej[r, 3], "->")))
    join_pair[gsub("\\s+", "_", trimws(ej[r, 3]))] <- pair(ej[r, 3], "+")
  }
  ne_row <- function(nm, idx) {
    if (grepl("^(t\\.)?Ne\\.anc_", nm)) {
      key  <- sub("^(t\\.)?Ne\\.anc_", "", nm)
      stem <- sub("^((t\\.)?Ne\\.anc).*$", "\\1", nm)
      if (key %in% names(join_pair)) return(paste(stem, join_pair[[key]]))
      return(paste(nm, lab(idx)))
    }
    stem <- sub("\\.pop[0-9]+$", "", nm)
    if (identical(stem, nm)) paste(nm, lab(idx)) else paste(stem, lab(idx))
  }
  if (!is.null(n))   for (r in seq_len(nrow(n)))   add(n[r, 1],   ne_row(n[r, 1],   n[r, 3]))
  if (!is.null(ens)) for (r in seq_len(nrow(ens))) add(ens[r, 1], ne_row(ens[r, 1], ens[r, 3]))
  if (!is.null(ent)) for (r in seq_len(nrow(ent))) add(ent[r, 1], ne_row(ent[r, 1], ent[r, 3]))

  # migration rates: "-m i j" = pop i receives from pop j  ->  "mig0 i<-j"
  mig_stem <- function(nm) { s <- sub("\\.[0-9]+_[0-9]+$", "", nm); if (identical(s, nm)) nm else s }
  if (!is.null(m))   for (r in seq_len(nrow(m)))   add(m[r, 1],   paste(mig_stem(m[r, 1]),   pair(m[r, 3], "<-")))
  if (!is.null(ems)) for (r in seq_len(nrow(ems))) add(ems[r, 1], paste(mig_stem(ems[r, 1]), pair(ems[r, 3], "<-")))
  # ancestral migration TIMES carry only the receiving pop; take the pair from the
  # matching rate row (t.mig1.i_j <-> mig1.i_j)
  if (!is.null(emt)) for (r in seq_len(nrow(emt))) {
    nm <- emt[r, 1]; rate <- sub("^t\\.", "", nm)
    idx <- if (!is.null(ems) && rate %in% ems[, 1]) ems[ems[, 1] == rate, 3][1] else emt[r, 3]
    add(nm, paste(mig_stem(nm), pair(idx, "<-")))
  }

  if (is.null(names)) return(out)
  res <- stats::setNames(as.character(names), names)
  hit <- names %in% names(out)
  res[hit] <- out[names[hit]]
  res
}

## Resolve display labels for a posterior-type object (internal).
## `labels` may be a Model, a named character vector, or NULL (then x$labels).
.resolve.param.labels <- function(x, labels, param_names) {
  if (is.null(labels)) labels <- x$labels
  if (is.null(labels)) return(NULL)
  if (inherits(labels, "Model")) labels <- param.labels(labels)
  if (!is.character(labels) || is.null(names(labels))) return(NULL)
  res <- stats::setNames(as.character(param_names), param_names)
  hit <- param_names %in% names(labels)
  res[hit] <- labels[param_names[hit]]
  res
}
