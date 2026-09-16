#' Posterior predictive check
#'
#' @description Simulates at the posterior draws and asks whether the observed
#'   statistics fall inside the predictive distribution.
#'
#' @details \strong{This is not optional, and no internal metric replaces it.}
#'   Validation R2 measures how well the network INVERTS the model; it says
#'   nothing about whether the model can PRODUCE the data. \code{OOD.pretrain} and
#'   \code{OOD.posttrain} have the same blind spot in reverse -- they test whether
#'   a model can produce the data, never whether it recovers the right parameter.
#'   A model can pass every OOD check, reach high R2 on every parameter, and still
#'   put most fitted statistics outside their 95% predictive interval.
#'
#'   Reported as PIT values: the fraction of simulations below the observation.
#'   0.5 is dead centre; near 0 or 1 means the observation sits outside. Grouped
#'   by statistic family, because failures are rarely uniform -- ratio statistics
#'   such as Fst survive misspecification that destroys count statistics such as
#'   private/shared.
#'
#'   \strong{Do not respond to a failure by dropping the offending statistics.}
#'   That makes the check unable to fail. Drop a statistic only for an a-priori
#'   reason -- being undefined for the sample configuration, say -- never because
#'   of its residual.
#'
#'   Every parameter of the model is fixed per draw, not only those carried by the
#'   posterior. Parameters tied by \code{"="} conditions are resolved to their tie
#'   group's representative first: \code{sample.w.cond} applies \code{"="} as
#'   \code{x[param2] <- x[param1]}, so a still-free tied parameter would otherwise
#'   overwrite a fixed one. Draws violating a \code{"<"} condition are dropped and
#'   counted rather than silently resampled.
#'
#' @param posterior a \code{posterior} object from \code{nn.predict()} or
#'   \code{abc.rejection()}.
#' @param model the \code{Model} object the posterior was estimated under.
#' @param observed named numeric vector of observed statistics.
#' @param draws number of POSTERIOR SAMPLES to simulate at. Default 1000.
#' @param reps simulations per draw. Total simulations are \code{draws * reps}.
#' @param method which posterior to draw from, \code{"abc"} or \code{"bootstrap"}.
#' @param sampling how parameter vectors are taken. \code{"density"} samples a
#'   variance-preserving kernel density over the accepted samples, giving distinct
#'   vectors that fill the space between them; \code{"points"} cycles the accepted
#'   samples themselves, which is what the posterior literally is but reuses each
#'   point when more draws are asked for than are available.
#' @param bounds clip sampled vectors to the prior bounds. A kernel places mass
#'   outside a bounded prior; without clipping a check can simulate parameters the
#'   prior forbids.
#' @param reftable optional reference table. When supplied the check also projects
#'   the prior-predictive cloud, the posterior-predictive cloud and the observation
#'   onto axes fitted on the reference table alone, and that projection becomes
#'   available as \code{plot(x, type = "pc")}.
#' @param engine \code{"scrm"} or \code{"msABC"}.
#' @param mu.rates,rec.rates rate specifications, as for \code{sim.scrm.sumstats()}.
#' @param stat.config statistics population structure, as for \code{sim.scrm.sumstats()}.
#' @param skip.zns,skip.sfs passed to the simulator.
#' @param ncores parallel worker processes (independent Rscript workers).
#' @param path working directory for simulation output.
#' @param seed RNG seed.
#' @param verbose progress messages.
#' @return an object of class \code{ppc}: simulated statistics, the observed
#'   vector, PIT values, and the parameter draw behind each simulation.
#' @author Marcelo Gehara
#' @export
posterior.predictive.check <- function(posterior, model, observed,
                                       draws = 1000L, reps = 1L,
                                       method = c("abc", "bootstrap"),
                                       sampling = c("density", "points"),
                                       bounds = TRUE, reftable = NULL,
                                       engine = c("scrm", "msABC"),
                                       mu.rates = NULL, rec.rates = NULL,
                                       stat.config = NULL,
                                       skip.zns = TRUE, skip.sfs = FALSE,
                                       ncores = 1L, path = tempdir(),
                                       seed = 42L, verbose = TRUE) {

  method <- match.arg(method); sampling <- match.arg(sampling)
  engine <- match.arg(engine)
  if (!inherits(posterior, "posterior"))
    stop("`posterior` must be a posterior object from nn.predict() or abc.rejection()")
  P <- posterior[[method]]
  if (is.null(P))
    stop(sprintf("the posterior carries no `%s` samples (available: %s)", method,
                 paste(intersect(c("bootstrap","abc"), names(posterior)), collapse = ", ")))
  P <- as.matrix(P); draws <- as.integer(draws); reps <- as.integer(reps)
  set.seed(seed)

  if (verbose)
    cat(sprintf("PipeMaster:: posterior.predictive.check -- %d draws x %d reps = %d simulations\n",
                draws, reps, draws * reps),
        sprintf("PipeMaster::   posterior '%s' (%d samples, %d parameters), sampling = '%s'\n",
                method, nrow(P), ncol(P), sampling), sep = "")

  theta <- .ppc.draw.theta(P, draws, sampling, bounds,
                           prior = posterior$prior, verbose = verbose)

  ## every model parameter -> its tie-group representative
  tie <- .ppc.tie.map(model, verbose = verbose)

  sims <- .ppc.simulate(model, theta, tie, reps = reps, engine = engine,
                        mu.rates = mu.rates, rec.rates = rec.rates,
                        stat.config = stat.config, skip.zns = skip.zns,
                        skip.sfs = skip.sfs, ncores = ncores, path = path,
                        verbose = verbose)
  S <- sims$stats; kept <- sims$kept

  sc <- intersect(names(observed), colnames(S))
  sc <- sc[vapply(sc, function(k) is.finite(observed[[k]]) &&
                    any(is.finite(S[, k])), logical(1))]
  S <- S[, sc, drop = FALSE]; obs <- observed[sc]

  pit <- vapply(sc, function(k) {
    v <- S[, k]; v <- v[is.finite(v)]
    if (!length(v)) return(NA_real_)
    (sum(v < obs[[k]]) + 0.5 * sum(v == obs[[k]])) / length(v) }, numeric(1))

  out <- list(sim = S, obs = obs, pit = pit, stats = sc,
              params = theta[kept, , drop = FALSE],
              n_sim = nrow(S), n_req = draws * reps,
              method = method, sampling = sampling)
  if (!is.null(reftable))
    out$pc <- .ppc.project(reftable, S, obs, sc, verbose = verbose)
  class(out) <- "ppc"
  if (verbose) print(out)
  out
}

## ---------------------------------------------------------------------------
##  internals
## ---------------------------------------------------------------------------

## Draw `draws` parameter vectors from the posterior.
.ppc.draw.theta <- function(P, draws, sampling, bounds, prior = NULL, verbose = TRUE) {
  if (sampling == "points") {
    idx <- ((seq_len(draws) - 1L) %% nrow(P)) + 1L
    th  <- P[idx, , drop = FALSE]
    if (verbose && draws > nrow(P))
      cat(sprintf("PipeMaster::   'points' reuses each sample ~%.1fx (%d draws from %d samples)\n",
                  draws / nrow(P), draws, nrow(P)))
  } else {
    ## variance-preserving kernel: shrink toward the mean, then add kernel noise,
    ## so the sampled cloud keeps the posterior's covariance instead of inflating it
    mu  <- colMeans(P)
    Sig <- stats::cov(P)
    d   <- ncol(P); n <- nrow(P)
    h2  <- n^(-2 / (d + 4))                       # Silverman-type bandwidth
    Sig <- Sig + diag(1e-12 * max(diag(Sig)) + 1e-300, d)
    L   <- tryCatch(chol(Sig * h2), error = function(e) {
      ev <- eigen(Sig * h2, symmetric = TRUE)
      ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), d) %*% t(ev$vectors) })
    i   <- sample(n, draws, replace = TRUE)
    Z   <- matrix(stats::rnorm(draws * d), draws, d)
    th  <- sweep(sqrt(1 - h2) * sweep(P[i, , drop = FALSE], 2, mu), 2, mu, "+") + Z %*% L
    colnames(th) <- colnames(P)
  }
  if (bounds && !is.null(prior)) {
    pr <- as.matrix(prior); nclip <- 0L
    for (k in colnames(th)) if (k %in% colnames(pr)) {
      lo <- min(pr[, k]); hi <- max(pr[, k])
      nclip <- nclip + sum(th[, k] < lo | th[, k] > hi)
      th[, k] <- pmax(lo, pmin(hi, th[, k]))
    }
    if (verbose && nclip > 0)
      cat(sprintf("PipeMaster::   clipped %d parameter values to the prior bounds\n", nclip))
  }
  th
}

## Map every model parameter to its tie-group representative.
## sample.w.cond applies "=" as x[param2] <- x[param1], so fixing only the free
## parameters lets a still-free tied parameter overwrite a fixed one.
.ppc.tie.map <- function(model, n = 200L, verbose = TRUE) {
  D <- do.call(rbind, lapply(seq_len(n), function(i) {
    cm <- msABC.commander(model, use.alpha = NULL, arg = "")[[2]]
    stats::setNames(as.numeric(cm[2, ]), cm[1, ]) }))
  key <- apply(D, 2, function(c) paste(format(c, digits = 15), collapse = "|"))
  grp <- split(names(key), key)
  map <- stats::setNames(rep(NA_character_, ncol(D)), colnames(D))
  for (g in grp) map[g] <- g[1]
  if (verbose) {
    nt <- sum(lengths(grp) > 1)
    cat(sprintf("PipeMaster::   %d model parameters in %d tie groups (%d groups are tied)\n",
                ncol(D), length(grp), nt))
  }
  map
}

## Write one parameter vector into the model as point priors.
.ppc.fix.model <- function(model, vals) {
  pin <- function(M) {
    if (is.null(M) || !nrow(M)) return(M)
    for (i in seq_len(nrow(M))) {
      nm <- M[i, 1]
      if (!is.na(vals[nm]) && is.finite(vals[nm])) {
        v <- format(vals[[nm]], scientific = TRUE)
        M[i, 4] <- v; M[i, 5] <- v
      }
    }
    M
  }
  model$flags$n       <- pin(model$flags$n)
  model$flags$m       <- pin(model$flags$m)
  model$flags$ej      <- pin(model$flags$ej)
  if (!is.null(model$flags$en)) {
    model$flags$en$size <- pin(model$flags$en$size)
    model$flags$en$time <- pin(model$flags$en$time)
  }
  if (!is.null(model$flags$em)) {
    model$flags$em$size <- pin(model$flags$em$size)
    model$flags$em$time <- pin(model$flags$em$time)
  }
  model
}

## Does a fixed vector satisfy the model's "<" conditions?
.ppc.conds.ok <- function(vals, model) {
  chk <- function(cl) {
    if (is.null(cl) || !length(cl)) return(TRUE)
    for (c in cl) {
      if (is.null(c$op) || c$op != "<") next
      a <- vals[c$param1]; b <- vals[c$param2]
      if (is.na(a) || is.na(b)) next
      if (!(a < b)) return(FALSE)
    }
    TRUE
  }
  chk(model$conds$time) && chk(model$conds$size) && chk(model$conds$mig)
}

## Simulate at each fixed parameter vector. Draws are chunked across independent
## Rscript workers (PM convention -- not mclapply), each looping its slice with
## ncores = 1.
.ppc.simulate <- function(model, theta, tie, reps, engine, mu.rates, rec.rates,
                          stat.config, skip.zns, skip.sfs, ncores, path, verbose) {

  dir.create(path, showWarnings = FALSE, recursive = TRUE)
  abs_path <- normalizePath(path, mustWork = TRUE)

  ## expand each draw over the tie map, and screen "<" conditions
  full <- vector("list", nrow(theta)); ok <- logical(nrow(theta))
  for (i in seq_len(nrow(theta))) {
    v <- stats::setNames(rep(NA_real_, length(tie)), names(tie))
    for (p in names(tie)) if (tie[[p]] %in% colnames(theta)) v[p] <- theta[i, tie[[p]]]
    full[[i]] <- v; ok[i] <- .ppc.conds.ok(v, model)
  }
  if (verbose && any(!ok))
    cat(sprintf("PipeMaster::   %d of %d draws violate a '<' condition and are DROPPED\n",
                sum(!ok), length(ok)))
  keep <- which(ok)
  if (!length(keep)) stop("every posterior draw violates a model condition")

  one <- function(i, wpath) {
    m <- .ppc.fix.model(model, full[[i]])
    nm <- sprintf("ppc_%06d", i)
    if (engine == "scrm") {
      sim.scrm.sumstats(model = m, nsims = reps, batch.size = reps,
                        mu.rates = mu.rates, rec.rates = rec.rates,
                        skip.zns = skip.zns, skip.sfs = skip.sfs,
                        stat.config = stat.config, ncores = 1,
                        path = wpath, output.name = nm, verbose = FALSE)
    } else {
      sim.sumstats(model = m, nsims = reps, batch.size = reps,
                   mu.rates = mu.rates, rec.rates = rec.rates,
                   path = wpath, output.name = nm, ncores = 1)
    }
    f <- file.path(wpath, paste0("SIMS_", nm, ".txt"))
    r <- utils::read.table(f, header = TRUE, sep = "\t", check.names = FALSE)
    unlink(f)
    r
  }

  if (ncores <= 1) {
    res <- vector("list", length(keep))
    for (j in seq_along(keep)) {
      res[[j]] <- tryCatch(one(keep[j], abs_path), error = function(e) NULL)
      if (verbose && (j %% max(1L, floor(length(keep) / 20)) == 0))
        cat(sprintf("PipeMaster::   %d/%d draws simulated\n", j, length(keep)))
    }
  } else {
    slices <- split(keep, rep(seq_len(ncores), length.out = length(keep)))
    save(model, full, reps, engine, mu.rates, rec.rates, stat.config,
         skip.zns, skip.sfs, slices,
         file = file.path(abs_path, ".PM_ppc_worker_params.RData"))
    ws <- paste(
      'args <- commandArgs(TRUE); w <- as.integer(args[1])',
      'suppressMessages(library(PipeMaster))',
      sprintf('bp <- "%s"', abs_path),
      'load(file.path(bp, ".PM_ppc_worker_params.RData"))',
      'wd <- file.path(bp, paste0(".ppc_worker_", w)); dir.create(wd, showWarnings = FALSE)',
      'out <- list()',
      'for (i in slices[[w]]) {',
      '  m <- PipeMaster:::.ppc.fix.model(model, full[[i]])',
      '  nm <- sprintf("ppc_%06d", i)',
      '  r <- tryCatch({',
      '    if (engine == "scrm")',
      '      sim.scrm.sumstats(model=m, nsims=reps, batch.size=reps, mu.rates=mu.rates,',
      '        rec.rates=rec.rates, skip.zns=skip.zns, skip.sfs=skip.sfs,',
      '        stat.config=stat.config, ncores=1, path=wd, output.name=nm, verbose=FALSE)',
      '    else sim.sumstats(model=m, nsims=reps, batch.size=reps, mu.rates=mu.rates,',
      '        rec.rates=rec.rates, path=wd, output.name=nm, ncores=1)',
      '    f <- file.path(wd, paste0("SIMS_", nm, ".txt"))',
      '    x <- utils::read.table(f, header=TRUE, sep="\\t", check.names=FALSE)',
      '    unlink(f); x }, error = function(e) NULL)',
      '  out[[as.character(i)]] <- r',
      '}',
      'saveRDS(out, file.path(bp, paste0(".ppc_worker_", w, ".rds")))',
      'write("done", file.path(bp, paste0(".ppc_worker_", w, ".done")))',
      'quit(save="no")', sep = "\n")
    writeLines(ws, file.path(abs_path, ".PM_ppc_worker.R"))
    for (w in seq_len(ncores)) {
      unlink(file.path(abs_path, paste0(".ppc_worker_", w, ".done")))
      system(paste("Rscript", file.path(abs_path, ".PM_ppc_worker.R"), w,
                   ">", file.path(abs_path, sprintf(".ppc_worker_%d.log", w)), "2>&1 &"))
    }
    if (verbose) cat(sprintf("PipeMaster::   launched %d worker processes\n", ncores))
    repeat {
      Sys.sleep(5)
      nd <- sum(file.exists(file.path(abs_path, paste0(".ppc_worker_", seq_len(ncores), ".done"))))
      if (verbose) cat(sprintf("PipeMaster::   %d/%d workers done\n", nd, ncores))
      if (nd == ncores) break
    }
    res <- unlist(lapply(seq_len(ncores), function(w) {
      f <- file.path(abs_path, paste0(".ppc_worker_", w, ".rds"))
      if (file.exists(f)) readRDS(f) else list() }), recursive = FALSE)
    res <- res[order(as.integer(names(res)))]
    unlink(file.path(abs_path, c(".PM_ppc_worker.R", ".PM_ppc_worker_params.RData")))
    unlink(file.path(abs_path, paste0(".ppc_worker_", seq_len(ncores))), recursive = TRUE)
    unlink(file.path(abs_path, paste0(".ppc_worker_", seq_len(ncores), c(".rds", ".done", ".log"))))
  }

  good <- !vapply(res, is.null, logical(1))
  if (!any(good)) stop("no posterior predictive simulation succeeded")
  if (verbose && any(!good))
    cat(sprintf("PipeMaster::   %d of %d simulations FAILED\n", sum(!good), length(good)))
  cols <- Reduce(intersect, lapply(res[good], colnames))
  S <- do.call(rbind, lapply(res[good], function(r) r[, cols, drop = FALSE]))
  list(stats = as.matrix(S), kept = rep(keep[good], each = reps))
}

## Project prior-predictive, posterior-predictive and observed onto axes fitted
## on the REFERENCE TABLE alone, so all three are measured against one basis.
.ppc.project <- function(reftable, SIM, obs, sc, k = 5L, verbose = TRUE) {
  R <- as.data.frame(reftable)
  sc <- intersect(sc, colnames(R))
  if (length(sc) < 3) return(NULL)
  X <- as.matrix(R[, sc, drop = FALSE])
  keep <- apply(X, 2, function(c) all(is.finite(c)) && stats::sd(c) > 0)
  sc <- sc[keep]; X <- X[, sc, drop = FALSE]
  mu <- colMeans(X); sdv <- apply(X, 2, stats::sd)
  z  <- function(M) sweep(sweep(M[, sc, drop = FALSE], 2, mu), 2, sdv, "/")
  pca <- stats::prcomp(z(X), center = FALSE, scale. = FALSE)
  nc  <- max(2L, min(20L, sum(cumsum(pca$sdev^2) / sum(pca$sdev^2) < 0.95) + 1L))
  Pp  <- z(X)            %*% pca$rotation[, seq_len(nc), drop = FALSE]
  Ps  <- z(SIM)          %*% pca$rotation[, seq_len(nc), drop = FALSE]
  Po  <- z(matrix(obs[sc], nrow = 1, dimnames = list(NULL, sc))) %*%
         pca$rotation[, seq_len(nc), drop = FALSE]
  knn <- function(Q, Ref, self) {
    if (requireNamespace("RANN", quietly = TRUE)) {
      r <- RANN::nn2(Ref, Q, k = min(k + self, nrow(Ref)))$nn.dists
      rowMeans(r[, (1 + self):(k + self), drop = FALSE])
    } else {
      apply(Q, 1, function(q) {
        d <- sqrt(colSums((t(Ref) - q)^2)); mean(sort(d)[(1 + self):(k + self)]) })
    }
  }
  if (verbose) cat(sprintf("PipeMaster::   projection: %d stats -> %d PCs (95%% var), k=%d\n",
                           length(sc), nc, k))
  list(pca = list(x.prior = Pp[, 1], x.post = Ps[, 1], x.obs = Po[, 1],
                  d.prior = knn(Pp, Pp, 1L), d.post = knn(Ps, Pp, 0L),
                  d.obs = knn(Po, Pp, 0L),
                  lab = sprintf("PC1 (%.1f%% var)", 100 * pca$sdev[1]^2 / sum(pca$sdev^2))),
       k = k, n.prior = nrow(Pp), n.post = nrow(Ps))
}

#' @export
print.ppc <- function(x, ...) {
  cat(sprintf("\nPosterior predictive check -- %d simulations at %d posterior draws ('%s', %s)\n",
              x$n_sim, length(unique(x$params[, 1])), x$method, x$sampling))
  o <- 100 * mean(x$pit < 0.025 | x$pit > 0.975, na.rm = TRUE)
  cat(sprintf("  %d statistics | %.1f%% outside their 95%% predictive interval\n",
              length(x$stats), o))
  fam <- sub("^(s_[a-z]+_[A-Za-z]+).*$", "\\1", x$stats)
  tab <- table(fam); tab <- tab[tab >= 3]
  if (length(tab)) {
    sev <- vapply(names(tab), function(f) {
      v <- x$pit[fam == f]; 100 * mean(v < 0.025 | v > 0.975, na.rm = TRUE) }, numeric(1))
    sev <- sort(sev, decreasing = TRUE)
    cat("  worst statistic families:\n")
    for (f in names(sev)[seq_len(min(8L, length(sev)))])
      cat(sprintf("    %-22s %5.1f%% outside  (n = %d)\n", sub("^s_", "", f),
                  sev[[f]], as.integer(tab[[f]])))
  }
  if (o > 20) cat("  NOTE: do NOT respond by dropping the offending statistics --\n",
                  "        that makes the check unable to fail.\n", sep = "")
  invisible(x)
}

#' Plot a posterior predictive check
#'
#' @description \code{type = "pit"} gives PIT histograms by statistic family;
#'   \code{type = "pc"} projects prior-predictive, posterior-predictive and the
#'   observation onto one basis fitted on the reference table.
#' @param x a \code{ppc} object.
#' @param type \code{"pit"} or \code{"pc"}.
#' @param families number of statistic families to show, or their names. Numeric
#'   picks the most severe.
#' @param col histogram fill.
#' @param main outer title.
#' @param ... ignored.
#' @author Marcelo Gehara
#' @export
plot.ppc <- function(x, type = c("pit", "pc"), families = 12L,
                     col = "#4292C6", main = NULL, ...) {
  type <- match.arg(type)
  if (type == "pc") {
    if (is.null(x$pc))
      stop("this check carries no projection. Re-run posterior.predictive.check() ",
           "with `reftable=`, or ask for the histograms with type = \"pit\".",
           call. = FALSE)
    return(.ppc.plot.pc(x, main = main))
  }
  fam <- sub("^(s_[a-z]+_[A-Za-z]+).*$", "\\1", x$stats)
  tab <- table(fam); tab <- tab[tab >= 3]; cand <- names(tab)
  if (!length(cand)) stop("no statistic family has enough values to plot.", call. = FALSE)
  if (is.character(families)) {
    fs <- intersect(families, cand)
    if (!length(fs)) stop("none of the named families are present. Available: ",
                          paste(cand, collapse = ", "), call. = FALSE)
  } else {
    sev <- vapply(cand, function(f) {
      v <- x$pit[fam == f]; mean(v < 0.025 | v > 0.975, na.rm = TRUE) }, numeric(1))
    fs <- cand[order(-sev, -as.integer(tab[cand]))][seq_len(min(families, length(cand)))]
  }
  nr <- ceiling(sqrt(length(fs))); nc <- ceiling(length(fs) / nr)
  op <- graphics::par(mfrow = c(nr, nc), mar = c(3.1, 3.1, 2.4, 0.8),
                      mgp = c(1.9, 0.6, 0),
                      oma = c(0, 0, if (is.null(main)) 0 else 2.2, 0))
  on.exit(graphics::par(op))
  for (f in fs) {
    v <- x$pit[fam == f]; v <- v[is.finite(v)]
    o <- 100 * mean(v < 0.025 | v > 0.975)
    h <- graphics::hist(v, breaks = seq(0, 1, length.out = 21), plot = FALSE)
    graphics::plot(h, col = col, border = "white", freq = FALSE,
                   xlab = "PIT", ylab = "density", main = "",
                   ylim = c(0, max(h$density, 1.6) * 1.08))
    graphics::rect(0, 0, 0.025, graphics::par("usr")[4], col = "#00000010", border = NA)
    graphics::rect(0.975, 0, 1, graphics::par("usr")[4], col = "#00000010", border = NA)
    graphics::abline(h = 1, lty = 2, col = "grey35")
    graphics::title(main = sprintf("%s  (n = %d)", sub("^s_", "", f), length(v)),
                    cex.main = 0.92)
    graphics::mtext(sprintf("%.1f%% outside", o), side = 3, line = -1.2, adj = 0.97,
                    cex = 0.62, col = if (o > 20) "#A50F15" else "grey30")
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, line = 0.6, cex = 1.15, font = 2)
  invisible(x)
}

.ppc.plot.pc <- function(x, main = NULL) {
  B <- x$pc$pca
  op <- graphics::par(c("mfrow", "mar", "oma")); on.exit(graphics::par(op))
  graphics::par(mfrow = c(1, 2), mar = c(4.4, 4.6, 3.4, 1),
                oma = c(0, 0, if (is.null(main)) 0 else 2.2, 0))
  q95 <- stats::quantile(B$d.prior, 0.95)
  yl <- range(c(B$d.prior, B$d.post, B$d.obs)); yl[1] <- max(yl[1], 1e-3)
  graphics::plot(B$x.prior, B$d.prior, pch = 16, cex = 0.30, col = "#9E9E9E44", log = "y",
       ylim = yl, xlim = range(c(B$x.prior, B$x.post, B$x.obs)), xlab = B$lab,
       ylab = sprintf("k-NN distance to the PRIOR cloud (full %d-comp space)", x$pc$k),
       main = "prior, posterior and observed on one scale")
  graphics::points(B$x.post, B$d.post, pch = 16, cex = 0.62, col = "#2166ACAA")
  graphics::points(B$x.obs, B$d.obs, pch = 24, cex = 2.3, bg = "#D7191C",
                   col = "white", lwd = 2)
  graphics::abline(h = q95, lty = 2, col = "grey45")
  graphics::legend("topleft", bty = "n", cex = 0.75, pch = c(16, 16, 24),
    pt.cex = c(0.9, 1.0, 1.6), col = c("#9E9E9E", "#2166AC", "white"),
    pt.bg = c(NA, NA, "#D7191C"),
    legend = c(sprintf("prior predictive (n=%d)", x$pc$n.prior),
               sprintf("posterior predictive (n=%d)", x$pc$n.post), "observed"))
  graphics::legend("bottomright", bty = "n", cex = 0.72, lty = 2, col = "grey45",
                   legend = "95% of prior LOO-NN")
  graphics::mtext("y measures distance to the SAME reference cloud for all three",
                  3, 0.2, cex = 0.58, col = "grey30")
  dp <- stats::density(log10(B$d.prior)); dq <- stats::density(log10(B$d.post))
  graphics::plot(dp, col = "#757575", lwd = 2.5,
       xlim = range(log10(c(B$d.prior, B$d.post, B$d.obs))),
       ylim = c(0, max(dp$y, dq$y)),
       xlab = expression(log[10] * " k-NN distance to the prior cloud"),
       ylab = "density", main = "the same three, as distributions")
  graphics::polygon(dp, col = "#9E9E9E44", border = NA)
  graphics::lines(dq, col = "#2166AC", lwd = 2.5)
  graphics::polygon(dq, col = "#2166AC33", border = NA)
  graphics::abline(v = log10(B$d.obs), col = "#D7191C", lwd = 3)
  graphics::legend("topright", bty = "n", cex = 0.75, lwd = c(2.5, 2.5, 3),
    col = c("#757575", "#2166AC", "#D7191C"),
    legend = c("prior predictive", "posterior predictive", "observed"))
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, line = 0.6, cex = 1.15, font = 2)
  invisible(x)
}
