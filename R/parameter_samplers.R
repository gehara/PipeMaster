# internal function of the ms.commander
sample.pars<-function(x){
  k<-sample(nrow(x),nrow(x))
  for(i in k){
    if(c(as.numeric(x[i,4])+as.numeric(x[i,5]))==0){
      next
    } else {
    if(x[i,6]=="rtnorm"){
      samp<-msm::rtnorm(1, as.numeric(x[i,4]), as.numeric(x[i,5]), lower=0)
    } else {
      samp<-do.call(x[i,6],args=list(1,as.numeric(x[i,4]),as.numeric(x[i,5])),quote=F)
      while(samp<=0){
        samp<-do.call(x[i,6],args=list(1,as.numeric(x[i,4]),as.numeric(x[i,5])),quote=F)
      }
    }
    }
  x[i,4:5]<-samp
  }
  return(x)
}

# internal function of the ms.commander
# Sample parameters with conditions
#
# cond.list: a list of conditions, each is list(param1, op, param2) where:
#   op = "<" : param1 < param2 (enforced via rejection sampling)
#   op = "=" : copy param1's value to param2 (param1 is source)
#
# Processing order: (1) sample all, (2) enforce "<" via rejection, (3) apply "="
#
sample.w.cond <- function(par.matrix, cond.list = NULL) {

  x <- sample.pars(par.matrix)

  if (is.null(cond.list) || length(cond.list) == 0) return(x)

  # Map parameter names to row indices in par.matrix
  find_row <- function(name) {
    idx <- which(par.matrix[, 1] == name)
    if (length(idx) == 1) idx else NA
  }

  # --- Step 1: Collect "<" and "=" conditions ---
  lt_pairs <- list()
  eq_pairs <- list()
  for (cond in cond.list) {
    r1 <- find_row(cond$param1)
    r2 <- find_row(cond$param2)
    if (is.na(r1) || is.na(r2)) next
    if (cond$op == "<") {
      lt_pairs[[length(lt_pairs) + 1]] <- c(r1, r2)
    } else if (cond$op == "=") {
      eq_pairs[[length(eq_pairs) + 1]] <- c(r1, r2)
    }
  }

  # --- Step 2: Enforce "<" via rejection sampling ---
  if (length(lt_pairs) > 0) {
    max_tries <- 10000
    tries <- 0
    while (tries < max_tries) {
      all_ok <- TRUE
      for (p in lt_pairs) {
        if (as.numeric(x[p[1], 4]) >= as.numeric(x[p[2], 4])) {
          all_ok <- FALSE
          break
        }
      }
      if (all_ok) break
      x <- sample.pars(par.matrix)
      tries <- tries + 1
    }
    if (tries == max_tries)
      warning("sample.w.cond: max iterations reached for '<' conditions")
  }

  # --- Step 3: Apply "=" (copy param1 value to param2) ---
  for (p in eq_pairs) {
    x[p[2], 4:5] <- x[p[1], 4:5]
  }

  return(x)
}

# Convert legacy condition matrix to condition list
# Used by convert.model.conds() to migrate old model objects
.matrix.to.cond.list <- function(cond.matrix) {
  if (is.null(cond.matrix)) return(NULL)
  nam <- rownames(cond.matrix)
  conds <- list()
  for (i in 1:nrow(cond.matrix)) {
    for (j in 1:ncol(cond.matrix)) {
      if (i == j) next
      val <- cond.matrix[i, j]
      if (is.na(val) || val == "0") next
      if (val == "<") {
        conds[[length(conds) + 1]] <- list(param1 = nam[i], op = "<", param2 = nam[j])
      } else if (val == "=" && i < j) {
        # Only upper triangle for "=" to avoid bidirectional
        conds[[length(conds) + 1]] <- list(param1 = nam[i], op = "=", param2 = nam[j])
      }
    }
  }
  conds
}

# Convert old model objects (matrix-based conditions) to new list format
# @param model A Model object with legacy matrix conditions
# @return The same model with list-based conditions (matrices removed)
convert.model.conds <- function(model) {
  if (!is.null(model$conds$size.matrix) && is.null(model$conds$size)) {
    model$conds$size <- .matrix.to.cond.list(model$conds$size.matrix)
  }
  if (!is.null(model$conds$time.matrix) && is.null(model$conds$time)) {
    model$conds$time <- .matrix.to.cond.list(model$conds$time.matrix)
  }
  if (!is.null(model$conds$mig.matrix) && is.null(model$conds$mig)) {
    model$conds$mig <- .matrix.to.cond.list(model$conds$mig.matrix)
  }
  model$conds$size.matrix <- NULL
  model$conds$mig.matrix  <- NULL
  model$conds$time.matrix <- NULL
  model
}

# internal function to generate the locus file
get.locfile<-function(model){
  nloci <- nrow(model$loci)
  npop <- as.numeric(model$I[1,3])
  nrows <- nloci * npop
  locfile <- matrix(NA_character_, nrow = nrows, ncol = 6)
  colnames(locfile) <- c("id","n","pop","length","mu","rec")
  idx <- 1L
  for(i in 1:nloci){
    for(j in 1:npop){
      locfile[idx, ] <- c(model$I[i,1], model$I[i,j+3], j,
                          formatC(as.integer(model$loci[i,2]), format = "d"),
                          model$loci[i,4], 0)
      idx <- idx + 1L
    }
  }
  return(locfile)
}

# Check if parent process is still alive (cross-platform)
#
# Used by simulation workers to detect when the parent process has been
# killed or crashed, so they can exit cleanly instead of running as orphans.
#
# @param pid_file Path to the .PM_parent.pid file written by the parent.
# @return TRUE if parent is alive (or pid_file doesn't exist, i.e. single-core
#   mode), FALSE if parent is confirmed dead.
#
.pm.parent.alive <- function(pid_file) {
  if (!file.exists(pid_file)) return(TRUE)  # no pid file = single-core, always OK
  pid <- tryCatch(
    as.integer(readLines(pid_file, n = 1L, warn = FALSE)),
    error = function(e) NA_integer_
  )
  if (is.na(pid)) return(TRUE)  # can't read pid, assume alive
  # signal 0 checks existence without killing (cross-platform via tools)
  tryCatch({ tools::pskill(pid, 0); TRUE },
           error = function(e) FALSE)
}

# Write parent PID file and register on.exit cleanup that removes it
# and kills all worker processes.
#
# @param pid_file Path to write the parent PID.
# @param worker_pids_env An environment to collect worker PIDs into.
#   The on.exit handler reads worker_pids_env$pids to kill them.
# @param envir The environment to register on.exit in (parent function).
#
.pm.register.parent <- function(pid_file, worker_pids_env, envir = parent.frame()) {
  writeLines(as.character(Sys.getpid()), pid_file)

  do.call(on.exit, list(substitute({
    # Remove PID file so workers detect parent death
    tryCatch(file.remove(pid_file), error = function(e) NULL)
    # Kill any remaining workers
    for (wpid in worker_pids_env$pids) {
      tryCatch({
        if (.Platform$OS.type == "unix") {
          system(sprintf("kill -9 -%d 2>/dev/null; kill -9 %d 2>/dev/null", wpid, wpid),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
        } else {
          system(sprintf("taskkill /F /T /PID %d 2>NUL", wpid),
                 ignore.stdout = TRUE, ignore.stderr = TRUE)
        }
      }, error = function(e) NULL)
    }
  }), add = TRUE), envir = envir)
}

# internal function to generate the locus file
sample.mu.rates<-function(model){
  MEAN <- runif(1, as.numeric(model$loci[1,4]), as.numeric(model$loci[1,5]))
  SD <- runif(1, as.numeric(model$loci[1,4]), as.numeric(model$loci[1,5]))
  rates<-rtnorm(nrow(model$loci), MEAN, SD, 1e-12)
  rates<-rep(rates, each=as.numeric(model$I[1,3]))
  return(list(rates,c(MEAN,SD)))
}
