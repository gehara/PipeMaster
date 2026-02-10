#' Simulate summary statistics using msABC
#' @param model A model object built by the main.menu function.
#' @param use.alpha Logical. If TRUE the most recent population size change will be exponential. If FALSE sudden demographic changes. Default is FALSE.
#' @param nsim.blocks Number of blocks to simulate. The total number of simulations is: nsim.blocks x block.size (x ncores when ncores > 1).
#' @param block.size Simulations are performed in blocks. This argument defines the size of the block in number of simulations, i.e. how many simulations to run per block.
#'                   A block of 1000 will work for most cases. Increase the total number of simulations with nsim.blocks argument.
#' @param path Path to write the output. By default outputs will be saved in the working directory.
#' @param output.name String. The prefix of the output names. Default is "model".
#' @param append.sims Logical. If TRUE simulations will be appended in the last output. Default is FALSE.
#' @param ncores Number of cores for parallel execution. When ncores > 1, separate R worker processes are spawned. Default is 1.
#' @return Writes simulations and parameters to the path directory. msABC outputs a bunch of summary stats by default. They need to be selected a posteriori.
#' @references Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338.
#' @references Pavlidis P., Laurent S., & Stephan W. (2010) msABC: A modification of Hudson’s ms to facilitate multi-locus ABC analysis. Molecular Ecology Resources, 10, 723–727.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
sim.msABC.sanger <- function(model, use.alpha=F, nsim.blocks=5, path=getwd(), append.sims=F, block.size=1000,
                    output.name="model", ncores=1) {

  # set working directory
  setwd(path)
  locfile <- PipeMaster:::get.locfile(model)

  ############### this section is just to get the names of the sumstats
  if(append.sims==F){
    com <- PipeMaster:::msABC.commander(model,use.alpha=use.alpha,arg=1)
    write.table(locfile,paste(".",1,"locfile.txt",sep=""),row.names = F,col.names = T,quote = F,sep=" ")
    x <- strsplit(run.msABC(com[[1]]),"\t")
    nam<-x[1][[1]]
    #TD_denom <- paste(nam[grep("pi",nam)],nam[grep("_w",nam)],sep="_")
    #nam<-nam[-grep("ZnS",nam)]
    #nam<-nam[-grep("thomson",nam)]
    #cols <- grep("fwh",nam)
    #cols <- grep("thomson",nam)
    #cols <- c(cols, grep("ZnS",nam))
    #cols <- c(cols,grep("_FayWuH",nam))
    #if(length(cols)!=0) nam <- nam[-cols]
    #nam <- c(nam, TD_denom)
    nam <- c(com[[2]][1,], model$loci[,1], nam)
    #t(paste(nam,"_skew",sep="")),t(paste(nam,"_var",sep="")))
    write.table(t(nam),file=paste("SIMS_",output.name,".txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  }
  #################################


  if(ncores > 1) {
    # === MULTI-CORE: spawn separate R processes ===
    abs.path <- normalizePath(getwd())
    save(model, nsim.blocks, block.size, use.alpha, output.name,
         file = file.path(abs.path, ".PM_worker_params.RData"))

    worker_script <- paste(
      'args <- commandArgs(TRUE)',
      'worker_id <- as.integer(args[1])',
      'suppressMessages(library(PipeMaster))',
      sprintf('base_path <- "%s"', abs.path),
      'load(file.path(base_path, ".PM_worker_params.RData"))',
      'worker_dir <- file.path(base_path, paste0(".worker_", worker_id))',
      'dir.create(worker_dir, showWarnings=FALSE)',
      'sim.msABC.sanger(model=model, nsim.blocks=nsim.blocks, block.size=block.size,',
      '                 path=worker_dir, use.alpha=use.alpha,',
      '                 append.sims=TRUE,',
      '                 output.name=output.name, ncores=1)',
      'write("done", file.path(base_path, paste0(".worker_", worker_id, ".done")))',
      'quit(save="no")',
      sep="\n")
    writeLines(worker_script, file.path(abs.path, ".PM_worker.R"))

    start.time <- Sys.time()
    for(w in 1:ncores) {
      system(paste("Rscript", file.path(abs.path, ".PM_worker.R"), w,
                    ">", file.path(abs.path, paste0(".worker_", w, ".log")), "2>&1"),
             wait = FALSE)
    }
    cat(paste("PipeMaster:: Launched", ncores, "worker processes"), "\n")

    total_expected <- nsim.blocks * block.size * ncores
    while(TRUE) {
      Sys.sleep(5)
      done_count <- sum(file.exists(file.path(abs.path, paste0(".worker_", 1:ncores, ".done"))))

      total_sims <- 0
      for(w in 1:ncores) {
        wf <- file.path(abs.path, paste0(".worker_", w), paste0("SIMS_", output.name, ".txt"))
        if(file.exists(wf)) {
          n <- length(readLines(wf))
          if(n > 0) total_sims <- total_sims + n
        }
      }

      elapsed_h <- as.numeric(difftime(Sys.time(), start.time, units="hours"))
      if(total_sims > 0 && elapsed_h > 0) {
        rate <- total_sims / elapsed_h
        remaining <- round(max(0, (total_expected - total_sims) / rate), 3)
        cat(paste0("PipeMaster:: ", total_sims, " of ", total_expected,
                   " (~", round(rate), " sims/h) | ~", remaining,
                   " hours remaining | ", done_count, "/", ncores, " workers done"), "\n")
      }
      if(done_count >= ncores) break
    }

    cat("Compiling results from workers", sep="\n")
    outfile <- file.path(abs.path, paste0("SIMS_", output.name, ".txt"))
    for(w in 1:ncores) {
      wf <- file.path(abs.path, paste0(".worker_", w), paste0("SIMS_", output.name, ".txt"))
      if(file.exists(wf)) {
        worker_data <- readLines(wf)
        if(length(worker_data) > 0) {
          cat(paste(worker_data, collapse="\n"), "\n", file=outfile, append=TRUE, sep="")
        }
      }
      unlink(file.path(abs.path, paste0(".worker_", w)), recursive=TRUE)
      f <- file.path(abs.path, paste0(".worker_", w, ".done"))
      if(file.exists(f)) file.remove(f)
    }
    file.remove(file.path(abs.path, ".PM_worker_params.RData"))
    file.remove(file.path(abs.path, ".PM_worker.R"))
    for(w in 1:ncores) {
      f <- file.path(abs.path, paste0(".worker_", w, ".log"))
      if(file.exists(f)) file.remove(f)
    }

    end.time <- Sys.time()
    elapsed_h <- as.numeric(difftime(end.time, start.time, units="hours"))
    cat(paste0("PipeMaster:: Done! ", total_expected, " simulations in ",
               round(elapsed_h, 3), " hours (~", round(total_expected / elapsed_h), " sims/h)"), "\n")

  } else {
    # === SINGLE CORE ===
    sanger.sim.func <- function(worker.id) {
      ss <- NULL
      param <- NULL
      for(i in 1:block.size) {
        com <- PipeMaster:::ms.commander2(model, use.alpha = use.alpha)
        SS <- list()
        suppressWarnings({
          for(u in 1:nrow(model$loci)) {
            SS[[u]] <- as.numeric(strsplit(
              run.msABC(paste(sum(as.numeric(model$I[u,4:ncol(model$I)])), 1, com[[u]]))[2],
              "\t")[[1]])
          }
        })
        SS <- do.call("rbind", SS)
        SS.means <- colMeans(SS, na.rm = TRUE)
        SS.vars <- diag(var(SS, na.rm = TRUE))
        SS <- as.vector(rbind(SS.means, SS.vars))
        ss <- rbind(ss, SS)
        par <- com[[nrow(model$loci) + 1]][2, ]
        param <- rbind(param, par[-length(par)])
      }
      data.frame(cbind(param, ss))
    }

    total.sims <- 0
    for(k in 1:nsim.blocks) {
      start.time <- Sys.time()
      simulations <- sanger.sim.func(1)

      cat("Writing simulations to file", sep="\n")
      write.table(simulations, file=paste("SIMS_",output.name,".txt",sep=""),
                  quote=F, row.names=F, col.names=F, append=T, sep="\t")

      end.time <- Sys.time()
      total.sims <- total.sims + block.size
      cycle.time <- (as.numeric(end.time) - as.numeric(start.time)) / 60 / 60
      total.time <- cycle.time * nsim.blocks
      passed.time <- cycle.time * k
      remaining.time <- round(total.time - passed.time, 3)
      cat(paste("PipeMaster:: ", total.sims, " (~", round(block.size / cycle.time),
                " sims/h) | ~", remaining.time, " hours remaining", sep=""), "\n")
    }
    print("Done!")
  }
}


