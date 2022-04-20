# Simulate summary statistics using msABC
# @param model A model object bult by the main.menu function.
# @param use.alpha Logical.If TRUE the most recent population size change will be exponential. If FALSE sudden demographic changes. Default is FALSE.
# @param nsim.blocks Number of blocks to simulate. The total number of simulations is: nsim.blocks x sim.block.size.
# @param sim.block.size Simulations are performed in blocks. This argument defines the size of the block in number of simulations, i.e. how many simulations to run per block.
#                       A block of 1000 will work for most cases. Increse the total number of simulations with nsim.block argument.
# @param path Path to write the output. By default outputs will be saved in the working directory.
# @param output.name String. The prefix of the output names. Defalt is "model"
# @param append.sims Logical. If TRUE simulations will be appended in the last output. Default is FALSE.
# @param get.moments Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. Defalt is FALSE.
# @param msABC.call String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.
#                   There is no need to change that unless you want to compile the program yourself and point the function to it.
# @return Writes simulations and parameters to the path directory. msABC outputs a bunch of summary stats by defalt. They need to be selectd a posteriori.
# @references Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338.
# @references Pavlidis P., Laurent S., & Stephan W. (2010) msABC: A modification of Hudson’s ms to facilitate multi-locus ABC analysis. Molecular Ecology Resources, 10, 723–727.
# @author Marcelo Gehara
# @note This function does not work on Windows systems.
#
#
sim.msABC.sanger <- function(model,use.alpha=F, nsim.blocks=5, path=getwd(), append.sims=F, block.size=1000,
                    msABC.call = get.msABC(), output.name="model", ncores=1) {

  # set working directory
  setwd(path)

  if(append.sims==F){
    options(warn=-1)
  com <- ms.commander2(model, use.alpha = use.alpha)
  nam <- strsplit(system(paste(msABC.call,sum(as.numeric(model$I[1,4:5])),1,com[[1]]), intern=T)[1],"\t")
  write.table(t(nam[[1]]),file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F, col.names = F, append=F,sep="\t")
  write.table(t(com[[nrow(model$loci)+1]][1,]),file=paste(output.name,"_param.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  options(warn=0)
   }

  dput(model, ".model")
  dput(list(msABC.call, block.size, use.alpha), ".objects")
  file.copy(paste(system.file(package="PipeMaster"),"/run.msABC.R", sep=""), paste(path,"/.msABC", sep=""), overwrite = T)

  total.sims<-0
  for(k in 1:nsim.blocks){

  start.time <- Sys.time()

  write(0,".log")
  for(j in 1:ncores){
    system(paste("Rscript .msABC",j), wait = F)
  }

  l<-"0"
  while(sum(as.numeric(unlist(strsplit(l, "")))) < ncores){
    Sys.sleep(3)
    l <- readLines(".log")
  }

  file.remove(".log")

  simulations <- NULL
  cat("Reading simulations from worker nodes", sep="\n")
  for(h in 1:ncores){
    res <- read.table(file = paste(".",h,"_stats",sep=""))
    simulations <- rbind(simulations, res)
  }

  parameters <- NULL
  for(h in 1:ncores){
    res <- read.table(file = paste(".",  h,"_param",sep=""))
    parameters <- rbind(parameters, res)
  }

  cat("Writing simulations to file", sep="\n")


  cat("Removing old simulations", sep="\n")

  for(t in 1:ncores){
    file.remove(paste(".",t,"_stats",sep=""))
    file.remove(paste(".",t,"_param",sep=""))
  }

  write.table(ss,file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  write.table(param,file=paste(output.name,"_param.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")

  end.time <- Sys.time()
  total.sims <- total.sims+(block.size*ncores)
  cycle.time <- (as.numeric(end.time)-as.numeric(start.time))/60/60
  total.time <- cycle.time*nsim.blocks
  passed.time <- cycle.time*k
  remaining.time <- round(total.time-passed.time,3)
  cat(paste("PipeMaster:: ",total.sims," (~",round((block.size*ncores)/cycle.time)," sims/h) | ~",remaining.time," hours remaining",sep=""),"\n")


  }
  print("Done!")
}


