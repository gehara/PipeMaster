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
sim.msABC <- function(model,use.alpha=F,nsim.blocks=1,path=getwd(),append.sims=F,block.size=1000,
                    msABC.call=get.msABC(),output.name="model") {

  # set working directory
  setwd(path)

  if(append.sims==F){
    options(warn=-1)
  com<-msABC.commander(model,use.alpha=use.alpha)
  nam<-strsplit(system(paste(msABC.call,com[[1]]), intern=T)[1],"\t")
  write.table(t(nam[[1]]),file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F, col.names = F, append=F,sep="\t")
  write.table(t(com[[nrow(model$loci)+1]][1,]),file=paste(output.name,"_param.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  options(warn=0)
   }

  thou<-0
  for(j in 1:nsim.blocks){
  ss<-NULL
  param<-NULL


  for(i in 1:block.size){

    com <- PipeMaster:::ms.commander2(model,use.alpha=use.alpha)

    SS<-list()
    options(warn=-1)
    for(u in 1:nrow(model$loci)){
    SS[[u]]<-as.numeric(strsplit(system(paste(msABC.call,com[[u]]), intern=T)[2],"\t")[[1]])
    }
    options(warn=0)

    #ss<-Reduce("+",SS)/nrow(model$loci)
    SS<-do.call("rbind", SS)

    #while(sum(as.numeric(is.na(colMeans(SS,na.rm = T))))>=1){
    # com<-msABC.commander(model,use.alpha=use.alpha)
    #  SS<-list()
    #  options(warn=-1)
    #  for(u in 1:nrow(model$loci)){
    #    SS[[u]]<-as.numeric(strsplit(system(paste(msABC.call,com[[u]]), intern=T)[2],"\t")[[1]])
    #  }
    #  options(warn=0)
    #  SS<-do.call("rbind", SS)
    #}
    ss<-rbind(ss,colMeans(SS,na.rm = T))
    param<-rbind(param,com[[nrow(model$loci)+1]][2,])
    print(thou+i)
    }

  write.table(ss,file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  write.table(param,file=paste(output.name,"_param.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")

  thou<-thou+block.size
  }
  print("Done!")
}
