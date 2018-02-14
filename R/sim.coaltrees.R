#' Simulate coalescent trees using ms
#' @description This function simulates coalescent trees and writes them into a file.
#' @param model A model object bult by the main.menu function.
#' @param use.alpha Logical.If TRUE the most recent population size change will be exponential. If FALSE sudden demographic changes. Default is FALSE.
#' @param nsim.blocks Number of blocks to simulate. The total number of simulations is: nsim.blocks x sim.block.size.
#' @param sim.block.size Simulations are performed in blocks. This argument defines the size of the block in number of simulations, i.e. how many simulations to run per block.
#'                       A block of 1000 will work for most cases. Increse the total number of simulations with nsim.block argument.
#' @param path Path to write the output. By default outputs will be saved in the working directory.
#' @param append.sims Logical. If TRUE simulations will be appended in the last output. Default is FALSE.
#' @return Writes trees and parameters to the path directory.
#' @references Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338.
#' @author Marcelo Gehara
#' @export
sim.coaltrees<-function(model,nsim.blocks,use.alpha=F,path=getwd(),append.sims=F,sim.block.size=1000){

  setwd(path)
  thou<-0

  for(j in 1:nsim.blocks){

    sims<-list(NULL)
    for(i in 1:nrow(model$I)){
      sims[[i]]<-vector()
    }

   if(append.sims==F){
      write.table(t(ms.commander2(model,use.alpha = use.alpha)[[nrow(model$loci)+1]][1,]),file="SampPars.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    }

    pars<-NULL
    for(k in 1:sim.block.size){
      com<-ms.commander2(model,use.alpha=F)
      pars<-rbind(pars,com[[nrow(model$I)+1]][2,])
      for(u in 1:nrow(model$I)){
        t<-ms(nreps=1,nsam=sum(as.numeric(model$I[u,4:ncol(model$I)])),opts=paste("-T",com[[u]]))
        sims[[u]]<-c(sims[[u]],t[3])
        }
      print(k+thou)
    }

    if(append.sims==F){
    for(u in 1:nrow(model$I)){
      write(sims[[u]],file=paste("locus",u,sep=""), append=F)
    }
    } else {
     for(u in 1:nrow(model$I)){
       write(sims[[u]],file=paste("locus",u,sep=""), append=T)
     }
    }

    write.table(pars,file="SampPars.txt",quote=F,row.names=F, col.names = F,sep="\t",append=T)
    }
    thou<-thou+sim.block.size
  }

