args = commandArgs(trailingOnly=TRUE)

print(paste("core",args,"running"))

library(msm)
model <- dget(".model")
obj <- dget(".objects")
block.size <- obj[[2]]
msABC.call <- obj[[1]]
use.alpha <- obj[[3]]

  ss <- NULL
  param <- NULL
  for(i in 1:block.size){

    com <- PipeMaster:::ms.commander2(model, use.alpha = use.alpha)

    SS<-list()
    options(warn=-1)
    for(u in 1:nrow(model$loci)){
      SS[[u]] <- as.numeric(strsplit(system(paste(msABC.call,sum(as.numeric(model$I[u,4:ncol(model$I)])),1,com[[u]]), intern=T)[2],"\t")[[1]])
    }
    options(warn=0)

    #ss<-Reduce("+",SS)/nrow(model$loci)
    SS <- do.call("rbind", SS)

    SS.means <- colMeans(SS, na.rm = T)
    SS.vars <- diag(var(SS, na.rm = T))

    SS <- as.vector(rbind(SS.means,SS.vars))

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
    ss <- rbind(ss, SS)
    par <- com[[nrow(model$loci)+1]][2,]
    param <- rbind(param, par[-length(par)])

  }

  write.table(data.frame(cbind(param,ss)), file=paste(".",args,"_stats",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")

  write(1, file=".log", append = T)

  print(paste("core", args, block.size, "sims done!"))
