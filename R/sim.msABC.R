
sim.msABC<-function(model,use.alpha=F,nsim.blocks,path=getwd(),append.sims=F,block.size=1000, msABC.call=get.msABC(),output.name="sims"){

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
    com<-msABC.commander(model,use.alpha=use.alpha)

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
