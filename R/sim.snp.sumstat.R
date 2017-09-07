
sim.snp.sumstat<-function(model,nsim.blocks,path=getwd(),use.alpha=F,
                          append.sims=F,block.size=100, msABC.call="./msABC",output.name){

  # set working directory
  setwd(path)

  if(append.sims==F){
  com<-ms.commander.snp(model,use.alpha=use.alpha,
                        msABC=msABC.call)
  x<-strsplit(system(com[[1]],intern=T),"\t")
  nam<-x[1]
  write.table(t(nam[[1]]),file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  write.table(t(com[[3]][1,]),file=paste(output.name,"_par.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  }

  thou<-0
  for(j in 1:nsim.blocks){
  SS<-NULL
  param<-NULL
  for(i in 1:block.size){
    com<-ms.commander.snp(model,msABC=msABC.call, use.alpha = use.alpha)

    system(paste(com[[1]],"> out.txt"))
    S1<-read.table("out.txt",header = T)
    if(model$loci[,1]=="0"){S2<-NULL
    }else{
    system(paste(com[[2]],"> out.txt"))
    S2<-read.table("out.txt",header = T)}
    sumstat<-rbind(S1,S2)
    sumstat<-colMeans(sumstat,na.rm = T)
    param<-rbind(param,com[[3]][2,])
    SS<-rbind(SS,sumstat)
    print(thou+i)

    }
  thou<-thou+block.size
  write.table(SS,file=paste(output.name,"_stats.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  write.table(param,file=paste(output.name,"_par.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  }
}
