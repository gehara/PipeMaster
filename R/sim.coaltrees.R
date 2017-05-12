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

