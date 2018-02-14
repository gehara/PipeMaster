
sim.ms<-function(model,nsim.blocks,path=getwd(),append.sims=F,sim.block.size=1000){

  setwd(path)
  thou<-0
  for(j in 1:nsim.blocks){

    sims<-list(NULL)

    for(i in 1:nrow(model$I)){
      sims[[i]]<-character()
      write(sims[[i]],file=paste("locus",i,sep=""), append=append.sims)
    }



    for(k in 1:sim.block.size){
      for(u in 1:nrow(model$I)){
        S<-ms(nreps=1,nsam=sum(as.numeric(model$I[u,4:ncol(model$I)])),opts=ms.commander(model)[[u]])
        sims[[u]]<-c(sims[[u]],S)
      }
      print(k+thou)
    }


    for(u in 1:nrow(model$I)){
      write(sims[[u]],file=paste("locus",u,sep=""), append=T)
    }

    thou<-thou+sim.block.size
  }
}



