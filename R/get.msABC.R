# internal function to find the simulator
get.msABC<-function(){
  pack<-find.package("PipeMaster")

  if(Sys.info()[1]=="Linux"){
    msABC.call<-paste(pack,"/msABClinux",sep="")
  }
  if(Sys.info()[1]=="Darwin"){
    msABC.call<-paste(pack,"/msABCmac",sep="")
  }
  Sys.chmod(paths=msABC.call, mode = 7777, use_umask = TRUE)
  return(msABC.call)
  }

# internal function to generate the locus file
get.locfile<-function(model){
  locfile<-NULL
  for(i in 1:nrow(model$loci)){
    for(j in 1:model$I[1,3]){
      locfile<-rbind(locfile,c(model$I[i,1],model$I[i,j+3],j,model$loci[i,2],model$loci[i,4],0))
    }
    colnames(locfile)<-c("id","n","pop","length","mu","rec")
  }
  return(locfile)
}

# internal function to generate the locus file
sample.mu.rates<-function(model){
  MEAN <- runif(1, as.numeric(model$loci[1,4]), as.numeric(model$loci[1,5]))
  SD <- runif(1, as.numeric(model$loci[1,4]), as.numeric(model$loci[1,5]))
  rates<-rtnorm(nrow(model$loci),MEAN,SD,0)
  rates<-rep(rates, each=as.numeric(model$I[1,3]))
  return(list(rates,c(MEAN,SD)))
    }
