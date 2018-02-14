#' internal function of the sim.sumstat
get.pops<-function(model){
  pops<-list(NULL)
  for(i in 1:nrow(model$I)){
    u<-1
    pop<-list()
    for(j in 4:ncol(model$I)){
      pop[[j-3]]<-c(u:(u-1+as.numeric(model$I[i,j])))
      u<-u+as.numeric(model$I[i,j])
    }
    pops[[i]]<-pop
  }
  return(pops)
}
