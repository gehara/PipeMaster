plot.update<-function(sim,obs){
  mylabels <- colnames(sim)
  par(mfrow=c(3,3))
  
  for (i in 1:5){
    plot(sim[,i], xlab=mylabels[i], main="",type="l")
  
  }
  
  for (i in c(6:ncol(sim))){
       hist(sim[,i],breaks=20, xlab=mylabels[i], main="")
       abline(v = obs[i-5], col = 2)
  }
  
}

