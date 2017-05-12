BN.MS<-function(MS.par,gene.prior){
  
  nspecies<-length(MS.par)
  nruns<-nrow(MS.par[[1]])
  sim<-list(NULL)
  
  for(ii in 1:nruns){
      nspecies<-nspecies
      for(xx in 1:nspecies){
        sims<-ms(gene.prior[xx,6],1,opts=paste("-t",MS.par[[xx]][ii,1],"-G",MS.par[[xx]][ii,4],"-eN",MS.par[[xx]][ii,2],1))
        sim[[xx]]<-sims 
      }
    }
 
 return(sim)
}