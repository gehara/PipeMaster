#' internal function
#' @description control ms simulations
#' @return ms simulations
coexp.MS<-function(MS.par,gene.prior,alpha){

  nspecies<-length(MS.par)
  nruns<-nrow(MS.par[[1]])
  sim<-list(NULL)

  if(alpha==TRUE){
    for(ii in 1:nruns){
      nspecies<-nspecies
      for(xx in 1:nspecies){
        sims<-ms(gene.prior[xx,6],1,opts=paste("-t",MS.par[[xx]][ii,1],"-G",MS.par[[xx]][ii,4],"-eG",MS.par[[xx]][ii,2],0,"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3]))
          while(as.numeric(strsplit(sims[3]," ")[[1]][2])==0){
            sims<-ms(gene.prior[xx,6],1,opts=paste("-t",MS.par[[xx]][ii,1],"-G",MS.par[[xx]][ii,4],"-eG",MS.par[[xx]][ii,2],0,"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3]))
        }
         sim[[xx]]<-sims
       }
     }
  }

  if(alpha==FALSE){
    for(ii in 1:nruns){
      nspecies<-nspecies
      for(xx in 1:nspecies){
        sims<-ms(gene.prior[xx,6],1,opts=paste("-t",MS.par[[xx]][ii,1],"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3]))
        while(as.numeric(strsplit(sims[3]," ")[[1]][2])==0){
          sims<-ms(gene.prior[xx,6],1,opts=paste("-t",MS.par[[xx]][ii,1],"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3]))
        }
        sim[[xx]]<-sims
      }
    }
  }

  return(sim)
}
