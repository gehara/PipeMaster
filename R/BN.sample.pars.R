BN.sample.pars<-function(nruns,
                          Ne.prior,
                          NeA.prior,
                          time.prior,
                          gene.prior){
  
  MS.par<-list(NULL)
  samppars<-matrix(nrow=nruns,ncol=4)
  nspecies<-nrow(Ne.prior)
  
  for(i in 1:nspecies){
    mat<-matrix(nrow=nruns,ncol=4)
    MS.par[[i]]<-mat
  }
  
  for(j in 1:nruns){
    
      Ne <- runif(1, Ne.prior[1,3], Ne.prior[1,4])
      Ne.EXP.t <- runif(1,time.prior[1,3],time.prior[1,4])/time.prior[1,5] # corrects by generations
      theta.A.ratio <- runif(1, NeA.prior[1,3], NeA.prior[1,4]) # thetaA (NeA) ratio  
      NeA <- Ne*theta.A.ratio
      mi <- runif(1,gene.prior[1,3],gene.prior[1,4])
      samppars<-cbind(Ne,NeA,Ne.EXP.t,mi)
      Ne <- Ne*gene.prior[1,7]
      
      theta=4*Ne*mi*gene.prior[1,5]
      scalar=4*Ne
      EXP.time=Ne.EXP.t/scalar
      
      g.rate=-log(NeA/Ne)/Ne.EXP.t
      
      ms.par<-cbind(theta,EXP.time,theta.A.ratio,g.rate)
      
      MS.par[[1]]<-ms.par
      
    }
  pars<-list(NULL,NULL)
  names(pars)<-c("samp.par","MS.par")
  pars$samp.par<-samppars
  pars$MS.par<-MS.par
  return(pars)
}
  

