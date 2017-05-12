mcmc.exp.pars<-function(gene.prior,
                        new.pars,
                        gentime){
  
      MS.par<-list(NULL)
      
      Ne <- new.pars[1]
      NeA <- new.pars[2]/Ne
      mi <- new.pars[4]
      Ne <- Ne*as.numeric(gene.prior[7])
      
      theta=4*Ne*mi*as.numeric(gene.prior[5])
      scalar=4*Ne
      Ne.EXP.t<-new.pars[3]/gentime
      EXP.time=Ne.EXP.t/scalar
      
      g.rate=-log(NeA/Ne)/Ne.EXP.t
      
      ms.par<-cbind(theta,EXP.time,NeA,g.rate)
      
      MS.par[[1]]<-ms.par
      
    
  pars<-list(NULL)
  names(pars)<-"MS.par"
  pars$MS.par<-MS.par
  return(pars)
}
