demog.sample.pars<-function(nruns,
                            var.zeta,
                            coexp.prior,
                            buffer,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior){

  MS.par<-list(NULL)
  pop.par<-list(NULL)
  coexp.par<-matrix(nrow=1,ncol=4)
  nspecies<-1
  MS.par[[1]]<-matrix(nrow=1,ncol=4)
  pop.par[[1]]<-matrix(nrow=1,ncol=4)

      ms.par<-NULL

      Ne <- runif(1, Ne.prior[1,3], Ne.prior[1,4])
      e.t<-runif(1, time.prior[1,3], time.prior[1,4])
      Ne.EXP.t <- e.t/time.prior[1,5] #corrects by generations
      theta.A.ratio <- runif(1, NeA.prior[1], NeA.prior[2])# thetaA (NeA) ratio
      NeA <- Ne*theta.A.ratio
      mi <- do.call(as.character(gene.prior[1,2]),args=list(1,gene.prior[1,3],gene.prior[1,4]),quote=F)
      while(mi<0){
        mi <- do.call(as.character(gene.prior[1,2]),args=list(1,gene.prior[1,3],gene.prior[1,4]),quote=F)
      }
      po.par<-c(Ne,e.t,NeA,mi)

      Ne <- Ne*gene.prior[1,7]
      theta=4*Ne*mi*gene.prior[1,5]
      scalar=4*Ne
      EXP.time=Ne.EXP.t/scalar

      g.rate=-log(NeA/Ne)/Ne.EXP.t

      ms.par<-cbind(theta,EXP.time,theta.A.ratio,g.rate)

      MS.par[[1]][1,]<-ms.par
      pop.par[[1]][1,]<-po.par

    pars<-list(NULL,NULL,NULL)
  names(pars)<-c("MS.par","pop.par")
  pars$MS.par<-MS.par
  pars$pop.par<-po.par
  return(pars)
}
