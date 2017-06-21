


coexp.sample.pars.pulses<-function(nruns,
                            max.pulses,
                            buffer,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior){

  dmode <- function(x) {
    den <- density(x, kernel=c("gaussian"))
    ( den$x[den$y==max(den$y)] )
  }

  MS.par<-list(NULL)
  pop.par<-list(NULL)
  coexp.par<-matrix(nrow=nruns,ncol=4)
  nspecies<-nrow(Ne.prior)

  for(i in 1:nspecies){
    mat<-matrix(nrow=nruns,ncol=4)
    MS.par[[i]]<-mat
    pop.par[[i]]<-mat
  }
  for(j in 1:nruns){
     time.prior.B<-time.prior

      max.pulse.space<-coexp.prior[2]/buffer

      range<-coexp.prior[2]-coexp.prior[1]
      x<-range/max.pulse.space
      priors<-(x*(1:max.pulse.space))+coexp.prior[1]
      priors<-c(coexp.prior[1],priors)

      e.t<-vector()
      for(i in 1:(length(priors)-1)){
        e.t[i]<-runif(1,priors[i],priors[i+1])
      }

      pulses<-sample(1:max.pulses,1)

      pulse.times<-sort(sample(e.t,pulses))


      time.prior.B[,c(3,4)]<-replicate(nrow(time.prior.B),sample(pulse.times,1))

    }

    Et<-mean(time.prior.B[,3])
    Disp.index<-var(time.prior.B[,3])/Et
    Ts<-dmode(time.prior.B[,3])
    coexp.par[j,]<-c(pulses,Ts,Et,Disp.index)

    for(i in 1:nspecies){
      ms.par<-NULL

      Ne <- runif(1, Ne.prior[i,3], Ne.prior[i,4])
      Ne.EXP.t <- time.prior.B[i,3]/time.prior[i,5] #corrects by generations
      theta.A.ratio <- runif(1, NeA.prior[i,3], NeA.prior[i,4])# thetaA (NeA) ratio
      NeA <- Ne*theta.A.ratio
      mi <- do.call(as.character(gene.prior[i,2]),args=list(1,gene.prior[i,3],gene.prior[i,4]),quote=F)
      while(mi<0){
        mi <- do.call(as.character(gene.prior[i,2]),args=list(1,gene.prior[i,3],gene.prior[i,4]),quote=F)
      }
      Ne <- Ne*gene.prior[i,7] # Ne times inheritance scalar: cerrects the Ne
      theta=4*Ne*mi*gene.prior[i,5] # 4N x mutation per pb x bp
      scalar=4*Ne # generates the coalescent scalar
      EXP.time=Ne.EXP.t/scalar # expansion time / 4N

      g.rate=-log(NeA/Ne)/Ne.EXP.t

      ms.par<-cbind(theta,EXP.time,theta.A.ratio,g.rate)
      po.par<-c(Ne,time.prior.B[i,3],NeA,mi)
      MS.par[[i]][j,]<-ms.par
      pop.par[[i]][j,]<-po.par
    }
    #print(j)
  }
  #write.table(coexp.par,file="sim.par.txt", quote=F,row.names=F, col.names=F, append=T, sep="\t")
  #return(MS.par)
  pars<-list(NULL,NULL,NULL)
  names(pars)<-c("coexp.par","MS.par","pop.par")
  pars$coexp.par<-coexp.par
  pars$MS.par<-MS.par
  pars$pop.par<-pop.par
  return(pars)
}
