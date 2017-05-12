coexp.sample.pars<-function(nruns,
                            var.zeta,
                            coexp.prior,
                            buffer,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior){
  
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
     if (var.zeta=="FREE") {
      zeta.space<-1/nspecies # creates prior for n coexpanding species
      zeta.space<-zeta.space*(1:nspecies) # creates prior for n coexpanding species
      zeta<-sample(zeta.space,1)
      zeta.b<-nspecies*zeta
    } else {
      zeta<-var.zeta
      zeta.b<-nspecies*zeta
    }
     
    if(zeta.b==1){
      for(u in 1:nrow(time.prior.B)){
      time.prior.B[u,3]<-runif(1,time.prior.B[u,3],time.prior.B[u,4])
      } 
      Ts <-runif(1,coexp.prior[1],coexp.prior[2])
    } else {
      
    Ts<-runif(1,coexp.prior[1],coexp.prior[2])
    coexp.sp<-sort(sample.int(nspecies,zeta.b))
    time.prior.B[c(coexp.sp),c(3,4)]<-Ts
   
    for(u in 1:nrow(time.prior.B)){
      if(!(u %in% coexp.sp)){
    time<-runif(1,time.prior.B[u,3],time.prior.B[u,4])
      while(abs(time-Ts)<buffer){
        time<-runif(1,time.prior.B[u,3],time.prior.B[u,4])
      }
    time.prior.B[u,3:4]<-time
      }
    }
    }
    Et<-mean(time.prior.B[,3])
    Disp.index<-var(time.prior.B[,3])/Et
      
    coexp.par[j,]<-c(zeta,Ts,Et,Disp.index)
    
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
