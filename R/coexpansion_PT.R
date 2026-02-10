#'  Simulation of codemographic models
#' @description Simulation of the PT model of Gehara et al. (2017)
#' @param nsims Total number of simulations
#' @param var.zeta Variation on zeta parameter. Can be "FREE" to vary or be set to a specific value (between 0-1).
#' @param coexp.prior Uniform prior for the coexpansion time. Vector of two numbers with the lower and upper boundary of the prior.
#' @param Ne.prior Data frame with the prior values for the Ne of each population.
#' @param NeA.prior Data frame with the prior values for the ancestral Ne of each population.
#' @param time.prior Data frame with parameter values for the priors of the time of demographic change of each population.
#' @param gene.prior Data frame with parameter values for the priors of the mutation rate of each species.
#' @param alpha logical. If TRUE all demographic changes are exponential. If FALSE sudden changes. Default is FALSE.
#' @param append.sims logical. If TRUE simulations are appended to the simulations file. Default is FALSE.
#' @param path Path to the directory to write the simulations. Default is the working directory.
#' @details To simulate the model of Chan et al. (2014), the Threshold model and the Narrow Coexpansion Time model use the sim.coexp function.
#' @details See references for more details.
#' @references Gehara M., Garda A.A., Werneck F.P. et al. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756–4771.
#' @references Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31, 2501–2515.
#' @export
sim.coexpPT<-function(nsims,
                    var.zeta,
                    coexp.prior,
                    Ne.prior,
                    NeA.prior,
                    time.prior,
                    gene.prior,
                    alpha=F,
                    append.sims=F,
                    path=getwd())
{

  setwd(path)


  if(append.sims==F){
    simulations<-matrix(nrow=1,ncol=20)
    simulations[1,]<-c("zeta","ts","E(t)","DI",
               "var.pi","mean.pi","skew.pi","kur.pi",
               "var.ss","mean.ss","skew.ss","kur.ss",
               "var.H","mean.H","skew.H","kur.H",
               "var.TD","mean.TD","skew.TD","kur.TD")

    write.table(simulations,file="simulations.txt", quote=F,row.names=F, col.names=F, sep="\t")

    #populations.par<-matrix(nrow=1,ncol=nrow(NeA.prior)*4)
    #pop.names<-NULL
    #for(i in 1:nrow(NeA.prior)){
    #pop.names<-cbind(pop.names,t(c(paste("Ne",i,sep=""),paste("Exp.time",i,sep=""),paste("NeA",i,sep=""),paste("mi",i,sep=""))))
    #}
    #populations.par[1,]<-pop.names
    #write.table(populations.par,file="pop_parameters.txt", quote=F,row.names=F, col.names=F, sep="\t")
  }


  TIME<-system.time(for (i in 1:nsims){

    x<-coexp.sample.pars2(nruns=1,var.zeta=var.zeta,coexp.prior=coexp.prior,Ne.prior=Ne.prior,
                         NeA.prior=NeA.prior,time.prior=time.prior,gene.prior=gene.prior)

    y<-coexp.MS(MS.par=x$MS.par, gene.prior = gene.prior,alpha=alpha)

    z<-coexp.sumstat(ms.output=y,gene.prior=gene.prior)

    simulations<-c(x$coexp.par,z)

    #populations.par<-unlist(x$pop.par)

    write.table(t(simulations),file="simulations.txt", quote=F,row.names=F, col.names=F, append=T, sep="\t")
    #write.table(t(populations.par),file="pop_parameters.txt", quote=F,row.names=F, col.names=F, append=T,sep="\t")
    print(paste(i,"sims of",nsims,"| zeta = ",x$coexp.par[,1]))
   })
  print(TIME)
}

# Internal function of the sim.coexp2 function
# @description sample the parameters of the models from prior distributions.
coexp.sample.pars2<-function(nruns,
                             var.zeta,
                             coexp.prior,
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

    range<-coexp.prior[2]-coexp.prior[1]
    x<-range/nspecies
    priors<-(x*(1:nspecies))+coexp.prior[1]
    priors<-c(coexp.prior[1],priors)

    e.t<-vector()
    for(i in 1:(length(priors)-1)){
      e.t[i]<-runif(1,priors[i],priors[i+1])
    }

    Ts <-sample(e.t,1)
    x<-match(Ts,e.t)
    e.t<-e.t[-x]

    if (var.zeta=="FREE") {
      zeta.space<-1/nspecies # creates prior for n coexpanding species
      zeta.space<-zeta.space*(1:nspecies) # creates prior for n coexpanding species
      zeta<-sample(zeta.space,1)
      zeta.b<-nspecies*zeta
    } else {
      zeta<-var.zeta
      zeta.b<-nspecies*zeta
    }

    if(zeta.b==nspecies){
      time.prior[1:nrow(time.prior),3]<-Ts
    } else {
      coexp.sp<-sort(sample.int(nspecies,zeta.b))
      time.prior[c(coexp.sp),c(3,4)]<-Ts
      time.prior[-c(coexp.sp),c(3,4)]<-sample(e.t,(nspecies-length(coexp.sp)))
    }

    Et<-mean(time.prior[,3])
    Disp.index<-var(time.prior[,3])/Et

    coexp.par[j,]<-c(zeta,Ts,Et,Disp.index)

    for(i in 1:nspecies){
      ms.par<-NULL

      Ne <- runif(1, Ne.prior[i,3], Ne.prior[i,4])
      Ne.EXP.t <- time.prior[i,3]/time.prior[i,5] #corrects by generations
      theta.A.ratio <- runif(1, NeA.prior[i,3], NeA.prior[i,4])# thetaA (NeA) ratio
      NeA <- Ne*theta.A.ratio
      mi <- do.call(as.character(gene.prior[i,2]),args=list(1,gene.prior[i,3],gene.prior[i,4]),quote=F)
      while(mi<0){
        mi <- do.call(as.character(gene.prior[i,2]),args=list(1,gene.prior[i,3],gene.prior[i,4]),quote=F)
      }
      Ne <- Ne*gene.prior[i,7]
      theta=4*Ne*mi*gene.prior[i,5]
      scalar=4*Ne
      EXP.time=Ne.EXP.t/scalar

      g.rate=-log(NeA/Ne)/Ne.EXP.t

      ms.par<-cbind(theta,EXP.time,theta.A.ratio,g.rate)
      po.par<-c(Ne,time.prior[i,3],NeA,mi)
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
