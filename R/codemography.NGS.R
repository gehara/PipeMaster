#' Simulation of codemographic models with ngs data
#' @description Simulation of codemographic models
#' @param nsims Total number of simulations
#' @param var.zeta Variation on zeta parameter. Can be "FREE" to vary or be set to a specific value (between 0-1).
#' @param coexp.prior Uniform prior for the coespansion time. Vector of two numbers with the lower and upper boudary of the prior.
#' @param th Threshold. Minimum time difference between Ts, time of simultaneous change and population specific times.
#' @param Ne.prior Data frame with the prior values for the Ne of each population.
#' @param NeA.prior Data frame with the prior values for the ancestral Ne of each population.
#' @param time.prior Data frame with parameter values for the priors of the time of demographic change of each population.
#' @param gene.prior Data frame with parameter values for the priors of the mutation rate of each species.
#' @param alpha logical. If TRUE all demographic chages are exponential. If FALSE sudden changes. Defaut is FALSE.
#' @param append.sims logical. If TRUE simulations are appended to the simulations file.  Defaut is FALSE.
#' @param mu.rates list. A list of parameters of mutation rate prior distribution. Fist parameter: the type of distribution. Secound element should be NA. Third and remainning elements are parameters of the distribution which should be givin in order. For instance if using the rnorm distribution the third parameter is the mean and the forth the SD.
#' @param path Path to the directiry to write the simulations. Defaut is the working directory.
#' @details To simulate the model of Chan et al. (2014) the th parameter should be set to zero and the time.prior should have the same value of the coexp.prior.
#' @details To simulate the Threshold model the th argument need to be higher than zero. To simulate the Narrow Coexpansion Time model the th argument need to be higher than zero and the boundaries of coexp.time shoud be narrower than the time.prior values.
#' See references for more details. Use the sim.coexp2 function to simulate the partitioned time model.
#' @references Gehara M., Garda A.A., Werneck F.P. et al. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756–4771.
#' @references Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31, 2501–2515.
#' @export
sim.coexp.ngs<-function(nsims,
                    var.zeta,
                    coexp.prior,
                    th,
                    Ne.prior,
                    NeA.prior,
                    time.prior,
                    gene.prior,
                    mu.rates,
                    alpha=F,
                    append.sims=F,
                    path=getwd())
{

  setwd(path)


  if(append.sims==F){
    simulations<-matrix(nrow=1,ncol=28)
    simulations[1,]<-c("zeta","ts","E(t)","DI",
                       "s_average_segs", "s_average_pi", "s_average_w"," s_average_tajd", "s_average_dvk", "s_average_dvh",
                       "s_variance_segs", "s_variance_pi", "s_variance_w", "s_variance_tajd", "s_variance_dvk", "s_variance_dvh",
                       "s_kurtosis_segs", "s_kurtosis_pi", "s_kurtosis_w", "s_kurtosis_tajd", "s_kurtosis_dvk", "s_kurtosis_dvh",
                       "s_skewness_segs", "s_skewness_pi", "s_skewness_w", "s_skewness_tajd", "s_skewness_dvk", "s_skewness_dvh")

    write.table(simulations,file="simulations.txt", quote=F,row.names=F, col.names=F, sep="\t")
    }


  TIME <- system.time(for (i in 1:nsims){

    x <- coexp.sample.pars.msABC(nruns=1,var.zeta = var.zeta, coexp.prior = coexp.prior, th = th, Ne.prior = Ne.prior,
                         NeA.prior = NeA.prior, time.prior = time.prior, gene.prior = gene.prior)

    y <- coexp.msABC(MS.par = x$MS.par, gene.prior = gene.prior, alpha = alpha, pop.par = x$pop.par, mu.rates = mu.rates)

    simulations <- c(x$coexp.par,y)

    write.table(t(simulations), file="simulations.txt", quote=F,row.names=F, col.names=F, append=T, sep="\t")
    print(paste(i,"sims of",nsims,"| zeta = ",x$coexp.par[,1]))
  })
  print(TIME)
}

# Internal function of the sim.coexp function
# @description sample the parameters of the models from prior distributions.
# @export
coexp.sample.pars.msABC<-function(nruns,
                            var.zeta,
                            coexp.prior,
                            th,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior) {

  MS.par<-list(NULL)
  pop.par<-list(NULL)
  coexp.par<-matrix(nrow=nruns,ncol=4)
  nspecies<-nrow(Ne.prior)

  for(i in 1:nspecies){
    mat<-matrix(nrow=nruns,ncol=4)
    MS.par[[i]]<-mat
    pop.par[[i]]<-mat
  }

  for(j in 1:nruns) {
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
          while(abs(time-Ts)<th){
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
      scalar = 4*Ne # generates the coalescent scalar
      EXP.time=Ne.EXP.t/scalar # expansion time / 4N

      g.rate=-log(NeA/Ne)/Ne.EXP.t

      ms.par<-cbind(0,EXP.time,theta.A.ratio,g.rate)
      po.par<-c(Ne,time.prior.B[i,3],NeA,0)
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

# internal function
# @description control ms simulations
# @return ms simulations
coexp.msABC<-function(MS.par, gene.prior, alpha, pop.par, mu.rates) {

  msABC <- PipeMaster:::get.msABC()
  nspecies <- length(MS.par)
  nruns<-nrow(MS.par[[1]])
  sim<-list(NULL)

  for(ii in 1:nruns){

    for(xx in 1:nspecies) {

      locfile <- gene.prior[[xx]]

      mu.rates[[2]] <- nrow(locfile)

      locfile[,5] <- do.call(mu.rates[[1]],args=mu.rates[2:length(mu.rates)])

      write.table(locfile,paste(".locfile.txt",sep=""),row.names = F,col.names = T,quote = F,sep=" ")

      sims <- read.table(text=system(paste(msABC,gene.prior[[xx]][1,2],1,"-eN",MS.par[[xx]][ii,2],MS.par[[xx]][ii,3],
                                           "--frag-begin --finp .locfile.txt --N",pop.par[[xx]][ii,1],"--frag-end"), intern = T), sep="\t", header=T)

      sims <- sims[,grep("average",colnames(sims))]
      sims <- sims[,-grep("thomson",colnames(sims))]
      sims <- sims[,-grep("ZnS",colnames(sims))]
      sims <- sims[,-grep("FayWuH",colnames(sims))]

      sim[[xx]] <- sims
    }

    sim <- matrix(unlist(sim), ncol = 6, byrow = TRUE)
    colnames(sim) <- colnames(sims)
  }
  average <- colMeans(sim)
  vari <- apply(sim, 2, var, na.rm = T)
  kur <- apply(sim, 2, kurtosis, na.rm = T)
  skew <- apply(sim, 2, skewness, na.rm = T)

  names(skew)<-gsub("average","skewness",names(skew))
  names(vari)<-gsub("average","variance",names(vari))
  names(kur)<-gsub("average","kurtosis",names(kur))

  return(c(average,vari,kur,skew))
}






