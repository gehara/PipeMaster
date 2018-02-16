#' @export
sim.coexp<-function(nsims,
                    var.zeta,
                    coexp.prior,
                    buffer,
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

    x<-coexp.sample.pars(nruns=1,var.zeta=var.zeta,coexp.prior=coexp.prior,buffer=buffer,Ne.prior=Ne.prior,
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

# Internal function of the sim.coexp function
# @description sample the parameters of the models from prior distributions.
# @export
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

# internal function
# @description control ms simulations
# @return ms simulations
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

# internal function of the coexpansion model
# @description calculates summary statistics on ms-like files.
coexp.sumstat<-function(ms.output,gene.prior){
  sum.stat<-NULL

  for (j in 1:length(ms.output)){

    ss<-as.numeric(strsplit(ms.output[[j]][3]," ")[[1]][2])
    x<-ms.output[[j]][5:length(ms.output[[j]])]
    x<-gsub("0","A",x)
    x<-gsub("1","C",x)
    se<-list(NULL,NULL,NULL,NULL)
    names(se)<-c("nb","seq","nam","com")
    se$nb<-length(x) # number of samples
    se$seq<-x # sequencies
    se$nam<-c(1:length(x)) # sequence names, just numbers
    se$com<-NA
    class(se)<-"alignment" # this is an alignment
    x<-as.DNAbin(se)

    pi<-nuc.div(x)*ss/gene.prior[j,5]
    H<-H.div(x)
    TD<-tajima.test(x)$D

    #sspec<-site.spectrum(x,folded=T)
    #sspec<-sspec/sum(na.omit(sspec))
    #sspec<-sspec[1:3]

    SS<-c(pi,ss,H[2],TD)
    sum.stat<-rbind(sum.stat,SS)
  }
  #SDs<-apply(sum.stat,2,sd)
  vari<-diag(var(sum.stat))
  means<-colMeans(sum.stat)
  skew<-NULL
  kur<-NULL
  for(u in 1:ncol(sum.stat)){
    s<-skewness(sum.stat[,u])
    skew<-c(skew,s)
    k<-kurtosis(sum.stat[,u])
    kur<-c(kur,k)
  }

  h.s<-c(vari[1],means[1],skew[1],kur[1],
         vari[2],means[2],skew[2],kur[2],
         vari[3],means[3],skew[3],kur[3],
         vari[4],means[4],skew[4],kur[4])

  return(h.s)
}

