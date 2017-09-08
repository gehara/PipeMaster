
sim.sumstat<-function(model,use.alpha=F,nsim.blocks,perpop.SS=T,overall.SS=T,path=getwd(),output.name,append.sims=F,sim.block.size=1000){

  # get population structure
  pops<-get.pops(model)
  # get sumstats names
  NAMES<-get.ss.pop.name(pops[[1]])
  # overall SS names
  overall.NAMES<-c("s.sites","pi","Hap.div","Taj.D","Fu.Li.D","Fu.Li.F")
  # get temp dir to save ms output
  swap<-path
  # set working directory
  setwd(path)
  # write output headings
  if(append.sims==F){
    if(perpop.SS==T){
    write.table(t(NAMES),file = paste(output.name, "_popstats_mean.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
      if(get.moments==T){
        write.table(t(NAMES),file = paste(output.name, "_popstats_var.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(NAMES),file = paste(output.name, "_popstats_kur.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(NAMES),file = paste(output.name, "_popstats_skew.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
      }
       }
    if(overall.SS==T){
    write.table(t(overall.NAMES),file = paste(output.name, "_overallstats_mean.txt",
                                                   sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
      if(get.moments==T){
        write.table(t(overall.NAMES),file = paste(output.name, "_overallstats_var.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(overall.NAMES),file = paste(output.name, "_overallstats_kur.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(overall.NAMES),file = paste(output.name, "_overallstats_skew.txt",sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
      }
    }
  write.table(t(ms.commander2(model,use.alpha = use.alpha)[[nrow(model$loci)+1]][1,]),file = paste(output.name, "_pars.txt",
                                                                                                        sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=F)
  }
  # beggin simulations
  thou<-0
  for(j in 1:nsim.blocks){
    # generate ms output vector
    sims<-list(NULL)
    for(i in 1:nrow(model$I)){
      sims[[i]]<-paste("ms",sum(as.numeric(model$I[i,4:ncol(model$I)])),sim.block.size)
    }
    sims[[nrow(model$loci)+1]]<-vector()
    # ms simulation
    for(k in 1:sim.block.size){
      com<-ms.commander2(model,use.alpha=use.alpha)
      for(u in 1:nrow(model$I)){
        S<-ms(nreps=1,nsam=sum(as.numeric(model$I[u,4:ncol(model$I)])),opts=com[[u]])
        while(as.numeric(strsplit(S[3]," ")[[1]][2])==0) S<-ms(nreps=1,nsam=sum(as.numeric(model$I[u,4:ncol(model$I)])),opts=com[[u]])
        sims[[u]]<-c(sims[[u]],S[2:length(S)])
      }
      sims[[nrow(model$loci)+1]]<-rbind(sims[[nrow(model$loci)+1]],com[[nrow(model$loci)+1]][2,])
      print(k+thou)
    }

    # write ms output to file
    for(u in 1:nrow(model$I)){
    write(sims[[u]],file=paste(swap,"/locus",u,sep=""), append=F)
    }

    # Start Summary stats calulation
    print("PopGenome!")

    ss<-list(NULL)
    OA.ss<-list(NULL)
    for(u in 1:nrow(model$I)){
      print(paste("reading simulations - locus",u))
      ss[[u]]<-readMS(paste(swap,"/locus",u,sep=""), big.data = F)
      if(overall.SS==T){
        print(paste("calculating overall sumary statistics - locus",u))
        ss[[u]]<-neutrality.stats(ss[[u]],FAST=T)
        ss[[u]]<-diversity.stats(ss[[u]])
        s.sites<-ss[[u]]@n.segregating.sites
        pi.within<-ss[[u]]@nuc.diversity.within/as.numeric(model$loci[u,2])
        Hap.div<-ss[[u]]@hap.diversity.within
        Taj.D<-ss[[u]]@Tajima.D
        Fu.Li.D<-ss[[u]]@Fu.Li.D
        Fu.Li.F<-ss[[u]]@Fu.Li.F
        OA.ss[[u]]<-cbind(s.sites,pi.within,Hap.div,Taj.D,Fu.Li.D,Fu.Li.F)
      }
      if(perpop.SS==T){
      ss[[u]]<-set.populations(ss[[u]],pops[[u]])
      print(paste("calculating per pop sumary statistics - locus",u))
      ss[[u]]<-neutrality.stats(ss[[u]],FAST=T)
      ss[[u]]<-F_ST.stats(ss[[u]],FAST=T)

      s.sites<-ss[[u]]@n.segregating.sites
      pi.within<-ss[[u]]@nuc.diversity.within/as.numeric(model$loci[u,2])
      Hap.div<-ss[[u]]@hap.diversity.within
      Taj.D<-ss[[u]]@Tajima.D
      Fu.Li.D<-ss[[u]]@Fu.Li.D
      Fu.Li.F<-ss[[u]]@Fu.Li.F
      Hap.Fst<-ss[[u]]@haplotype.F_ST
      nuc.Fst<-ss[[u]]@nucleotide.F_ST
      ss[[u]]<-cbind(s.sites,pi.within,Hap.div,Taj.D,Fu.Li.D,Fu.Li.F,Hap.Fst,nuc.Fst)
      }

      }
    # mean sumstats
    if(perpop.SS==T){
      SS<-NULL
      kur<-NULL
      vari<-NULL
      skew<-NULL
      for(j in 1:sim.block.size){
        x<-NULL
        for(i in 1:nrow(model$loci)){
          x<-rbind(x,ss[[i]][j,])
          }
        SS<-rbind(SS,colMeans(x, na.rm=T))
        if(get.moments==T){
          vari<-rbind(vari,diag(var(x, na.rm=T)))
          kk<-NULL
          sk<-NULL
          for(u in 1:ncol(x)){
            s<-skewness(x[,u],na.rm=T)
            sk<-c(sk,s)
            k<-kurtosis(x[,u],na.rm=T)
            kk<-c(kk,k)
          }
          kur<-rbind(kur,kk)
          skew<-rbind(skew,sk)
        }
    }
    # write outputs
    write.table(SS,file=paste(output.name, "_popstats_mean.txt",sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
    if(get.moments==T){
      write.table(vari,file=paste(output.name, "_popstats_var.txt",sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(kur,file=paste(output.name, "_popstats_kur.txt",sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(skew,file=paste(output.name, "_popstats_skew.txt",sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
      }
    }
    if(overall.SS==T){
      SS<-NULL
      kur<-NULL
      vari<-NULL
      skew<-NULL
      for(j in 1:sim.block.size){
        x<-NULL
        for(i in 1:nrow(model$loci)){
          x<-rbind(x,OA.ss[[i]][j,])
        }
        SS<-rbind(SS,colMeans(x, na.rm=T))
        if(get.moments==T){
          vari<-rbind(vari,diag(var(x, na.rm=T)))
          kk<-NULL
          sk<-NULL
          for(u in 1:ncol(x)){
            s<-skewness(x[,u],na.rm=T)
            sk<-c(sk,s)
            k<-kurtosis(x[,u],na.rm=T)
            kk<-c(kk,k)
          }
          kur<-rbind(kur,kk)
          skew<-rbind(skew,sk)
        }
      }


    write.table(SS,file=paste(output.name, "_overallstats_mean.txt",
                                 sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
    if(get.moments==T){
    write.table(vari,file=paste(output.name, "_overallstats_var.txt",
                              sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
    write.table(kur,file=paste(output.name, "_overallstats_kur.txt",
                              sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
    write.table(skew,file=paste(output.name, "_overallstats_skew.txt",
                              sep = ""),quote=F,row.names=F, col.names = F, append=T,sep="\t")
    }

    }
    write.table(sims[[nrow(model$loci)+1]],file=paste(output.name, "_pars.txt",
                                                      sep = ""),quote=F,row.names=F, col.names = F,sep="\t",append=T)
    # report job
    print(paste(j,"000 sims done!"))
    list = ls()
    x<-grep("thou",list)
    list<-list[-x]
    rm(list)
    thou<-thou+sim.block.size
  }
}


