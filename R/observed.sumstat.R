observed.sumstat<-function(model,path.to.fasta,fasta.files=list.files(),overall.SS=T,perpop.SS=T,get.moments=T){

  fasta.files<-fasta.files[grep(".fas",fasta.files)]

  fasta2ms(path.to.fasta,fasta.files,write.file=T)
  nrow(model$loci)# get population structure
  if(perpop.SS==T){
  pops<-get.pops(model)
  # get sumstats names
  NAMES<-get.ss.pop.name(pops[[1]])
  }
  # overall SS names
  overall.NAMES<-c("s.sites","pi","Hap.div","Taj.D","Fu.Li.D","Fu.Li.F")

  setwd(path.to.fasta)
  if(perpop.SS==T){
      write.table(t(NAMES),file ="observed_popstats_mean.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      if(get.moments==T){
        write.table(t(NAMES),file = "observed_popstats_var.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(NAMES),file = "observed_popstats_kur.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(NAMES),file = "observed_popstats_skew.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      }
    }
    if(overall.SS==T){
      write.table(t(overall.NAMES),file = "observed_overallstats_mean.txt",
                                                quote=F,row.names=F, col.names = F,sep="\t",append=F)
      if(get.moments==T){
        write.table(t(overall.NAMES),file = "observed_overallstats_var.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(overall.NAMES),file = "observed_overallstats_kur.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
        write.table(t(overall.NAMES),file = "observed_overallstats_skew.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      }
    }


  ss<-list(NULL)
  OA.ss<-list(NULL)
  for(u in 1:length(fasta.files)){
    ss[[u]]<-readMS(paste(fasta.files[u],".ms",sep=""), big.data = F)
    if(overall.SS==T){
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

  if(perpop.SS==T){
    SS<-NULL
    kur<-NULL
    vari<-NULL
    skew<-NULL
    x<-NULL
    for(jj in 1:nrow(model$loci)){
      x<-rbind(x,ss[[jj]][1,])
      }
      SS<-rbind(SS,colMeans(x, na.rm=T))
      if(get.moments==T){
        vari<-rbind(vari,diag(var(x, na.rm=T)))
        kk<-NULL
        sk<-NULL
        for(uu in 1:ncol(x)){
          s<-skewness(x[,uu],na.rm=T)
          sk<-c(sk,s)
          k<-kurtosis(x[,uu],na.rm=T)
          kk<-c(kk,k)
        }
        kur<-rbind(kur,kk)
        skew<-rbind(skew,sk)
      }
    
    # write outputs
    write.table(SS,file="observed_popstats_mean.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
    if(get.moments==T){
      write.table(vari,file="observed_popstats_var.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(kur,file="observed_popstats_kur.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(skew,file="observed_popstats_skew.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
    }
}
  if(overall.SS==T){
    SS<-NULL
    kur<-NULL
    vari<-NULL
    skew<-NULL
    x<-NULL
    for(jj in 1:nrow(model$loci)){
      x<-rbind(x,OA.ss[[jj]][1,])
      }
      SS<-rbind(SS,colMeans(x, na.rm=T))
      if(get.moments==T){
        vari<-rbind(vari,diag(var(x, na.rm=T)))
        kk<-NULL
        sk<-NULL
        for(uu in 1:ncol(x)){
          s<-skewness(x[,uu],na.rm=T)
          sk<-c(sk,s)
          k<-kurtosis(x[,uu],na.rm=T)
          kk<-c(kk,k)
        }
        kur<-rbind(kur,kk)
        skew<-rbind(skew,sk)
      }
    


    write.table(SS,file="observed_overallstats_mean.txt",
                              quote=F,row.names=F, col.names = F, append=T,sep="\t")
    if(get.moments==T){
      write.table(vari,file="observed_overallstats_var.txt",
                                  quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(kur,file="observed_overallstats_kur.txt",
                                 quote=F,row.names=F, col.names = F, append=T,sep="\t")
      write.table(skew,file="observed_overallstats_skew.txt",
                                 quote=F,row.names=F, col.names = F, append=T,sep="\t")
    }
  }
}
