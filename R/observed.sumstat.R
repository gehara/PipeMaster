observed.sumstat<-function(model,path.to.fasta,fasta.files=list.files(),overall.SS=T,perpop.SS=T){
  
  # get population structure
  pops<-get.pops(model)
  # get sumstats names
  NAMES<-get.ss.pop.name(pops[[1]])
  # overall SS names
  overall.NAMES<-c("s.sites","pi","Hap.div","Taj.D","Fu.Li.D","Fu.Li.F")
  
  sewd(path.to.fasta)
  if(perpop.SS==T){
    write.table(t(NAMES),file="ObservedSS.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    }
  if(overall.SS==T){
    write.table(t(overall.NAMES),file="ObservedoverallSS.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    }
  
  
  fasta2ms(path.to.fasta,fasta.files,write.file=T)
  
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
      ss[[u]]@Tajima.D[is.na(ss[[u]]@Tajima.D)]<-0
      ss[[u]]@Fu.Li.D[is.na(ss[[u]]@Fu.Li.D)]<-0
      ss[[u]]@Fu.Li.F[is.na(ss[[u]]@Fu.Li.F)]<-0
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
      ss[[u]]@Tajima.D[is.na(ss[[u]]@Tajima.D)]<-0
      ss[[u]]@Fu.Li.D[is.na(ss[[u]]@Fu.Li.D)]<-0
      ss[[u]]@Fu.Li.F[is.na(ss[[u]]@Fu.Li.F)]<-0
      Taj.D<-ss[[u]]@Tajima.D
      Fu.Li.D<-ss[[u]]@Fu.Li.D
      Fu.Li.F<-ss[[u]]@Fu.Li.F
      Hap.Fst<-ss[[u]]@haplotype.F_ST
      nuc.Fst<-ss[[u]]@nucleotide.F_ST
      ss[[u]]<-cbind(s.sites,pi.within,Hap.div,Taj.D,Fu.Li.D,Fu.Li.F,Hap.Fst,nuc.Fst)
    }
  }
  if(perpop.SS==T){
    ss<-Reduce("+",ss)/nrow(model$I)
    write.table(ss,file="ObservedSS.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
  }
  if(overall.SS==T){
    OA.ss<-Reduce("+",OA.ss)/nrow(model$I)
    write.table(OA.ss,file="ObservedoverallSS.txt",quote=F,row.names=F, col.names = F, append=T,sep="\t")
  }
}