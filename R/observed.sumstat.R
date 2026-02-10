#' Observed summary statistics
#' @description This function calculates the observed summary statistics from an empirical data.
#'              This summary statistics are the same as those simulated by the sim.sumstat function.
#'              It is optimized for sanger data.
#' @param model A model object bult by the main.menu function.
#' @param path.to.fasta Path to the folder containing all fastas to be included in the calculation.
#' @param list.files A list of fasta files to be included in the calculation. By default all fastas present in the folder are included.
#' @param moments Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. Defalt is FALSE.
#' @param pop.assign A two-column data frame with sample names in the first column and the corresponding population membership as numbers in the second column. The numbers should match the population number in the model object.
#' @param msABC.call String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.
#'                   There is no need to change that unless you want to compile the program yourself and point the function to it.
#' @return A list of vectors containing the observed summary stats.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
observed.sumstat<-function(model,path.to.fasta,fasta.files=list.files(),overall.SS=T,perpop.SS=T,get.moments=T){
  .Deprecated("obs.sumstat.ngs", package = "PipeMaster",
              msg = "observed.sumstat() is deprecated. Use obs.sumstat.ngs() instead, which uses the built-in msABC C code and does not require PopGenome.")
  if (!requireNamespace("PopGenome", quietly = TRUE)) {
    stop("observed.sumstat() requires the PopGenome package. Please use obs.sumstat.ngs() instead.")
  }

  wd.path <- getwd()

  setwd(path.to.fasta)

  fasta.files <- list.files(pattern = ".fas")

  fasta2ms(path.to.fasta, fasta.files, write.file=T)
  # get population structure
  if(perpop.SS==T){
    pops<-get.pops(model)
    # get sumstats names
    NAMES<-get.ss.pop.name(pops[[1]])
  }
  # overall SS names
  overall.NAMES<-c("s.sites","pi","Hap.div","Taj.D","Fu.Li.D","Fu.Li.F")

  setwd(path.to.fasta)
  if(perpop.SS==T){
    write.table(t(paste(NAMES,"_mean",sep="")),file ="observed_popstats_mean.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    if(get.moments==T){
      write.table(t(paste(NAMES,"_var",sep="")),file = "observed_popstats_var.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      write.table(t(paste(NAMES,"_kur",sep="")),file = "observed_popstats_kur.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      write.table(t(paste(NAMES,"_skew",sep="")),file = "observed_popstats_skew.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    }
  }
  if(overall.SS==T){
    write.table(t(paste(overall.NAMES,"_mean",sep="")),file = "observed_overallstats_mean.txt",
                quote=F,row.names=F, col.names = F,sep="\t",append=F)
    if(get.moments==T){
      write.table(t(paste(overall.NAMES,"_var",sep="")),file = "observed_overallstats_var.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      write.table(t(paste(overall.NAMES,"_kur",sep="")),file = "observed_overallstats_kur.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
      write.table(t(paste(overall.NAMES,"_skew",sep="")),file = "observed_overallstats_skew.txt",quote=F,row.names=F, col.names = F,sep="\t",append=F)
    }
  }


  ss<-list(NULL)
  OA.ss<-list(NULL)
  for(u in 1:length(fasta.files)){
    ss[[u]]<-readMS(paste(strsplit(fasta.files[u],".fas")[[1]],".ms",sep=""), big.data = F)

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
  setwd(wd.path)
}
