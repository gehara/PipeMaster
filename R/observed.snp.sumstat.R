#' Observed summary statistics over thousends of loci
#' @description his function calculates the observed summary statistics from an empirical data. This summary statistics are the same as those simulated by the sim.snp.sumstat function. It is optimized for nexgen data.
#' @param model A model object bult by the main.menu function.
#' @param path.to.fasta Path to the folder containing all fastas to be included in the calculation.
#' @param moments Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. The defalt is FALSE.
#' @param pop.assign A two-column data frame with sample names in the first column and the corresponding population membership as numbers in the second column. The numbers should match the population number in the model object.
#' @param msABC.call String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.
#'                   There is no need to change that unless you want to compile the program yourself and point the function to it.
#' @return A list of vectors containing the observed summary stats.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
observed.snp.sumstat<-function(model,path.to.fasta,pop.assign,moments=F,msABC.call=get.msABC()){

  com<-ms.commander.snp(model,use.alpha=F,
                        msABC=msABC.call)

  x<-strsplit(system(com[[1]],intern=T),"\t")
  nam<-x[1][[1]]

  TD_denom<-paste(nam[grep("pi",nam)],nam[grep("theta_w",nam)],sep="_")
  nam<-nam[-grep("ZnS",nam)]
  nam<-nam[-grep("thomson",nam)]
  nam<-c(nam, TD_denom)
  mean_nam<-paste(nam,"_mean",sep="")
  if(moments==T){
    kur_nam<-paste(nam,"_kur",sep="")
    skew_nam<-paste(nam,"_skew",sep="")
    var_nam<-paste(nam,"_var",sep="")
  }

  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fa",fasta.files,fixed=T)]
  observed<-list()
  for(i in 1:length(fasta.files)){
    ms.output<-fasta.snp.2ms(path.to.fasta,fasta.files[i],write.file=T,pop.assign)
    locus.name<-strsplit(fasta.files[i],".",fixed=T)[[1]][1]
    xx<-strsplit(ms.output[[1]][1]," ")
    xx<-xx[[1]][2:length(xx[[1]])]
    xx<-paste(xx, collapse=" ")
    system(paste(msABC.call," ",xx," --obs ",locus.name,".ms > ",locus.name,".out",sep=""),wait=T)
    observed[[i]]<-read.table(file=paste(locus.name,".out",sep=""),header=T)
    snps<-grep("segs",names(observed[[i]]))
    cat(paste(i,"   ",sum(observed[[i]][snps[1:(length(snps)-1)]]),"SNPs"))
  }
  observed<-matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

  colnames(observed)<-x[1][[1]]

  TD_denom<-data.frame(observed[,grep("pi",colnames(observed))]-observed[,grep("theta_w",colnames(observed))])

  observed<-observed[,-grep("ZnS",colnames(observed))]
  observed<-observed[,-grep("thomson",colnames(observed))]
  observed<-cbind(observed, TD_denom)

  mean<-colMeans(observed,na.rm = T)
  names(mean)<-mean_nam
  if(moments==T){
    kur<-apply(observed,2,kurtosis, na.rm=T)
    names(kur)<-kur_nam
    skew<-apply(observed,2,skewness, na.rm=T)
    names(skew)<-skew_nam
    var<-apply(observed,2,var, na.rm=T)
    names(var)<-var_nam
  }

  observed<-list(mean,NULL,NULL,NULL)
  names(observed)<-c("mean","variance","kurtosis","skewness")
  if(moments==T){
    observed[[2]]<-var
    observed[[3]]<-kur
    observed[[4]]<-skew
  }
  return(observed)
}
