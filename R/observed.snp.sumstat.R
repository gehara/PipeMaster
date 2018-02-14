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

  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fa",fasta.files,fixed=T)]
  observed<-list()
  pop.assign<-read.table(pop.assign, header=T)
  for(i in 1:length(fasta.files)){
    ms.output<-fasta.snp.2ms(path.to.fasta,fasta.files[i],write.file=T,pop.assign)
    locus.name<-strsplit(fasta.files[i],".",fixed=T)[[1]][1]
    xx<-strsplit(ms.output[[1]][1]," ")
    xx<-xx[[1]][2:length(xx[[1]])]
    xx<-paste(xx, collapse=" ")
    system(paste(msABC.call," ",xx," --obs ",locus.name,".ms > ",locus.name,".out",sep=""),wait=T)
    observed[[i]]<-read.table(file=paste(locus.name,".out",sep=""),header=T)
    snps<-grep("segs",names(observed[[i]]))
    print(paste(i,"   ",sum(observed[[i]][snps[1:(length(snps)-1)]]),"SNPs"))
  }
  observed<-matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

  com<-ms.commander.snp(model,use.alpha=F,
                        msABC=msABC.call)

  nam<-strsplit(system(com[[1]],intern=T),"\t")[1][[1]]

  if(ncol(observed)!=length(nam)){
    stop(cat("your model pop structure does not match your assignment file pop structure.",
              paste("you have",length(unique(pop.assign[,2])),"populations in your assignment file"),
              paste("you have",model$I[,3],"populations in your model"),sep="\n"))
  }
  colnames(observed)<-nam

  TD_denom<-data.frame(observed[,grep("pi",colnames(observed))]-observed[,grep("theta_w",colnames(observed))])
  colnames(TD_denom)<-paste(nam[grep("pi",nam)],nam[grep("theta_w",nam)],sep="_")

  observed<-observed[,-grep("ZnS",colnames(observed))]
  observed<-observed[,-grep("thomson",colnames(observed))]
  observed<-cbind(observed, TD_denom)

  ## get names of sumstats
  nam<-colnames(observed)

  mean<-colMeans(observed,na.rm = T)
  names(mean)<-paste(nam,"_mean",sep="")
  if(moments==T){
    kur<-apply(observed,2,kurtosis, na.rm=T)
    names(kur)<-paste(nam,"_kur",sep="")
    skew<-apply(observed,2,skewness, na.rm=T)
    names(skew)<-paste(nam,"_skew",sep="")
    var<-apply(observed,2,var, na.rm=T)
    names(var)<-paste(nam,"_var",sep="")
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
