#' Observed summary statistics for nexgen data
#' @description This function calculates the observed summary statistics from an empirical data. This summary statistics are the same as those simulated by the sim.snp.sumstat function. It is optimized for nexgen data.
#' @param model A model object bult by the main.menu function. Any model with the same number of populations of your empirical data will do the job. This is just to build the sumstats names correctly.
#' @param path.to.fasta Path to the folder containing all fastas to be included in the calculation.
#' @param moments Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. The defalt is FALSE.
#' @param pop.assign A two-column data frame with sample names in the first column and the corresponding population membership preferably as numbers in the second column (If you have a single population the numbers wont matter). The numbers should match the population number in the model object.
#' @param msABC.call String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.
#'                   There is no need to change that unless you want to compile the program yourself and point the function to it.
#' @return A list of vectors containing the observed summary stats.
#' @author Marcelo Gehara
#' @note This function does not work on Windows systems.
#' @export
obs.sumstat.ngs<-function(model,path.to.fasta,pop.assign,moments=F,msABC.call=get.msABC()){

  if(is.na(get0("model",ifnotfound=NA))){
    stop("You need to include a model to calculate sumstats.")
  }

  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fa",fasta.files,fixed=T)]
  observed<-list()
  #pop.assign<-read.table(pop.assign, header=T)
  for(i in 1:length(fasta.files)){
    ms.output <- fasta.snp.2ms(path.to.fasta,fasta.files[i],write.file=T,pop.assign)
    locus.name <- strsplit(fasta.files[i],".",fixed=T)[[1]][1]
    xx<-strsplit(ms.output[[1]][1]," ")
    xx<-xx[[1]][2:length(xx[[1]])]
    xx<-paste(xx, collapse=" ")
    system(paste(msABC.call," ",xx," --obs ",locus.name,".ms > ",locus.name,".out",sep=""),wait=T)
    observed[[i]]<-read.table(file=paste(locus.name,".out",sep=""),header=T)
    snps<-grep("segs",names(observed[[i]]))
    print(paste(i,"   ",observed[[i]][snps[(length(snps))]],"SNPs"))
  }
  observed<-matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

  com<-ms.commander2(model,use.alpha=F)
  options(warn=-1)
  x<-strsplit(system2(msABC.call, args=paste(sum(as.numeric(model$I[1,4:ncol(model)])),1,com[[1]]), stdout = T,stderr=T,wait=T),"\t")
  options(warn=0)
  nam<-x[1][[1]]

  if(ncol(observed)!=length(nam)){
  stop(cat("your model pop structure does not match your assignment file pop structure.",
            paste("you have",length(unique(pop.assign[,2])),"populations in your assignment file"),
            paste("you have",model$I[,3],"populations in your model"),sep="\n"))
  }
  colnames(observed)<-nam

  TD_denom<-data.frame(observed[,grep("pi",colnames(observed))]-observed[,grep("_w",colnames(observed))])
  colnames(TD_denom)<-paste(nam[grep("pi",nam)],nam[grep("_w",nam)],sep="_")

  #observed<-observed[,-grep("ZnS",colnames(observed))]
  #observed<-observed[,-grep("thomson",colnames(observed))]
  observed <- cbind(observed, TD_denom)

  ## get names of sumstats
  nam <- colnames(observed)
  #observed[observed==0]<-NA
  mean <- colMeans(observed,na.rm = T)
  var <- apply(observed,2,var, na.rm=T)
  observed <- cbind(mean,var)

  observed<-as.vector(t(observed))

  com<-msABC.commander(model,use.alpha=F,arg=1)
  options(warn=-1)
  x<-strsplit(system2(msABC.call, args=com[[1]], stdout = T,stderr=T,wait=T),"\t")
  options(warn=0)
  nam<-x[1][[1]]

  names(observed) <- c(nam,paste(nam[grep("pi",nam)],nam[grep("_w",nam)],sep="_"))

  return(t(data.frame(observed)))
}
