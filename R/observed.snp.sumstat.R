
observed.snp.sumstat<-function(model,path.to.fasta,pop.assign,moments=F,msABC.call){


  com<-ms.commander.snp(model,use.alpha=F,
                        msABC=msABC.call)
  x<-strsplit(system(com[[1]],intern=T),"\t")
  nam<-x[1][[1]]
  x<-strsplit(system(com[[1]],intern=T),"\t")
  nam<-x[1][[1]]
  TD_denom<-paste(nam[grep("pi",nam)],nam[grep("theta_w",nam)],sep=" - ")
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
    print(paste(i,"   ",sum(observed[[i]][snps[1:(length(snps)-1)]]),"SNPs"))
  }
  observed<-matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

  colnames(observed)<-x[1][[1]]

  TD_denom<-observed[,grep("pi",colnames(observed))]-observed[,grep("theta_w",colnames(observed))]
  colnames(TD_denom)<-paste(colnames(observed[,grep("pi",colnames(observed))]),
                            colnames(observed[,grep("theta_w",colnames(observed))]),sep=" - ")
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
