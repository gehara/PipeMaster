
observed.snp.sumstat<-function(model,path.to.fasta,pop.assign,moments=F,msABC.call){


    com<-ms.commander.snp(model,use.alpha=F,
                          msABC=msABC.call)
    x<-strsplit(system(com[[1]],intern=T),"\t")
    nam<-x[1]
    mean_nam<-paste(nam[[1]],"_mean",sep="")
    if(moments==T){
      kur_nam<-paste(nam[[1]],"_kur",sep="")
      skew_nam<-paste(nam[[1]],"_skew",sep="")
      var_nam<-paste(nam[[1]],"_var",sep="")
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
  print(paste(i,"   ",observed[[i]][1],"SNPs"))
  }
  observed<-matrix(unlist(observed), ncol = length(observed[[1]]), byrow = TRUE)

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
