
observed.snp.sumstat<-function(path.to.fasta,pop.assign,msABC.call){
  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fas",fasta.files,fixed=T)]
  observed<-NULL
  for(i in 1:length(fasta.files)){
  ms.output<-fasta.snp.2ms(path.to.fasta,fasta.files[i],write.file=T,pop.assign)
  locus.name<-strsplit(fasta.files[i],".",fixed=T)[[1]][1]
  xx<-strsplit(ms.output[[1]][1]," ")
  xx<-xx[[1]][2:length(xx[[1]])]
  xx<-paste(xx, collapse=" ")
  system(paste(msABC.call," ",xx," --obs ",locus.name,".ms > ",locus.name,".out",sep=""))
  ss<-read.table(file=paste(locus.name,".out",sep=""),header=T)
  observed<-rbind(observed,ss)
  print(i)
  }
  observed<-colMeans(observed,na.rm = T)
  return(observed)
  }  
    


