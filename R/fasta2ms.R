fasta2ms<-function(path.to.fasta,fasta.files,write.file=T){

  setwd(path.to.fasta)
  ms.out<-list()
  for(u in 1:length(fasta.files)){
    fas<-read.dna(file=fasta.files[u], format="fasta")
    fas<-as.character(fas)
    bin<-NULL
    pos<-NULL

    if(ncol(fas)==0){
      stop(paste(fasta.files[u],"has 0 base pairs! Check the alignment"))
    }

    for(i in 1:ncol(fas)){
      a<-length(grep("a",fas[,i]))
      c<-length(grep("c",fas[,i]))
      g<-length(grep("g",fas[,i]))
      t<-length(grep("t",fas[,i]))
      n<-length(grep("n",fas[,i]))
      gap<-length(grep("-",fas[,i]))
      if(nrow(fas) %in% c(a,c,g,t)){
      } else if (gap>0){
        } else if(n>0){
        } else {bin<-cbind(bin,fas[,i])
        pos<-c(pos,i)
        }
    }

    pos<-pos/ncol(fas)
    for(j in 1:ncol(bin)){
      a<-length(grep("a",bin[,j]))/nrow(bin)
      c<-length(grep("c",bin[,j]))/nrow(bin)
      g<-length(grep("g",bin[,j]))/nrow(bin)
      t<-length(grep("t",bin[,j]))/nrow(bin)
      bases<-c(a,c,g,t)
      names(bases)<-c("a","c","g","t")
      base<-which.max(bases)
      bin[,j]<-gsub(names(base),0,bin[,j])
      bases<-bases[-base]
      for(i in 1:3){
        bin[,j]<-gsub(names(bases[i]),1,bin[,j])
      }
      }
    seqs<-NULL
    for(i in 1:nrow(bin)){
      seqs<-c(seqs,paste(bin[i,],collapse=""))
    }

    if(write.file==T){
      write(file=paste(fasta.files[u],".ms",sep=""),paste("ms",nrow(fas),1))
      write(file=paste(fasta.files[u],".ms",sep=""),"//",append=T)
      write(file=paste(fasta.files[u],".ms",sep=""),paste("segsites:",ncol(bin)),append=T)
      write(file=paste(fasta.files[u],".ms",sep=""),paste("positions:   ",paste(pos,collapse="    ")),append=T)
      write(file=paste(fasta.files[u],".ms",sep=""),seqs,sep="\n",append=T)
    }

    if(write.file==F){
      ms.out[[u]]<-paste("ms",nrow(fas),1)
      ms.out[[u]]<-c(ms.out[[u]],paste("//"))
      ms.out[[u]]<-c(ms.out[[u]],paste("segsites:",ncol(bin)))
      ms.out[[u]]<-c(ms.out[[u]],paste("positions:   ",paste(pos,collapse="    ")))
      ms.out[[u]]<-c(ms.out[[u]],paste(seqs,sep="\n"))
    }
  }
  if(write.file==F){
    return(ms.out)
  }
}
