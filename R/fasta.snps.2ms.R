fasta.snp.2ms<-function(path.to.fasta,fasta.files,write.file=T,pop.assign){

  ms.out<-list()
  for(u in 1:length(fasta.files)){
    fas<-read.dna(file=fasta.files[1], format="fasta")
    fas<-as.character(fas)
    pops<-read.table(pop.assign, header=T)
    pops<-pops[with(pops, order(pops[,2])), ]
    fasta<-NULL
    p<-NULL
    for(j in 1:nrow(pops)){
      x<-match(pops[j,1],rownames(fas))
      if (is.na(x)==F) {
        p <- rbind(p, pops[j, ])
        fasta <- rbind(fasta, fas[x, ])
      }
    }
    fas<-fasta
    pops<-p
    npops<-length(unique(pops[,2]))
    pop.list<-list()
    for(i in 1:npops){
      pop.list[[i]]<-length(grep(i,pops[,2]))
    }

    string<-paste("-I",npops,paste(unlist(pop.list),collapse=" "))

    bin<-NULL
    pos<-NULL
    for(i in 1:ncol(fas)){
      a<-length(grep("a",fas[,i]))
      c<-length(grep("c",fas[,i]))
      g<-length(grep("g",fas[,i]))
      t<-length(grep("t",fas[,i]))

      r<-length(grep("r",fas[,i]))
      y<-length(grep("y",fas[,i]))
      m<-length(grep("m",fas[,i]))
      k<-length(grep("k",fas[,i]))
      s<-length(grep("s",fas[,i]))
      w<-length(grep("w",fas[,i]))

      h<-length(grep("h",fas[,i]))
      b<-length(grep("b",fas[,i]))
      v<-length(grep("v",fas[,i]))
      d<-length(grep("d",fas[,i]))
      n<-length(grep("n",fas[,i]))
      gap<-length(grep("-",fas[,i]))
      if(gap>0){next
      } else if (n>0){next}
      if(!nrow(fas) %in% c(a,c,g,t,r,y,m,k,s,w,h,b,v,d)){
        bin<-cbind(bin,fas[,i])
        pos<-c(pos,i)
        }
      }
    pos<-pos/ncol(fas)
    if(!(is.null(bin))){
    for(j in 1:ncol(bin)){
      g<-length(grep("g",bin[,j]))/nrow(bin)
      a<-length(grep("a",bin[,j]))/nrow(bin)
      t<-length(grep("t",bin[,j]))/nrow(bin)
      c<-length(grep("c",bin[,j]))/nrow(bin)

      r<-length(grep("r",bin[,j]))/nrow(bin)
      y<-length(grep("y",bin[,j]))/nrow(bin)
      m<-length(grep("m",bin[,j]))/nrow(bin)
      k<-length(grep("k",bin[,j]))/nrow(bin)
      s<-length(grep("s",bin[,j]))/nrow(bin)
      w<-length(grep("w",bin[,j]))/nrow(bin)

      h<-length(grep("h",bin[,j]))/nrow(bin)
      b<-length(grep("b",bin[,j]))/nrow(bin)
      v<-length(grep("v",bin[,j]))/nrow(bin)
      d<-length(grep("d",bin[,j]))/nrow(bin)
      bases<-c(a,c,g,t,r,y,m,s,k,w,h,b,v,d)
      names(bases)<-c("a","c","g","t","r","y","m","s","k","w","h","b","v","d")
      base<-which.max(bases[1:4])
      bin[,j]<-gsub(names(base),0,bin[,j])
      bases<-bases[-base]
      for(i in 1:length(bases)){
        bin[,j]<-gsub(names(bases[i]),1,bin[,j])
      }
      }
    seqs<-NULL
    for(i in 1:nrow(bin)){
      seqs<-c(seqs,paste(bin[i,],collapse=""))
    }
    ss<-ncol(bin)
    }else{seqs<-NULL
    ss<-0}

    if(write.file==T){
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("ms",nrow(fas),1,string))
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),"//",append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("segsites:",ss),append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),paste("positions:   ",paste(pos,collapse="    ")),append=T)
      write(file=paste(strsplit(fasta.files[u],".",fixed=T)[[1]][1],".ms",sep=""),seqs,sep="\n",append=T)
    }
      if(npops>1){
      ms.out[[u]]<-paste("ms",nrow(fas),1,string)
      } else { ms.out[[u]]<-paste("ms",nrow(fas),1)}
      ms.out[[u]]<-c(ms.out[[u]],paste("//"))
      ms.out[[u]]<-c(ms.out[[u]],paste("segsites:",ss))
      ms.out[[u]]<-c(ms.out[[u]],paste("positions:   ",paste(pos,collapse="    ")))
      ms.out[[u]]<-c(ms.out[[u]],paste(seqs,sep="\n"))


  }
  return(ms.out)
}
