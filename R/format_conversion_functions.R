# internal function
# @description transforms fasta alignments in ms-like.
# @export
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

    if(is.null(pos)){
      if(write.file==T){
        write(file=paste(fasta.files[u],".ms",sep=""),paste("ms",nrow(fas),1))
        write(file=paste(fasta.files[u],".ms",sep=""),"//",append=T)
        write(file=paste(fasta.files[u],".ms",sep=""),paste("segsites:",0),append=T)
        write(file=paste(fasta.files[u],".ms",sep=""),paste("positions:"),append=T)
      }
      next
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

# internal function of the obs.snp.sumstat function
# @description transforms fasta alignments in ms-like.
# @export
fasta.snp.2ms<-function(path.to.fasta,fasta.files,write.file=T,pop.assign){

  ms.out<-list()
  for(u in 1:length(fasta.files)){
    fas<-read.dna(file=fasta.files[u], format="fasta")
    fas<-as.character(fas)
    if(is.matrix(fas)==F){
      stop(cat(paste("Something is wrong with alignment",fasta.files[u]),
               paste("Some potential problems:"),
               paste("1) sequences are not aligned"),
               paste("2) uknown character in the alignemt (only IUPAC nucleotide codes allowed, no question marks)"),
               paste("3) something else..."),sep="\n"))
    }

    pops<-pop.assign
    pops<-pops[with(pops, order(pops[,2])), ]


    if(length(grep(-9,match(rownames(fas),pops[, 1],nomatch=-9)))>0){
      rownames(fas)[grep(-9,match(rownames(fas),pops[, 1],nomatch=-9))]
      stop(cat(paste("There is at least one sequence in the alignment",fasta.files[u],"without assignment."),
               paste("The sequence name", rownames(fas)[grep(-9,match(rownames(fas),pops[, 1],nomatch=-9))]),
               paste("has no match in the assignment file."),sep="\n"))
    }

    fasta<-NULL
    p<-NULL
    for (j in 1:nrow(pops)) {
      x <- match(rownames(fas),pops[j, 1])
      x <- which(x==1)
      if (length(x)!=0) {
        p <- rbind(p, pops[rep(j,length(x)), ])
        fasta <- rbind(fasta, fas[x, ])
      }
    }
    fas<-fasta
    #ape:::write.dna(as.DNAbin(fas),fasta.files[u],format = "fasta", colw = 10000)
    pops<-p

    npops<-length(unique(pops[,2]))
    pop.list<-list()
    for(i in 1:npops){
      pop.list[[i]]<-length(which(pops[,2]==i))
    }

    string<-paste("-I",npops,paste(unlist(pop.list),collapse=" "))

    # keep only variable sites and get their position in the alignment
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

      x<-as.vector(bin)
      if(length(c(grep("0",x),grep("1",x)))!=length(x)){
        stop(paste("Something is wrong with alignment",fasta.files[u]))
      }

      ### if there is no variation
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

#' Converts ms output to DNAbin file
#' @description This function will take a ms-like output and convert it to DNAbin format.
#' @param ms.output A list of strings representing a ms output.
#' @param bp.length The number of base pairs used to calculate theta for the ms simulation.
#' @return a DNAbin object
#' @note This function is used internally by all the main functions used to simulate coexpantion models.
#' One can read in an ms output with readLines().
#' @author Marcelo Gehara
#' @examples # theta = 4Ne x mi x bp
#' Ne = 100000 # effective pop size
#' mi = 1e-8  # per base pair per generation mutation rate
#' bp = 1000  # number of base pairs
#' theta<-4*Ne*mi*bp
#' x<-ms(nsam=10, nrep=1, opts = paste("-t",theta))
#' y<-ms.to.DNAbin(x,bp=1000)
#' nuc.div(y)
#' @export
ms.to.DNAbin<-function(ms.output, bp.length){
  ss<-as.numeric(strsplit(ms.output[3]," ")[[1]][2]) # get seg sites

  if(ss>0){
    x<-ms.output[5:length(ms.output)]
    x<-gsub("0","A",x)
    x<-gsub("1","C",x)
  } else {
    x<-vector(mode="character",length=as.numeric(strsplit(ms.output[1]," ")[[1]][2]))
  }

  if(bp.length>0){
    for(i in 1:length(x)){
      x[i]<-paste(x[i],paste(rep("G",(bp.length-ss)),collapse=""),sep="")
    }
  }
  se<-list(NULL,NULL,NULL,NULL)
  names(se)<-c("nb","seq","nam","com")
  se$nb<-length(x) # number of samples
  se$seq<-x # sequencies
  se$nam<-c(1:length(x)) # sequence names, just numbers
  se$com<-NA
  class(se)<-"alignment" # this is alignment
  x<-as.DNAbin(se) # convert to DNAbin
  return(x)
}

