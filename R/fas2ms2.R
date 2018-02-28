fas2ms2<-function (fas) 
{
  
  fas <- as.character(fas)
  
  bin<-NULL
  pos<-NULL
  gaps<-0
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
    if(gap>0){
      gaps<-gaps+1
      next
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
  
  
  ms.out <- paste("ms", nrow(fas), 1)
  ms.out <- c(ms.out, paste("//"))
  ms.out <- c(ms.out, paste("segsites:", 
                            ncol(bin)))
  ms.out <- c(ms.out, paste("positions:   ", 
                            paste(pos, collapse = "    ")))
  ms.out <- c(ms.out, paste(seqs, sep = "\n"))
  
  return(list(ms.out,gaps))
  
}