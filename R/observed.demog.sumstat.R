observed.demog.sumstat<-function(path.to.fasta,fasta.files){

  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fas",fasta.files)]

  ms.output<-fasta2ms(path.to.fasta,fasta.files,write.file=F)
  bp.length<-list()
  ss<-list()
  for(i in 1:length(ms.output)){
    fas<-read.dna(file=fasta.files[i], format="fasta")
    bp.length[[i]]<-ncol(fas)
    ss[[i]]<-as.numeric(strsplit(ms.output[[i]][3]," ")[[1]][2])
    if(is.null(ss[[i]])){
      stop(paste("sequence",fasta.files[i],"has zero variation"))
    }

  }

  sum.stat<-NULL
  for (j in 1:length(ms.output)){
    x<-ms.to.DNAbin(ms.output = ms.output[[j]],bp.length = bp.length[[j]])

    pi<-nuc.div(x)
    H<-H.div(x)[2]
    TD<-tajima.test(x)$D
    #R2<-R2.test(x,B=0,plot = F,quiet = T)
    spec<-site.spectrum(x)[1:3]
    SS<-c(pi,ss[[j]],H,TD,spec)

    sum.stat<-rbind(sum.stat,SS)
    }
rownames(sum.stat)<-fasta.files
colnames(sum.stat)<-c("pi","ss","H","TD","ss1","ss2","ss3")
return(sum.stat)
}
