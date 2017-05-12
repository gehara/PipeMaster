observed.coexp.sumstat<-function(path.to.fasta){
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
  }

  sum.stat<-NULL
  for (j in 1:length(ms.output)){
    x<-ms.to.DNAbin(ms.output = ms.output[[j]],bp.length = bp.length[[j]])

    pi<-nuc.div(x)
    H<-H.div(x)[2]
    TD<-tajima.test(x)$D

    #sfs<-site.spectrum(x,folded=T)
    #sfs<-sfs/sum(na.omit(sfs))
    #sfs<-sfs[1:3]

    SS<-c(pi[[1]],ss[[j]],H,TD[1])
    sum.stat<-rbind(sum.stat,SS)
    }

vari<-diag(var(sum.stat))
means<-colMeans(sum.stat,na.rm=T)
skew<-NULL
kur<-NULL
for(u in 1:ncol(sum.stat)){
  s<-skewness(sum.stat[,u])
  skew<-c(skew,s)
  k<-kurtosis(sum.stat[,u])
  kur<-c(kur,k)
}

h.s<-c(vari[1],means[1],skew[1],kur[1],
       vari[2],means[2],skew[2],kur[2],
       vari[3],means[3],skew[3],kur[3],
       vari[4],means[4],skew[4],kur[4])

names(h.s)<-c("var.pi","mean.pi","skew.pi","kur.pi",
              "var.ss","mean.ss","skew.ss","kur.ss",
              "var.H","mean.H","skew.H","kur.H",
              "var.TD","mean.TD","skew.TD","kur.TD")
#write.table(t(h.s),file="h.sum.stat.txt",col.names = F, row.names=F, append=T,sep="\t")
return(h.s)
}
