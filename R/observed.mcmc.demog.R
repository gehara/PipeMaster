observed.mcmc.demog<-function(path.to.fasta){
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
    
    SS<-c(pi[[1]],ss[[j]],H,TD[1])
    sum.stat<-rbind(sum.stat,SS)
    }
    
#write.table(t(h.s),file="h.sum.stat.txt",col.names = F, row.names=F, append=T,sep="\t")  
return(sum.stat)
}
