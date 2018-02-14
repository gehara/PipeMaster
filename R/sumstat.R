
sumstat<-function(ms.output,gene.prior){
sum.stat<-NULL

for (j in 1:length(ms.output)){

  ss<-as.numeric(strsplit(ms.output[[j]][3]," ")[[1]][2])
  x<-ms.output[[j]][5:length(ms.output[[j]])]
  x<-gsub("0","A",x)
  x<-gsub("1","C",x)
  se<-list(NULL,NULL,NULL,NULL)
  names(se)<-c("nb","seq","nam","com")
  se$nb<-length(x) # number of samples
  se$seq<-x # sequencies
  se$nam<-c(1:length(x)) # sequence names, just numbers
  se$com<-NA
  class(se)<-"alignment" # this is an alignment
  x<-as.DNAbin(se)

  pi<-nuc.div(x)*ss/gene.prior[j,5]
  H<-H.div(x)
  TD<-tajima.test(x)$D
  #R2<-R2.test(x,B=0,plot = F,quiet = T)
  spec<-site.spectrum(x)[1:3]
  SS<-c(pi,ss,H[2],TD,spec)
  sum.stat<-rbind(sum.stat,SS)
}
return(sum.stat)
}
