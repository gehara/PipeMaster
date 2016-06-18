coexp.sumstat<-function(ms.output,gene.prior){
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
  
  sspec<-site.spectrum(x,folded=T)
  sspec<-sspec/sum(na.omit(sspec))
  sspec<-sspec[1:3]
  
  SS<-c(pi,ss,H[2],TD,sspec)
  sum.stat<-rbind(sum.stat,SS)
}
#SDs<-apply(sum.stat,2,sd)
vari<-diag(var(sum.stat))
means<-colMeans(sum.stat)
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
       vari[4],means[4],skew[4],kur[4],
       vari[5],means[5],skew[5],kur[5],
       vari[6],means[6],skew[6],kur[6],
       vari[7],means[7],skew[7],kur[7])

return(h.s)
}
