mcmc.sumstat<-function(ms.output,gene.prior){
sum.stat<-NULL

for (j in 1:length(ms.output)){
  
  ss<-as.numeric(strsplit(ms.output[[j]][3]," ")[[1]][2])
  if(ss==0){
  pi=0
  H=c(0,0)
  TD<-NA
  } else {
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
  H<-H.div(x)[2]
  TD<-tajima.test(x)$D
  
  }
  SS<-c(pi[[1]],ss,H,TD[1])
  sum.stat<-rbind(sum.stat,SS)
}
return(sum.stat)
}
