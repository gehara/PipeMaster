# read ms output from Phyclust to DNAbin file

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

