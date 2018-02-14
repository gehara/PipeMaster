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
#' 
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

