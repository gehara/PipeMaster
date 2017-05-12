# calculates haplotipe diversity from DNAbin sequence file 
# needs pegas package
# x = a DNAbin object

H.div<-function(x){
  h<-haplotype(x)
  hap<-attr(h, "index")
  n.hap<-length(hap)
  
  h.freqs<-NULL
  for(i in 1:n.hap){
    freq<-length(hap[[i]])/nrow(x)
    h.freqs<-c(h.freqs,freq)
    }
  
  H.d = (nrow(x)/(nrow(x)-1))*(1 - sum(h.freqs^2))
  return(c(n.hap,H.d))
}







