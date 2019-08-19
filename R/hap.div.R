#' Haplotype Diversity
#' @description This function calculates haplotipe diversity from DNAbin sequence file
#' @param x a DNAbin object
#' @return Number of haplotypes and haplotype diversity and of x.
#' @author Marcelo Gehara
#' @references Nei, M., & Tajima, F. (1981). DNA polymorphism detectable by restriction endonucleases. Genetics, 97, 145–163.
#' @note requires Pegas package
#' @export
H.div<-function(x){
  h<-pegas::haplotype(x)
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

