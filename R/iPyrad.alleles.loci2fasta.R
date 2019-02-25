#' iPyrad alleles.loci file to fasta alignments conversion
#' @description Converts iPyrad alleles.loci file to fasta alignmets.
#'              The alleles.loci is a specific iPyrad format that contains the phased alleles for all loci.
#'              To get this file use the * option in the output_formats parameter of iPyrad.
#' @param alleles.loci The name of the alleles.loci file as a character string. The alleles.loci file must be in the working directory.
#' @param output.dir The path to the output directory to write the fasta alignments.
#' @return Writes one fasta alignment per locus in the output.dir
#' @author Marcelo Gehara (Thanks Gustavo Cabanne for the help)
#' @export
iPyrad.alleles.loci2fasta <- function(alleles.loci, output.dir){

  x <- readLines(alleles.loci)
  setwd(output.dir)
  breaks <- grep("//",x)
    z <- 1
  for(i in 1:length(breaks)){
      fas <- x[z:(breaks[i]-1)]
      fas <-sapply(fas,strsplit," ")
      fas2 <- NULL
      for(j in 1:length(fas)){
        y <- paste(">",fas[[j]][1],sep="")
        y <- c(y, fas[[j]][length(fas[[j]])])
        fas2 <- c(fas2,y)
      }
      nam <- strsplit(x[breaks[i]],"|")[[1]]
      nam <- nam[grep("|",nam, fixed=T)[1]:grep("|",nam, fixed=T)[2]]
      nam <- gsub("|","",nam, fixed=T)
      nam <- paste(nam, collapse="")
      write(paste(fas2, sep="\n"), paste("locus_",nam,".fas",sep=""))
      z <- 1+breaks[i]
      print(paste("locus_",nam,".fas",sep=""))
  }

}
