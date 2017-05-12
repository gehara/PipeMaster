print.gene.menu<-function()
  {
  
  if(.e$loci[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$loci[1,6]=="uniform")
    dist.par<-"Min - Max"
  
  cat(paste("M > Mutation rate prior distribution:  ",.e$loci[1,6]),
      paste("P > priors                          ",dist.par),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$loci[,1],"  ",.e$loci[,4]," ",.e$loci[,5]),
      paste(" "),
      paste("1 > number of base pairs"),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$I[,1],":",.e$loci[,2]),
      paste(" "),
      paste("2 > locus inheritance scalar"),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$I[,1],":",.e$loci[,3]),
      paste(" "),
      paste("3 > Sample structure"),
      paste("                    ",.e$I[,1],"pop:",c(1:as.numeric(.e$npops))," samples:",.e$I[,4:(3+as.numeric(.e$npops))]),
      paste(" "),
      paste("B > Back to main menu"),
      sep="\n")
}
