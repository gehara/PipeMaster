print.genomic.menu<-function()
  {

  if(.e$loci[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$loci[1,6]=="uniform")
    dist.par<-"Min - Max"

  cat(paste("M > Mutation rate prior distribution:  ",.e$loci[1,6]),
      paste("P > priors                          ",dist.par),
      paste("                              ",.e$loci[,4]," ",.e$loci[,5]),
      paste(" "),
      paste("1 > percentage of missing data"),
      paste("                              ",.e$loci[,1],"%"),
      paste(" "),
      paste("2 > number of base pairs"),
      paste("                              ",.e$loci[,2]),
      paste(" "),
      paste("3 > number of loci"),
      paste("                              ",.e$loci[,3]),
      paste(" "),
      paste("B > Back to main menu"),
      sep="\n")
}
