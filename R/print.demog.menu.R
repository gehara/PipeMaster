print.demog.menu<-function()
  {
  if(.e$n[1,6]=="normal")
    dist.par<-"Mean, SD"
  if(.e$n[1,6]=="uniform")
    dist.par<-"min, max"

  cat(paste("N > Ne prior distribution:              ",.e$n[1,6]),
      paste("D > Different ancestral Ne?             ",exists("en",envir=.e)),
      paste("C > current Ne priors                      ",dist.par),
      paste("                ",c(1:nrow(.e$n))," ",.e$n[,1],"      ",.e$n[,4],"     ",.e$n[,5]),
      paste(" "),
      if(exists("en",envir=.e))
      paste("A > ancestral Ne priors"),
      if(exists("en",envir=.e))
      paste("                ",c(1:nrow(.e$en$size))," ",.e$en$size[,1],"      ",.e$en$size[,4],"    ",.e$en$size[,5]),
      paste("B > Back to main menu"),
      sep="\n")
}
