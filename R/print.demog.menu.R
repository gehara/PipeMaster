print.demog.menu<-function()
  {
  if(.e$n[1,6]=="normal")
    dist.par<-"Mean, SD"
  if(.e$n[1,6]=="uniform")
    dist.par<-"min, max"
  
  cat(paste("1 > Ne prior distribution:              ",.e$n[1,6]),
      paste("2 > Different ancestral Ne?             ",exists("en",envir=.e)),
      paste("3 > current Ne priors                      ",dist.par),
      paste("                ",c(1:nrow(.e$n))," ",.e$n[,1],"      ",.e$n[,4],"     ",.e$n[,5]),
      paste(" "),
      if(exists("en",envir=.e))
      paste("4 > ancestral Ne priors"),
      if(exists("en",envir=.e))
      paste("                ",c(1:nrow(.e$en$size))," ",.e$en$size[,1],"      ",.e$en$size[,4],"    ",.e$en$size[,5]),
      paste("B > Back to main menu"),
      sep="\n")
}
