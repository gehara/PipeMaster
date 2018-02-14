#' internal function of the Model Builder
print.mig.menu<-function()
  {

  if(.e$m[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$m[1,6]=="uniform")
    dist.par<-"Min - Max"

  cat(paste("M > Migration prior distribution:       ",.e$m[1,6]),
      paste("D > Different ancestral migration?         ",exists("em",envir=.e)),
      paste("P > priors                          ",dist.par),
      paste("                    ",c(1:nrow(.e$m))," ",.e$m[,1],"  ",.e$m[,4]," ",.e$m[,5]),
      paste(" "),
      if(exists("em", envir=.e))
        paste("A >  ancestral migrations:         ",dist.par),
      if(exists("em", envir=.e))
      paste("                    ",c(1:nrow(.e$em$size))," ",.e$em$size[,1],"  ",.e$em$size[,4]," ",.e$em$size[,5]),
      paste("B > Back to main menu"),
      sep="\n")
}
