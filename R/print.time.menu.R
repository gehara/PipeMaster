print.time.menu<-function()
  {
  if(is.null(.e$ej)==F){
    if(.e$ej[1,6]=="normal")
      dist.par<-"Mean - SD"
    if(.e$ej[1,6]=="uniform")
      dist.par<-"Min - Max"
  }

  cat(if(is.null(.e$ej)==F)
      paste("P > Time prior distribution:    ",.e$ej[1,6]),
      if(is.null(.e$ej)==F)
      paste("    Time priors                  ",dist.par),
      if(is.null(.e$ej)==F)
      paste("   J >  time of junctions: "),
      if(is.null(.e$ej)==F)
      paste("                    ",c(1:nrow(.e$ej)),"  ",.e$ej[,1],"  ",.e$ej[,4]," ",.e$ej[,5]),
      paste(" "),
      if(exists("en", envir=.e))
      paste("   N >  time of ancestral Ne change: "),
      if(exists("en", envir=.e))
      paste("            ",c(1:nrow(.e$en$time)),"  ",.e$en$time[,1],"  ",.e$en$time[,4]," ",.e$en$time[,5]),
      paste(" "),
      if(exists("em", envir=.e))
      paste("   M >  time of migration change: "),
      if(exists("em", envir=.e))
      paste("                ",c(1:nrow(.e$em$time)),"  ",.e$em$time[,1],"  ",.e$em$time[,4]," ",.e$em$time[,5]),
      paste("B > Back to Ne menu"),
      sep="\n")
}
