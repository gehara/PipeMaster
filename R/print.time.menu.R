print.time.menu<-function()
  {
  if(exists("ej",envir=.e)==T){
    if(.e$ej[1,6]=="normal")
      dist.par<-"Mean - SD"
    if(.e$ej[1,6]=="uniform")
      dist.par<-"Min - Max"
  } else {
    .e$ej<-NULL
    dist.par<-NULL}

  cat(paste("A > Time prior distribution:    ",.e$ej[1,6]),
      paste("    Time priors                  ",dist.par),
      paste("   j >  time of junctions: "),
      paste("                    ",c(1:nrow(.e$ej)),"  ",.e$ej[,1],"  ",.e$ej[,4]," ",.e$ej[,5]),
      paste(" "),
      if(exists("en", envir=.e))
      paste("   n >  time of ancestral Ne change: "),
      if(exists("en", envir=.e))
      paste("            ",c(1:nrow(.e$en$time)),"  ",.e$en$time[,1],"  ",.e$en$time[,4]," ",.e$en$time[,5]),
      paste(" "),
      if(exists("em", envir=.e))
      paste("   m >  time of migration change: "),
      if(exists("em", envir=.e))
      paste("                ",c(1:nrow(.e$em$time)),"  ",.e$em$time[,1],"  ",.e$em$time[,4]," ",.e$em$time[,5]),
      paste("B > Back to Ne menu"),
      sep="\n")
}
