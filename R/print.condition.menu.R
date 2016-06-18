
print.condition.menu<-function(){
   
{cat(paste("size parameter               -- ",colnames(.e$size.matrix)),
     paste("   "),
     if(exists("mig.matrix",envir=.e))
     paste("mig parameter                -- ",colnames(.e$mig.matrix)),
     paste("   "),
     paste("time parameter               -- ",colnames(.e$time.matrix)),
     paste("   "),
     paste("s > place a size condition"),
     if(exists("mig.matrix",envir=.e))
     paste("m > place a mig condition"),
     paste("t > place a time condition"),
     paste(""),
     paste("1 > see size matrix"),
     if(exists("mig.matrix",envir=.e))
     paste("2 > see mig  matrix"),
     paste("3 > see time matrix"),
      
     paste("B > back to main menu"),
     sep="\n")
}
}
