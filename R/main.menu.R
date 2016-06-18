#setwd("/media/marcelo/HD2/Dropbox/Programacao/Rcodes/ABC.package/menu")
#setwd("C:/Users/Guilherme/Dropbox/ABC.package/menu")
library(ape)
.e<-new.env()

main.menu<-function(input=NULL)
    {
  if(is.null(input)==T){} else{
    read.model.input(input)
  }
  if(exists("tree", envir=.e)){
  }else{join.par()}
  if(exists("n", envir=.e)){
  }else{cur.Ne.par()}
  
  print.main.menu()  

    letter<<-readline(">>>>")
  
  switch.main.menu()

  }



