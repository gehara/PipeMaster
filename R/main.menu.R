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
    while(letter %in% c("A","B","C","D","E","F","G","H","I","Q")==F){
     cat(paste("Choose a valid letter. You typed:",letter))
     letter<<-readline(">>>>")
    }
  switch.main.menu()

  }



