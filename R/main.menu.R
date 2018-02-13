#' Starts the menu for model building.
#'
#' @param input a model object to be used as template. Can be left blank.
#' @return An object representing a diversification model.
#' @examples
#' my.model<-main.menu()
#' sim.sumstat(my.model)

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



