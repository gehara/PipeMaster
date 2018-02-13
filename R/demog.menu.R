#' internal function of the Model Builder
#'
demog.menu<-function(){

  print.demog.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("N","D","C","A","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.demog.menu()

}

