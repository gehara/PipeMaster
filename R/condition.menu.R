#' internal function of the Model Builder
condition.menu<-function(){

  print.condition.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("S","M","T","1","2","3","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.condition.menu()

}

