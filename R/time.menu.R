
time.menu<-function(){

  print.time.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("P","J","M","N","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.time.menu()

}

