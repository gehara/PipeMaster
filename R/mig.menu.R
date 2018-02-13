
mig.menu<-function(){

  print.mig.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("M","D","P","A","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.mig.menu()

}

