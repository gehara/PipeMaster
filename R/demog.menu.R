
demog.menu<-function(){

  print.demog.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("N","D","C","A","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    sys.call(which = 0)
    print.demog.menu()
    letter<<-readline(">>>>")
  }

  switch.demog.menu()

}

