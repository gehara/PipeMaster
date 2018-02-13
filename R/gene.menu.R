gene.menu<-function(){

print.gene.menu()

letter<<-readline(">>>>")
while(letter %in% c("M","P","1","2","3","B")==F){
  cat(paste("Choose a valid letter. You typed",letter))
  letter<<-readline(">>>>")
}

switch.gene.menu()

}

