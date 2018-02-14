#' internal function of the Model Builder
genomic.menu<-function(){

print.genomic.menu()

letter<<-readline(">>>>")
while(letter %in% c("M","P","1","2","3","B")==F){
  cat(paste("Choose a valid letter. You typed",letter))
  letter<<-readline(">>>>")
}

switch.genomic.menu()

}

