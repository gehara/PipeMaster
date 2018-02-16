# internal function of the Model Builder
gene.menu<-function(){

print.gene.menu()

letter<<-readline(">>>>")
while(letter %in% c("M","P","1","2","3","B")==F){
  cat(paste("Choose a valid letter. You typed",letter))
  letter<<-readline(">>>>")
}

switch.gene.menu()

}
# internal function of the Model Builder
print.gene.menu<-function(){

  if(.e$loci[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$loci[1,6]=="uniform")
    dist.par<-"Min - Max"

  cat(paste("M > Mutation rate prior distribution:  ",.e$loci[1,6]),
      paste("P > priors                          ",dist.par),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$loci[,1],"  ",.e$loci[,4]," ",.e$loci[,5]),
      paste(" "),
      paste("1 > number of base pairs"),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$I[,1],":",.e$loci[,2]),
      paste(" "),
      paste("2 > locus inheritance scalar"),
      paste("                    ",c(1:nrow(.e$loci))," ",.e$I[,1],":",.e$loci[,3]),
      paste(" "),
      paste("3 > Sample structure"),
      paste("                    ",.e$I[,1],"pop:",c(1:as.numeric(.e$npops))," samples:",.e$I[,4:(3+as.numeric(.e$npops))]),
      paste(" "),
      paste("B > Back to main menu"),
      sep="\n")
}
# internal function of the Model Builder
switch.gene.menu<-function(){

  switch(letter,

         "M" = {prior.dist.mut<-readline("Mutation rate prior distribution (normal or uniform?): ")
         while (prior.dist.mut %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.mut<-readline("Migration prior distribution: ")
         }
         .e$loci[,6]<-prior.dist.mig
         sys.call(which = -1)
         gene.menu()},

         "P" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$loci))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         if(.e$loci[1,6]=="normal"){
           .e$loci[xrow,4]<-readline(paste("per base pair mutation prior",.e$I[xrow,1],"mean: "))
           .e$loci[xrow,5]<-readline(paste("per base pair mutation prior",.e$I[xrow,1],"Standard Deviation: "))
         }
         if(.e$loci[1,6]=="uniform"){
           .e$loci[xrow,4]<-readline(paste("per base pair mutation prior",.e$I[xrow,1],"min: "))
           .e$loci[xrow,5]<-readline(paste("per base pair mutation prior",.e$I[xrow,1],"max: "))
         }
         sys.call(which = -1)
         gene.menu()},

         "1" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$loci))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         .e$loci[xrow,2]<-readline(paste("number of base pairs for",.e$I[xrow,1],":"))
         sys.call(which = -1)
         gene.menu()},

         "2" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$loci))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         .e$loci[xrow,3]<-readline(paste("inheritance scalar for",.e$I[xrow,1],":"))
         sys.call(which = -1)
         gene.menu()},

         "3" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$loci))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         for(i in 1:as.numeric(.e$npops)){
           .e$loci[xrow,3+i]<-readline(paste("number of samples for",.e$I[xrow,1],"pop",i,":"))
         }
         sys.call(which = -1)
         gene.menu()},

         "B" = {sys.call(which = -1)
           main.menu()})

}
# internal function of the Model Builder
loci.par<-function(){
  ## get topology and number of nodes
  .e$ngenes<-as.numeric(readline("how many loci to simulate?: "))

  tot.gene.par<-NULL
  for (i in 1:.e$ngenes){
    gene.par<-paste("rate",i,sep="")
    tot.gene.par<-c(tot.gene.par,gene.par)
  }
  .e$loci<-matrix(nrow=.e$ngenes,ncol=6)
  .e$loci[,1]<-tot.gene.par
  .e$loci[,2]<-500
  .e$loci[,3]<-1
  .e$loci[,4]<-5e-9
  .e$loci[,5]<-1.5e-8
  .e$loci[,6]<-"uniform"

}
