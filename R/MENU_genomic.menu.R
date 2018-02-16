# internal function of the Model Builder
genomic.menu<-function(){

print.genomic.menu()

letter<<-readline(">>>>")
while(letter %in% c("M","P","1","2","3","B")==F){
  cat(paste("Choose a valid letter. You typed",letter))
  letter<<-readline(">>>>")
}

switch.genomic.menu()

}
# internal function of the Model Builder
print.genomic.menu<-function(){

  if(.e$loci[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$loci[1,6]=="uniform")
    dist.par<-"Min - Max"

  cat(paste("M > Mutation rate prior distribution:  ",.e$loci[1,6]),
      paste("P > priors                          ",dist.par),
      paste("                                 ",.e$loci[,4]," ",.e$loci[,5]),
      paste(" "),
      paste("1 > percentage of missing data"),
      paste("                              ",.e$loci[,1],"%"),
      paste(" "),
      paste("2 > number of base pairs"),
      paste("                              ",.e$loci[,2]),
      paste(" "),
      paste("3 > number of loci"),
      paste("                              ",.e$loci[,3]),
      paste(" "),
      paste("B > Back to main menu"),
      sep="\n")
}
# internal function of the Model Builder
switch.genomic.menu<-function(){

  switch(letter,

         "M" = {prior.dist.mut<-readline("Mutation rate prior distribution (normal or uniform?): ")
         while (prior.dist.mut %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.mut<-readline("Mutation rate prior distribution: ")
         }
         .e$loci[,6]<-prior.dist.mut
         sys.call(which = -1)
         genomic.menu()},

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
         genomic.menu()},

         "1" = {.missdata<-readline("percentage of missing data:")
         while(is.na(as.numeric(.missdata))){
           .missdata<-readline("percentage of missing data:")
         }
         .e$loci[1,1]<-.missdata
         sys.call(which = -1)
         genomic.menu()},

         "2" = {.bp<-readline("average number of base pairs:")
         while(is.na(as.numeric(.bp))){
           .bp<-readline("average number of base pairs:")
         }
         .e$loci[1,2]<-.bp
         sys.call(which = -1)
         genomic.menu()},

         "3" = {.lo<-readline("number of loci to simulate:")
         while(is.na(as.numeric(.lo))){
           .lo<-readline("average number of base pairs:")
         }
         .e$loci[1,3]<-.lo
         sys.call(which = -1)
         genomic.menu()},

         "B" = {sys.call(which = -1)
           main.menu()})

}
# internal function of the Model Builder
genomic.samples.par<-function(){

  tot.gene.par<-NULL
  for (i in 1:.e$ngenes){
    gene.par<-paste("genomic",sep="")
    tot.gene.par<-c(tot.gene.par,gene.par)
  }
  .e$I<-matrix(nr=.e$ngenes,nc=3+.e$npops)
  .e$I[,1]<-tot.gene.par
  .e$I[,2]<-"-I"
  .e$I[,3]<-.e$npops

  for(j in 1:.e$ngenes){
    for(i in 1:.e$npops){
      .e$I[j,i+3]<-readline(paste("number of samples for pop",i,":"))
    }
  }
}
# internal function of the Model Builder
genomic.par<-function(){
  ## get topology and number of nodes
  .e$ngenes<-1
  .e$nloci<-as.numeric(readline("how many loci to simulate?: "))
  bp<-as.numeric(readline("average number of base pairs: "))
  .e$loci<-matrix(nrow=1,ncol=6)
  .e$loci[,1]<-0
  .e$loci[,2]<-bp
  .e$loci[,3]<-.e$nloci
  .e$loci[,4]<-1e-11
  .e$loci[,5]<-1e-9
  .e$loci[,6]<-"uniform"

}

