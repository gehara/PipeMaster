
switch.genomic.menu<-function(){

  switch(letter,

         M = {prior.dist.mut<-readline("Mutation rate prior distribution (normal or uniform?): ")
         while (prior.dist.mut %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.mut<-readline("Mutation rate prior distribution: ")
         }
         .e$loci[,6]<-prior.dist.mut
         sys.call(which = -1)
         genomic.menu()},

         P = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
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

          "1" = {.e$loci[1,1]<-readline("percentage of missing data:")
                  sys.call(which = -1)
                  genomic.menu()},

         "2" = {.e$loci[1,2]<-readline("average number of base pairs:")
                  sys.call(which = -1)
                  genomic.menu()},

          "3" = {.e$loci[1,3]<-readline("number of loci to simulate:"))
                  sys.call(which = -1)
                  gene.menu()},

          B = {sys.call(which = -1)
            main.menu()})

}
