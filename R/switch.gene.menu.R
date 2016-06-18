
switch.gene.menu<-function(){
  
  switch(letter,
         
         M = {prior.dist.mut<-readline("Mutation rate prior distribution (normal or uniform?): ")
         while (prior.dist.mut %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.mut<-readline("Migration prior distribution: ")
         }
         .e$loci[,6]<-prior.dist.mig
         sys.call(which = -1)
         gene.menu()},
         
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
              gene.menu()},
              
          "1" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
                .e$loci[xrow,2]<-readline(paste("number of base pairs for",.e$I[xrow,1],":"))
                sys.call(which = -1)
                gene.menu()},
              
          "2" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
                  .e$loci[xrow,3]<-readline(paste("inheritance scalar for",.e$I[xrow,1],":"))
                  sys.call(which = -1)
                  gene.menu()},
         
          "3" = {xrow<-as.numeric(readline("Which locus do you want to set up? (write the reference number from the menu): "))
                for(i in 1:as.numeric(.e$npops)){
                .e$loci[xrow,3+i]<-readline(paste("number of samples for",.e$I[xrow,1],"pop",i,":"))
                }
                sys.call(which = -1)
                gene.menu()},
         
          B = {sys.call(which = -1)
            main.menu()})
  
}
