
switch.mig.menu<-function(){
  
  switch(letter,
         
         M = {prior.dist.mig<-readline("Migration prior distribution (normal or uniform?): ")
         while (prior.dist.mig %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.mig<-readline("Migration prior distribution: ")
         }
         .e$m[,6]<-prior.dist.mig
         if(exists("em",envir=.e)){
           .e$em$size[,6]<-prior.dist.mig
         }
         sys.call(which=-1)
         mig.menu()},
         
         D = {anc.mig<-readline("Different ancestral migration (YES or NO?): ")
         if(anc.mig %in% .e$YES){
           anc.mig.par()
           } else if (exists("em",envir=.e)){
           rm(em, envir=.e)
          } 
         sys.call(which=-1)
         mig.menu()},
         
         P = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
              if(.e$m[1,6]=="normal"){
              .e$m[xrow,4]<-readline(paste("migration prior (4Nm)",.e$m[xrow,1],"mean: "))
              .e$m[xrow,5]<-readline(paste("migration prior (4Nm)",.e$m[xrow,1],"Standard Deviation: "))
              }
              if(.e$m[1,6]=="uniform"){
              .e$m[xrow,4]<-readline(paste("migration prior (4Nm)",.e$m[xrow,1],"min: "))
              .e$m[xrow,5]<-readline(paste("migration prior (4Nm)",.e$m[xrow,1],"max: "))
              }
         sys.call(which=-1)
         mig.menu()},
              
          A = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
              if(.e$m[1,6]=="normal"){
                .e$em$size[xrow,4]<-readline(paste("migration prior (4Nm)",.e$em$size[xrow,1],"mean: "))
                .e$em$size[xrow,5]<-readline(paste("migration prior (4Nm)",.e$em$size[xrow,1],"Standard Deviation: "))
                }
                if(.e$m[1,6]=="uniform"){
                  .e$em$size[xrow,4]<-readline(paste("migration prior (4Nm)",.e$em$size[xrow,1],"min: "))
                  .e$em$size[xrow,5]<-readline(paste("migration prior (4Nm)",.e$em$size[xrow,1],"max: "))
                }
          sys.call(which=-1)
          mig.menu()},
              
         B = {sys.call(which=-1)
           main.menu()})
  
}
