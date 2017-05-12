
switch.demog.menu<-function(){
  
  switch(letter,
         
         "1" = {prior.dist.Ne<-readline("Ne prior distribution (normal or uniform): ")
         while (prior.dist.Ne %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.Ne<-readline("Ne prior distribution: ")
         }
         .e$n[,6]<-prior.dist.Ne
         if(exists("en",envir=.e)){
           .e$en$size[,6]<-prior.dist.Ne
         }
         sys.call(which=-1)
         demog.menu()},
         
         "2" = {anc.Ne<-readline("Different ancestral NE (YES or NO?): ")
         if(anc.Ne %in% .e$YES){
           anc.Ne.par()
           sys.parent(n=1)
           } else if (anc.Ne %in% .e$NO){
             if (exists("en",envir=.e))
               rm(en, envir=.e)
           }
         sys.call(which=-1)
         demog.menu()
         },
         
         "3" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         if(.e$n[1,6]=="normal"){
           .e$n[xrow,4]<-readline(paste("Ne prior (4Nm)",.e$n[xrow,1],"mean: "))
           .e$n[xrow,5]<-readline(paste("Ne prior (4Nm)",.e$n[xrow,1],"Standard Deviation: "))
         }
         if(.e$n[1,6]=="uniform"){
           .e$n[xrow,4]<-readline(paste("Ne prior (4Nm)",.e$n[xrow,1],"min: "))
           .e$n[xrow,5]<-readline(paste("Ne prior (4Nm)",.e$n[xrow,1],"max: "))
         }
         sys.call(which=-1)
         demog.menu()
           },
           
           
         "4" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         if(.e$en$size[1,6]=="normal"){
           .e$en$size[xrow,4]<-readline(paste("Ancestral Ne prior (4Nm)",.e$en$size[xrow,1],"mean: "))
           .e$en$size[xrow,5]<-readline(paste("Ancestral Ne prior (4Nm)",.e$en$size[xrow,1],"Standard Deviation: "))
         }
         if(.e$en$size[1,6]=="uniform"){
           .e$en$size[xrow,4]<-readline(paste("Ancestral Ne prior (4Nm)",.e$en$size[xrow,1],"min: "))
           .e$en$size[xrow,5]<-readline(paste("Ancestral Ne prior (4Nm)",.e$en$size[xrow,1],"max: "))
         }
         sys.call(which=-1)
         demog.menu()
         },
         
         B = {sys.call(which=-1)
           main.menu()})
  
}

