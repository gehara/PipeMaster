
switch.time.menu<-function(){
  switch(letter,
         "P" = {prior.dist.Ne<-readline("time prior distribution (normal or uniform): ")
              while (prior.dist.Ne %in% c("normal","uniform")==F){
           print("Possible distributions are normal or uniform!")
           prior.dist.Ne<-readline("time prior distribution: ")
         }
         .e$ej[,6]<-prior.dist.Ne
         if(exists("en",envir=.e)){
           .e$en$time[,6]<-prior.dist.Ne
         }
         if(exists("em",envir=.e)){
           .e$em$time[,6]<-prior.dist.Ne
         }
         sys.call(which = -1)
        time.menu()},

         "J" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$ej))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         if(.e$ej[1,6]=="normal"){
           .e$ej[xrow,4]<-readline(paste("Time of junction in generations",.e$ej[xrow,1],"mean: "))
           .e$ej[xrow,5]<-readline(paste("Time of junction in generations",.e$ej[xrow,1],"Standard Deviation: "))
         }
         if(.e$n[1,6]=="uniform"){
           .e$ej[xrow,4]<-readline(paste("Time of junction in generations",.e$ej[xrow,1],"min: "))
           .e$ej[xrow,5]<-readline(paste("Time of junction in generations",.e$ej[xrow,1],"max: "))
         }
         sys.call(which = -1)
           time.menu()},

          "N" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
          while(xrow %in% c(1:nrow(.e$en$time))==F){
            cat(paste("Type a valid number. You typed:",xrow))
            xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
          }
          if(.e$en$time[1,6]=="normal"){
            .e$en$time[xrow,4]<-readline(paste("Time of Ne change in generations",.e$en$time[xrow,1],"mean: "))
            .e$en$time[xrow,5]<-readline(paste("Time of Ne change in generations",.e$en$time[xrow,1],"Standard Deviation: "))
          }
          if(.e$en$time[1,6]=="uniform"){
            .e$en$time[xrow,4]<-readline(paste("Time of Ne change in generations",.e$en$time[xrow,1],"min: "))
            .e$en$time[xrow,5]<-readline(paste("Time of Ne change in generations",.e$en$time[xrow,1],"max: "))
          }
          sys.call(which = -1)
          time.menu()
           },


         "M" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$em$time))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
         if(.e$em$time[1,6]=="normal"){
           .e$em$time[xrow,4]<-readline(paste("Time of mig change in generations",.e$em$time[xrow,1],"mean: "))
           .e$em$time[xrow,5]<-readline(paste("Time of mig change in generations",.e$em$time[xrow,1],"Standard Deviation: "))
         }
         if(.e$em$time[1,6]=="uniform"){
           .e$em$time[xrow,4]<-readline(paste("Time of mig change in generations",.e$em$time[xrow,1],"min: "))
           .e$em$time[xrow,5]<-readline(paste("Time of mig change in generations",.e$em$time[xrow,1],"max: "))
         }
         sys.call(which = -1)
         time.menu()},

         "B" = {sys.call(which = -1)
           main.menu()})

}

