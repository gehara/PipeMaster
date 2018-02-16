#' internal function of the Model Builder
mig.menu<-function(){

  print.mig.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("M","D","P","A","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.mig.menu()

}

#' internal function of the Model Builder
print.mig.menu<-function(){

  if(.e$m[1,6]=="normal")
    dist.par<-"Mean - SD"
  if(.e$m[1,6]=="uniform")
    dist.par<-"Min - Max"

  cat(paste("M > Migration prior distribution:       ",.e$m[1,6]),
      paste("D > Different ancestral migration?         ",exists("em",envir=.e)),
      paste("P > priors                          ",dist.par),
      paste("                    ",c(1:nrow(.e$m))," ",.e$m[,1],"  ",.e$m[,4]," ",.e$m[,5]),
      paste(" "),
      if(exists("em", envir=.e))
        paste("A >  ancestral migrations:         ",dist.par),
      if(exists("em", envir=.e))
        paste("                    ",c(1:nrow(.e$em$size))," ",.e$em$size[,1],"  ",.e$em$size[,4]," ",.e$em$size[,5]),
      paste("B > Back to main menu"),
      sep="\n")
}

#' internal function of the Model Builder
switch.mig.menu<-function(){

  switch(letter,

         "M" = {prior.dist.mig<-readline("Migration prior distribution (normal or uniform?): ")
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

         "D" = {anc.mig<-readline("Different ancestral migration (YES or NO?): ")
         while(anc.mig %in% c(.e$YES,.e$NO)==F){
           cat(paste("Type a valid letter. You typed:",anc.mig))
           anc.mig<-readline("Different ancestral migration (YES or NO?): ")
         }
         if(anc.mig %in% .e$YES){
           anc.mig.par()
         } else if (exists("em",envir=.e)){
           rm(em, envir=.e)
         }
         sys.call(which=-1)
         mig.menu()},

         "P" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$m))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
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

         "A" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$em$size))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
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

         "B" = {sys.call(which=-1)
           main.menu()})

}

#' internal function of the Model Builder
mig.par<-function(){
  mig.par<-NULL
  pops<-NULL
  for(i in 1:nrow(.e$n)){
    for(j in 1:nrow(.e$n)){
      if(i==j){}
      else{
        m0<-paste("mig0.",i,"_",j,sep="")
        mig.par<-c(mig.par,m0)
        pops<-c(pops,paste(i,j))
      }
    }
  }
  .e$m<-matrix(nrow=length(mig.par),ncol=6)
  .e$m[,1]<-mig.par
  .e$m[,2]<-"-m"
  .e$m[,3]<-pops
  .e$m[,4]<-0.1
  .e$m[,5]<-1
  .e$m[,6]<-'uniform'

}

#' internal function of the Model Builder
anc.mig.par<-function(){
  mig.par<-NULL
  pops<-NULL
  for(i in 1:nrow(.e$n)){
    for(j in 1:nrow(.e$n)){
      if(i==j){
      }else{x<-readline(paste("How many changes in migration for mig",i,"_",j,": ",sep=""))
      if(x==0){
      }else{
        for(l in 1:x){
          m<-paste("mig",l,".",i,"_",j,sep="")
          mig.par<-c(mig.par,m)
          pops<-c(pops,paste(i,j))
        }
      }
      }
    }
  }

  .e$em$size<-matrix(nrow=length(mig.par),ncol=6)
  .e$em$size[,1]<-mig.par
  .e$em$size[,2]<-"-em"
  .e$em$size[,3]<-pops
  .e$em$size[,4]<-0
  .e$em$size[,5]<-0
  .e$em$size[,6]<-'uniform'

  t.mig.par<-mig.par
  for(i in 1:length(t.mig.par)){
    t.mig.par[i]<-paste("t.",mig.par[i],sep="")
  }

  .e$em$time<-matrix(nrow=length(t.mig.par),ncol=6)
  .e$em$time[,1]<-t.mig.par
  .e$em$time[,2]<-"-em"
  .e$em$time[,3]<-pops
  .e$em$time[,4]<-10000
  .e$em$time[,5]<-20000
  .e$em$time[,6]<-'uniform'

}

