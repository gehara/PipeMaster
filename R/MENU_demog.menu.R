# internal function of the Model Builder
demog.menu<-function(){

  print.demog.menu()

  letter<<-toupper(readline(">>>>"))
  while(letter %in% c("N","D","C","A","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-toupper(readline(">>>>"))
  }

  switch.demog.menu()

}
# internal function of the Model Builder
print.demog.menu<-function(){
  if(.e$n[1,6]=="normal")
    dist.par<-"Mean, SD"
  if(.e$n[1,6]=="uniform")
    dist.par<-"min, max"

  cat(paste("N > Ne prior distribution:              ",.e$n[1,6]),
      paste("D > Different ancestral Ne?             ",exists("en",envir=.e)),
      paste("C > current Ne priors                      ",dist.par),
      paste("                ",c(1:nrow(.e$n))," ",.e$n[,1],"      ",.e$n[,4],"     ",.e$n[,5]),
      paste(" "),
      if(exists("en",envir=.e))
        paste("A > ancestral Ne priors"),
      if(exists("en",envir=.e))
        paste("                ",c(1:nrow(.e$en$size))," ",.e$en$size[,1],"      ",.e$en$size[,4],"    ",.e$en$size[,5]),
      paste("B > Back to main menu"),
      sep="\n")
}
# internal function of the Model Builder
switch.demog.menu<-function(){

  switch(letter,

         "N" = {prior.dist.Ne<-readline("Ne prior distribution (normal or uniform): ")
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

         "D" = {anc.Ne<-readline("Different ancestral Ne (YES or NO?): ")
         while(anc.Ne %in% c(.e$YES,.e$NO)==F){
           cat(paste("Type a valid letter. You typed:",anc.Ne))
           anc.Ne<-readline("Different ancestral migration (YES or NO?): ")
         }
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

         "C" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$n))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
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


         "A" = {xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         while(xrow %in% c(1:nrow(.e$en$size))==F){
           cat(paste("Type a valid number. You typed:",xrow))
           xrow<-as.numeric(readline("Which parameter do you want to set up? (write the reference number from the menu): "))
         }
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

         "B" = {sys.call(which=-1)
           main.menu()})

}
# internal function of the Model Builder
anc.Ne.par<-function(){
  anc.Ne.par<-NULL
  time.anc.Ne.par<-NULL
  pop<-NULL
  for (i in 1:nrow(.e$n)){
    n.anc.pop<-readline(paste("How many changes in Ne for pop ",i,"?",sep=""))
    if (n.anc.pop==0){
    } else (
      for (j in 1:n.anc.pop){
        Ne.par<-paste("Ne",j,".pop",i,sep="")
        anc.Ne.par<-c(anc.Ne.par,Ne.par)
        time.Ne.par<-paste("t.Ne",j,".pop",i,sep="")
        time.anc.Ne.par<-c(time.anc.Ne.par,time.Ne.par)
        pop<-c(pop,i)
      }
    )
  }

  .e$en$size<-matrix(nrow=length(anc.Ne.par),ncol=6)
  .e$en$size[,1]<-anc.Ne.par
  .e$en$size[,2]<-'-en'
  .e$en$size[,3]<-pop
  .e$en$size[,4]<-1000
  .e$en$size[,5]<-10000
  .e$en$size[,6]<-'uniform'

  .e$en$time<-matrix(nrow=length(anc.Ne.par),ncol=6)
  .e$en$time[,1]<-time.anc.Ne.par
  .e$en$time[,2]<-'-en'
  .e$en$time[,3]<-pop
  .e$en$time[,4]<-10000
  .e$en$time[,5]<-100000
  .e$en$time[,6]<-'uniform'
}
# internal function of the Model Builder
cur.Ne.par<-function(){
  list.Ne.pars<-NULL
  for (i in 1:.e$npops){
    Ne0.par<-paste("Ne0.pop",i,sep="")
    list.Ne.pars<-c(list.Ne.pars,Ne0.par)
  }
  .e$n<-matrix(nrow=length(list.Ne.pars), ncol=6)
  .e$n[,1]<-list.Ne.pars
  .e$n[,2]<-'-n'
  .e$n[,3]<-c(1:length(list.Ne.pars))
  .e$n[,4]<-100000
  .e$n[,5]<-500000
  .e$n[,6]<-"uniform"

}

