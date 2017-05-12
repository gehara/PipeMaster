
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

