#' internal function of the Model Builder
#'

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

