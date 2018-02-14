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
