#' internal function of the ms.commander
#'
sample.pars<-function(x){
  k<-sample(nrow(x),nrow(x))
  for(i in k){
    x[i,4:5]<-do.call(x[i,6],args=list(1,as.numeric(x[i,4]),as.numeric(x[i,5])),quote=F)
  }
  return(x)
}
