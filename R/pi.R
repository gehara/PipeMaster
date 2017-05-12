pi<-function(x){
  nuc<-NULL
  for(i in 1:length(x[1,])){
  y<-H.div(x[,i])[2]
  nuc<-c(nuc,y)
}
pi<-mean(nuc)
return(pi)
}