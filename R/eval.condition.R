eval.condition<-function(x,y){
  value<-NULL
  for(i in 1:length(y)){
    value[i]<-as.numeric(x[y[[i]][1],4])>as.numeric(x[y[[i]][2],4])
  }
  return(sum(value))
}