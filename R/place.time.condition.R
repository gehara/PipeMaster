#' internal function of the Model Builder
place.time.condition<-function(){
  print(.e$time.matrix)
cond<-readline("Write the name of 2 parameters with a logic sign inbetween ( >  or < or = ) separated by a space.
                Ex: Ne0.pop1 > Ne0.pop2 :   ")
cond<-strsplit(cond," ")
y<-grep(cond[[1]][1],rownames(.e$time.matrix))
x<-grep(cond[[1]][3],colnames(.e$time.matrix))
.e$time.matrix[y,x]<-cond[[1]][2]
.e$time.matrix<-inv.mirror.lower(.e$time.matrix)

}

