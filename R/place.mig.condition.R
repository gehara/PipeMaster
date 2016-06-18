place.mig.condition<-function(){
  print(.e$mig.matrix)
  cond<-readline("Write the name of 2 parameters with a logic sign inbetween ( >  or < or = ) separated by a space.
                Ex: Ne0.pop1 > Ne0.pop2 :   ")
  cond<-strsplit(cond," ")
  y<-grep(cond[[1]][1],rownames(.e$mig.matrix))
  x<-grep(cond[[1]][3],colnames(.e$mig.matrix))
  .e$mig.matrix[y,x]<-cond[[1]][2]
  .e$mig.matrix<-inv.mirror.lower(.e$mig.matrix)
  
}