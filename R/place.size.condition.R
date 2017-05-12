place.size.condition<-function(){
  print(.e$size.matrix)
cond<-readline("Write the name of 2 parameters with a logic sign inbetween ( >  or < or = ) separated by a space.
               Ex.: Ne0.pop1 < Ne0.pop2")
cond<-strsplit(cond," ")
y<-grep(cond[[1]][1],rownames(.e$size.matrix))
x<-grep(cond[[1]][3],colnames(.e$size.matrix))
.e$size.matrix[y,x]<-cond[[1]][2]
.e$size.matrix<-inv.mirror.lower(.e$size.matrix)

}
