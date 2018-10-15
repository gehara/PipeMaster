# internal function of the Model Builder
condition.menu<-function(){

  print.condition.menu()

  letter<<-readline(">>>>")
  while(letter %in% c("S","M","T","1","2","3","B")==F){
    cat(paste("Choose a valid letter. You typed",letter))
    letter<<-readline(">>>>")
  }

  switch.condition.menu()

}
# internal function of the Model Builder
switch.condition.menu<-function(){
  switch(letter,

         "S" = {place.size.condition()
           sys.call(which = -1)
           condition.menu()},

         "M" = {place.mig.condition()
           sys.call(which = -1)
           condition.menu()},

         "T" = {place.time.condition()
           sys.call(which = -1)
           condition.menu()},


         "1" = {print(.e$size.matrix)
           print("-----------------")
           sys.call(which = -1)
           condition.menu()},

         "2" = {print(.e$mig.matrix)
           print("-----------------")
           sys.call(which = -1)
           condition.menu()},

         "3" = {print(.e$time.matrix)
           print("-----------------")
           sys.call(which = -1)
           condition.menu()},

         "B" = {sys.call(which = -1)
           main.menu()}

  )}
# internal function of the Model Builder
print.condition.menu<-function(){

  {cat(paste("size parameter               -- ",colnames(.e$size.matrix)),
       paste("   "),
       if(exists("mig.matrix",envir=.e))
         paste("mig parameter                -- ",colnames(.e$mig.matrix)),
       paste("   "),
       paste("time parameter               -- ",colnames(.e$time.matrix)),
       paste("   "),
       paste("S > place a size condition"),
       if(exists("mig.matrix",envir=.e))
         paste("M > place a mig condition"),
       paste("T > place a time condition"),
       paste(""),
       paste("1 > see size matrix"),
       if(exists("mig.matrix",envir=.e))
         paste("2 > see mig  matrix"),
       paste("3 > see time matrix"),

       paste("B > back to main menu"),
       sep="\n")
  }
}
# internal function of the Model Builder
condition.matrix<-function(){
  size <- .e$n[,1]

  if(exists("ej", envir=.e)){
    time<-.e$ej[,1]
  }

  if(exists("m",envir=.e)){
    mig<-.e$m[,1]
  }

  if(exists("en", envir=.e)){
    size<-c(size,.e$en$size[,1])
    time<-c(time,.e$en$time[,1])
  }

  if(exists("em", envir=.e)){
    mig<-c(mig,.e$em$size[,1])
    time<-c(time,.e$em$time[,1])
  }

  .e$size.matrix<-matrix(nrow=length(size),ncol=length(size))
  colnames(.e$size.matrix)<-size
  rownames(.e$size.matrix)<-size
  diag(.e$size.matrix)<-0

  if(exists("mig")){
    .e$mig.matrix<-matrix(nrow=length(mig),ncol=length(mig))
    colnames(.e$mig.matrix)<-mig
    rownames(.e$mig.matrix)<-mig
    diag(.e$size.matrix)<-0
  }

  if(exists("time")){
    .e$time.matrix<-matrix(nrow=length(time),ncol=length(time))
    colnames(.e$time.matrix)<-time
    rownames(.e$time.matrix)<-time
    diag(.e$time.matrix)<-0
  }
}
# internal function of the Model Builder
inv.mirror.lower<-function(x) {
  x1<-t(x)[lower.tri(x, diag=F)]
  for(i in 1:length(x1)) {
    if(!x1[i] %in% c("<",">")) next
    if(x1[i]=="<") {
      x1[i]<-">"
      next}
    else x1[i]<-"<"
  }
  x[lower.tri(x, diag=F)]<-x1
  return(x)
}
# internal function of the Model Builder
inv.mirror.upper<-function(x) {
  x1<-t(x)[upper.tri(x, diag=F)]
  for(i in 1:length(x1)) {
    if(!x1[i] %in% c("<",">")) next
    if(x1[i]=="<") {
      x1[i]<-">"
      next}
    else x1[i]<-"<"
  }
  x[upper.tri(x, diag=F)]<-x1
  return(x)
}
# internal function of the Model Builder
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
# internal function of the Model Builder
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
# internal function of the Model Builder
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

