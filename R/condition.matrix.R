condition.matrix<-function(){
  size<-.e$n[,1]
  time<-.e$ej[,1]
  
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

  .e$time.matrix<-matrix(nrow=length(time),ncol=length(time))
  colnames(.e$time.matrix)<-time
  rownames(.e$time.matrix)<-time
  diag(.e$time.matrix)<-0
}

