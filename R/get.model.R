

get.model<-function(){

    
    .e$ej<-gsub("normal","rtnorm",.e$ej)
    .e$n<-gsub("normal","rtnorm",.e$n)
    .e$m<-gsub("normal","rtnorm",.e$m)
    .e$en$size<-gsub("normal","rtnorm",.e$en$size)
    .e$em$size<-gsub("normal","rtnorm",.e$em$size)
    .e$en$time<-gsub("normal","rtnorm",.e$en$time)
    .e$em$time<-gsub("normal","rtnorm",.e$em$time)
    .e$loci<-gsub("normal","rtnorm",.e$loci)
    
    .e$ej<-gsub("uniform","runif",.e$ej)
    .e$n<-gsub("uniform","runif",.e$n)
    .e$m<-gsub("uniform","runif",.e$m)
    .e$en$size<-gsub("uniform","runif",.e$en$size)
    .e$em$size<-gsub("uniform","runif",.e$em$size)
    .e$en$time<-gsub("uniform","runif",.e$en$time)
    .e$em$time<-gsub("uniform","runif",.e$em$time)
    .e$loci<-gsub("uniform","runif",.e$loci)
    
    if(is.null(nrow(.e$em$size))){rm("em",envir=.e)}
    if(is.null(nrow(.e$en$size))){rm("en",envir=.e)}
    if(is.null(nrow(.e$m))){rm("m",envir=.e)}
    
      
    model<-list(NULL,NULL,NULL,NULL,NULL)
    names(model)<-c("loci","I","flags","conds","tree")
    model$loci<-.e$loci
    model$I<-.e$I

    flags<-list(NULL,NULL,NULL,NULL,NULL)
    names(flags)<-c("n","m","en","em","ej")
    flags$n <- .e$n
    flags$m <- .e$m
    flags$en <- .e$en
    flags$em <- .e$em
    flags$ej <- .e$ej

    model$flags<-flags

    conds<-list(NULL,NULL,NULL)
    names(conds)<-c("size.matrix","mig.matrix","time.matrix")
    conds$size.matrix<-.e$size.matrix
    conds$mig.matrix<-.e$mig.matrix
    conds$time.matrix<-.e$time.matrix

    model$conds<-conds
    model$tree<-.e$tree

    rm(list=ls(envir=.e),envir=.e)
    return(model)

}

