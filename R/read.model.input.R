#' internal function of the Model Builder
#'
read.model.input<-function(input){

  .e$ej<-input$flags$ej
  .e$n<-input$flags$n
  .e$en<-input$flags$en
  .e$m<-input$flags$m
  .e$em<-input$flags$em

  .e$ej<-gsub("rtnorm","normal",.e$ej)
  .e$n<-gsub("rtnorm","normal",.e$n)
  .e$m<-gsub("rtnorm","normal",.e$m)
  .e$en$size<-gsub("rtnorm","normal",.e$en$size)
  .e$em$size<-gsub("rtnorm","normal",.e$em$size)
  .e$en$time<-gsub("rtnorm","normal",.e$en$time)
  .e$em$time<-gsub("rtnorm","normal",.e$em$time)
  .e$loci<-gsub("rtnorm","normal",.e$loci)

  .e$ej<-gsub("runif","uniform",.e$ej)
  .e$n<-gsub("runif","uniform",.e$n)
  .e$m<-gsub("runif","uniform",.e$m)
  .e$en$size<-gsub("runif","uniform",.e$en$size)
  .e$em$size<-gsub("runif","uniform",.e$em$size)
  .e$en$time<-gsub("runif","uniform",.e$en$time)
  .e$em$time<-gsub("runif","uniform",.e$em$time)
  .e$loci<-gsub("runif","uniform",.e$loci)

  if(is.null(nrow(.e$em$size))){rm("em",envir=.e)}
  if(is.null(nrow(.e$en$size))){rm("en",envir=.e)}
  if(is.null(nrow(.e$m))){rm("m",envir=.e)}
  if(is.null(nrow(.e$ej))){rm("ej",envir=.e)}

  .e$loci<-input$loci
  .e$I<-input$I

  .e$size.matrix<-input$conds$size.matrix
  .e$mig.matrix<-input$conds$mig.matrix
  .e$time.matrix<-input$conds$time.matrix
  if(is.null(nrow(.e$mig.matrix))){rm("mig.matrix",envir=.e)}
  .e$npops<-nrow(.e$n)
  .e$tree<-input$tree

}

