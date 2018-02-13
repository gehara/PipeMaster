#' internal function of the Model Builder
#'
remove.all.par<-function(){
  options(warn=-1)
  rm(list=c("loci","I","n","en","m","em","ej","conds","tree","npops"), envir=.e)
  options(warn=0)
}
