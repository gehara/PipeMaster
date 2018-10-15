PopPlanner<-function(){
  pack <- find.package("PipeMaster")
  Sys.chmod(paths=paste(pack,"/PopPlanner.jar"), mode = 7777, use_umask = TRUE)
  system(paste(pack,"/PopPlanner.jar"))
}
