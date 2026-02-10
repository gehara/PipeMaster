#' Calls PopPlanner, a java application with GUI to build ms population models generating the corresponding ms string.
#' @description When you run this function PopPlanner will pop up. PopPlanner has an easy to use interface to build diversification models.
#'            The resulting ms string will be generated in the bottom of the PopPlanner window. This can be used as an input to the main.menu function.
#' @export
PopPlanner<-function(){
  x<-getwd()
  pack <- find.package("PipeMaster")
  Sys.chmod(paths = pack, mode = "7777", use_umask = T)
  setwd(pack)
  system(paste(pack,"/PopPlanner.jar", sep=""))
  setwd(x)
}
