#' Calls PopPlanner, a java aplication with GUI to build ms population models generating the coresponding ms string.
#' @description When you run this function PopPlanner will pop up. PopPlanner has a easy to use interface to build a diversification models.
#'            The resulting ms string will be generated in the botton of the PopPlanner window. This can be used as an input to the main.menu function.
#' @export
PopPlanner<-function(){
  pack <- find.package("PipeMaster")
  Sys.chmod(paths = pack, mode = "7777", use_umask = T)
  system(paste(pack,"/PopPlanner.jar", sep=""))
}
