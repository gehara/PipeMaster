# internal function of the Model Builder
get.msABC<-function(){
  pack<-find.package("PipeMaster")

  if(Sys.info()[1]=="Linux"){
    msABC.call<-paste(pack,"/msABClinux",sep="")
  }
  if(Sys.info()[1]=="Darwin"){
    msABC.call<-paste(pack,"/msABCmac",sep="")
  }
  Sys.chmod(paths=msABC.call, mode = 777, use_umask = TRUE)
  return(msABC.call)
  }
