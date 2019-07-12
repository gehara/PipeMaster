# internal function to find the simulator
get.msABC<-function(){
  pack<-find.package("PipeMaster")

  if(Sys.info()[1]=="Linux"){
    msABC.call<-paste(pack,"/msABClinux",sep="")
  }
  if(Sys.info()[1]=="Darwin"){
    msABC.call<-paste(pack,"/msABCmac",sep="")
  }
  Sys.chmod(paths=msABC.call, mode = 7777, use_umask = TRUE)
  return(msABC.call)
  }

