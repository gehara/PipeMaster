get.msABC<-function(){
  pack<-find.package("PipeMaster")

  if(Sys.info()[1]=="Linux"){
    msABC.call<-paste(pack,"/msABClinux",sep="")
    }
  return(msABC.call)
  }
