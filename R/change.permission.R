
change.permission<-function(mode="0777"){

pack<-find.package("PipeMaster")
if(Sys.info()[1]=="Darwin"){
  msABC.call<-paste(pack,"/msABCmac",sep="")
}
if(Sys.info()[1]=="Linux"){
  msABC.call<-paste(pack,"/msABClinux",sep="")
}
Sys.chmod(paths=msABC.call, mode = mode, use_umask = TRUE)
}
