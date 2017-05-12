get.ss.pop.name<-function(pops){
  
  S.s<-NULL
  pi<-NULL
  H.pop<-NULL
  TD<-NULL
  FuLiD<-NULL
  FuLiF<-NULL

  for(i in 1:length(pops)){
    x<-paste("ss.pop",i,sep="")
    S.s<-c(S.s,x)
    }

  for(i in 1:length(pops)){
    x<-paste("pi.pop",i,sep="")
    pi<-c(pi,x)
    }

  for(i in 1:length(pops)){
    x<-paste("H.pop",i,sep="")
    H.pop<-c(H.pop,x)
    }

  for(i in 1:length(pops)){
    x<-paste("TajD.pop",i,sep="")
    TD<-c(TD,x)
    }

  for(i in 1:length(pops)){
    x<-paste("FuLiD.pop",i,sep="")
    FuLiD<-c(FuLiD,x)
    }

  for(i in 1:length(pops)){
    x<-paste("FuLiF.pop",i,sep="")
    FuLiF<-c(FuLiF,x)
    }

  HapFst<-paste("Hap.Fst")
  nucFst<-paste("nuc.Fst")

  name<-c(S.s,pi,H.pop,TD,FuLiD,FuLiF,HapFst,nucFst)
  
  return(name)
}