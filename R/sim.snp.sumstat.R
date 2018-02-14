#' Simulate summary statistics using msABC
#'
#' @param model A model object bult by the main.menu function.
#' @param use.alpha Logical.If TRUE the most recent population size change will be exponential. If FALSE sudden demographic changes. Default is FALSE.
#' @param nsim.blocks Number of blocks to simulate. The total number of simulations is: nsim.blocks x sim.block.size.
#' @param sim.block.size Simulations are performed in blocks. This argument defines the size of the block in number of simulations, i.e. how many simulations to run per block.
#'                       A block of 1000 will work for most cases. Increse the total number of simulations with nsim.block argument.
#' @param path Path to write the output. By default outputs will be saved in the working directory.
#' @param output.name String. The prefix of the output names. Defalt is "model"
#' @param append.sims Logical. If TRUE simulations will be appended in the last output. Default is FALSE.
#' @param get.moments Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. Defalt is FALSE.
#' @param msABC.call String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.
#'                   There is no need to change that unless you want to compile the program yourself and point the function to it.
#' @return Writes simulations and parameters to the path directory.
#' @references Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338.
#' @references Pavlidis P., Laurent S., & Stephan W. (2010) msABC: A modification of Hudson’s ms to facilitate multi-locus ABC analysis. Molecular Ecology Resources, 10, 723–727.
#' @author Marcelo Gehara
#'

sim.snp.sumstat<-function(model,nsim.blocks,path=getwd(),use.alpha=F,moments=F,
                          append.sims=F,block.size=100, msABC.call=get.msABC(),output.name){

  # set working directory
  setwd(path)

  if(append.sims==F){
  com<-ms.commander.snp(model,use.alpha=use.alpha,
                        msABC=msABC.call)

  x<-strsplit(system(com[[1]],intern=T),"\t")
  nam<-x[1][[1]]
  TD_denom<-paste(nam[grep("pi",nam)],nam[grep("theta_w",nam)],sep="_")
  nam<-nam[-grep("ZnS",nam)]
  nam<-nam[-grep("thomson",nam)]
  nam<-c(nam, TD_denom)



  write.table(t(paste(nam,"_mean",sep="")),file=paste(output.name,"_mean.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  write.table(t(com[[3]][1,]),file=paste(output.name,"_par.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  if(moments==T){
    write.table(t(paste(nam,"_kur",sep="")),file=paste(output.name,"_kur.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
    write.table(t(paste(nam,"_skew",sep="")),file=paste(output.name,"_skew.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
    write.table(t(paste(nam,"_var",sep="")),file=paste(output.name,"_var.txt",sep=""),quote=F,row.names = F,col.names = F, append=F,sep="\t")
  }
  }

  thou<-0
  for(j in 1:nsim.blocks){
  SS<-NULL
  param<-NULL

  if(moments==T){
    VAR<-NULL
    KUR<-NULL
    SKEW<-NULL
  }

  for(i in 1:block.size){
    com<-ms.commander.snp(model,msABC=msABC.call, use.alpha = use.alpha)

    system(paste(com[[1]],"> out.txt"))
    S1<-read.table("out.txt",header = T)
    if(model$loci[,1]=="0"){S2<-NULL
    }else{
    system(paste(com[[2]],"> out.txt"))
    S2<-read.table("out.txt",header = T)}
    sumstat<-rbind(S1,S2)
    TD_denom<-data.frame(sumstat[,grep("pi",colnames(sumstat))]-sumstat[,grep("theta_w",colnames(sumstat))])
    colnames(TD_denom)<-paste(colnames(sumstat)[grep("pi",colnames(sumstat))],
                              colnames(sumstat)[grep("theta_w",colnames(sumstat))],sep="_")
    sumstat<-sumstat[,-grep("ZnS",colnames(sumstat))]
    sumstat<-sumstat[,-grep("thomson",colnames(sumstat))]
    sumstat<-cbind(sumstat, TD_denom)

    mean<-colMeans(sumstat,na.rm = T)
    if(moments==T){
      kur<-apply(sumstat,2,kurtosis, na.rm=T)
      skew<-apply(sumstat,2,skewness, na.rm=T)
      var<-apply(sumstat,2,var, na.rm=T)
      VAR<-rbind(VAR,var)
      SKEW<-rbind(SKEW,skew)
      KUR<-rbind(KUR,kur)
    }
    param<-rbind(param,com[[3]][2,])
    SS<-rbind(SS,mean)

    print(thou+i)

    }
  thou<-thou+block.size
  write.table(SS,file=paste(output.name,"_mean.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  write.table(param,file=paste(output.name,"_par.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  if(moments==T){
    write.table(KUR,file=paste(output.name,"_kur.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
    write.table(SKEW,file=paste(output.name,"_skew.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
    write.table(VAR,file=paste(output.name,"_var.txt",sep=""),quote=F,row.names = F,col.names = F, append=T,sep="\t")
  }

  }
}
