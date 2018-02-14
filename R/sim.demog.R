#' @export
sim.demog<-function(nsims,
                    coexp.prior,
                    Ne.prior,
                    NeA.prior,
                    time.prior,
                    gene.prior,
                    alpha=alpha,
                    append.sims=F,
                    path=getwd())
{

  setwd(path)


  if(append.sims==F){
    simulations<-matrix(nrow=1,ncol=7)
    simulations[1,]<-c("pi","ss","H","TD","ss1","ss2","ss3")
    write.table(simulations,file=paste(Ne.prior[,1],"demog_sim.txt",sep=""), quote=F,row.names=F, col.names=F, sep="\t")
    populations.par<-matrix(c("Ne","Exp.time","NeA","mi"),nrow=1,ncol=4)
    write.table(populations.par,file=paste(Ne.prior[,1],"pop_parameters.txt",sep=""), quote=F,row.names=F, col.names=F, sep="\t")
  }


  TIME<-system.time(for (i in 1:nsims){

    x<-demog.sample.pars(nruns=1,coexp.prior=coexp.prior,Ne.prior=Ne.prior,
                         NeA.prior=NeA.prior,time.prior=time.prior,gene.prior=gene.prior)

    y<-coexp.MS(MS.par=x$MS.par, gene.prior = gene.prior,alpha=alpha)

    z<-sumstat(ms.output=y,gene.prior=gene.prior)

    populations.par<-unlist(x$pop.par)

    write.table(z,file=paste(Ne.prior[,1],"demog_sim.txt",sep=""), quote=F,row.names=F, col.names=F, append=T, sep="\t")
    write.table(t(populations.par),file=paste(Ne.prior[,1],"pop_parameters.txt",sep=""), quote=F,row.names=F, col.names=F, append=T,sep="\t")
    print(paste(i,"sims of",nsims,"| single pop demogragphic test"))
   })
  print(TIME)
}
