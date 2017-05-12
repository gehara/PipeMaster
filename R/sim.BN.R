sim.BN<-function(nsims,
                    Ne.prior,
                    NeA.prior,
                    time.prior,
                    gene.prior,
                    append.sims=F,
                    path=getwd())
{

  setwd(path)
  
  
  if(append.sims==F){
    simulations<-matrix(nrow=1,ncol=8)
    simulations[1,]<-c("Ne","NeA","tNeA","mi",
                       "pi","nH","H","TD")
    
    write.table(simulations,file="simulations.txt", quote=F,row.names=F, col.names=F, sep="\t")
    
  }
  
  
  TIME<-system.time(for (i in 1:nsims){
    
    x<-BN.sample.pars(nruns=1,Ne.prior=Ne.prior,
                         NeA.prior=NeA.prior,time.prior=time.prior,gene.prior=gene.prior)
    
    y<-BN.MS(MS.par=x$MS.par, gene.prior = gene.prior)
    
    z<-BN.sumstat(ms.output=y,gene.prior=gene.prior)
    
    simulations<-c(x$samp.par,z)
    
    write.table(t(simulations),file="simulations.txt", quote=F,row.names=F, col.names=F, append=T, sep="\t")
    print(paste(i,"sims of",nsims,"| Ne = ",x$coexp.par[,1]))
    
  })
  
  print(TIME)
}
