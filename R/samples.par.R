
samples.par<-function()
{

  tot.gene.par<-NULL
  for (i in 1:.e$ngenes){
    gene.par<-paste("locus",i,sep="")
    tot.gene.par<-c(tot.gene.par,gene.par)
  }
  .e$I<-matrix(nr=.e$ngenes,nc=3+.e$npops)
  .e$I[,1]<-tot.gene.par
  .e$I[,2]<-"-I"
  .e$I[,3]<-.e$npops

  for(j in 1:.e$ngenes){
    for(i in 1:.e$npops){
      .e$I[j,i+3]<-readline(paste("number of samples for pop",i,"locus",j,":"))
    }
  }
}
