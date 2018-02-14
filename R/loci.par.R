#' internal function of the Model Builder
loci.par<-function()
  {
    ## get topology and number of nodes
    .e$ngenes<-as.numeric(readline("how many loci to simulate?: "))

    tot.gene.par<-NULL
    for (i in 1:.e$ngenes){
      gene.par<-paste("rate",i,sep="")
      tot.gene.par<-c(tot.gene.par,gene.par)
    }
    .e$loci<-matrix(nrow=.e$ngenes,ncol=6)
    .e$loci[,1]<-tot.gene.par
    .e$loci[,2]<-500
    .e$loci[,3]<-1
    .e$loci[,4]<-5e-9
    .e$loci[,5]<-1.5e-8
    .e$loci[,6]<-"uniform"

}


