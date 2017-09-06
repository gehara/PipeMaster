
genomic.par<-function()
  {
    ## get topology and number of nodes
    .e$ngenes<-1
    .e$nloci<-as.numeric(readline("how many loci to simulate?: "))
     bp<-as.numeric(readline("average number of base pairs: "))
    .e$loci<-matrix(nrow=1,ncol=6)
    .e$loci[,1]<-0
    .e$loci[,2]<-bp
    .e$loci[,3]<-.e$nloci
    .e$loci[,4]<-1e-9
    .e$loci[,5]<-1e-11
    .e$loci[,6]<-"uniform"

}
