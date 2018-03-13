get.NGSdata.structure<-function(model,path.to.fasta,pop.assign){

  setwd(path.to.fasta)
  pops<-pop.assign
  pops<-pops[with(pops, order(pops[,2])), ]
  n_pops<-unique(pops[,2])

  # get all fasta files in WD
  fasta<-list.files(pattern = ".fa")

  # all samples per pop in different files
  pops_samples<-list()
  for(i in n_pops){
    pops_samples[[i]]<-pops[which(pops[,2]==i),]
  }

  #### get the population structure for every locus
  base_pairs<-NULL
  pop_str<-NULL
  for(i in 1:length(fasta)){
    seq<-read.dna(fasta[i], format="fasta")

    bp<-nrow(seq)

    npop<-list()
    for(j in 1:length(pops_samples)){
      npop[[j]]<-0
      for(u in 1:nrow(pops_samples[[j]])){
        x<-grep(pops_samples[[j]][u,1],rownames(seq))
        if(length(x)>0){
          npop[[j]]<-npop[[j]]+length(x)
        }
      }
    }
    pop_str<-rbind(pop_str,unlist(npop))
    base_pairs<-c(base_pairs,bp)
    cat(paste(i,"loci"),"\n")
  }

  LOCI<-cbind(rep("rate",as.numeric(model$loci[,3])),
              base_pairs,
              rep(1,as.numeric(model$loci[,3])),
              rep(model$loci[,4],as.numeric(model$loci[,3])),
              rep(model$loci[,5],as.numeric(model$loci[,3])),
              rep(model$loci[,6],as.numeric(model$loci[,3])))
  colnames(LOCI)<-NULL

  I<-cbind(paste("locus",1:as.numeric(model$loci[,3]),sep=""),
           rep("-I",as.numeric(model$loci[,3])),
           rep(model$I[,3],as.numeric(model$loci[,3])),
           pop_str)

  model$loci<-LOCI
  model$I<-I

  return(model)
}
