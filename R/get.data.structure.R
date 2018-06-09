#' @export
get.data.structure<-function(model,path.to.fasta,pop.assign){

  setwd(path.to.fasta)
  pops<-pop.assign
  pops<-pops[with(pops, order(pops[,2])), ]
  n_pops<-unique(pops[,2])

  # get all fasta files in WD
  fasta<-list.files(pattern = ".fa")
  fasta<-fasta[grep(".fa",fasta,fixed=T)]

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

    bp<-ncol(seq)

    npop<-list()
    for(j in 1:length(pops_samples)){
      npop[[j]]<-length(na.omit(match(rownames(seq),as.character(pops_samples[[j]][,1]))))
      }
    if(sum(unlist(npop))!=nrow(seq)){
      stop(paste("one or more samples in locus",fasta[i],"have no assignment in your pop.assign file"))
    }

    pop_str<-rbind(pop_str,unlist(npop))
    base_pairs<-c(base_pairs,bp)
    cat(paste(i,"loci"),"\n")
  }
  if(model$I[1,1]=="genomic"){
  LOCI<-cbind(rep("rate",as.numeric(model$loci[,3])),
              base_pairs,
              rep(1,as.numeric(model$loci[,3])),
              rep(model$loci[,4], as.numeric(model$loci[,3])),
              rep(model$loci[,5], as.numeric(model$loci[,3])),
              rep(model$loci[,6], as.numeric(model$loci[,3])))
  colnames(LOCI)<-NULL

  I<-cbind(paste("locus",1:as.numeric(model$loci[,3]),sep=""),
           rep("-I",as.numeric(model$loci[,3])),
           rep(model$I[,3],as.numeric(model$loci[,3])),
           pop_str)

  model$loci<-LOCI
  model$I<-I
  } else {
    model$loci[,2] <- base_pairs
    model$I[,4:5] <- pop_str
    }
  return(model)
}
