#' Read the observed data to get the simulation parameters; base pairs and number of individuals per population.
#' @param model A model object generated with main.menu function.
#' @param path.to.fasta path to the fasta alignments directory.
#' @param pop.assign a two column data frame with poppulation assignemt of individuals. First column must have the samples name, the secound column the population in numbers, 1 for pop1, 2 for pop2 etc.
#' @param sanger logical. If TRUE the inheritance scalar and mutation rates set up in the main.menu are kept. If FALSE the mutation rate is propagated to all loci. Defaut is FALSE.
#' @author Marcelo Gehara
#' @return Model object with updated gene parameters.
#' @export
get.data.structure <- function(model, path.to.fasta, pop.assign, sanger=F)
  {
  orig.path <- getwd()
  setwd(path.to.fasta)
  pops <- pop.assign
  pops <- pops[with(pops, order(pops[,2])), ]
  n_pops <- unique(pops[,2])

  # get all fasta files in WD
  fasta <- list.files(pattern = ".fa")
  fasta <- fasta[grep(".fa",fasta,fixed=T)]

  # all samples per pop in different files
  pops_samples <- list()
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

  LOCI<-cbind(paste("locus",1:length(base_pairs),sep=""),
              base_pairs,
              rep(1,length(base_pairs)),
              rep(model$loci[1,4], length(base_pairs)),
              rep(model$loci[1,5], length(base_pairs)),
              rep(model$loci[1,6], length(base_pairs)))
  colnames(LOCI) <- NULL

  I <- cbind(paste("locus", 1:length(base_pairs), sep=""),
           rep("-I",length(base_pairs)),
           rep(model$I[1,3],length(base_pairs)),
           pop_str)

  if(sanger==T){
    model$loci[,2] <- LOCI[,2]
    model$I[,4:5] <- I[,4:5]
  } else {
  model$loci <- LOCI
  model$I <- I
  }
  setwd(orig.path)
  return(model)
}
