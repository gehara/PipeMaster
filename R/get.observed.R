#' Observed pairwise distances between tips of tree
#' @description Calculation of paiwise distance between the tips of the tree.
#' @param tree The species tree topology in newick format.
#' @param path path to the directory containing fasta files to be included in the calculation.
#' @return Pairwise distances betwen tips of species tree.
#' @export
observed.pw.distances<-function(tree, path){

  moments=T
  setwd(path)
  x<-list.files()
  x<-x[grep(".fas",x)]
  tab.tab<-NULL
  basepairs<-NULL
  for(i in 1:length(x)){

    fas<-read.dna(x[i],format="fasta")
    fas<-fas[match(tree[[1]],rownames(fas)),]
    length<-ncol(fas)
    fas<-fas2ms2(fas)
    basepairs<-c(basepairs,length-fas[[2]])
    fas<-ms.to.DNAbin(fas[[1]],bp.length=length-fas[[2]])

    d<-dist.dna(fas, model="raw")
    sim.t<-as.vector(d)
    print(i)
  }
  if(moments==T){
    observed<-c(apply(tab.tab,2,mean))#,apply(tab.tab,2,var),apply(tab.tab,2,kurtosis),apply(tab.tab,2,skewness))
  } else {observed<-tab.tab}

  basepairs<-c(mean(basepairs),sd(basepairs))
  names(basepairs)<-c("basepairs mean and SD")
  return(list(t(data.frame(observed)),basepairs))
}
