#' Observed pairwise distances between tips of tree
#' @description Calculation of paiwise distance between the tips of the tree.
#' @param tree The species tree topology in newick format.
#' @param path path to the directory containing fasta files to be included in the calculation.
#' @return Pairwise distances betwen tips of species tree.
#' @export
observed.pw.distances<-function(tree, path){
  tree<-get.tree.info(tree)
  setwd(path)
  x<-list.files()
  x<-x[grep(".fas",x)]

  print("Calculating distances. It may take a while.")
  tab.tab <- foreach(i = 1:length(x), .combine="rbind", .verbose=F) %dopar% {

    fas<-read.dna(x[i],format="fasta")
    fas<-fas[match(tree[[1]],rownames(fas)),]
    length<-ncol(fas)
    fas<-fas2ms2(fas)
    fas<-ms.to.DNAbin(fas[[1]],bp.length=length-fas[[2]])

    d<-dist.dna(fas, model="raw")
    sim.t<-c(as.vector(d),length)
    setTxtProgressBar(pb, i)
    #print(i)
    return(sim.t)
  }
  print("Done!")

  nam<-t(combn(attr(d,"Labels"),2))
  nam<-apply(nam,1,paste,collapse="_")
  nam<-c(nam,"bp_lengh")
  mean_dist<-colMeans(tab.tab,na.rm=T)#,apply(tab.tab,2,var),apply(tab.tab,2,kurtosis),apply(tab.tab,2,skewness))
  SD_dist<-apply(tab.tab,2,sd)
  names(mean_dist)<-paste(nam,"mean",sep="_")
  names(SD_dist)<-paste(nam,"SD",sep="_")
  return(list(mean_dist,SD_dist))
}
