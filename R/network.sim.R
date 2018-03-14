## get a node in a phylogeny
get.node<-function(tree){
  output<-list()
  tree<-strsplit(tree,"")

  for(i in 1:length(tree[[1]])){
    if(tree[[1]][i]==")"){
      brack1<-i
      break}
  }
  for(i in brack1:1){
    if(tree[[1]][i]=="("){
      brack2<-i
      break}
  }

  for(i in brack1:1){
    if(tree[[1]][i]==","){
      comma<-i
      break}
  }

  # for(i in brack1:length(tree[[1]])){
  # if(tree[[1]][i]==","){
  #    comma2<-i
  #    break}
  #}

  junction<-tree[[1]][brack2:brack1]
  junction<-paste(junction,collapse="")
  t<-paste(tree[[1]],collapse="")
  output[[1]]<-gsub(junction,paste(tree[[1]][(comma+1):(brack1-1)],collapse=""),t,fixed=T)
  junction<-gsub("(","",junction,fixed=T)
  junction<-gsub(")","",junction,fixed=T)
  junction<-gsub(","," ",junction,fixed=T)
  #join.time<-paste(tree[[1]][(brack1+2):(comma2-1)],collapse="")
  output[[2]]<-junction
  #output[[3]]<-join.time
  return(output)
}

## get rid of branch lengths
getrid<-function(x){
  z<-c(0:9,".",":")
  for(i in z){
    x<-gsub(i,"",x,fixed=T)
  }
  return(x)
}


### get divergence time in the same order as nodes above
get.time<-function(tree,tree2){
  output<-list()

  node.matrix<-read.tree.nodes(paste(tree2[[1]],collapse=""))

  tree<-strsplit(tree,"")

  for(i in 1:length(tree[[1]])){

    if(tree[[1]][i]==")"){
      brack1<-i
      break}
  }
  for(i in brack1:1){
    if(tree[[1]][i]=="("){
      brack2<-i
      break}
  }

  for(i in brack1:1){
    if(tree[[1]][i]==","){
      comma<-i
      break}
  }

  for(i in comma:1){
    if(tree[[1]][i]==":"){
      comma2<-i
      break}
  }

  junction<-tree[[1]][brack2:brack1]
  junction<-paste(junction,collapse="")
  t<-paste(tree[[1]],collapse="")
  output[[1]]<-gsub(junction,paste(tree[[1]][(comma+1):(brack1-1)],collapse=""),t,fixed=T)
  junction<-gsub("(","",junction,fixed=T)
  junction<-gsub(")","",junction,fixed=T)
  nodes<-strsplit(junction,",")[[1]]
  junction<-gsub(","," ",junction,fixed=T)
  output[[2]]<-junction

  output[[3]]<-coaltime(which(node.matrix$names==getrid(nodes[1])),which(node.matrix$names==getrid(nodes[2])),
                        nodematrix = node.matrix$nodes,nspecies=length(node.matrix$names))

  return(output)
}


## change taxon names to numbers
taxons.to.numbers<-function(tree,namelist=F){

  if(namelist[1]==F){
    x<-tree
    x<-gsub(")","",x,fixed=T)
    x<-gsub("(","",x,fixed=T)
    x<-gsub("#H","",x,fixed=T)
    x<-gsub(";","",x,fixed=T)
    x<-strsplit(x,",",fixed=T)[[1]]
    x<-x[x!=""]
  } else { x<-namelist}

  for(i in 1:length(x)){
    tree<-gsub(x[i],i,tree)
  }
  output<-list(x,tree)
  return(output)
}


#'  Get tree information
#' @description Get information of tip names, associated numbers and node heights.
#' @param tree The bifurcating tree topology in newick format.
#' @return A list containing the tip names, the newick tree with numbers and a matrix with the node hwighs.
#' @export
get.tree.info<-function(tree){

  tree.b<-getrid(tree)
  tree.b<-taxons.to.numbers(tree.b)

  ### get all joints of the tree
  input<-tree.b[[2]]
  join<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.node(input)
    join<-rbind(join,xx[[2]])
    input<-xx[[1]]
  }

  ### get all joint times of the tree
  input<-tree
  tree2<-tree
  join.t<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.time(input,tree2)
    join.t<-rbind(join.t,xx[[3]])
    input<-xx[[1]]
  }
  join<-cbind(join,join.t)
  tree.b[[3]]<-join
  return(tree.b)
}


#'  Simulation of phylogenetic networks under the coalescent
#' @description Simulation of bifurcating and nonbifurcating species tree
#' @param tree The bifurcating tree topology in newick format.
#' @param Ne.prior A vector of two values representing the min and max boundaries of a uniform prior for Ne.
#' @param bifurcating Logical. If TRUE the bifurcating topology is simulated.
#' @param migration Logical. If TRUE gene flow between two branches is allowed. Migrarion prior need to be specified in the 'mig' argument. hib.clade, hib.priors, major.sister and minor.sister need to be specified.
#' @param admixture Logical. If TRUE an admixture event between two branches is allowed. hib.clade, hib.priors, major.sister and minor.sister need to be specified.
#' @param mig A vector of two values representing the min and max boundaries of a uniform prior for gene flow in number of migrant copies (4Nm).
#' @param hib.clade a vector of numbers indicating the hybrid clade. Number associated with terminals can be checked with the get.tree.info function.
#' @param hib.priors a vector of 5 numbers representing the lower, upper boundaries of the hybridization time; the upper boundary for the connection with the minor sister; lower and upper boundaries of the admixture proportion (between 0-1).
#' @param major.sister vector of numbers indicating the the major sister clade. Number associated with terminals can be checked with the get.tree.info function.
#' @param minor.sister vector of numbers indicating the minor sister clade. Number associated with terminals can be checked with the get.tree.info function.
#' @param bp a vector of two numbers indicating the mean and SD of base pairs across all loci.
#' @param mi a vector of two numbers indicating the min and max mutation rate across all loci.
#' @param nsims Total number of simulations.
#' @param nloci Number of loci to be simulated in each iteration.
#' @param gen.time Generation time in years.
#' @param time.modif A time modifier to alter the age of the nodes in the newick tree. This is a uniform prior for node heights and it is provided as a multiplyer. A vector of two numbers, min and max, need to be specified. For instace, if you like to simulate node heights that are min 1/2 the heights of the input tree and max 2x the heights of the input tree you should use: c(0.5,2)
#' @param time.scalar multiplier to scale the three heights in the newick tree to years. For instance, if the node heights in the input tree are in Mya, the time.scalar should be 1000000. If time is in years, the time.scalar should be 1.
#' @export
sim.sp.tree<-function(tree,
                      Ne.prior,
                      bifurcating=T,
                      migration=F,
                      admixture=F,
                      mig,
                      hib.clade,
                      hib.priors,
                      major.sister,
                      minor.sister,
                      bp,
                      mi,
                      nsims,
                      nloci,
                      gen.time,
                      time.modif,
                      time.scalar=1)
{

  #### exclude branch lengths and change taxon to numbers
  tree.b<-getrid(tree)
  tree.b<-taxons.to.numbers(tree.b)

  ### get all joints of the tree
  input<-tree.b[[2]]
  join<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.node(input)
    join<-rbind(join,xx[[2]])
    input<-xx[[1]]
  }

  ### get all joint times of the tree
  input<-tree
  tree2<-tree
  join.t<-NULL
  while(length(grep("(",input,fixed=T))>=1){
    xx<-get.time(input,tree2)
    join.t<-rbind(join.t,xx[[3]])
    input<-xx[[1]]
  }

  ####### generate ej (nodes) flag string. Nodes ages are rescaled.
  ej<-cbind(join.t,join)
  ej<-cbind(rep("-ej",nrow(ej)),ej)
  ej[,2]<-as.numeric(ej[,2])*time.scalar

  ###### generate en (ancestral Ne) flag. Sample random Nes. Ne values will change in the simulation cicle
  en<-cbind(rep("-en",nrow(ej)),ej[,2],ej[,3],runif(nrow(ej),Ne.prior[1],Ne.prior[2]))
  for(i in 1:nrow(en)){
    en[i,3]<-strsplit(ej[i,3]," ")[[1]][2]
  }

  #######
  ms.string<-list()
  simulated<-NULL
  ## master Ne
  master.Ne<-mean(Ne.prior)
  ## get nodes
  nodes<-as.numeric(ej[,2])

  #######################
  #######################
  # simulations #########
  for(j in 1:nsims){
    sim.t<-NULL
    trees<-list()

    # put the original dates back on time strings
    ej[,2]<-nodes
    en[,2]<-nodes
    # rescale to coalescent (Ne proportion). Rescale to number of generations instead of years
    ej[,2]<-as.numeric(ej[,2])/(4*master.Ne)/gen.time
    en[,2]<-as.numeric(ej[,2])
    # sample time modifier
    time.mod<-runif(1,time.modif[1],time.modif[2])
    # modify time according to sampled time.mod
    ej[,2]<-as.numeric(ej[,2])*time.mod
    en[,2]<-as.numeric(en[,2])*time.mod

    # sample an Ne mean
    Ne.mean<-runif(1,Ne.prior[1],Ne.prior[2])
    # sample Ne Standart deviation
    Ne.SD<-Ne.mean*runif(1,0.1,1)
    # sample contemporary Nes (truncated to min 100 individuals)
    Nes<-rtnorm((nrow(ej)+1),Ne.mean,Ne.SD,lower=100)
    # Sample ancestral Nes (truncated to min 100 individuals)
    anc.Nes<-rtnorm(nrow(ej),Ne.mean,Ne.SD,lower=100)
    # Sample migration rate
    if(migration==T){
      Mig.rate<-runif(1,mig[1],mig[2])
      }else{
        Mig.rate<-0
      }
    # Sample admixture proportions
    if(admixture==T){
      Admix.prob.minor<-runif(1,hib.priors[4],hib.priors[5])
    }else{
      Admix.prob.minor<-0
    }

    mi.mean<-runif(1,mi[1],mi[2])
    mi.SD<-mi.mean*runif(1,0.1,1)

    ### simulate n-loci #################################################################################
    for(i in 1:nloci) {
      master.theta<-0
      while(master.theta<0.000001){
      # sample mutation rate per site per year
      rate<-rnorm(1,mi.mean,mi.SD)
      while(rate<=0){
        rate<-rnorm(1,mi.mean,mi.SD)
      }
      # sample sequence length
      seq.length<-rnorm(1,bp[1],bp[2])
      # generate theta
      master.theta<-master.Ne*4*seq.length*(rate*gen.time)
      }
      # theta string
      ms.string[[1]]<-paste("-t",master.theta)
      ### pop structure
      ms.string[[2]]<-paste(c("-I",(nrow(ej)+1),(rep(1,nrow(ej)+1))),collapse=" ")
      ### current popsize strings
      n<-cbind(rep("-n",(nrow(ej)+1)),c(1:(nrow(ej)+1)),rep(0,(nrow(ej)+1)),Nes)
      n[,3]<-as.numeric(n[,4])/master.Ne
      ms.string[[3]]<-paste(apply(n[,1:3],1,paste,collapse=" "),collapse=" ")
      ### ancestral pop sizes string
      en[,4]<-anc.Nes
      en[,4]<-as.numeric(en[,4])/master.Ne
      ms.string[[4]]<-paste(apply(en,1,paste,collapse=" "),collapse=" ")
      ###  node string
      ms.string[[5]]<-paste(apply(ej,1,paste,collapse=" "),collapse=" ")


      #### network string

      if(bifurcating==F){
        if(admixture==T){
        split.time<-runif(1,hib.priors[1],hib.priors[2])*time.scalar
        ej.hib<-runif(1,split.time,(hib.priors[2]*time.scalar))
        split.time<-((split.time/gen.time)*time.mod)/(4*master.Ne)
        ej.hib<-((ej.hib/gen.time)*time.mod)/(4*master.Ne)
        es<-paste("-es",split.time,max(hib.clade),(1-Admix.prob.minor))
        ejh<-paste("-ej",ej.hib,(length(tree.b[[1]])+1),max(minor.sister))
        ms.string[[6]]<-paste(es,ejh)
      }
      if(migration==T){
        mig.time<-((hib.priors[1]/gen.time)*time.mod)/(4*master.Ne)
        em<-paste("-em",mig.time,max(hib.clade),max(minor.sister),Mig.rate)
        ms.string[[6]]<-em
      }
      }

      ######
      # combining all string peaces in one ms string
      ms.string.final<-paste(unlist(ms.string),collapse=" ")
      #print(ms.string.final)
      ########################################
      ########################################
      ########################################
      #### simulated segregating sites
        fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)

        while(length(fas)==0){
          fas<-ms.to.DNAbin(ms(nreps = 1, nsam=(nrow(ej)+1),opts=ms.string.final),bp.length = 0)
        }

        d<-dist.dna(fas, model="N")/seq.length

        sim.t<-rbind(sim.t, as.vector(d))

      rm(fas)
    }

          if(nrow(sim.t)>1){
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-c(apply(sim.t,2,mean))#,apply(sim.t,2,var),apply(sim.t,2,kurtosis),apply(sim.t,2,skewness))
            mean<-paste("mean",names(sim[1:length(nam)]),sep="_")

            #var<-paste("var",names(sim[(length(nam)+1):(2*length(nam))]),sep="_")
            #kur<-paste("kur",names(sim[((2*length(nam))+1):(3*length(nam))]),sep="_")
            #skew<-paste("skew",names(sim[((3*length(nam))+1):(4*length(nam))]),sep="_")
            names(sim)<-mean#,var,kur,skew)
          } else {
            nam<-t(combn(attr(d,"Labels"),2))
            nam<-apply(nam,1,paste,collapse="_")
            colnames(sim.t)<-nam
            sim<-t(sim.t)
            }


    #print(i)
    sim<-cbind(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.rate,Admix.prob.minor,t(sim))
    simulated<-rbind(simulated,sim)
    rm(time.mod,Ne.mean,Ne.SD,mi.mean,mi.SD,Mig.rate,Admix.prob.minor,sim,sim.t)
    print(j)
    }
  return(simulated)
}


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
    #setTxtProgressBar(pb, i)
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


