#' Model Builder
#' @description This function starts the menu for model building.
#' @param  input A model object to be used as template. Can be left blank, defalt NULL.
#' @return An object representing a diversification model.
#' @examples
#' my.model<-main.menu()
#' sim.sumstat(my.model)
#
#' @export
main.menu<-function(input=NULL){

  if(is.null(input)==T){} else{
    read.model.input(input)
  }
  if(exists("tree", envir=.e)){
  }else{join.par()}
  if(exists("n", envir=.e)){
  }else{cur.Ne.par()}

  print.main.menu()

    letter<<-readline("Model Builder >>>>")
    while(letter %in% c("A","B","C","D","E","F","G","H","I","Q")==F){
     cat(paste("Choose a valid letter. You typed:",letter))
     letter<<-readline("Model Builder >>>>")
    }
  switch.main.menu()
  }
### hidden environment
.e<-new.env()
# internal function of the Model Builder
print.main.menu<-function(){

  {cat(paste("A > Number of populations to simulate     ",nrow(.e$n)),
       paste("B > Tree (go here to remove nodes)        ",.e$tree),
       paste("    Number of nodes (junctions)           ",nrow(.e$ej)),
       paste("C > Migration                             ",exists("m",envir=.e)),
       paste("D > Pop size change through time          ",exists("en",envir=.e)),
       paste("E > Setup Ne priors"),
       paste("    Population size parameters (total)    ",sum(nrow(.e$n),nrow(.e$en$size))),
       paste("    current Ne parameters                 ",nrow(.e$n)),
       paste("    ancestral Ne parameters               ",sum(nrow(.e$en$size))),
       paste("F > Setup Migration priors"),
       paste("    current migration                     ",sum(nrow(.e$m))),
       paste("    ancestral migration parameters        ",sum(nrow(.e$em$size))),
       paste("G > Setup time priors "),
       paste("    time of join parameters               ",sum(nrow(.e$ej))),
       paste("    time of Ne change parameters          ",sum(nrow(.e$en$time))),
       paste("    time of Migration change parameters   ",sum(nrow(.e$em$time))),
       paste("H > Conditions"),
       paste("I > Gene setup"),
       paste("Q > Quit, my model is Ready!"),
       sep="\n")
  }
}
# internal function of the Model Builder
switch.main.menu<-function(){
  .e$YES<-c("Y","y","yes","YES","Yes")
  .e$NO<-c("N","n","No","NO","no")

  switch(letter,

         A = {remove.all.par()
           join.par()
           cur.Ne.par()
           sys.call(which = 0)
           main.menu()},

         B = {.e$island.presence<-readline("would you like to remove a node (YES or NO)?: ")
         if(.e$island.presence %in% .e$YES){
           print(.e$ej)
           node<-readline("which node to remove (write all node rows separated by spaces)?: ")
           node<-as.numeric(strsplit(node," ")[[1]])
           if(length(node)==nrow(.e$ej)){
             .e$ej<-NULL
           } else {
             .e$ej<-.e$ej[,-node]
           }
           .e$tree<-"non tree-like model"
           mig.par()
         }
         sys.call(which = 0)
         main.menu()},

         C = {.e$mig.presence<-readline("Migration among populations (YES or NO)?: ")
         while(.e$mig.presence %in% c(.e$YES,.e$NO)==F){
           cat(paste("Type a valid letter. You typed:",.e$mig.presence))
           .e$mig.presence<-readline("Migration among populations (YES or NO)?: ")
         }
         if(.e$mig.presence %in% .e$YES){
           mig.par()
         } else if (.e$mig.presence %in% .e$NO){
           options(warn=-1)
           rm(m,em,envir=.e)
           options(warn=0)
         }
         sys.call(which = 0)
         main.menu()},


         D = {hist.demog<-readline("Ne change throgh time (YES or NO?): ")
         if (hist.demog %in% .e$YES){
           anc.Ne.par()
         } else if (hist.demog %in% .e$NO){
           options(warn=-1)
           rm(en,envir=.e)
           options(warn=0)
         }
         sys.call(which = 0)
         main.menu()},

         E = {sys.call(which = 0)
           demog.menu()},

         "F" = {if(exists("m",envir=.e)){
         }else{
           .e$mig.presence<-readline("Migration among populations (YES or NO)?: " )
           if(.e$mig.presence=="NO"){
             sys.call(which = 0)
             main.menu()
           } else {mig.par()}
         }
           sys.call(which = 0)
           mig.menu()},

         G = {time.menu()},

         H = {condition.matrix()
           sys.call(which = 0)
           condition.menu()},

         I = {.e$data.type<-readline("What type of data to simulate (sanger or genomic)?:" )

         if(.e$data.type %in% c("sanger","Sanger","S","s")){
           sys.call(which = 0)
           loci.par()
           samples.par()
           gene.menu()
         } else {
           sys.call(which = 0)
           genomic.par()
           genomic.samples.par()
           genomic.menu()
         }
         },



         Q={if(exists("size.matrix",envir=.e)){
         } else {condition.matrix()}
           if(exists("m",envir=.e) & exists("mig.matrix",envir=.e)==F){
             condition.matrix()
           }
           get.model()}

  )}
# internal function of the Model Builder
join.par<-function(){

  check.tree<-function(){
    x<-strsplit(.e$tree,"")
    y<-as.numeric(length(grep("(",x[[1]],fixed=T)))
    z<-as.numeric(length(grep(")",x[[1]],fixed=T)))
    w<-as.numeric(length(grep(",",x[[1]],fixed=T)))
    while(z!=y){print("a parenthesis is missing, write the tree correctly!")
      .e$tree<-readline("write bifurcating topology in newick format: ")
      x<-strsplit(.e$tree,"")
      y<-as.numeric(length(grep("(",x[[1]],fixed=T)))
      z<-as.numeric(length(grep(")",x[[1]],fixed=T)))
      w<-as.numeric(length(grep(",",x[[1]],fixed=T)))
    }
    while (z!=w){print("either a comma is missing or this is a nonbifurcating tree!")
      .e$tree<-readline("write bifurcating topology in newick format: ")
      x<-strsplit(.e$tree,"")
      y<-as.numeric(length(grep("(",x[[1]],fixed=T)))
      z<-as.numeric(length(grep(")",x[[1]],fixed=T)))
      w<-as.numeric(length(grep(",",x[[1]],fixed=T)))
      while(z!=y){print("a parenthesis is missing, write the tree correctly!")
        .e$tree<-readline("write bifurcating topology in newick format: ")
        x<-strsplit(.e$tree,"")
        y<-as.numeric(length(grep("(",x[[1]],fixed=T)))
        z<-as.numeric(length(grep(")",x[[1]],fixed=T)))
        w<-as.numeric(length(grep(",",x[[1]],fixed=T)))
      }
    }
  }

  ## get topology and number of nodes
  .e$tree<-readline("write bifurcating topology in newick format or 1 for single population: ")
  if (.e$tree=="1"){
    .e$npops<-1
    .e$ej<-NULL
  } else {
    check.tree()

    .e$npops<-as.numeric(nchar(gsub("(","",gsub(")","",gsub(",","",gsub(";","",.e$tree,fixed=T),fixed=T),fixed=T),fixed=T)))


    .e$t<-.e$tree
    .e$joints<-NULL

    get.joint<-function(){
      tree<-strsplit(.e$t,"")
      for(i in 1:length(tree[[1]])){
        if(tree[[1]][i]==")"){
          junction<-tree[[1]][(i-4):(i)]
          junction<-paste(junction,collapse="")
          t<-paste(tree[[1]],collapse="")
          t<-gsub(junction,tree[[1]][i-1],t, fixed=T)
          .e$t<-t
          junction<-gsub("(","",junction,fixed=T)
          junction<-gsub(")","",junction,fixed=T)
          junction<-gsub(","," ",junction,fixed=T)
          .e$joints<-c(.e$joints,junction)
          break
        }
      }
    }

    while(length(strsplit(.e$t,"")[[1]])>2){
      get.joint()
    }

    ## generate parameters of time of join of populations

    tot.join.par<-NULL
    for (i in 1:length(.e$joints)){
      join.par<-paste("join",paste(strsplit(.e$joints[i]," ")[[1]],collapse="_"),sep="")
      tot.join.par<-c(tot.join.par,join.par)
    }
    .e$ej<-matrix(nrow=length(.e$joints),ncol=6)
    .e$ej[,1]<-tot.join.par
    .e$ej[,2]<-'-ej'
    .e$ej[,3]<-.e$joints
    .e$ej[,4]<-500000
    .e$ej[,5]<-1500000
    .e$ej[,6]<-"uniform"

  }
}
# internal function of the Model Builder
get.model<-function(){


  .e$ej<-gsub("normal","rtnorm",.e$ej)
  .e$n<-gsub("normal","rtnorm",.e$n)
  .e$m<-gsub("normal","rtnorm",.e$m)
  .e$en$size<-gsub("normal","rtnorm",.e$en$size)
  .e$em$size<-gsub("normal","rtnorm",.e$em$size)
  .e$en$time<-gsub("normal","rtnorm",.e$en$time)
  .e$em$time<-gsub("normal","rtnorm",.e$em$time)
  .e$loci<-gsub("normal","rtnorm",.e$loci)

  .e$ej<-gsub("uniform","runif",.e$ej)
  .e$n<-gsub("uniform","runif",.e$n)
  .e$m<-gsub("uniform","runif",.e$m)
  .e$en$size<-gsub("uniform","runif",.e$en$size)
  .e$em$size<-gsub("uniform","runif",.e$em$size)
  .e$en$time<-gsub("uniform","runif",.e$en$time)
  .e$em$time<-gsub("uniform","runif",.e$em$time)
  .e$loci<-gsub("uniform","runif",.e$loci)

  if(is.null(nrow(.e$em$size))){rm("em",envir=.e)}
  if(is.null(nrow(.e$en$size))){rm("en",envir=.e)}
  if(is.null(nrow(.e$m))){rm("m",envir=.e)}
  if(is.null(nrow(.e$ej))){rm("ej",envir=.e)}

  model<-list(NULL,NULL,NULL,NULL,NULL)
  names(model)<-c("loci","I","flags","conds","tree")
  model$loci<-.e$loci
  model$I<-.e$I

  flags<-list(NULL,NULL,NULL,NULL,NULL)
  names(flags)<-c("n","m","en","em","ej")
  flags$n <- .e$n
  flags$m <- .e$m
  flags$en <- .e$en
  flags$em <- .e$em
  flags$ej <- .e$ej

  model$flags<-flags

  conds<-list(NULL,NULL,NULL)
  names(conds)<-c("size.matrix","mig.matrix","time.matrix")
  conds$size.matrix<-.e$size.matrix
  conds$mig.matrix<-.e$mig.matrix
  conds$time.matrix<-.e$time.matrix

  model$conds<-conds
  model$tree<-.e$tree

  rm(list=ls(envir=.e),envir=.e)
  class(model)<-"Model"
  return(model)

}
# internal function of the Model Builder
remove.all.par<-function(){
  options(warn=-1)
  rm(list=c("loci","I","n","en","m","em","ej","conds","tree","npops"), envir=.e)
  options(warn=0)
}

# internal function of the Model Builder
read.model.input<-function(input){

  .e$ej<-input$flags$ej
  .e$n<-input$flags$n
  .e$en<-input$flags$en
  .e$m<-input$flags$m
  .e$em<-input$flags$em

  .e$ej<-gsub("rtnorm","normal",.e$ej)
  .e$n<-gsub("rtnorm","normal",.e$n)
  .e$m<-gsub("rtnorm","normal",.e$m)
  .e$en$size<-gsub("rtnorm","normal",.e$en$size)
  .e$em$size<-gsub("rtnorm","normal",.e$em$size)
  .e$en$time<-gsub("rtnorm","normal",.e$en$time)
  .e$em$time<-gsub("rtnorm","normal",.e$em$time)
  .e$loci<-gsub("rtnorm","normal",.e$loci)

  .e$ej<-gsub("runif","uniform",.e$ej)
  .e$n<-gsub("runif","uniform",.e$n)
  .e$m<-gsub("runif","uniform",.e$m)
  .e$en$size<-gsub("runif","uniform",.e$en$size)
  .e$em$size<-gsub("runif","uniform",.e$em$size)
  .e$en$time<-gsub("runif","uniform",.e$en$time)
  .e$em$time<-gsub("runif","uniform",.e$em$time)
  .e$loci<-gsub("runif","uniform",.e$loci)

  if(is.null(nrow(.e$em$size))){rm("em",envir=.e)}
  if(is.null(nrow(.e$en$size))){rm("en",envir=.e)}
  if(is.null(nrow(.e$m))){rm("m",envir=.e)}
  if(is.null(nrow(.e$ej))){rm("ej",envir=.e)}

  .e$loci<-input$loci
  .e$I<-input$I

  .e$size.matrix<-input$conds$size.matrix
  .e$mig.matrix<-input$conds$mig.matrix
  .e$time.matrix<-input$conds$time.matrix
  if(is.null(nrow(.e$mig.matrix))){rm("mig.matrix",envir=.e)}
  .e$npops<-nrow(.e$n)
  .e$tree<-input$tree

}

# internal function of the Model Builder
samples.par<-function(){

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

