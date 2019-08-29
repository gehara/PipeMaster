#' Model Builder
#' @description This function starts the menu for model building.
#' @param  input A model object to be used as template. Can be left blank, defalt NULL.
#' @param  ms.string a character string describing the model. This can be obtained by running
#' @return An object representing a diversification model.
#' @examples
#' my.model<-main.menu()
#' sim.sumstat(my.model)
#
#' @export
main.menu<-function(input = NULL, ms.string = NULL){


  if(is.null(ms.string)==F){
    read.ms.string(ms.string = ms.string)
  } else {

  if(is.null(input)==F){
    read.model.input(input)
  }
  if(!(exists("tree", envir=.e))) join.par()
  if(!(exists("n", envir=.e))) cur.Ne.par()
}
  print.main.menu()

    letter<<-readline("Model Builder >>>>")
    while(letter %in% c("A","B","C","D","E","F","G","H","I","Q","P")==F){
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
       paste("P > Plot Model"),
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

         H = {if(exists("size.matrix", envir=.e)){
                x<-readline("Would you like to overwrite the existing matrix?  ")
                if(x %in% .e$YES) condition.matrix()
                } else {
                condition.matrix()
              }
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

         P = {if(exists("size.matrix",envir=.e)){
         } else {sys.call(which = 0)
                   condition.matrix()}

           if(exists("m",envir=.e) & exists("mig.matrix",envir=.e)==F){
             sys.call(which = 0)
             condition.matrix()
           }

           if(exists("loci",envir=.e)){
           } else {print("you need to go to gene menu first!")
               sys.call(which = 0)
               main.menu()
           }

           alpha <- as.logical(readline("exponential size change (TRUE or FALSE)? "))
           if(alpha==T){
             .e$exp.pops <- readline("Indicate pop numbers separated by comma for exponential change ")
             .e$exp.pops <- as.numeric(strsplit(.e$exp.pops,",")[[1]])
             alpha <- c(T,.e$exp.pops)
           }
           model <- get.model()

          PipeMaster::PlotModel(model = model, use.alpha = alpha)

          sys.call(which = 0)
          main.menu(model)
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
      while(z!=y){print("a parenthesis is missing, please write the tree correctly!")
        .e$tree<-readline("write bifurcating topology in newick format: ")
        x<-strsplit(.e$tree,"")
        y<-as.numeric(length(grep("(",x[[1]],fixed=T)))
        z<-as.numeric(length(grep(")",x[[1]],fixed=T)))
        w<-as.numeric(length(grep(",",x[[1]],fixed=T)))
      }
    }
  }

  ## get topology and number of nodes
  .e$tree<-readline("write bifurcating topology in newick format or 1 for a single population >>> ")
  while(.e$tree==""){
    .e$tree<-readline("write bifurcating topology in newick format or 1 for a single population >>> ")
  }
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
### read ms string from popplanner
read.ms.string<-function(ms.string){

  ms.string <- strsplit(ms.string,"-")[[1]]

  if(length(ms.string[grep("g", ms.string)])!=0)
   stop("Exponential demographich change is not allowed in the ms.string option. It will be added in a letter stage in the simulation function.")

  if(length(ms.string[-grep("ms", ms.string)])!=0) ms.string <- ms.string[-grep("ms", ms.string)]

  .e$npops <- ms.string[grep("I", ms.string)]
  .e$npops <- as.numeric(strsplit(.e$npops," ")[[1]][2])

  ## get -n parameters
  list.Ne.pars<-NULL
  for (i in 1:.e$npops){
    Ne0.par<-paste("Ne0.pop",i,sep="")
    list.Ne.pars<-c(list.Ne.pars,Ne0.par)
  }
  .e$n<-matrix(nrow=length(list.Ne.pars), ncol=6)
  .e$n[,1]<-list.Ne.pars
  .e$n[,2]<-'-n'
  .e$n[,3]<-c(1:length(list.Ne.pars))
  .e$n[,4]<-100000
  .e$n[,5]<-500000
  .e$n[,6]<-"uniform"



  #### get -ej nodes
  if(length(ms.string[grep("ej", ms.string)])>0){
    nodes <- ms.string[grep("ej", ms.string)]
    nodes <- sapply(nodes, strsplit," ")
    .e$joints <- NULL
    for(i in 1:length(nodes)){
      .e$joints<-c(.e$joints, paste(nodes[[i]][3:4],collapse = " "))
    }

    .e$tree<-paste("ms string [",paste(.e$joints, collapse="]["),"]",sep="")

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

  # get ancestral pop sizes
  if(length(ms.string[grep("en", ms.string)])>0){
    .e$en<-ms.string[grep("en", ms.string)]

    if(length(.e$en)!=0){
      x<-NULL
      for(i in 1:length(.e$en)){
        x <- c(x,strsplit(.e$en[i]," ")[[1]][3])
      }
      .e$en<-x

      anc.Ne.par<-NULL
      time.anc.Ne.par<-NULL
      pop<-NULL
      for (i in 1:.e$npops){
        n.anc.pop <- length(grep(i,.e$en))
        if (n.anc.pop==0){
        } else (
          for (j in 1:n.anc.pop){
            Ne.par<-paste("Ne",j,".pop",i,sep="")
            anc.Ne.par<-c(anc.Ne.par,Ne.par)
            time.Ne.par<-paste("t.Ne",j,".pop",i,sep="")
            time.anc.Ne.par<-c(time.anc.Ne.par,time.Ne.par)
            pop<-c(pop,i)
          }
        )
      }

      .e$en$size<-matrix(nrow=length(anc.Ne.par),ncol=6)
      .e$en$size[,1]<-anc.Ne.par
      .e$en$size[,2]<-'-en'
      .e$en$size[,3]<-pop
      .e$en$size[,4]<-1000
      .e$en$size[,5]<-10000
      .e$en$size[,6]<-'uniform'

      .e$en$time<-matrix(nrow=length(anc.Ne.par),ncol=6)
      .e$en$time[,1]<-time.anc.Ne.par
      .e$en$time[,2]<-'-en'
      .e$en$time[,3]<-pop
      .e$en$time[,4]<-10000
      .e$en$time[,5]<-100000
      .e$en$time[,6]<-'uniform'
    }
  }

  #### get -m string
  if(length(ms.string[grep("m", ms.string)])>0){
    .e$m <- ms.string[grep("m", ms.string)]

    if(length(.e$m[grep("em", .e$m)])!=0) .e$m <-.e$m[-grep("em", .e$m)]

    mig.par<-NULL
    pops<-NULL
    for(i in 1:length(.e$m)){
      x <- strsplit(.e$m[i]," ")[[1]]
      m0<-paste("mig0.",x[2],"_",x[3],sep="")
      mig.par<-c(mig.par,m0)
      pops<-c(pops,paste(x[2],x[3]))
    }

    .e$m<-matrix(nrow=length(mig.par),ncol=6)
    .e$m[,1]<-mig.par
    .e$m[,2]<-"-m"
    .e$m[,3]<-pops
    .e$m[,4]<-0.1
    .e$m[,5]<-1
    .e$m[,6]<-'uniform'
  }
  ### get -em string
  if(length(ms.string[grep("em", ms.string)])>0){
    .e$em <- ms.string[grep("em", ms.string)]

    x<-NULL
    for(i in 1:length(.e$em)){
      x<-c(x,paste(strsplit(.e$em[i]," ")[[1]][3:4],collapse ="_"))
    }
    migs<-unique(x)

    mig.par<-NULL
    pops<-NULL
    for(i in 1:length(migs)){
      pop.migs<-grep(migs[i],x)

      for(j in 1:length(pop.migs)){
        mig.par<-c(mig.par,paste("mig",j,".",x[pop.migs[j]], sep=""))
        pops<-c(pops,x[pop.migs[j]])
      }

    }
    .e$em<-NULL
    .e$em$size<-matrix(nrow=length(mig.par),ncol=6)
    .e$em$size[,1]<-mig.par
    .e$em$size[,2]<-"-em"
    .e$em$size[,3]<-pops
    .e$em$size[,4]<-0
    .e$em$size[,5]<-0
    .e$em$size[,6]<-'uniform'


    t.mig.par<-mig.par
    for(i in 1:length(t.mig.par)){
      t.mig.par[i]<-paste("t.",mig.par[i],sep="")
    }

    .e$em$time<-matrix(nrow=length(t.mig.par),ncol=6)
    .e$em$time[,1]<-t.mig.par
    .e$em$time[,2]<-"-em"
    .e$em$time[,3]<-pops
    .e$em$time[,4]<-10000
    .e$em$time[,5]<-20000
    .e$em$time[,6]<-'uniform'

  }
  condition.matrix()

  if(length(ms.string[grep("ej", ms.string)])>0) if(nrow(.e$ej)>1) update.matrix(nodes=nodes)

}
# update join matrix conditions according to the ms string from popplanner
update.matrix<-function(nodes){


  joints <- NULL
  for(i in 1:length(nodes)){
    joints<-rbind(joints, c(nodes[[i]][3],nodes[[i]][4]))
  }

  x=NULL
  for(i in 1:.e$npops) if(length(grep(i,joints))>1) x<-c(x,i)


  for(i in 1:length(x)){
    y<-which(joints == x[i], arr.ind=TRUE)
    w<-nodes[y[,1]]
    for(j in 1:length(w)){
      w[[j]]<-as.numeric(w[[j]][2:4])
    }
    w <- matrix(unlist(w), ncol = 3, byrow = TRUE)
    w <- w[match(sort(w[,1]),w[,1]),]
    cb <- combn(1:nrow(w),2)
    for(u in 1: ncol(cb)){
      cond<-paste("join",w[cb[1,u],2],"_",w[cb[1,u],3]," < ","join",w[cb[2,u],2],"_",w[cb[2,u],3],sep="")
      cond<-strsplit(cond," ")
      yy<-grep(cond[[1]][1],rownames(.e$time.matrix))
      xx<-grep(cond[[1]][3],colnames(.e$time.matrix))
      .e$time.matrix[yy,xx]<-cond[[1]][2]
      if(yy>xx) .e$time.matrix<-inv.mirror.upper(.e$time.matrix)
      if(xx>yy) .e$time.matrix<-inv.mirror.lower(.e$time.matrix)
    }
  }


}


