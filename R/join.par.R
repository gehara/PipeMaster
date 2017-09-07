
join.par<-function()
  {

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
      join.par<-paste("join_",.e$joints[i],sep="")
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

