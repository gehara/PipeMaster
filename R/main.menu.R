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

    letter<<-readline(">>>>")
    while(letter %in% c("A","B","C","D","E","F","G","H","I","Q")==F){
     cat(paste("Choose a valid letter. You typed:",letter))
     letter<<-readline(">>>>")
    }
  switch.main.menu()
  }

### hidden environment
.e<-new.env()

#' internal function of the Model Builder
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


#' internal function of the Model Builder
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




