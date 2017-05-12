
print.main.menu<-function(){
    
{cat(paste("A > Number of populations to simulate     ",nrow(.e$n)),
     paste("B > Tree                                  ",.e$tree),
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

