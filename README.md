#PipeMaster

Pipemaster can help you build coalescent models, add prior information on model parameters and simulate data sampling parameter values from those priors.

!!This is a beta version of the package. I have not tested very complex models with more than 4 populations.!!

You can use PipeMaster to simulate summay statistics and coalescent trees. You can also calculate the same summary statistics on your empirical data. You can then perform an abc analysis using the "abc" R-package or use a machine learning algorithm to do model and/or parameter inference.

Instaling the package:
dowload the .zip file from this page unzip it, open R and tipe the following code adding the path to where the package was downloaded.

> install.packages("path to PipeMaster", repos=NULL)

Using the menu to build your model:
You start by setting up your model throgh the menu, you run the main.menu() function and directs the output to an R oject that will store your model.

> model<-main.menu()

You can use a previous setup model as a template:

> model2<-main.menu(model1)

The function will ask you to tipe a topology in newick format and the menu will pop up with the following options.



R version 3.3.0 (2016-05-03) -- "Supposedly Educational"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> library(devtools)
> library(help=devtools)
> :create
Error: unexpected ':' in ":"
> ?create
> getwd()
[1] "/home/marcelo/Github/codes"
> list.files
function (path = ".", pattern = NULL, all.files = FALSE, full.names = FALSE, 
    recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, 
    no.. = FALSE) 
.Internal(list.files(path, pattern, all.files, full.names, recursive, 
    ignore.case, include.dirs, no..))
<bytecode: 0x37612e0>
<environment: namespace:base>
> list.files()
[1] "DESCRIPTION"      "man"              "NAMESPACE"        "PipeMaster.Rproj" "R"               
> create(path="../codes")
Error: Directory exists and is not empty
> build
function (pkg = ".", path = NULL, binary = FALSE, vignettes = TRUE, 
    manual = FALSE, args = NULL, quiet = FALSE) 
{
    pkg <- as.package(pkg)
    if (is.null(path)) {
        path <- dirname(pkg$path)
    }
    check_build_tools(pkg)
    compile_rcpp_attributes(pkg)
    if (binary) {
        args <- c("--build", args)
        cmd <- paste0("CMD INSTALL ", shQuote(pkg$path), " ", 
            paste0(args, collapse = " "))
        if (.Platform$OS.type == "windows") {
            ext <- ".zip"
        }
        else if (grepl("darwin", R.version$os)) {
            ext <- ".tgz"
        }
        else {
            ext <- paste0("_R_", Sys.getenv("R_PLATFORM"), ".tar.gz")
        }
    }
    else {
        args <- c(args, "--no-resave-data")
        if (manual && !has_latex(verbose = TRUE)) {
            manual <- FALSE
        }
        if (!manual) {
            args <- c(args, "--no-manual")
        }
        if (!vignettes) {
            args <- c(args, "--no-build-vignettes")
        }
        cmd <- paste0("CMD build ", shQuote(pkg$path), " ", paste0(args, 
            collapse = " "))
        ext <- ".tar.gz"
    }
    withr::with_temp_libpaths(R(cmd, path, quiet = quiet))
    targz <- paste0(pkg$package, "_", pkg$version, ext)
    file.path(path, targz)
}
<environment: namespace:devtools>
> ?bbuild
No documentation for ‘bbuild’ in specified packages and libraries:
you could try ‘??bbuild’
> ?build
> build()
'/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore --quiet CMD build '/home/marcelo/Github/codes'  \
  --no-resave-data --no-manual 

* checking for file ‘/home/marcelo/Github/codes/DESCRIPTION’ ... OK
* preparing ‘PipeMaster’:
* checking DESCRIPTION meta-information ... OK
* excluding invalid files
Subdirectory 'R' contains invalid file names:
  ‘.RData’ ‘README.md’
* checking for LF line-endings in source and make files
* checking for empty or unneeded directories
* building ‘PipeMaster_0.1.tar.gz’

[1] "/home/marcelo/Github/PipeMaster_0.1.tar.gz"

> main.menu()  
write bifurcating topology in newick format: (1,2);  
A > Number of populations to simulate      2  
B > Tree                                   (1,2);  
    Number of nodes (junctions)            1  
C > Migration                              FALSE  
D > Pop size change through time           FALSE  
E > Setup Ne priors  
    Population size parameters (total)     2  
    current Ne parameters                  2  
    ancestral Ne parameters                0  
F > Setup Migration priors  
    current migration                      0  
    ancestral migration parameters         0  
G > Setup time priors   
    time of join parameters                1  
    time of Ne change parameters           0  
    time of Migration change parameters    0  
H > Conditions  
I > Gene setup  
Q > Quit, my model is Ready!  

This will give you a simple isolation model. You can add complexity to the model by using the interactive menu.
You must go to the "H > Conditions" and the "I > Gene setup" before quitting.

Simulating summary stats:

To simulate summary statistics you run the following function.

> sim.sumstat(model,path="",overall.SS=F,perpop.SS=T,nsim.blocks=1,use.alpha = T, append.sims=F, sim.block.size = 1000)

model = your model object.  
path = path to write the output.  
overall.SS = if TRUE calculates the summary stats across all your populations. Defaut is FALSE.  
perpop.SS = if TRUE calculates the summary stats per population. Defalt is TRUE.  
use.alpha = if TRUE the most recent population size chage will be exponential. Sudden change if FALSE. Defaut is FALSE.  
sim.block.size = simulations are perforemed in blocks. This argument the fine the size of the block. How many simmulations per block.  
nsim.blocks = that setup how many blocks you what to simulate. The total number of simulations is: nsim.blocks x sim.block.size.  
append.sims = if TRUE simmulations will be appended in the last simulation output. Defaut is FALSE.  

To simmulate coalescent trees:

> sim.coaltrees(model,path="",use.alpha=F,nsim.blocks=1, append.sims=F, sim.block.size = 1000)


