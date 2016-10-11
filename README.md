#PipeMaster

Pipemaster can help you build coalescent models, add prior information on model parameters and simulate data sampling parameter values from those priors.

**This is a beta version of the package. I have not tested very complex models with more than 4 populations!**

You can use PipeMaster to simulate summary statistics and coalescent trees. You can also calculate the same summary statistics on your empirical data. You can then perform an abc analysis using the "abc" R-package or use a machine learning algorithm to do model and/or parameter inference.

####Installing the package

Install all dependencies.  

> install.packages("ape")  
> install.packages("pegas")  
> install.packages("phyclust")  
> install.packages("e1071")  
> install.packages("PopGenome")  

download the .zip file from this page and unzip it. Open R and type the following code adding the path to where the package was downloaded.
> install.packages("**_path_**/PipeMaster_0.0.3.tar.gz", repos=NULL)

####Using the menu to build your model
You start by setting up your model through the menu, you run the main.menu() function and directs the output to an R object that will store your model.

> model<-main.menu()

_You can use a previously build model as a template_:

> model2<-main.menu(model1)

The function will ask you to write a topology in newick format and the menu will pop up with the following options.

> main.menu()  
write bifurcating topology in newick format: (1,2);  
```
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
Q > Quit, my model is ready!
 ```

This will give you a simple isolation model. You can add complexity to the model by using the interactive menu.
You must go to the "_H > Conditions_" and the "_I > Gene setup_" before quitting.

####Simulating summary stats

To simulate summary statistics you run the following function.

> sim.sumstat(model,path="",overall.SS=F,perpop.SS=T,nsim.blocks=1,use.alpha = T, append.sims=F, sim.block.size = 1000)

* model = your model object.  
* path = path to write the output.  
* overall.SS = if TRUE calculates the summary stats across all your populations. Defaut is FALSE.  
* perpop.SS = if TRUE calculates the summary stats per population. Default is TRUE.  
* use.alpha = if TRUE the most recent population size change will be exponential. Sudden change if FALSE. Default is FALSE.  
* sim.block.size = simulations are performed in blocks. This argument defines the size of the block. How many simulations you what to run per block.
* nsim.blocks = that setup how many blocks you what to simulate. The total number of simulations is: nsim.blocks x sim.block.size.  
* append.sims = if TRUE simulations will be appended in the last simulation output. Default is FALSE.  

####Simulating coalescent trees

> sim.coaltrees(model,path="",use.alpha=F,nsim.blocks=1, append.sims=F, sim.block.size = 1000)


