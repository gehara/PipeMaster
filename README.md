# PipeMaster

PipeMaster is an R-package to build and simulate diversification models under the coalescent. Current implementation can simulate diversification models for single species or species complex diversification studies, hierarchical  demographic models for comparative studies and species trees with one horizontal connection.

**This is a beta version of the package!**

PipeMaster simulates summary statistics and coalescent trees. It calculates the same summary statistics on an empirical data. The user can use these sumary statistics to perform aproximate Bayesian computation (ABC) or supervized machine learning (SML) for model and parameter inference.

#### Installing the package  

> install.packages("devtools")  
> install_github("gehara/PipeMaster@developing")  

### Hierarchical codemographic model

PipeMaster can simulate a hierarchical demographic model for comparative analys of populations/species. The hABC method used in the package was first described in Chan et al 2014 and improved in Gehara et al 2017. These two papers show good examples of this analysis. If you use this tool please cite these references.

	Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community
	assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31,
	2501–2515.

	Gehara M, Garda AA, Werneck FP, et al. Estimating synchronous demographic changes across populations 
	using hABC and its application for a herpetological community from northeastern Brazil.
	Mol Ecol. 2017;00:1–16. https://doi.org/10.1111/mec.14239

[Instructions to run the hABC](hABC_manual.md)


### Nonhierarchical models
#### Using the menu to build your model
You start by setting up your model through the menu, you run the main.menu() function and directs the output to an R object that will store your model.

> model<-main.menu()

_You can use a previously build model as a template_:

> model2<-main.menu(model1)

The function will ask you to write a topology in newick format or "1" if you want to simulate a single pop model. 
After typing your option the menu will pop up.

> main.menu()  
write bifurcating topology in newick format or 1 for single population: (1,2);
```
A > Number of populations to simulate      2  
B > Tree (go here to remove nodes)         (1,2);  
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

#### Simulating summary stats

To simulate summary statistics you run the following function.

> sim.sumstat(model,path="",overall.SS=F,perpop.SS=T,nsim.blocks=1,use.alpha = T, append.sims=F, sim.block.size = 1000)

* model = your model object.  
* path = path to write the output.  
* overall.SS = if TRUE calculates the summary stats across all your populations. Defaut is FALSE.  
* perpop.SS = if TRUE calculates the summary stats per population. Default is TRUE.  
* use.alpha = if TRUE the most recent population size change will be exponential. Sudden change if FALSE. Default is FALSE.  
* sim.block.size = simulations are performed in blocks. This argument defines the size of the block, i.e. how many simulations you what to run per block.
* nsim.blocks = number of blocks you want to simulate. The total number of simulations is: nsim.blocks x sim.block.size.  
* append.sims = if TRUE simulations will be appended in the last simulation output. Default is FALSE.  

####Simulating coalescent trees

> sim.coaltrees(model,path="",use.alpha=F,nsim.blocks=1, append.sims=F, sim.block.size = 1000)


