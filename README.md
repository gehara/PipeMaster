# PipeMaster
** under development**
  
  PipeMaster is an R-package to build demographic models and simulate data under the coalescent. Current implementation can simulate sanger-type and nexgen data for single species or complex of species; single-locus data for hierarchical demographic models of comparative phylogeography; and species trees with one horizontal connection.

**This is a beta version of the package!**
  
  PipeMaster simulates summary statistics and coalescent trees. It calculates the same summary statistics on an empirical data. The user can use these sumary statistics to perform aproximate Bayesian computation (ABC) or supervized machine learning (SML) for model and parameter inference.

### Nonhierarchical models tutorial with instalation instructions.

[Main Tutorial](PipeMaster_tutorial.md)

[Tutorial em Português](PipeMaster_tutorial_PT.md)

[CompPhylo Workshop Tutorial](https://compphylo.github.io/Oslo2019/PM_files/Dermatonotus_example.html)

[Example of analysis of nexgen data](Dermatonotus_example.md)

[Example of analysis of sanger data](Agkistrodon_example.md)


### Hierarchical codemographic model

PipeMaster can simulate a hierarchical demographic model for comparative analys of populations/species. The hABC method used in the package was first described in Chan et al 2014 and improved in Gehara et al 2017. These two papers show good examples of this analysis. If you use this tool please cite these references.

Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community
assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31,
2501–2515.

Gehara M, Garda AA, Werneck FP, et al. Estimating synchronous demographic changes across populations 
using hABC and its application for a herpetological community from northeastern Brazil.
Mol Ecol. 2017;00:1–16. https://doi.org/10.1111/mec.14239

[Instructions to run the hABC](hABC_manual.md)



