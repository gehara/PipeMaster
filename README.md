# PipeMaster

PipeMaster is an R package for building demographic models and simulating data under the coalescent. It supports:

- **Site Frequency Spectrum (SFS)** simulation for 1-pop, 2-pop, and 3-pop models
- **Summary statistics** simulation for single and multi-population models
- **Shiny GUI** (`main.menu.gui()`) for interactive model building
- Approximate Bayesian computation (ABC) and supervised machine learning for model selection and parameter estimation

PipeMaster uses a built-in msABC C implementation for fast coalescent simulation and summary statistic computation.

### Tutorials

[Main Tutorial (stdpopsim models)](PipeMaster_tutorial.md) - Covers SFS and summary statistics workflows using well-characterized stdpopsim demographic models

[Tutorial em Portugues](PipeMaster_tutorial_PT.md)

[CompPhylo Workshop Tutorial](https://compphylo.github.io/Oslo2019/PM_files/Dermatonotus_example.html)

[Legacy Tutorial (Dermatonotus example)](PipeMaster_tutorial_old.md)

### Hierarchical codemographic model

PipeMaster can simulate a hierarchical demographic model for comparative analysis of populations/species. The hABC method used in the package was first described in Chan et al 2014 and improved in Gehara et al 2017.

Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community
assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31,
2501-2515.

Gehara M, Garda AA, Werneck FP, et al. Estimating synchronous demographic changes across populations
using hABC and its application for a herpetological community from northeastern Brazil.
Mol Ecol. 2017;00:1-16. https://doi.org/10.1111/mec.14239

[Instructions to run the hABC](hABC_manual.md)
