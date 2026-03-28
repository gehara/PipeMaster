# PipeMaster

> **PipeMaster v2 is under active development.** Models built with previous
> versions may not work out of the box. See [Breaking Changes](#breaking-changes)
> below for migration instructions.

PipeMaster is an R package for building demographic models and simulating data under the coalescent. It supports:

- **Site Frequency Spectrum (SFS)** simulation for 1-pop, 2-pop, and 3-pop models
- **Summary statistics** simulation for single and multi-population models
- **Two simulation engines**: built-in msABC (RADseq) and vendored scrm 1.7.5 (WGS-scale with recombination)
- **Shiny GUI** (`main.menu.gui()`) for interactive model building and visualization
- **Neural network SML** (`tune.nn()`, `nn.predict()`) with keras and torch backends
- Approximate Bayesian computation (`abc.rejection()`) for parameter estimation
- **Out-of-distribution diagnostics** (`OOD.diagnose()`) for model validation

## Breaking Changes

### Migration parameterization (v2)

Migration priors are now specified as **Nm** (number of migrants per generation)
instead of 4Nm. This is more intuitive: Nm = 1 means one migrant per generation.

| | Old (v1) | New (v2) |
|---|---|---|
| **Default scale** | 4Nm | Nm |
| **Value for 1 migrant/gen** | 4.0 | 1.0 |
| **Model field** | (none) | `model$mig_scale = "Nm"` |

**To update old models:** divide migration prior bounds by 4, or add
`model$mig_scale <- "Nm"` and adjust values accordingly. Old models without a
`mig_scale` field will use the legacy 4Nm convention.

A new **per-generation rate** option (`model$mig_scale = "m"`) is also available
for specifying symmetric migration. In this mode, both directions of a migration
pair share the same prior and are constrained equal via conditions. This is useful
for models like Gutenkunst et al. (2009) where migration is symmetric.

### Conditions format (v2)

Conditions (`model$conds`) are now specified as a **list of lists** instead of
matrices. Each condition is a list with `param1`, `op`, and `param2` fields:

```r
# Old format (matrix):
conds$size.matrix <- matrix(c("Ne1.pop1", "<", "Ne0.pop1"), nrow = 1)

# New format (list of lists):
conds$size <- list(
  list(param1 = "Ne1.pop1", op = "<", param2 = "Ne0.pop1")
)
```

Old models using matrix conditions will need to be converted to the list format.

### Labels (v2)

Models now support optional `model$labels` with a model name and population names
for GUI visualization:

```r
model$labels <- list(
  name = "MyModel",
  pops = setNames(c("PopA", "PopB"), c("1", "2"))
)
```

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
