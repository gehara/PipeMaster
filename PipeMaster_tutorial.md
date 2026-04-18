---
title: "PipeMaster v2 Tutorial"
subtitle: "Demographic Inference from RAD-seq and WGS Data"
output:
  html_document:
    toc: true
    toc_depth: 3
  pdf_document:
    toc: true
    toc_depth: '3'
---

This tutorial demonstrates PipeMaster's inference pipeline using three
demographic models of increasing complexity. All data files are included
with the package — no external downloads are needed.

## Contents

1. [Installation](#1-installation)
2. [Example Models](#2-example-models)
3. [Computing Observed Statistics](#3-computing-observed-statistics)
4. [Simulating Reference Tables](#4-simulating-reference-tables)
5. [ABC Rejection](#5-abc-rejection)
6. [Supervised Machine Learning (SML)](#6-supervised-machine-learning)
7. [Out-of-Distribution Diagnostics](#7-out-of-distribution-diagnostics)

---

## 1. Installation

Install PipeMaster from GitHub:

```r
install.packages("devtools")
devtools::install_github("gehara/PipeMaster")
```

Load the package:

```r
library(PipeMaster)
```

PipeMaster requires `ape`, `e1071`, and `msm` (installed automatically).
For the Shiny GUI, you also need `shiny`, `shinydashboard`, `shinyjs`, and
`DT`. For the torch neural network backend, install the `torch` package.

### Simulation engines

PipeMaster includes two coalescent simulation engines:

- **msABC** (built-in C) — a modified version of msABC (Pavlidis et al.
  2010), itself based on Hudson's ms (Hudson 2002). Fast for short loci
  (RAD-seq, 100 bp). Uses `sim.sumstats()`.
- **scrm** (vendored C++) — a modified version of scrm (Staab et al.
  2015), which uses the Sequential Markov Coalescent approximation to
  efficiently simulate long loci with recombination (WGS, 100 kb). Uses
  `sim.scrm.sumstats()`.

Both engines compute summary statistics and the site frequency spectrum (SFS)
in C during simulation — no intermediate files or parsing overhead.

### Parameter scaling

You specify parameters in **natural units** (Ne, generations, per-bp mutation
rate). PipeMaster converts them to coalescent units internally:

| Natural unit | Coalescent scale | Formula |
|---|---|---|
| Population size (Ne) | Relative size | Ne / Ne0 |
| Time (generations) | Coalescent time | t / (4 * Ne0) |
| Mutation rate (per bp per gen) | Theta per locus | 4 * Ne0 * mu * L |
| Migration (Nm) | Scaled migration | 4 * Ne0 * m |

Where Ne0 = 100,000 is the reference population size used for scaling.

---

## 2. Example Models

PipeMaster ships three example models based on well-characterized
[stdpopsim](https://popsim-consortium.github.io/stdpopsim-docs/stable/)
demographic histories. Each comes with pseudo-observed RAD-seq data (VCF)
and pre-computed WGS observed statistics.

### 2.1 Vaquita two-epoch bottleneck

A single-population bottleneck model for the vaquita porpoise (*Phocoena
sinus*) (Robinson et al. 2022). Three parameters.

| Parameter | True value | Description |
|-----------|------------|-------------|
| Ne0.pop1 | 2,807 | Current Ne |
| Ne1.pop1 | 4,485 | Ancestral Ne |
| t.Ne1.pop1 | 2,162 gen | Time of size change |

Tree: `(1);` | Constraint: Ne0 < Ne1 (bottleneck)

### 2.2 Orangutan isolation with migration

Two-population model for Bornean and Sumatran orangutans (Locke et al.
2011). Eight parameters with asymmetric migration in Nm scale.

| Parameter | True value | Description |
|-----------|------------|-------------|
| Ne0.pop1 | 8,805 | Current Bornean Ne |
| Ne0.pop2 | 37,661 | Current Sumatran Ne |
| Ne1.pop1 | 10,617 | Ancestral Bornean Ne |
| Ne1.pop2 | 7,317 | Ancestral Sumatran Ne |
| Ne.anc_1_2 | 17,934 | Ancestral Ne at split |
| join1 | 20,157 gen | Divergence time |
| mig0.1_2 | 0.06 Nm | Migration into Bornean |
| mig0.2_1 | 0.415 Nm | Migration into Sumatran |

Tree: `(1,2);` | Migration scale: Nm

### 2.3 Human Out-of-Africa

Three-population model for YRI, CEU, and CHB (Gutenkunst et al. 2009).
Fourteen free parameters with symmetric migration in per-generation rate (m).

| Parameter | True value | Description |
|-----------|------------|-------------|
| Ne0.pop1 | 12,300 | Current YRI Ne |
| Ne0.pop2 | 29,725 | Current CEU Ne |
| Ne0.pop3 | 54,090 | Current CHB Ne |
| Ne1.pop1 | 7,300 | Ancestral Ne |
| Ne1.pop2 | 1,000 | Initial CEU Ne at EU/AS split (N_EU0) |
| Ne1.pop3 | 510 | Initial CHB Ne at EU/AS split (N_AS0) |
| Ne.anc_3_2 | 2,100 | OoA bottleneck Ne (between EU/AS split and OoA) |
| join3_2 | 848 gen | CEU-CHB split time |
| join2_1 | 5,600 gen | OoA split time |
| t.Ne1.pop1 | 8,800 gen | YRI Ne change time |
| m_YRI_CEU | 2.47e-5 | Per-gen rate |
| m_YRI_CHB | 7.67e-6 | Per-gen rate |
| m_CEU_CHB | 6.38e-5 | Per-gen rate |
| m_ancestral | 1.71e-3 | Per-gen rate (YRI-OoA) |

Tree: `((3,2),1);` | Migration scale: m (symmetric via conditions)

### Loading example data

```r
library(PipeMaster)

# File paths for Vaquita example
extdir <- system.file("extdata", "Vaquita2Epoch", package = "PipeMaster")
list.files(extdir)
# [1] "model_diffuse.RData"      "model_informative.RData"
# [2] "obs_allstats_wgs.txt"     "pop_assign.txt"
# [3] "Vaquita2Epoch_RAD.vcf"

# Load model
load(file.path(extdir, "model_informative.RData"))
model <- get(ls(pattern = "Vaquita"))

# File paths
vcf_file <- file.path(extdir, "Vaquita2Epoch_RAD.vcf")
pop_file <- file.path(extdir, "pop_assign.txt")
```

---

## 3. Computing Observed Statistics

`observed.sumstats()` computes summary statistics and the folded SFS from
a VCF or PHYLIP file. It uses the same C stat functions as the simulation
engines, ensuring consistency between observed and simulated data.

### From a VCF file

```r
# Compute observed stats from the shipped RAD VCF
obs <- observed.sumstats(
  model       = model,
  path.to.vcf = vcf_file,
  pop.assign  = pop_file
)
cat(sprintf("Observed: %d summary statistics + SFS bins\n", ncol(obs)))
```

### Multi-population example (PonAbe)

```r
extdir_pa <- system.file("extdata", "PonAbe_TwoSpecies", package = "PipeMaster")
load(file.path(extdir_pa, "model_informative.RData"))
model_pa <- get(ls(pattern = "PonAbe"))

obs_pa <- observed.sumstats(
  model       = model_pa,
  path.to.vcf = file.path(extdir_pa, "PonAbe_TwoSpecies_RAD.vcf"),
  pop.assign  = file.path(extdir_pa, "pop_assign.txt")
)
```

### Using pre-computed WGS observed stats

WGS VCFs are too large to ship with the package. Pre-computed observed
statistics are provided instead:

```r
obs_wgs <- read.table(
  file.path(extdir, "obs_allstats_wgs.txt"),
  header = TRUE
)
```

---

## 4. Simulating Reference Tables

A reference table contains many simulations under the model with parameters
sampled from the prior distribution. Each row contains the sampled parameters
and the resulting summary statistics.

### RAD-seq reference table (msABC engine)

For short loci (e.g., 100 bp RAD-seq), use `sim.sumstats()`:

```r
sim.sumstats(
  model       = model,
  nsim.blocks = 10,       # blocks per worker
  block.size  = 10,       # simulations per block
  mu.rates    = 5.83e-9,  # mutation rate per bp per gen
  ncores      = 4,        # parallel workers
  path        = "output/vaquita",
  output.name = "reftable"
)
# Output: output/vaquita/SIMS_reftable.txt
```

Total simulations = nsim.blocks x block.size x ncores. For the Vaquita
model with 3 parameters, 10,000-50,000 simulations are sufficient.

### WGS reference table (scrm engine)

For long loci with recombination (e.g., 100 kb WGS), use
`sim.scrm.sumstats()`:

```r
sim.scrm.sumstats(
  model       = model,
  nsim.blocks = 10,
  block.size  = 10,
  mu.rates    = 5.83e-9,  # mutation rate
  rec.rates   = 1e-8,     # recombination rate
  skip.zns    = TRUE,     # skip ZnS (O(segsites^2))
  ncores      = 4,
  path        = "output/vaquita_wgs",
  output.name = "reftable"
)
```

### Loading a reference table

```r
reftable <- read.table("output/vaquita/SIMS_reftable.txt", header = TRUE)

# Identify parameter and statistic columns
param_cols <- c("Ne0.pop1", "Ne1.pop1", "t.Ne1.pop1")
nuisance   <- c("mean.rate", "sd.rate")
stat_cols  <- setdiff(colnames(reftable), c(param_cols, nuisance))
```

### Simulation budget guidelines

| Model complexity | Parameters | Recommended sims |
|------------------|------------|------------------|
| 1-pop (Vaquita) | 3 | 10,000-50,000 |
| 2-pop IM (PonAbe) | 8 | 50,000-100,000 |
| 3-pop (OoA) | 14-19 | 100,000-300,000 |

---

## 5. ABC Rejection

`abc.rejection()` performs ABC with optional PLS distance projection.
It accepts the simulated reference table and observed statistics, and
returns posterior samples at multiple tolerance levels.

```r
# Prepare observed vector (matching stat columns)
obs_raw <- as.numeric(obs[1, stat_cols])

# Run ABC rejection at three tolerance levels
result <- abc.rejection(
  reftable   = reftable,
  observed   = obs_raw,
  param.cols = param_cols,
  tol        = 0.01,       # accept closest 1%
  distance   = "sd",       # standardized Euclidean distance
  pls        = TRUE,       # PLS distance projection
  verbose    = TRUE
)

# Point estimate (weighted median of accepted simulations)
result$point_estimate

# Posterior samples
head(result$abc)
```

### Plotting ABC posteriors

```r
true_vals <- c(Ne0.pop1 = 2807, Ne1.pop1 = 4485, t.Ne1.pop1 = 2162)
plot(result, true_values = true_vals, show_prior = TRUE)
```

---

## 6. Supervised Machine Learning

PipeMaster's SML pipeline trains a neural network to learn the mapping from
summary statistics to parameters, then predicts parameters for the observed
data with uncertainty quantification.

### Step 1: Hyperparameter search with `tune.nn()`

```r
tune_result <- tune.nn(
  reftable,
  param.cols   = param_cols,
  type         = "sumstat",
  exclude.cols = nuisance,
  max_epochs   = 100,      # max training epochs per model
  eta          = 3,        # Hyperband reduction factor
  n_searches   = 20,       # parallel Hyperband searches
  top_k        = 3,        # keep top 3 models as ensemble
  cores        = 4,        # parallel workers
  seed         = 42,
  backend      = "torch",  # or "keras"
  verbose      = TRUE
)

# Save for later use
save.tune.result(tune_result, "output/vaquita/tune_result")
```

The `backend` parameter selects the neural network framework:
- `"torch"` — pure R, no Python dependency, lower memory usage
- `"keras"` — requires Python + TensorFlow, slightly faster training

### Step 2: Prediction with `nn.predict()`

```r
pred <- nn.predict(
  tune_result,
  observed   = obs_raw,
  reftable   = reftable,
  param.cols = param_cols,
  method     = c("bootstrap", "ABC_NN_regression"),
  n_boot     = 5000,       # bootstrap resamples
  tolerance  = 0.05,       # ABC-NN acceptance tolerance
  seed       = 42,
  verbose    = TRUE
)

# Point estimate (weighted ensemble prediction)
pred$point_estimate

# Posterior summary
summary(pred)
```

Two uncertainty methods are available:

- **Residual bootstrap**: Locally-weighted resampling of NN prediction
  residuals. Captures prediction uncertainty without reference to the prior.
- **ABC-NN regression**: Beaumont et al. (2002) regression adjustment.
  Proper Bayesian posterior incorporating prior and likelihood.

### Step 3: Plot posteriors

```r
plot(pred, true_values = true_vals, bw.adjust = 2)
```

The `bw.adjust` parameter controls density smoothing (default 1, increase
for smoother curves).

### Loading saved results

```r
tune_result <- load.tune.result("output/vaquita/tune_result")
```

---

## 7. Out-of-Distribution Diagnostics

`OOD.diagnose()` checks whether the observed data is compatible with the
simulated reference table. This is critical — if the observed data falls
outside the range of simulations, inference results are unreliable.

### Without a trained model (stat-level checks only)

```r
ood <- OOD.diagnose(
  trained.nn = NULL,
  observed   = obs_raw,
  reftable   = reftable[, stat_cols],
  stat.cols  = stat_cols,
  plot       = TRUE
)
cat(sprintf("Overall: %s\n", ood$overall))
```

Three checks are performed:
1. **Mahalanobis distance** — is the observed point far from the centroid?
2. **PCA projection** — does the observed point fall within the PCA cloud?
3. **Per-stat percentiles** — are individual statistics within range?

### With a trained model (adds model disagreement check)

```r
ood <- OOD.diagnose(
  trained.nn = tune_result,
  observed   = obs_raw,
  reftable   = reftable[, stat_cols],
  stat.cols  = stat_cols,
  plot       = TRUE
)
```

Check 4 (model disagreement) compares predictions across ensemble models.
Large disagreement suggests the observed data is in a region where the
neural network is uncertain.

### Using OOD to filter statistics

```r
# Drop outlier statistics before inference
ood_pass <- !ood$percentiles$outlier
stat_cols_filtered <- stat_cols[ood_pass]
cat(sprintf("Kept %d/%d statistics\n",
            length(stat_cols_filtered), length(stat_cols)))
```

---

## Migration Parameterization

PipeMaster v2 supports two migration scales, selected via
`model$mig_scale`:

### Nm (number of migrants) — default

Migration priors specify Nm, the number of migrants per generation into the
receiving population. This is population-size-dependent: the same Nm implies
different per-generation rates for populations with different Ne.

```r
model$mig_scale  # "Nm"
model$flags$m    # migration priors in Nm units
```

### m (per-generation rate) — symmetric migration

For models with symmetric migration (e.g., Gutenkunst 2009), use
`model$mig_scale = "m"`. Both directions share the same prior and are
constrained equal:

```r
model$mig_scale <- "m"
# In model$conds$mig:
#   mig0.1_2 = mig0.2_1  (symmetric YRI-CEU)
#   mig0.1_3 = mig0.3_1  (symmetric YRI-CHB)
```

---

## Interactive Model Building

PipeMaster includes a Shiny GUI for interactive model building:

```r
main.menu.gui()
```

The GUI provides:
- Population tree builder with drag-and-drop
- Prior distribution specification with real-time validation
- Migration rate configuration (Nm or m scale)
- Condition constraints (parameter ordering, equality)
- Model visualization with `PlotModel()`
- Export to R object or file

---

## References

- Gutenkunst R.N. et al. (2009) Inferring the joint demographic history of
  multiple populations from multidimensional SNP frequency data. *PLoS
  Genetics* 5: e1000695.
- Locke D.P. et al. (2011) Comparative and demographic analysis of
  orang-utan genomes. *Nature* 469: 529-533.
- Robinson J.A. et al. (2022) The critically endangered vaquita is not
  doomed to extinction by inbreeding depression. *Science* 376: 635-639.
- Beaumont M.A. et al. (2002) Approximate Bayesian computation in
  population genetics. *Genetics* 162: 2025-2035.
- Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral
  model of genetic variation. *Bioinformatics* 18: 337-338.
- Pavlidis P., Laurent S. & Stephan W. (2010) msABC: a modification of
  Hudson's ms to facilitate multi-locus ABC analysis. *Molecular Ecology
  Resources* 10: 723-727.
- Staab P.R. et al. (2015) scrm: efficiently simulating long sequences
  using the approximated coalescent with recombination. *Bioinformatics*
  31: 1680-1682.
