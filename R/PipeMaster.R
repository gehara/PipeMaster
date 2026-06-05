#' PipeMaster: Coalescent Simulation and Demographic Inference
#'
#' PipeMaster builds demographic models, simulates genetic data under the
#' coalescent, and infers demographic parameters via approximate Bayesian
#' computation (ABC) or neural-network supervised machine learning (SML).
#' It targets both RADseq-scale and whole-genome-scale data through two
#' native simulation engines and provides a Shiny graphical interface for
#' model construction.
#'
#' @section Model building:
#' [main.menu.gui()] launches the Shiny GUI for constructing
#' \code{Model} objects with priors, population structure, demographic
#' events, migration, and locus structure. [optimize.sfs.model()] solves
#' the uniform-sampling problem required for SFS-based inference.
#'
#' @section Empirical data:
#' [observed.sumstats()] computes the full summary-statistic panel and
#' (optionally) the joint site frequency spectrum from PHYLIP or VCF
#' files; column layout matches the simulators so observed and simulated
#' rows are reftable-compatible. [observed.diagnose()] flags non-modelled
#' signals (background selection, recombination-rate heterogeneity,
#' chromosome-landscape effects) before reference-table generation.
#' Format converters: [alleles2phylip()], [filter.phylip()],
#' [read.phylip.loci()], [one.snp.per.locus()],
#' [get.data.structure()].
#'
#' @section Coalescent simulation:
#' [sim.sumstats()] drives the native msABC engine (RADseq-scale, many
#' short loci). [sim.scrm.sumstats()] drives the vendored scrm SMC' engine
#' (whole-genome-scale, with optional per-locus per-simulation mutation-
#' and recombination-rate heterogeneity). Both compute summary statistics
#' and the joint SFS directly in C/C++. [run.msABC()] exposes the raw
#' msABC engine for advanced use.
#'
#' @section Dimensionality reduction:
#' [pls.fit()] and [pls.project()] implement partial least-squares
#' projection via NIPALS; [reduce.sfs()] marginalises or coarse-bins
#' high-dimensional joint SFS for tractable inference; [plot.2D.sfs()]
#' visualises 2D joint spectra in the dadi/moments convention.
#'
#' @section Neural-network inference:
#' [tune.nn()] tunes ResNet (summary statistics) or 1D/2D CNN (SFS)
#' regression architectures via Hyperband, optionally in parallel across
#' GPUs. [nn.predict()] returns point estimates plus residual-bootstrap
#' and ABC-NN-regression posteriors. [tune.nn.classify()] trains a
#' multi-class classifier for demographic model selection; [nn.predict.classify()]
#' provides softmax, Mahalanobis, and latent-NN output rules.
#' [tune.nn.sensitivity()] computes the inverse Jacobian and elasticity
#' at observed data. [save.tune.result()] / [load.tune.result()] and
#' [save.classifier()] / [load.classifier()] persist trained models.
#'
#' @section Out-of-distribution diagnostics:
#' Two tiers, each in regression and classification flavours. Pre-training
#' ([OOD.pretrain()], [OOD.pretrain.classify()]) operates on the
#' reference table alone and checks marginal support plus joint density in
#' PCA and PLS / PLS-DA space. Post-training ([OOD.posttrain()],
#' [OOD.posttrain.classify()]) adds the network's latent-space density
#' and ensemble-disagreement signals. [OOD.attribute.classify()] uses
#' integrated gradients to attribute classifier-rejected observations to
#' specific summary statistics. [OOD.projection.diagnose()] forensically
#' inspects projections; [OOD.outliers()] inspects per-statistic outlier
#' distributions; [OOD.priors.bestfit()] and [OOD.priors.outliers()]
#' translate diagnostic results into actionable prior-widening suggestions.
#'
#' @docType package
#' @name PipeMaster
NULL

