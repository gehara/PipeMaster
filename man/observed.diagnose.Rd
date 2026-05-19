% Manually written; regenerate via devtools::document() if you have roxygen2.
\name{observed.diagnose}
\alias{observed.diagnose}
\alias{print.observed_diagnostic}
\alias{summary.observed_diagnostic}
\alias{plot.observed_diagnostic}
\title{Diagnose non-modeled signals in observed per-locus statistics}
\usage{
observed.diagnose(per_locus, locus_table, cds_bed = NULL,
                  reftable = NULL, pop = 1, alpha = 0.01,
                  subsample = NULL, plot = TRUE, verbose = TRUE)

\method{print}{observed_diagnostic}(x, \dots)
\method{summary}{observed_diagnostic}(object, \dots)
\method{plot}{observed_diagnostic}(x, subsample = NULL, \dots)
}
\arguments{
\item{per_locus}{data.frame of per-locus statistics. Required columns:
  \code{locus_id}, \code{length}, \code{S}, \code{pi_per_site}, \code{D}.
  Recommended: \code{thetaW_per_site}.}

\item{locus_table}{data.frame mapping loci to genome coordinates.
  Required columns: \code{locus_id}, \code{chrom}, \code{start},
  \code{end}.}

\item{cds_bed}{Optional path to a BED file of CDS (or any selection-
  target) intervals. If \code{NULL}, the BGS-gradient test is skipped.}

\item{reftable}{Optional PipeMaster simulation reference table as
  produced by \code{sim.scrm.sumstats()}. If \code{NULL}, the
  over-dispersion test is skipped.}

\item{pop}{Integer population index (default 1) selecting which
  per-population moment columns to compare against.}

\item{alpha}{Significance level for the spatial-autocorrelation flag.}

\item{subsample}{Optional number of loci to subsample for plotting.}

\item{plot}{Logical. Draw the diagnostic figure.}

\item{verbose}{Logical. Print progress.}

\item{x, object}{An \code{observed_diagnostic} object.}

\item{\dots}{Ignored.}
}
\value{
An object of class \code{c("observed_diagnostic", "OOD_diagnostic")}
with fields \code{test_A}, \code{test_B}, \code{test_C}, \code{test_D},
\code{verdict}, and \code{verdict_reasons}.
}
\description{
Pre-training data diagnostic that flags departures from the assumptions
of a neutral, uniform-rate coalescent before generating a reference
table or fitting a model. Runs up to four tests on per-locus
statistics: joint D-pi distribution (Test A), per-locus moment
over-dispersion vs a PipeMaster reference table (Test B), per-chromosome
landscape with lag-1 spatial autocorrelation (Test C), and a
background-selection gradient with respect to distance from CDS or
conserved-region annotation (Test D).
}
\details{
Verdict combines Tests B, C, and D into \code{"pass"} / \code{"warn"}
/ \code{"fail"}. Test A is reported as informational only --- a strong
positive correlation between per-locus D and pi can arise from
background selection, mutation-rate heterogeneity, or as a statistical
property of the Tajima's D estimator on short loci under non-stationary
demography. Use \code{cds_bed} (Test D) to distinguish background
selection. A future release will add a neutral-coalescent simulation
companion to definitively distinguish the third case.

Per-locus statistics can be computed externally (e.g., bcftools query
plus the standard Tajima's D formula) and supplied as a data.frame
matching the schema in the \code{per_locus} argument description.
}
\seealso{
\code{\link{observed.sumstats}} for moment-level observed statistics;
\code{\link{OOD.pretrain}} for the post-reftable, pre-training OOD
diagnostic on the simulation cloud.
}
