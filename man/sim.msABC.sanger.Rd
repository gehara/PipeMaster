% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.msABC.R
\name{sim.msABC.sanger}
\alias{sim.msABC.sanger}
\title{Simulate summary statistics using msABC}
\usage{
sim.msABC.sanger(
  model,
  use.alpha = F,
  nsim.blocks = 5,
  path = getwd(),
  append.sims = F,
  block.size = 1000,
  msABC.call = get.msABC(),
  output.name = "model",
  ncores = 1
)
}
\arguments{
\item{model}{A model object bult by the main.menu function.}

\item{use.alpha}{Logical.If TRUE the most recent population size change will be exponential. If FALSE sudden demographic changes. Default is FALSE.}

\item{nsim.blocks}{Number of blocks to simulate. The total number of simulations is: nsim.blocks x sim.block.size.}

\item{path}{Path to write the output. By default outputs will be saved in the working directory.}

\item{append.sims}{Logical. If TRUE simulations will be appended in the last output. Default is FALSE.}

\item{msABC.call}{String. Path to the msABC executable. msABC binaries for Mac's and Linux are included in the package and should work in most computers.}

\item{output.name}{String. The prefix of the output names. Defalt is "model"}

\item{sim.block.size}{Simulations are performed in blocks. This argument defines the size of the block in number of simulations, i.e. how many simulations to run per block.}

\item{get.moments}{Logical. If TRUE computes the four moments (mean, variance, kurtosis, skewness) of each summary statistics across loci. If False only mean is computed. Defalt is FALSE.}
}
\value{
Writes simulations and parameters to the path directory. msABC outputs a bunch of summary stats by defalt. They need to be selectd a posteriori.
}
\description{
Simulate summary statistics using msABC
}
\note{
This function does not work on Windows systems.
}
\references{
Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338.

Pavlidis P., Laurent S., & Stephan W. (2010) msABC: A modification of Hudson’s ms to facilitate multi-locus ABC analysis. Molecular Ecology Resources, 10, 723–727.
}
\author{
Marcelo Gehara
}
