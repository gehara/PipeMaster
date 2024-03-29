% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coexpansion_PT.R
\name{sim.coexpPT}
\alias{sim.coexpPT}
\title{Simulation of codemographic models}
\usage{
sim.coexpPT(
  nsims,
  var.zeta,
  coexp.prior,
  Ne.prior,
  NeA.prior,
  time.prior,
  gene.prior,
  alpha = F,
  append.sims = F,
  path = getwd()
)
}
\arguments{
\item{nsims}{Total number of simulations}

\item{var.zeta}{Variation on zeta parameter. Can be "FREE" to vary or be set to a specific value (between 0-1).}

\item{coexp.prior}{Uniform prior for the coespansion time. Vector of two numbers with the lower and upper boudary of the prior.}

\item{Ne.prior}{Data frame with the prior values for the Ne of each population.}

\item{NeA.prior}{Data frame with the prior values for the ancestral Ne of each population.}

\item{time.prior}{Data frame with parameter values for the priors of the time of demographic change of each population.}

\item{gene.prior}{Data frame with parameter values for the priors of the mutation rate of each species.}

\item{alpha}{logical. If TRUE all demographic chages are exponential. If FALSE sudden changes. Defaut is FALSE.}

\item{append.sims}{logical. If TRUE simulations are appended to the simulations file.  Defaut is FALSE.}

\item{path}{Path to the directiry to write the simulations. Defaut is the working directory.}
}
\description{
Simulation of the PT model of Gehara et al. (2017)
}
\details{
To simulate the model of Chan et al. (2014), the Threshold model and the Narrow Coexpansion Time model use the sim.coexp function.

See references for more details.
}
\references{
Gehara M., Garda A.A., Werneck F.P. et al. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756–4771.

Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31, 2501–2515.
}
