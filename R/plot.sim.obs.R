#' Plot simulated data against observed
#' @description This function plots an histogram of the simulated summary statistics with the corresponding observed value as a red line for comparison.
#' @param sim A data frame with simulated data.
#' @param obs A vector of observed summary stats corresponding to the sim object.
#' @return Graphic
#' @author Marcelo Gehara
#'
plot.sim.obs<-function (sim, obs)
  {
    mylabels <- colnames(sim)
    for (i in 1:ncol(sim)) {
      hist(sim[, i], breaks = 20, xlab = mylabels[i], main = "")
      abline(v = obs[i], col = 2)
    }
  }
