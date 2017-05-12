plot.sim.obs<-function (sim, obs) 
{
  mylabels <- colnames(sim)
  for (i in 1:ncol(sim)) {
    hist(sim[, i], breaks = 20, xlab = mylabels[i], main = "")
    abline(v = obs[i], col = 2)
  }
}
