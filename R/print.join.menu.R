print.join.menu<-function()
  {
  if(eJt[[1]][3]=="normal")
    dist.par<-"Mean, SD"
  if(eJt[[1]][3]=="uniform")
    dist.par<-"min, max"

  cat(paste("A > Junction prior distribution:       ",eJt[[1]][3]),
      paste("P > Priors                          c(",dist.par,")"),
      paste("                    ",names(eJt),"     ",eJt),
      paste(" "),
      paste("T > Topology                          ",tree),
      paste("B > Back to main menu"),
      sep="\n")
}
