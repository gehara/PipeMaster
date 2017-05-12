
mcmc.demog<-function(samples,
                     thin,
                     burnin,
                     Ne.prior,
                     NeA.prior,
                     time.prior,
                     mi.prior,
                     gene.prior,
                     alpha,
                     obs,
                     factor=c(20,30,20,10),
                     mcmc.par=1.05,
                     update){


  
priors<-rbind(Ne.prior,NeA.prior,time.prior,mi.prior)

pars<-numeric()
for(i in 1:4){
pars[i]<-runif(1,priors[i,1],priors[i,2])
}

distance<-10
PAR<-matrix(ncol=9,nrow=samples)
colnames(PAR)<-c("Distace","Ne","NeA","Exp.time","mi",
                "Pi","Seg.Sites","H","TD")
k<-1
n<-0

jitter.fac<-factor

for(j in 1:(samples+burnin)){

  for(i in 1:thin){

    
    new.pars<-pars
    
    new.pars[k] <- jitter(new.pars[k], jitter.fac[k])
    
    tes<-as.numeric(new.pars[k]<priors[k,1])
    tes<-tes+as.numeric(new.pars[k]>priors[k,2])
    
    while(tes>=1) {
      new.pars[k] <- jitter(pars[k], factor[k])
      tes<-as.numeric(new.pars[k]<priors[k,1])
      tes<-tes+as.numeric(new.pars[k]>priors[k,2])
      }

    x<-mcmc.exp.pars(gene.prior=as.matrix(gene.prior),new.pars = new.pars,gentime=1)

    y<-coexp.MS(MS.par=x$MS.par, gene.prior = gene.prior,alpha=alpha)

    ss<-mcmc.sumstat(ms.output=y,gene.prior=gene.prior)

    while(ss[1,1]==0) {
      
      x<-mcmc.exp.pars(gene.prior=as.matrix(gene.prior), new.pars = new.pars,gentime=1)
      
      y<-coexp.MS(MS.par=x$MS.par, gene.prior = gene.prior,alpha=alpha)
      
      ss<-mcmc.sumstat(ms.output=y,gene.prior=gene.prior)
      
    }

    new.dist<-mean(abs((obs-ss)/obs))

    if(new.dist<distance){
      pars<-new.pars
      distance<-new.dist
      xx<-sample(length(pars),1)
      while(xx==k) xx<-sample(length(pars),1)
      k<-xx
      }else{distance<-distance*mcmc.par}
    
    if(distance==0){
      distance<-0.5
    }
    
    }
  
  if(j<=burnin) print(paste(j,"burnin"))
  
  if(j>burnin){
    print(paste(j-burnin,"samples",round(distance,4),floor(pars[1]),signif(pars[2],digits=4),floor(pars[3]),
                paste(jitter.fac,collapse = " ")))
    PAR[(j-burnin),]<-c(distance,pars,ss)
    n<-n+1
  }
  
  if(n==update){
    plot.update(PAR,obs)
    n<-0
  }
  
}
return(PAR)
}

#obs<-observed.mcmc.demog("~/Dropbox/AMNH/Cophylogeography/fastas")[3,]
#obs[1]<-obs[1]*917
#posterior<-mcmc.demog(samples=10000,thin=20,burnin=1000,factor=c(20,30,20,20),
#                      Ne.prior = c(1000,20000000),
#                      NeA.prior = c(100,10000000),
#                      time.prior = c(20000,5000000),
#                      mi.prior = c(9e-9,1.1e-8),
#                      gene.prior=gene.prior[3,],
#                      obs=obs,
#                      mcmc.par = 1.01,
#                      alpha=F,
#                      update=100)

#summary(posterior)


