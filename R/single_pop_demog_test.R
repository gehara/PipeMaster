#' @export
single.pop.demog<-function(nsims,
                     Ne.prior,
                     time.prior,
                     gene.prior,
                     observed,
                     alpha=F,
                     tol=0.01,
                     nval=100,
                     method="rejection",
                     do.ABC=F,
                     do.PCA=F,
                     CV=F,
                     mod=cbind(c(1,0.001,2),c(1,0.1,20)),
                     path=path){
  tabela<-NULL


  for(i in 1:nrow(Ne.prior)){

    parameters<-NULL
    models<-NULL
    index<-NULL
    mod<-mod
    rownames(mod)<-c("CS","Exp","BOTT")

    for(j in 1:nrow(mod)){
      sim.demog(nsims=nsims,coexp.prior=coexp.prior,Ne.prior=Ne.prior[i,], alpha=alpha,
                NeA.prior=mod[j,],time.prior=time.prior[i,],gene.prior=gene.prior[i,],append.sims = F, path=path)
      result<-read.table(file=paste(Ne.prior[i,1],"demog_sim.txt",sep=""), header=T)
      pars<-read.table(file=paste(Ne.prior[i,1],"pop_parameters.txt",sep=""), header=T)
      models<-rbind(models,result)
      parameters<-rbind(parameters,pars)
      index<-c(index,rep(rownames(mod)[j],nrow(result)))
    }

    if(do.ABC==T){
      if(method=="rejection"){
      prob<-postpr(observed[i,],index,models,method=method, tol=tol)
      prob<-summary(prob)
      tabela<-rbind(tabela,prob$Prob)
      }

      if(method=="neuralnet"){
      prob<-postpr(observed[i,],index,models,method=method, tol=tol)
      prob<-summary(prob)
      tabela<-rbind(tabela,prob$neuralnet$Prob)
      }
    }

    if(CV==T){
      cv<-cv4postpr(index,models,nval=nval,tols=tol,method=method)
      pdf(paste("CV_",rownames(observed)[i],".pdf",sep=""), paper="a4r", width=10, pointsize=10)
      plot(cv)
      graphics.off()
    }
    if(do.PCA==T){

      theme_set(theme_grey(base_size = 30))
      data<-c(index,"observed")
      x<-rbind(models,observed[i,])

      remove.zero.var<-function(x){
        tab<-NULL
        for(u in 1:ncol(x)){
          if(length(unique(x[,u]))>1)
            {
            tab<-cbind(tab,x[,u])
            }
        }
        return(tab)
      }

      x<-remove.zero.var(x)

      PCA<-prcomp(x, center = T, scale. = T, retx=T)
      scores <- data.frame(PCA$x[,1:2])

    PCA.plot<-ggplot(scores, aes(x=PC1, y=PC2))+
      theme(legend.text = element_text(size = 30, face = "bold"))+
      geom_point(aes(colour=data, size=data, shape=data))+
      scale_size_manual(values=c(3,3,3,10))+
      scale_colour_brewer(palette="Spectral")
    pdf(paste("PCA12",rownames(observed)[i],".pdf",sep=""), paper="a4r", width=10, pointsize=10)
    plot(PCA.plot)
    graphics.off()
}
    write.table(cbind(index,models),file=paste(Ne.prior[i,1],"demog_sim.txt",sep=""), quote=F,row.names=F, col.names=T, append=F, sep="\t")
    write.table(parameters,file=paste(Ne.prior[i,1],"pop_parameters.txt",sep=""), quote=F,row.names=F, col.names=T, append=F,sep="\t")

  }

  if(do.ABC==T){
    rownames(tabela)<-rownames(observed)
    return(tabela)
  }

}

#' @export
sim.demog<-function(nsims,
                    coexp.prior,
                    Ne.prior,
                    NeA.prior,
                    time.prior,
                    gene.prior,
                    alpha=alpha,
                    append.sims=F,
                    path=getwd())
{

  setwd(path)


  if(append.sims==F){
    simulations<-matrix(nrow=1,ncol=7)
    simulations[1,]<-c("pi","ss","H","TD","ss1","ss2","ss3")
    write.table(simulations,file=paste(Ne.prior[,1],"demog_sim.txt",sep=""), quote=F,row.names=F, col.names=F, sep="\t")
    populations.par<-matrix(c("Ne","Exp.time","NeA","mi"),nrow=1,ncol=4)
    write.table(populations.par,file=paste(Ne.prior[,1],"pop_parameters.txt",sep=""), quote=F,row.names=F, col.names=F, sep="\t")
  }


  TIME<-system.time(for (i in 1:nsims){

    x<-demog.sample.pars(nruns=1,coexp.prior=coexp.prior,Ne.prior=Ne.prior,
                         NeA.prior=NeA.prior,time.prior=time.prior,gene.prior=gene.prior)

    y<-coexp.MS(MS.par=x$MS.par, gene.prior = gene.prior,alpha=alpha)

    z<-sumstat(ms.output=y,gene.prior=gene.prior)

    populations.par<-unlist(x$pop.par)

    write.table(z,file=paste(Ne.prior[,1],"demog_sim.txt",sep=""), quote=F,row.names=F, col.names=F, append=T, sep="\t")
    write.table(t(populations.par),file=paste(Ne.prior[,1],"pop_parameters.txt",sep=""), quote=F,row.names=F, col.names=F, append=T,sep="\t")
    print(paste(i,"sims of",nsims,"| single pop demogragphic test"))
  })
  print(TIME)
}

# internal function of the test.demog function
demog.sample.pars<-function(nruns,
                            var.zeta,
                            coexp.prior,
                            buffer,
                            Ne.prior,
                            NeA.prior,
                            time.prior,
                            gene.prior){

  MS.par<-list(NULL)
  pop.par<-list(NULL)
  coexp.par<-matrix(nrow=1,ncol=4)
  nspecies<-1
  MS.par[[1]]<-matrix(nrow=1,ncol=4)
  pop.par[[1]]<-matrix(nrow=1,ncol=4)

  ms.par<-NULL

  Ne <- runif(1, Ne.prior[1,3], Ne.prior[1,4])
  e.t<-runif(1, time.prior[1,3], time.prior[1,4])
  Ne.EXP.t <- e.t/time.prior[1,5] #corrects by generations
  theta.A.ratio <- runif(1, NeA.prior[1], NeA.prior[2])# thetaA (NeA) ratio
  NeA <- Ne*theta.A.ratio
  mi <- do.call(as.character(gene.prior[1,2]), args=list(1, gene.prior[1,3], gene.prior[1,4]), quote=F)
  while(mi<0){
    mi <- do.call(as.character(gene.prior[1,2]), args=list(1, gene.prior[1,3], gene.prior[1,4]), quote=F)
  }
  po.par<-c(Ne, e.t, NeA,mi)

  Ne <- Ne*gene.prior[1,7]
  theta=4*Ne*mi*gene.prior[1,5]
  scalar=4*Ne
  EXP.time=Ne.EXP.t/scalar

  g.rate=-log(NeA/Ne)/Ne.EXP.t

  ms.par<-cbind(theta,EXP.time, theta.A.ratio,g.rate)

  MS.par[[1]][1,]<-ms.par
  pop.par[[1]][1,]<-po.par

  pars<-list(NULL,NULL,NULL)
  names(pars)<-c("MS.par", "pop.par")
  pars$MS.par<-MS.par
  pars$pop.par<-po.par
  return(pars)
}

# internal function of the test.demog function
sumstat<-function(ms.output, gene.prior){
  sum.stat<-NULL
  for (j in 1:length(ms.output)){

    ss<-as.numeric(strsplit(ms.output[[j]][3]," ")[[1]][2])
    x<-ms.output[[j]][5:length(ms.output[[j]])]
    x<-gsub("0","A",x)
    x<-gsub("1","C",x)
    se<-list(NULL,NULL,NULL,NULL)
    names(se)<-c("nb","seq","nam","com")
    se$nb<-length(x) # number of samples
    se$seq<-x # sequencies
    se$nam<-c(1:length(x)) # sequence names, just numbers
    se$com<-NA
    class(se)<-"alignment" # this is an alignment
    x<-as.DNAbin(se)

    pi<-nuc.div(x)*ss/gene.prior[j,5]
    H<-H.div(x)
    TD<-tajima.test(x)$D
    #R2<-R2.test(x,B=0,plot = F,quiet = T)
    spec<-site.spectrum(x)[1:3]
    SS<-c(pi,ss,H[2],TD,spec)
    sum.stat<-rbind(sum.stat,SS)
  }
  return(sum.stat)
}

