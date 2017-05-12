
ms.commander<-function(model){
  
  # bing Ne, mig and Time priors
  size.pars<-rbind(model$flags$n,model$flags$en$size)
  mig.pars<-rbind(model$flags$m,model$flags$em$size)
  time.pars<-rbind(model$flags$ej,model$flags$en$time,model$flags$em$time)
  # sample Ne, div.time and mutation rate
  size.pars<-sample.w.cond(par.matrix=size.pars,cond.matrix=model$conds$size.matrix)
  time.pars<-sample.w.cond(par.matrix=time.pars,cond.matrix=model$conds$time.matrix)
  loci<-sample.pars(model$loci)
  # empty parameter vector
  parameters<-vector()
  # bind Ne and time sampled parameters
  parameters<-rbind(parameters,size.pars[,c(1,4)])
  parameters<-rbind(parameters,time.pars[,c(1,4)])
  # generate coalescent scalar. Arbitrary value 
  if(model$flags$n[1,6]=="runif") {Ne0<-min(as.numeric(model$flags$n[,4]))
  } else {Ne0<-mean(as.numeric(model$flags$n[,4]))}
  ms.scalar<-4*Ne0
  # transform parameters to fit the scalar
  size.pars[,4:5]<-as.numeric(size.pars[,4])/ms.scalar
  time.pars[,4:5]<-as.numeric(time.pars[,4])/ms.scalar
  # empty string, ms flags
  string<-list()
  #
  if(is.null(mig.pars)==T){
  } else {
    mig.pars<-sample.w.cond(par.matrix=mig.pars,cond.matrix=model$conds$mig.matrix)
    #bind sampled migration parameters
    parameters<-rbind(parameters,mig.pars[,c(1,4)])
    #transform mig parameters
    mig.pars[,4:5]<-as.numeric(mig.pars[,4])*as.numeric(size.pars[,4])
    #
    emt<-subset(time.pars, time.pars[,2]=="-em")
    em<-subset(mig.pars, mig.pars[,2]=="-em")
    mig.pars<-subset(mig.pars, mig.pars[,2]=="-m")
    ## generate migration string
    m<-NULL
    for(i in 1:nrow(mig.pars)){
      m<-c(m,paste(mig.pars[i,c(2:4)],collapse=" "))
    }
    string[[3]]<-paste(m, collapse=" ")
    ## generate ancestral migration string
    if(nrow(em)!=0){
    m<-NULL
    for(i in 1:nrow(emt)){
      m<-c(m,paste(emt[i,2],emt[i,4],em[i,3],em[i,4],collapse=" "))  
      }
    string[[4]]<-paste(m, collapse=" ")
  }
    #
  }
  #
  #### bind sampled mutation rate
  parameters<-rbind(parameters,loci[,c(1,4)])
  #### bind scaled theta per gene (4N*m*pb)
  loci<-cbind(loci,4*Ne0*as.numeric(loci[,4])*as.numeric(loci[,2])*as.numeric(loci[,3]))
  # generate Ne string
  l<-NULL
  for(i in 1:as.numeric(max(size.pars[,3]))){
    l<-c(l,paste0(size.pars[i,c(2:4)],collapse=" "))
    }
  string[[1]]<-paste0(l,collapse = " ")
  # generate ancestral Ne string   
  ent<-subset(time.pars, time.pars[,2]=="-en")
  en<-subset(size.pars, size.pars[,2]=="-en")
  n<-NULL
  if(nrow(en)!=0){
  for(i in 1:nrow(en)){
    n<-c(n,paste(ent[i,2],ent[i,4],en[i,3],en[i,4],collapse=" "))
    }
  string[[2]]<-paste(n, collapse=" ")
  }
  # generate ej string
  ej<-subset(time.pars, time.pars[,2]=="-ej")
  j<-NULL
  for(i in 1:nrow(ej)){
    j<-c(j,paste(ej[i,2],ej[i,4],ej[i,3],collapse=" "))  
    }
  string[[5]]<-paste(j, collapse=" ")
  # paste strings
  string<-paste(unlist(string),collapse=" ")
  # generate -t and -I part of the command
  commands<-list(NULL)
    for(i in 1:nrow(loci)){
    y<-paste("-t",loci[i,7],paste(model$I[i,2:ncol(model$I)],collapse=" "),collapse=" ")
    commands[[i]]<-paste(y,string, collapse=" ")
    }
  commands[[nrow(loci)+1]]<-t(parameters)
  return(commands)
}
