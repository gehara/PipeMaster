# internal function of the ms.commander
sample.pars<-function(x){
  k<-sample(nrow(x),nrow(x))
  for(i in k){
    if(c(as.numeric(x[i,4])+as.numeric(x[i,5]))==0){
      next
    } else {
    samp<-do.call(x[i,6],args=list(1,as.numeric(x[i,4]),as.numeric(x[i,5])),quote=F)
    while(samp<=0){
      samp<-do.call(x[i,6],args=list(1,as.numeric(x[i,4]),as.numeric(x[i,5])),quote=F)
    }
    }
  x[i,4:5]<-samp
  }
  return(x)
}

# internal function of the ms.commander
sample.w.cond<-function(par.matrix,cond.matrix){

  nam<-rownames(cond.matrix)
  x<-sample.pars(par.matrix)

  y<-which(cond.matrix=="<", arr.ind=T)

  if(nrow(y)!=0){
    maior<-list(NULL)
    for(i in 1:nrow(y)){
      mm<-NULL
      for(j in 1:2){
        m<-which(par.matrix==nam[y[i,j]])
        mm<-c(mm,m)
      }
      maior[[i]]<-mm
    }

    while(eval.condition(x,y=maior)>0){
      for(j in 1:length(maior)){
        x[c(maior[[j]][1],maior[[j]][2]),]<-sample.pars(par.matrix[c(maior[[j]][1],maior[[j]][2]),])
      }
    }

  }

  z<-which(cond.matrix=="=", arr.ind=T)
  z<-z[order(z[,1]),]
  if(nrow(z)!=0){
    for(i in 1:nrow(z)){
      equal<-NULL
      for(j in 1:2){
        eq<-which(par.matrix==nam[z[i,j]])
        equal<-c(equal,eq)
      }
      x[equal[1],4:5]<-x[equal[2],4:5]
    }
  }
  return(x)
}

# internal function of the Model Builder
eval.condition<-function(x,y){
  value<-NULL
  for(i in 1:length(y)){
    value[i]<-as.numeric(x[y[[i]][1],4])>as.numeric(x[y[[i]][2],4])
  }
  return(sum(value))
}
