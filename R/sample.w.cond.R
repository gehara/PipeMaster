#' internal function of the ms.commander
#'
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
