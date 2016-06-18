cur.Ne.par<-function(){
 list.Ne.pars<-NULL
  for (i in 1:.e$npops){
    Ne0.par<-paste("Ne0.pop",i,sep="")
    list.Ne.pars<-c(list.Ne.pars,Ne0.par)
  }
  .e$n<-matrix(nrow=length(list.Ne.pars), ncol=6) 
  .e$n[,1]<-list.Ne.pars
  .e$n[,2]<-'-n'
  .e$n[,3]<-c(1:length(list.Ne.pars))
  .e$n[,4]<-100000
  .e$n[,5]<-500000
  .e$n[,6]<-"uniform"
  
}
