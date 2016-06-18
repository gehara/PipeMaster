inv.mirror.lower<-function(x) {
  x1<-t(x)[lower.tri(x, diag=F)]
  for(i in 1:length(x1)) {
    if(!x1[i] %in% c("<",">")) next
    if(x1[i]=="<") {
      x1[i]<-">"
      next}
    else x1[i]<-"<"
  }
  x[lower.tri(x, diag=F)]<-x1
  return(x)
}

inv.mirror.upper<-function(x) {
  x1<-t(x)[upper.tri(x, diag=F)]
  for(i in 1:length(x1)) {
    if(!x1[i] %in% c("<",">")) next
    if(x1[i]=="<") {
      x1[i]<-">"
      next}
    else x1[i]<-"<"
  }
  x[upper.tri(x, diag=F)]<-x1
  return(x)
}
