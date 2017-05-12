BN.sumstat<-function(ms.output,gene.prior){
sumstat<-NULL

  for (j in 1:length(ms.output)){
    
    x<-ms.to.DNAbin(ms.output = ms.output[[j]],bp.length = as.numeric(gene.prior[1,5]))
    pi<-nuc.div(x)
    TD<-tajima.test(x)$D
    H<-H.div(x)

    sumstat<-c(pi[[1]],H,TD[1])
    }
return(sumstat)
}
