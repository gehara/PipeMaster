#' Observed summary statistics calculation for single pop demography
#' @description This function calculates the summary statistics for each fasta file present in a particular folder. This calculates the observed for simulations with the single.pop.demog function.
#' @param path.to.fasta The path to the directory containing the fasta files for all species/populations.
#' @author Marcelo Gehara
#' @return A table of summary statistics.
#' @references Gehara M., Garda A.A., Werneck F.P., Oliveira E.F., da Fonseca E.M., Camurugi F., Magalhães F. de M., Lanna F.M., Sites J.W., Marques R., Silveira-Filho R., São Pedro V.A., Colli G.R., Costa G.C., & Burbrink F.T. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756–4771.
#' @export
observed.singlepop.sumstat <- function(path.to.fasta){

  setwd(path.to.fasta)
  fasta.files<-list.files()
  fasta.files<-fasta.files[grep(".fas",fasta.files)]

  ms.output<-fasta2ms(path.to.fasta,fasta.files,write.file=F)
  bp.length<-list()
  ss<-list()
  for(i in 1:length(ms.output)){
    fas<-read.dna(file=fasta.files[i], format="fasta")
    bp.length[[i]]<-ncol(fas)
    ss[[i]]<-as.numeric(strsplit(ms.output[[i]][3]," ")[[1]][2])
    if(is.null(ss[[i]])){
      stop(paste("sequence",fasta.files[i],"has zero variation"))
    }

  }

  sum.stat<-NULL
  for (j in 1:length(ms.output)){
    x<-ms.to.DNAbin(ms.output = ms.output[[j]],bp.length = bp.length[[j]])

    pi<-nuc.div(x)
    H<-H.div(x)[2]
    TD<-tajima.test(x)$D
    #R2<-R2.test(x,B=0,plot = F,quiet = T)
    spec<-site.spectrum(x)[1:3]
    SS<-c(pi,ss[[j]],H,TD,spec)

    sum.stat<-rbind(sum.stat,SS)
    }
rownames(sum.stat)<-fasta.files
colnames(sum.stat)<-c("pi","ss","H","TD","ss1","ss2","ss3")
return(sum.stat)
}
