
parallel.sims<-function(ncores,code,model)
{

  model.name<-strsplit(code,",")[[1]]
  model.name<-model.name[grep("output.name",model.name)]
  model.name<-gsub(")","",model.name)
  model.name<-gsub(" ","",model.name)
  model.name<-gsub("\"","",model.name)
  model.name<-strsplit(model.name,"=")[[1]][2]

  if(length(list.files()[grep(paste("SIMS_",model.name,".txt",sep=""),list.files())])>0){
    resp<-readline("Simulation file alread exists. Do you want to overwrite the simulations?")

    while (!(resp %in% c("Y","y","Yes","YES","n","N","No","NO"))){
      resp<-readline("Simulation file alread exists. Do you want to overwrite the simulations?")
    }

    if(resp %in% c("Y","y","Yes","YES")){
      file.remove(list.files()[grep(paste("SIMS_",model.name,".txt",sep=""),list.files())])
    } else if (resp %in% c("n","N","No","NO")){
      stop()
    }
  }


  for(i in 1:ncores){

    setwd(getwd())
    #dir.create(paste(tempdir(),"/",i,sep=""),showWarnings = F)
    #setwd(paste("./",i,sep=""))
    dput(model,".model.txt")
    write("suppressMessages(library(PipeMaster))",".script.R")
    write("model<-dget('.model.txt')",".script.R",append = T)
    write(code,".script.R",append=T)
    write('quit("no")',".script.R",append=T)
    system("Rscript .script.R",wait=F)
    cat(paste("Parallel PipeMaster:: starting",i,"core"),sep="\n")
    #setwd("../")
  }
  Sys.sleep(10)
  S<-list.files()[grep("SIMS",list.files())]
  while(length(S)<1){
    S<-list.files()[grep("SIMS",list.files())]
  }
  rows<-nrow(read.table(S,header=T,sep="\t"))

  tot<-strsplit(code,",")[[1]]
  a<-tot[grep("nsim",tot)]
  a<-gsub(")","",a)
  a<-gsub(" ","",a)
  a<-sum(as.numeric(strsplit(a,"=")[[1]][2]))

  b<-tot[grep("size",tot)]
  b<-gsub(")","",b)
  b<-gsub(" ","",b)
  b<-sum(as.numeric(strsplit(b,"=")[[1]][2]))

  tot<-a*b*ncores


  while(rows!=tot){
    Sys.sleep(10)
    if(rows!=nrow(read.table(S,header=T,sep="\t"))){
      rows<-nrow(read.table(S,header=T,sep="\t"))
      cat(paste("Parallel PipeMaster::",rows,"simulations of",tot),sep="\n")
    }
  }
file.remove(c(".script.R",".model.txt"))
}
