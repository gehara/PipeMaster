
parallel.sims<-function(ncores,code,model)
{

  model.name<-strsplit(code,",")[[1]]
  model.name<-model.name[grep("output.name",model.name)]
  if(length(model.name)==0){
    stop("Argument model.name needs to be specified in the code.")
  }
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

  tot<-strsplit(code,",")[[1]]
  nsim.blocks<-tot[grep("nsim.blocks",tot)]
  if(length(nsim.blocks)==0){
    stop("Argument nsim.blocks needs to be specified in the code.")
  }
  nsim.blocks<-gsub(")","",nsim.blocks)
  nsim.blocks<-gsub(" ","",nsim.blocks)
  nsim.blocks<-sum(as.numeric(strsplit(nsim.blocks,"=")[[1]][2]))

  block.size<-tot[grep("block.size",tot)]
  if(length(block.size)==0){
    block.size<-100
  } else {
    block.size<-gsub(")","",block.size)
    block.size<-gsub(" ","",block.size)
    block.size<-sum(as.numeric(strsplit(block.size,"=")[[1]][2]))
  }
  tot<-nsim.blocks*block.size*ncores

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
  S<-list.files()[grep(paste("SIMS_",model.name,sep=""),list.files())]
  while(length(S)<1){
    S<-list.files()[grep(paste("SIMS_",model.name,sep=""),list.files())]
  }
  rows<-nrow(read.table(S,header=T,sep="\t"))




  while(rows!=tot){
    Sys.sleep(10)
    if(rows!=nrow(read.table(S,header=T,sep="\t"))){
      rows<-nrow(read.table(S,header=T,sep="\t"))
      cat(paste("Parallel PipeMaster::",rows,"simulations of",tot),sep="\n")
    }
  }
file.remove(c(".script.R",".model.txt"))
}
