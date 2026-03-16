ms.string.generator<-function(model,size.pars,mig.pars,time.pars,use.alpha,scalar=1){

  # rescale to inheritance scalar and transform size parameters to relative to Ne0
  size.pars[,4:5] <- as.numeric(size.pars[,4]) * scalar
  # rescale migration to inheritance scalar
  if(!is.null(mig.pars)) {
    mig.pars[,4:5] <- as.numeric(mig.pars[,4]) * scalar
  }

  # empty string for ms flags
  string<-list()

################### size parameters ############################
################################################################

  curr.Ne <- subset(size.pars, size.pars[,2]=="-n")
  if(!is.null(time.pars)) {
    ent <- subset(time.pars, time.pars[,2]=="-en")
  } else {
    ent <- time.pars
  }
  en <- subset(size.pars, size.pars[,2]=="-en")


  # generate Ne string

    if(nrow(curr.Ne)==1){
      string[[1]]<-paste(curr.Ne[2:4],collapse = " ")
      } else {
        l<-apply(curr.Ne[,c(2:4)],1,paste,collapse=" ")
        string[[1]]<-paste(l,collapse = " ")
      }


  # generate alpha string (exponential growth)
  # Only use regular Ne changes (skip Ne.anc entries)
  if(use.alpha[1]==T){
    en_reg <- en[!grepl("^Ne\\.anc_", en[,1]), , drop = FALSE]
    ent_reg <- if(!is.null(ent)) ent[!grepl("^t\\.Ne\\.anc_", ent[,1]), , drop = FALSE] else ent
    alpha<-NULL
    if(nrow(en_reg) > 0){
      for(i in as.numeric(unique(en_reg[,3]))){
        ne0_row <- which(curr.Ne[,3] == as.character(i))
        ne1_row <- which(en_reg[,3] == as.character(i))[1]
        t_row   <- which(ent_reg[,3] == as.character(i))[1]
        if(length(ne0_row) > 0 && length(ne1_row) > 0 && length(t_row) > 0){
          g_rate <- -(1/as.numeric(ent_reg[t_row,4])) *
                     log(as.numeric(en_reg[ne1_row,4]) / as.numeric(curr.Ne[ne0_row,4]))
          alpha <- c(alpha, paste("-g", i, g_rate))
        }
      }
    }
    if(length(alpha) > 0)
      string[[2]]<-paste(alpha[use.alpha[2:length(use.alpha)]], collapse=" ")
    }

# generate ancestral Ne string
  if(nrow(en)!=0){
    if(nrow(en)>1){
      n<-apply(cbind(ent[,c(2,4)],en[,3:4]),1,paste,collapse=" ")
      string[[3]]<-paste(n, collapse=" ")
    } else {string[[3]]<-paste(c(ent[c(2,4)],en[3:4]), collapse=" ")}
  }

######### migration parameters #########################
####################################################
  if(is.null(mig.pars)==F){

  ###### transform current mig parameters
    curr.mig<-subset(mig.pars, mig.pars[,2]=="-m")
    for(i in 1:nrow(curr.mig)){
      curr.mig[i,3]<-strsplit(curr.mig[i,3]," ")[[1]][1]
      }
    curr.mig[,4] <- as.numeric(curr.mig[,4]) / as.numeric(curr.Ne[match(curr.mig[,3],curr.Ne[,3]),4])
    curr.mig[,3] <- mig.pars[1:nrow(curr.mig),3]

  ###### generate current migration string
    m<-apply(curr.mig[,c(2:4)],1,paste,collapse=" ")
    string[[4]]<-paste(m, collapse=" ")

  ########################################
  ###### ancestral migration conversion ##
    emt<-subset(time.pars, time.pars[,2]=="-em")
    em<-subset(mig.pars, mig.pars[,2]=="-em")

    if(nrow(em)!=0){
      for(i in 1:nrow(emt)){
        emt[i,3]<-strsplit(emt[i,3]," ")[[1]][1]
        }

      if(nrow(en)==0){
        em[,4]<-as.numeric(em[,4])/as.numeric(curr.Ne[match(emt[,3],curr.Ne[,3]),4])
        } else {
          if(sum(as.numeric(em[,4]))>0){
            for(j in 1:nrow(em)){
              x<-which(ent[,3]==emt[j,3])
              if(length(x)==0){
                em[j,4]<-as.numeric(em[j,4])/as.numeric(curr.Ne[match(emt[j,3],curr.Ne[,3]),4])
              } else {
                y<-which(as.numeric(ent[x,4])<=as.numeric(emt[j,4]))
                if(length(y)==0){
                em[j,4]<-as.numeric(em[j,4])/as.numeric(curr.Ne[match(emt[j,3],curr.Ne[,3]),4])
                } else {
                  y<-which(as.numeric(ent[x,4])==max(as.numeric(ent[x[y],4])))[1]
                  em[j,4]<-as.numeric(em[j,4])/as.numeric(en[x[y],4])
                }
            }
          }
        }
      }
    ## generate ancestral migration string
      m<-apply(cbind(emt[,c(2,4)],em[,3:4]),1,paste,collapse=" ")
      string[[5]]<-paste(m, collapse=" ")
    }
  }

############### joint parameters ##################################
###############################################################
########
  if(is.null(time.pars)==F){
    ej<-subset(time.pars, time.pars[,2]=="-ej")
    if(nrow(ej)==1){
      string[[6]]<-paste(ej[c(2,4,3)], collapse=" ")
    } else {
      j<-apply(ej[,c(2,4,3)],1,paste,collapse=" ")
      string[[6]]<-paste(j, collapse=" ")}
    }

# paste strings
  string<-paste(unlist(string),collapse=" ")
  return(string)
}
