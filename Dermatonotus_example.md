---
output: html_document
---
# Simulations and analyses with PipeMaster (nextgen example)

This is a script showing how to simulate data, test model and estimate parameters using [PipeMaster](https://github/gehara/PipeMaster), abc and caret packages.
The data set used is the same as the one used in Gehara et al (in prep), and represents 2177 UCE loci for the neotropical frog Dermatonotus muelleri. For more information about this species see Gehara et al (in prep). and [Oliveira et al. 2018](https://www.researchgate.net/profile/Adrian_Garda/publication/327624820_Phylogeography_of_Muller%27s_termite_frog_suggests_the_vicariant_role_of_the_Central_Brazilian_Plateau/links/5c40f99f92851c22a37d572c/Phylogeography-of-Mullers-termite-frog-suggests-the-vicariant-role-of-the-Central-Brazilian-Plateau.pdf)


### 1) Install and load all necessary packages

#### 1.1 instalation

> install.packages(c("devtools","ggplot2","abc","caret","doMC"))

> library(devtools) # devtools: necessary for PipeMaster instalation

> install_github("gehara/PipeMaster@developing")

#### 1.2 load packages

> library(devtools)

> library(PipeMaster) # PipeMaster: used to simulate data and some additional tools

> library(abc) # abc: used to perform approximate Bayesian computation (ABC)

> library(caret) # caret: used to perform the superevised machine-learning (SML)

> library(doMC) # doMC: necessary to run the SML in parallel

> library(ggplot2) # ggplot2: used to plot PCA

### 2) Create a working directory to save results

##### get the working directory
> path <- getwd()

##### create a new directory to save outputs
> dir.create(paste(path,"/PM_example",sep=""))

##### set working directory
> setwd(paste(path,"/PM_example",sep=""))

### 3) Load example data
This example data is based on Gehara et al. 

#### observed summary statistics
> data("observed_Dermatonotus", package = "PipeMaster")
  
##### models used in Gehara et al

> data("models", package="PipeMaster")

### 4) Usefull tips and tools

##### *dput* is usefull to save the loaded models in the working directory as .txt files

> dput(Is,"Is.txt")

##### You can use *dget* to retrieve models from .txt files saved with *dput*

> Is <- dget("Is.txt")

#### With the function bellow you can see the parameters and priors of the models

> tab <- get.prior.table(model=Is)

> tab

```
   Parameter prior.1 prior.2 distribution
1   Ne0.pop1  100000 5000000        runif
2   Ne0.pop2  100000 5000000        runif
3   Ne1.pop2  100000 5000000        runif
4 t.Ne1.pop2 1500000 8000000        runif
5      join1 1500000 8000000        runif
```

#### Visualize prior distributions

```
plot.priors(Is)
```
![Prior distribution]https://github.com/gehara/PipeMaster/priors.png

####4.5) it is possible to update the priors using the table

```
tab <- get.prior.table(model=Is)
tab[,2:3]<-1000
tab
Is2 <- PipeMaster:::update.priors(tab = tab, model = Is)
get.prior.table(model=Is2)
```

###5) Simulate data

```
sim.msABC.sumstat(Is, nsim.blocks = 1, use.alpha = F, output.name = "Is", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IM, nsim.blocks = 1, use.alpha = F, output.name = "IM", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IsD, nsim.blocks = 1, use.alpha = c(T,1,2), output.name = "IsD", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IMD, nsim.blocks = 1, use.alpha = c(T,1,2), output.name = "IMD", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IsBot, nsim.blocks = 1, use.alpha = c(T,1,2), output.name = "IsBot", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IMBot, nsim.blocks = 1, use.alpha = c(T,1,2), output.name = "IMBot", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IsBot2, nsim.blocks = 1, use.alpha = c(T,1), output.name = "IsBot2", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IMBot2, nsim.blocks = 1, use.alpha = c(T,1), output.name = "IMBot2", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IsD2, nsim.blocks = 1, use.alpha = c(T,1), output.name = "IsD2", append.sims = F, block.size = 500, ncores = 2)
sim.msABC.sumstat(IMD2, nsim.blocks = 1, use.alpha = c(T,1), output.name = "IMD2", append.sims = F, block.size = 500, ncores = 2)
```
###6) Read simulations
```{r}
Is.sim <- read.table("SIMS_Is.txt", header=T)
IM.sim <- read.table("SIMS_IM.txt", header=T)
IsD.sim <- read.table("SIMS_IsD.txt", header=T)
IMD.sim <- read.table("SIMS_IMD.txt", header=T)
IsD2.sim <- read.table("SIMS_IsD2.txt", header=T)
IMD2.sim <- read.table("SIMS_IMD2.txt", header=T)
IsBot.sim <- read.table("SIMS_IsBot.txt", header=T)
IMBot.sim <- read.table("SIMS_IMBot.txt", header=T)
IsBot2.sim <- read.table("SIMS_IsBot2.txt", header=T)
IMBot2.sim <- read.table("SIMS_IMBot2.txt", header=T)
```
###7) See and select summary stats
```{r}
# see sumstats names
colnames(observed)
# select summary stats to be excluded using grep
cols <- c(grep("thomson", names(observed)),
          grep("pairwise_fst", names(observed)),
          grep("Fay", names(observed)),
          grep("fwh", names(observed)),
          grep("_dv", names(observed)),
          grep("_s_", names(observed)))
# exclude
observed <- observed[,-cols]
# visualize sumary stats names aggain
colnames(observed)
```
###8) combine simulations in a single matrix matching by summary stats names in the observed and run a Principal Components Analyses to visualize model-fit
```{r, fig.width = 10, fig.height = 10}
# models without migration
models <- rbind(Is.sim[,colnames(Is.sim) %in% colnames(observed)],
                IsD.sim[,colnames(IsD.sim) %in% colnames(observed)],
                IsD2.sim[,colnames(IsD2.sim) %in% colnames(observed)],
                IsBot.sim[,colnames(IsBot.sim) %in% colnames(observed)],
                IsBot2.sim[,colnames(IsBot2.sim) %in% colnames(observed)])

# create an index for the models
data <- c(rep("Is", nrow(Is.sim)),
         rep("IsD", nrow(IsD.sim)),
         rep("IsD2", nrow(IsD2.sim)),
         rep("IsBot", nrow(IsBot.sim)),
         rep("IsBot2", nrow(IsBot2.sim)))


# exclude missing data for pca plot
data.PCA <- data[complete.cases(models)]
models.PCA <- models[complete.cases(models),]

 subsample 1000 simulations for PCA
x <- sample(1:length(data.PCA),1000)
data.PCA <- data.PCA[x]
models.PCA <- models.PCA[x,]

# run PCA
PCA <- prcomp(rbind(models.PCA,observed), center = T, scale. = T, retx=T)
# get scores
scores <- data.frame(PCA$x[,1:ncol(PCA$x)])
# multiplot
multiplot(ggplot(scores, aes(x=PC1, y=PC2))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"),   shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            theme(legend.position="top", legend.direction="horizontal",legend.title = element_blank())+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC3))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC4))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC5))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC6))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC7))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC8))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC9))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC10))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),cols = 3)

# models with migration
models<-rbind(IM.sim[,colnames(IM.sim) %in% colnames(observed)],
              IMD.sim[,colnames(IMD.sim) %in% colnames(observed)],
              IMD2.sim[,colnames(IMD2.sim) %in% colnames(observed)],
              IMBot.sim[,colnames(IMBot.sim) %in% colnames(observed)],
              IMBot2.sim[,colnames(IMBot2.sim) %in% colnames(observed)])

data<-c(rep("IM",nrow(IM.sim)),
        rep("IMD",nrow(IMD.sim)),
        rep("IMD2",nrow(IMD2.sim)),
        rep("IMBot",nrow(IMBot.sim)),
        rep("IMBot2",nrow(IMBot2.sim)))



data.PCA<-data[complete.cases(models)]
models.PCA<-models[complete.cases(models),]

x<-sample(1:length(data.PCA),1000)
data.PCA <- data.PCA[x]
models.PCA <- models.PCA[x,]

PCA<-prcomp(rbind(models.PCA,observed), center = T, scale. = T, retx=T)
scores <- data.frame(PCA$x[,1:ncol(PCA$x)])
#pdf("PCA13.pdf", paper="a4r", width=10, pointsize=10)
multiplot(ggplot(scores, aes(x=PC1, y=PC2))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            theme(legend.position="top", legend.direction="horizontal",legend.title = element_blank())+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC3))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC4))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC5))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC6))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC7))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC8))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC9))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),
          ggplot(scores, aes(x=PC1, y=PC10))+
            theme(legend.position = "none")+
            geom_point(aes(colour=c(data.PCA,"observed"), size=c(data.PCA,"observed"), shape=c(data.PCA,"observed")))+
            scale_size_manual(values=c(2,2,2,2,2,10))+
            scale_colour_brewer(palette="Set1"),cols = 3)

```
###9) Run a supervised machine learning analysis for model classification
```{R}
# combine all models for analysis
models<-rbind(Is.sim[,colnames(Is.sim) %in% colnames(observed)],
              IsD.sim[,colnames(IsD.sim) %in% colnames(observed)],
              IsD2.sim[,colnames(IsD2.sim) %in% colnames(observed)],
              IsBot.sim[,colnames(IsBot.sim) %in% colnames(observed)],
              IsBot2.sim[,colnames(IsBot2.sim) %in% colnames(observed)],
              IM.sim[,colnames(IM.sim) %in% colnames(observed)],
              IMD.sim[,colnames(IMD.sim) %in% colnames(observed)],
              IMD2.sim[,colnames(IMD2.sim) %in% colnames(observed)],
              IMBot.sim[,colnames(IMBot.sim) %in% colnames(observed)],
              IMBot2.sim[,colnames(IMBot2.sim) %in% colnames(observed)])

# generate index
data<-c(rep("Is", nrow(Is.sim)),
        rep("IsD", nrow(IsD.sim)),
        rep("IsD2", nrow(IsD2.sim)),
        rep("IsBot", nrow(IsBot.sim)),
        rep("IsBot2", nrow(IsBot2.sim)),
        rep("IM",nrow(IM.sim)),
        rep("IMD",nrow(IMD.sim)),
        rep("IMD2",nrow(IMD2.sim)),
        rep("IMBot",nrow(IMBot.sim)),
        rep("IMBot2",nrow(IMBot2.sim)))


# set up number of cores for SML
registerDoMC(2)
## combine simulations and index
models <- cbind(models,data)
## setup the outcome (name of the models, cathegories)
outcomeName <- 'data'
## set up predictors (summary statistics)
predictorsNames <- names(models)[names(models) != outcomeName]
## split the data into trainning and testing sets; 75% for trainning, 25% testing
splitIndex <- createDataPartition(models[,outcomeName], p = 0.75, list = FALSE, times = 1)
train <- models[ splitIndex,]
test  <- models[-splitIndex,]
## bootstraps and other controls
objControl <- trainControl(method='boot', number=10, returnResamp='final',
                           classProbs = TRUE)
## train the algoritm
nnetModel_select <- train(train[,predictorsNames], train[,outcomeName],
                          method="nnet",maxit=5000,
                          trControl=objControl,
                          metric = "Accuracy",
                          preProc = c("center", "scale"))
## predict outcome for testing data, classification
predictions <- predict(object=nnetModel_select, test[,predictorsNames], type='raw')
## calculate accuracy in model classification
accu <- postResample(pred=predictions, obs=as.factor(test[,outcomeName]))
## predict probabilities of the models for the observe data
pred <- predict(object=nnetModel_select, observed, type='prob')

# visualize results
t(c(pred,accu))
# write results to file
write.table(c(pred,accu),"results.selection.txt")
```
###10) estimate parameters using abc and neural networks
```{r}
# read selected model
IsBot2.sim <- read.table("SIMS_IsBot2.txt", header=T)

# separate summary statistics from parameters
sims <- IsBot2.sim[,colnames(IsBot2.sim) %in% colnames(observed)]
param <- IsBot2.sim[,1:11]

# estimate posterior distribution of parameters
post <- abc(target = observed,
            param = param,
            sumstat = sims,
            sizenet = 30,
            method = "neuralnet",
            MaxNWts = 5000,
            tol=0.1) # adjust tolerance level according to the total number of simulations.

# write results to file
write.table(summary(post), "parameters_est.txt")

# plot posterior probabilities
#plot(post, param = param)

# cross-validation for parameter estimates
cv3 <- cv4abc(param = param,
              sumstat = sims,
              nval=10,
              sizenet = 20,
              method = "neuralnet",
              MaxNWts = 5000,
              tol = 0.1)
plot(cv3)



