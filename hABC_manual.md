# hABC
The hierarchical model implemented here simulates single locus data for a group of populations/species undergoing a demographic change.

### 1) inputs

To simulate the hierarchical model, four input data frames are necessary. These inputs inform the priors to draw the paramenters of the simulartions. You can buld them in a spreed sheet editor, save them as a tab delimited table of csv file and read them into R using the <read.table> function (more instructions bellow). It is important that these inputs keep the same population order and they must have the same number of rows. These inputs are:

1) Ne.prior
2) NaA.prior
3) time.prior
4) gene.prior

#### 1.1) Ne.prior
A four column data frame (table) containing information about the effective population size prior. The first column indicates the name of the species/population. The secound column indicates the distribution of the prior. The third column indicates the first parameter of the prior and the fourth column indicates the secound parameter of the prior distribution. For a uniform prior these will be minumum and maximum values for the Ne. For a normal prior these will be the mean and standard deviation of the normal distribution.

To see an example of the Ne.prior input that was used in Gehara et al 2017 run:
> library(PipeMaster)

> data("hABC.priors")

> hABC.priors$Ne.prior


```
          species distribution   minNe   maxNe
1        C_ocel_1        runif 1500000 4000000
2        C_ocel_4        runif  900000 3500000
3        D_muel_1        runif  120000 1200000
4        D_muel_2        runif  700000 1800000
5      L_klugei_1        runif 1300000 3700000
6     P_nattereri        runif 1900000 4100000
7  P_pollicaris_1        runif  140000 1050000
8  P_pollicaris_7        runif  460000 1480000
9  Micra_concat_1        runif  870000 2400000
10 V_rubri_cytb_2        runif  480000 1900000
11 P_nordestina_2        runif  800000 4850000
12 P_nordestina_4        runif  860000 4470000
13       Rhinella        runif  290000 1600000
14  P_diplolister        runif 4000000 7020000
```

#### 1.2) NeA.prior

A four column data frame (table) containing information about the ancestral effective population size prior. The first column indicates the name of the species/population. The secound column indicates the distribution of the prior. The third column indicates the first parameter of the prior and the fourth column indicates the secound parameter of the prior distribution. The ancestral population is implemented as a multiplyer of the current population size. For instance if I define the ancestral Ne prior to be uniform, between 0.01 - 0.1, this will mean that the current Ne will be multiplyed by a value sampled from this distrubution. If the current Ne is 100,000 and the ancestral value sampled was 0.1, the ancestral Ne will be: 100,000 x 0.1 = 10,000.

To see an example of the NeA.prior input that was used in Gehara et al 2017 run:
> library(PipeMaster)

> data("hABC.priors")

> hABC.priors$NeA.prior

```
          species distribution maxexp minexp
1        C_ocel_1        runif  0.001    0.1
2        C_ocel_4        runif  0.001    0.1
3        D_muel_1        runif  0.001    0.1
4        D_muel_2        runif  0.001    0.1
5      L_klugei_1        runif  0.001    0.1
6     P_nattereri        runif  0.001    0.1
7  P_pollicaris_1        runif  0.001    0.1
8  P_pollicaris_7        runif  0.001    0.1
9  Micra_concat_1        runif  0.001    0.1
10 V_rubri_cytb_2        runif  0.001    0.1
11 P_nordestina_2        runif  0.001    0.1
12 P_nordestina_4        runif  0.001    0.1
13       Rhinella        runif  0.001    0.1
14  P_diplolister        runif  0.001    0.1
```

#### 1.3) time.priors

A five column data frame (table) containing information about the expansion time priors for each population/species. The first column indicates the name of the species/population. The secound column indicates the distribution of the prior. The third column indicates the first parameter of the prior and the fourth column indicates the secound parameter of the prior distribution. For uniform distributions these will be the min and max expantion times in years. The fith column specifies the generation length of the species. 

To see an example of the time.prior input that was used in Gehara et al 2017 run:
> library(PipeMaster)

> data("hABC.priors")

> hABC.priors$time.prior

```
     species distribution   min    max generationtime
1   A_ocel_1        runif 20000 1000000             2
2   A_ocel_2        runif 20000 1000000             2
3   D_muel_1        runif 20000 1000000             1
4   D_muel_2        runif 20000 1000000             1
5   L_klugei        runif 20000 1000000             1
6  P_nattere        runif 20000 1000000             2
7  P_polli_1        runif 20000 1000000             1
8  P_polli_2        runif 20000 1000000             1
9  M_maximil        runif 20000 1000000             1
10 V_multisc        runif 20000 1000000             1
11 P_norde_1        runif 20000 1000000             1
12 P_norde_2        runif 20000 1000000             1
13    R_jimi        runif 20000 1000000             1
14 P_diploli        runif 20000 1000000             1
```

#### 1.4) gene.prior

A seven column data frame (table) containing information about the loci to be simulated for each population, including mutation rate, number of base pairs, number of samples and inheritance scalar. The first column indicates the name of the species/population. The secound column indicates the distribution of the prior. The third column indicates the first parameter of the prior, the fourth column indicates the secound parameter of the prior distribution. For a uniform distribution the first parameter respresents the mean and the secound paramter the standard deviation. The fith column indicates the number of base pairs to simulate. The sixth column indicates the number of samples to simulate. And the Seventh column indicates the inheritance scalar. Mitochondrial DNA have an inheritance factor of 0.25 since it is maternaly inherited and haployd. Autosomal loci have an inheritance of 1.

To see an example of the gene.prior input that was used in Gehara et al 2017 run:
> library(PipeMaster)

> data("hABC.priors")

> hABC.priors$gene.prior

```
    species distribution   min     max   bp samples inheritance
1   A_ocel_1        rnorm 1e-08 1.5e-09  370      58        0.25
2   A_ocel_2        rnorm 1e-08 1.5e-09  370      31        0.25
3   D_muel_1        rnorm 1e-08 1.5e-09  514      16        0.25
4   D_muel_2        rnorm 1e-08 1.5e-09  514      32        0.25
5   L_klugei        rnorm 1e-08 1.5e-09  679      45        0.25
6  P_nattere        rnorm 1e-08 1.5e-09  686      58        0.25
7  P_polli_1        rnorm 1e-08 1.5e-09  942      17        0.25
8  P_polli_2        rnorm 1e-08 1.5e-09  942      46        0.25
9  M_maximil        rnorm 1e-08 1.5e-09 1524      31        0.25
10 V_multisc        rnorm 1e-08 1.5e-09  735      33        0.25
11 P_norde_1        rnorm 1e-08 1.5e-09  527      12        0.25
12 P_norde_2        rnorm 1e-08 1.5e-09  527      14        0.25
13    R_jimi        rnorm 1e-08 1.5e-09  533      31        0.25
14 P_diploli        rnorm 1e-08 1.5e-09  601     166        0.25
```

You can create all prior tables using Excel, Libre office or a text editor. They can also be created within R and saved as a text file. If you use Excel you should save the table as text. To read the table in R you should use the read.table function. Your .txt file should be in your working directory or you should specify the path to your table.

> Ne_prior <- read.table("myNepriortable.txt", header = T, sep = "\t")

or

> Ne_prior <- read.table("path/myNepriortable.txt", header = T, sep = "\t")


### 2) Hierarchical model simulation

After setting up the inputs you can run simulations using two functions: sim.coexp and sim.coexpPT. The sim.coexp function simulates the implementation of Chan et al 2014, Chan et al with threshold and the narow coexpantion time used in Gehara et al 2017. Explanations about the specifics of each implementation can be found in Gehara et al 2017. The sim.coexpPT simulates the partitioned time model used in Gehara et al 2017. Here is a link to the paper: https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14239

To simulate data used in Gehara et al 2017 run:

> library(PipeMaster)

> data("hABC.priors")

> sim.coexp(nsims = 1000, var.zeta="FREE", th=50000, coexp.prior = c(20000,1000000), 
          Ne.prior = hABC.priors$Ne.prior, alpha=F, NeA.prior = hABC.priors$NeA.prior, 
          time.prior = hABC.priors$time.prior, gene.prior = hABC.priors$gene.prior,
          append.sims = F, path = getwd())

The arguments of this function are:
```
# Arguments

nsims        Total number of simulations

var.zeta     Variation on zeta parameter. Can be "FREE" to vary or be set to a specific value (between 0-1).

coexp.prior  Uniform prior for the coespansion time. Vector of two numbers with the lower and upper boudary of the prior.

th           Threshold. Minimum time difference between Ts, time of simultaneous change and population specific times.

Ne.prior     Data frame with the prior values for the Ne of each population.

NeA.prior	   Data frame with the prior values for the ancestral Ne of each population.

time.prior   Data frame with parameter values for the priors of the time of demographic change of each population.

gene.prior   Data frame with parameter values for the priors of the mutation rate of each species.

alpha        logical. If TRUE all demographic chages are exponential. If FALSE sudden changes. Defaut is FALSE.

append.sims  logical. If TRUE simulations are appended to the simulations file. Defaut is FALSE.

path         Path to the directiry to write the simulations. Defaut is the working directory.
```
run

> ?sim.coexp 

for more details about the function.

The output of this function is a tab delimited table, saved in your working directory or specified path, containing the parameters in the first four columns and the hiper summary statistics in the remainning 16 columns.

While running the simulations the R console shows the number of simulations finished, the total number of simulations and the zeta value being simulated.

```
[1] "1 sims of 1000 | zeta =  0.357142857142857"
[1] "2 sims of 1000 | zeta =  1"
[1] "3 sims of 1000 | zeta =  0.214285714285714"
[1] "4 sims of 1000 | zeta =  0.571428571428571"
[1] "5 sims of 1000 | zeta =  0.928571428571428"
[1] "6 sims of 1000 | zeta =  0.714285714285714"
[1] "7 sims of 1000 | zeta =  0.285714285714286"
...
```

### 3) Calculating the the observed hipersummary statistics (hss) and estimating parameters.

The most simple ABC algorithm, the rejection algorithm, calculates the euclidian distances between each simulated data and your empirical data and rejects all simulated datasets that are too distant from the empirical data. Hopefully, the retained simulations will have information on the posterior probability of the model parameters.
Check the first 60 secounds of this video for a lightning explanation of ABC: https://www.youtube.com/watch?v=EUCl4v_NIRs

The euclidian distance threshold (or tolerace) is arbitrary but a cross-validation experiment should help you decide what tolerance level to use. See the abc package vignette for a good example of how to perform an ABC analisys in R. 
https://cran.r-project.org/web/packages/abc/vignettes/abcvignette.pdf

#### These references have good examples and explanations of ABC:
Csilléry K., Blum M.G.B., Gaggiotti O.E., & François O. (2010) Approximate Bayesian Computation (ABC) in practice. Trends in Ecology and Evolution, 25, 410–418. 

Fagundes N.J., Ray N., Beaumont M., Neuenschwander S., Salzano F.M., Bonatto S.L., & Excoffier L. (2007) Statistical evaluation of alternative models of human evolution. Proc Natl Acad Sci U S A, 104, 17614–17619.

Beaumont M.A. (2011) Approximate Bayesian Computation in Evolution and Ecology. Annual Review of Ecology, Evolution, and Systematics, 41, 379–406.

#### Another emerging method in population genetics is supervised machine learning (SML):

Schrider D.R. & Kern A.D. (2018) Supervised Machine Learning for Population Genetics: A New Paradigm. Trends in Genetics, xx, 1–12.

Sheehan S. & Song Y.S. (2016) Deep Learning for Population Genetic Inference. PLoS Computational Biology, 12,


To perform the ABC analysis or SML you need to calculate the observed summary stats.
To calculate the hss for your empirical data you should place aligned fasta files in a folder and rund the following function. 

> observed <- observed.coexp.sumstat(path.to.fasta = "path to the folder where the fastas are")

You will also need to read the simulations table generated by the sim.coexp function back to R. If you have many simulations it is advisable to use the bigmemory R package, which helps handling large amounts of data in R.

> install.packages("bigmemory")

> library(bigmemory)

> simulated <- read.big.matrix(file="simulations.txt", header=T, type="float", sep="\t")


### To estimate the parameters with abc you can use the abc R-package

> install.packages("abc")

> library(abc)

> abc.rej <- abc(target = observed, param = simulated[,1:4], sumstat = simulated[,5:20],
               tol = 0.001, method = "rejection")

> summary(abc.rej)
               
Check the abc R-package vignette for more details on this function.
For SML check the caret R-package. It has a vast number of machine learning algorithms and a good online manual.

### 4) Citation

#### If you use the hABC functions please cite:
Gehara M., Garda A.A., Werneck F.P., Oliveira E.F., da Fonseca E.M., Camurugi F., Magalhães F. de M., Lanna F.M., Sites J.W., Marques R., Silveira-Filho R., São Pedro V.A., Colli G.R., Costa G.C., & Burbrink F.T. (2017) Estimating synchronous demographic changes across populations using hABC and its application for a herpetological community from northeastern Brazil. Molecular Ecology, 26, 4756–4771.

Chan Y.L., Schanzenbach D., & Hickerson M.J. (2014) Detecting concerted demographic response across community assemblages using hierarchical approximate Bayesian computation. Molecular Biology and Evolution, 31, 2501–2515. 

#### PipeMaster uses ms to simulate the data and pegas to calculate the summaries: 
Hudson R.R. (2002) Generating samples under a Wright-Fisher neutral model of genetic variation. Bioinformatics, 18, 337–338. 

Paradis E. 2010. pegas: an R package for population genetics with an integrated-modular approach. Bioinformatics 26: 419-420. 

#### If you use the abc package to estimate parameters you should also cite:
Csilléry K., François O., & Blum M.G.B. (2012) Abc: An R package for approximate Bayesian computation (ABC). Methods in Ecology and Evolution, 3, 475–479. 
