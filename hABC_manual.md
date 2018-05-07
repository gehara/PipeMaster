# hABC
The hierarchical model implemented here simulates single locus data for a group of populations/species undergoing a demographic change.

### 1) inputs

To simulate the hierarchical model four input data frames are necessary. These inputs inform the priors to draw the paramenters of the simulartions. You can buld them in a spreed sheet editor, save them as a tab delimited table of csv file and read them into R using the <read.table> function (more instructions bellow). The inputs are:

1) Ne.prior
2) NaA.prior
3) time.prior
4) gene.prior

#### 1.1) Ne.prior
A four column data frame (table) containing information about the effective population size prior. The first column indicates the name of the species/population. The secound column indicates the distribution of the prior. The third column indicates the first parameter of the prior and the fourth column indicates the secound parameter of the prior distribution. For a uniform prior these will be minumum and maximum values for the Ne. For a normal prior these will be the mean and standard deviation of the normal distribution.

To see an example of the input which was used in Gehara et al 2017 run:
> library(PipeMaster)
>
> data("Ne.prior")
>
>Ne.prior

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

To see an example of the input which was used in Gehara et al 2017 run:
> library(PipeMaster)
>
> data("NeA.prior")
>
> NeA.prior

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
