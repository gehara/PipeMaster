-e 
#5

	 *** objective: Simulate 10 times, 6 (independent) loci from a sample of 12 lines.

The length of the loci is variable as well.


The file fragments.txt describes the configuration of the loci:
id	n	length
f1	12	1000
f2	12	500
f3	12	100
f4	12	2200
f5	12	1000
f6	12	100

Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt.

Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.



Parameters:
------------
p_theta : the value of parameter theta used in the simulations. It is constant.
p_rho : the value of parameter rho. It is drawn from U(10; 20)
p_npop : the number of sub-populations
p_isl_totmig : the total migration rate 
p_pop_change_time : time at which (total) population size changes
p_pop_change_newpop : the ratio of the changed population to the population at the present

Summaries
------------
s_average_segs_1 : the average number of segregating sites for sub-population 1
s_variance_segs_1 : the variance of the segregating sites for sub-population 1
s_average_segs_2 : the average of segregating sites for sub-population 2
s_variance_segs_2 : the variance of segregating sites for sub-population 2
s_average_segs : the average of the segregating sites for the total sample
s_variance_segs : the variance of the segregating sites for the total sample
s_average_pi_1 : the average pi for sub-population 1
s_variance_pi_1 : the variance of pi for sub-population 1
s_average_pi_2 : the average pi for sub-population 2
s_variance_pi_2 : the variance of pi for sub-population 2
s_average_pi : the average of pi for the total sample
s_variance_pi : the variance of pi for the total sample
s_average_w_1 : the average of theta for sub-population 1
s_variance_w_1 : the variance of theta for sub-population 1
s_average_w_2 : the average of theta for sub-population 2
s_variance_w_2 : the variance of the theta for sub-population 2
s_average_w : the average of theta for the total sample
s_variance_w : the variance of theta for the total sample
s_average_tajd_1 : the average Tajima's D for sub-population 1
s_variance_tajd_1 : the variance of Tajima's D for sub-population 1
s_average_tajd_2 : the average of Tajima's D for sub-population 2
s_variance_tajd_2 : the variance of Tajima's D for sub-population 2
s_average_tajd : the average of Tajima's D for the total sample
s_variance_tajd : the variance of Tajima's D for the total sample
s_average_ZnS_1 : the average of Zns for sub-population 1
s_variance_ZnS_1 : the variance of ZnS for sub-population 1
s_average_ZnS_2 : the average of ZnS for sub-population 2
s_variance_ZnS_2 : the variance of ZnS for sub-population 2
s_average_ZnS : the average of ZnS for the total sample
s_variance_ZnS : the variance of ZnS for the total sample
s_average_Fst : the average Fst (total sample, hbk calculation)
s_variance_Fst : the variance of Fst (total sample, hbk calculation)
s_average_shared_1_2 : the average percentage of shared polymorphisms between sub-populations 1 and 2
s_variance_shared_1_2 : the variance of the percentage of shared polymorphisms between sub-populations 1 and 2
s_average_private_1_2 : the average percentage of private polymorphisms between sub-populations 1 and 2
s_variance_private_1_2 : the variance of percentage of private polymorphisms between sub-populations 1 and 2
s_average_fixed_dif_1_2 : the average percentage of fixed differences between sub-populations 1 and 2
s_variance_fixed_dif_1_2 : the variance of percentage of fixed differences between sub-populations 1 and 2
s_average_pairwise_fst_1_2 : the average Fst between sub-populations 1 and 2
s_variance_pairwise_fst_1_2 : the variance of Fst between sub-populations 1 and 2
s_average_fwh_1 : the average H in sub-population 1
s_variance_fwh_1 : the variance of H in sub-population 1
s_average_fwh_2 : the average H in sub-population 2
s_variance_fwh_2 : the variance of H in sub-population 2
s_average_FayWuH : the average H in the total sample
s_variance_FayWuH : the variance of H in the total sample

***********************************************
IMPORTANT: theta AND rho values in the command line refer to a locus of length 1000 (as it is specified by the nsites).
----------
This number is rescaled for each of the loci (since their length may vary).
***********************************************

NOTICE!! 

1. When pop == 0, then the n should be equal to the total n given in the 
command line. Also, the sample configuration is the one given in the command
line.

2. When pop != 0, then information for ALL subpopulations is required. 
The order is important. It is supposed that first information is given
for sub-population 1, then for sub-population 2 etc. See the next example.
Then, the total sample size is the sum of the samples of the demes. 


