-e 
#6 
	 *** objective: Simulate 10 times, 2 (independent) loci. Sample size
is variable among loci. 

The samples originate from two populations. For the first locus 24 lines are from each population.
For the second locus 24 lines are derived from the first population and 23 lines from the second.

Data includes incomplete information, i.e. in the observed fasta alignment there are 'N's.

workflow
--------
1. use the script missing_ms_report.pl to obtain the information about the missing data. 
	 this script will transform in continuous form the missing positions.
	 this script should run for each of the loci:

	 missing_ms_report.pl in=locus1.fa > locus1.missing
	 missing_ms_report.pl in=locus2.fa > locus2.missing

Then, we create a list that contains the files that describe the missing information:
	 echo "locus1.missing" > missing.list
	 echo "locus2.missing" >> missing.list

the file missing.list will be used in the command line of msABC after the flag --missing


Polymorphism data are stored in tempoutput_ms.txt.
NOTICE!!!! : the polymorphism data contain except 0 and 1, the state 2, which is the missing information.


The parameter values with the summary statistics are stored in out.txt.

Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.

Parameters
----------
p_npop : the number of sub-populations
p_isl_totmig : the total migration rate
p_theta : the value of parameter theta used in the simulations. It is constant.
p_rho : the value of parameter rho. It is constant.
p_en_time_pop_2 : the time at which the population 2 changes size.
p_en_size_pop_2 : the new size of population 2.

Summaries
---------
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


