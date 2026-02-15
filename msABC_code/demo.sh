for i in `seq -w 1 8`;do
    if [ ! -d example$i ] ; then mkdir example$i; fi;
done;

echo "example 1"
echo -e " Examples
# 1
\t *** objective: Simulate 10 times one fragment from a sample of 12 lines, all the parameters fixed. 

This is like the normal ms function. 
Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt.
Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.

Parameters:
p_theta	: the value parameter theta used in the simulations. In this example it's constant.
p_rho	: the value parameter rho used in the simulations. In this example it's constat.
p_pop_change_time : the value of time that the size of population changes.
p_pop_change_newpop : the ratio of the changed population size to the present population size.

Summaries:
s_segs	: the number of segregating sites.
s_theta_pi : the pi estimator of theta.
s_theta_w	: the Watterson's estimator of theta.
s_tajimasD	: Tajima's D.
s_ZnS	: ZnS. This is a summary of LD.
s_FayWuH : The Fay and Wu's H.

" > example1/readme.txt
./msABC 12 10 -t 10 -r 10 100 -eN 0.01 0.3 --verbose -seeds 1 2 3 > example1/out.txt; 
mv log.txt example1; # just mv the header because it will be overwritten 
mv tempoutput_ms.txt example1/
mv seedms example1/

echo "example 2"
echo -e "
#2
\t *** objective: Simulate 10 times one fragment from a sample of 12 lines,
\t    theta is drawn from a prior U(10; 20). 

Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt. 
Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.


Parameters:
p_theta	: the value of parameter theta used in the simulations. It is drawn from U(10; 20).
p_rho	: the value of parameter rho used in the simulations. In this example it's constat.
p_pop_change_time : the value of time that the size of population changes.
p_pop_change_newpop : the ratio of the changed population size to the present population size.

Summaries:
s_segs	: the number of segregating sites.
s_theta_pi : the pi estimator of theta.
s_theta_w	: the Watterson's estimator of theta.
s_tajimasD	: Tajima's D.
s_ZnS	: ZnS. This is a summary of LD.
s_FayWuH : The Fay and Wu's H.
 
" > example2/readme.txt
./msABC 12 10 -t -U 10 20 -r 10 100 -eN 0.01 0.3 --verbose -seeds 1 2 3 > example2/out.txt; 
mv log.txt example2; # just mv the header because it will be overwritten 
mv tempoutput_ms.txt example2/
mv seedms example2/

echo "example 3"
echo -e "
#3
\t *** objective: Simulate 10 times 4 (independent) loci from a sample of 12 lines.

The length of loci 1,2,3 is 500 and the length of the locus 4 is 100.

The file fragments.txt describes the configuration of the loci:
id      n       length
f1      12      500     
f2      12      500
f3      12      500
f4      12      100

Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt.
Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.




Parameters:
p_theta	: the value of parameter theta used in the simulations. It is constant.
p_rho	: the value of parameter rho used in the simulations. It is constant.
p_pop_change_time : the value of time that the size of population changes.
p_pop_change_newpop : the ratio of the changed population size to the present population size.

***********************************************
IMPORTANT: theta AND rho values in the command line refer to a locus of length 100 (as it is specified by the nsites).
This number is rescaled for each of the locus (since their length may vary).
***********************************************


Summaries (averages and variances are calculated among loci).
s_average_segs  : average number of segregating sites.
s_variance_segs : variance of the number of segregating sites.
s_average_pi    : average pi.
s_variance_pi       : variance of pi.
s_average_w     : average Watterson's theta.
s_variance_w    : variance of Watterson's theta.
s_average_tajd  : average Tajima's D.
s_variance_tajd : variance of Tajima's D.
s_average_ZnS   : average ZnS.
s_variance_ZnS  : variance of ZnS.
s_average_FayWuH   : average Fay and Wu's H.
s_variance_FayWuH  : variance of Fay and Wu's H.
 
" > example3/readme.txt
./msABC 12 10 -t 10 -r 10 100 -eN 0.01 0.3 --frag-begin --finp fragments.txt --frag-end --verbose -seeds 1 2 3 > example3/out.txt; 
mv log.txt example3; # just mv the header because it will be overwritten
mv tempoutput_ms.txt example3/
mv seedms example3/


echo "example 4"
echo -e "
#4
\t *** objective: Simulate 10 times, 4 (independent) loci from a sample of 12 lines.

The length of loci 1,2,3 is 500 and the length of the locus 4 is 100.

The file fragments.txt describes the configuration of the loci:
id      n       length
f1      12      500     
f2      12      500
f3      12      500
f4      12      100

Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt.

Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.

Parameters:
p_theta	: the value of parameter theta used in the simulations. It is constant.
p_rho	: the value of parameter rho used in the simulations. It is drawn from a prior U(10; 20).
p_pop_change_time : the value of time that the size of population changes. It is drawn from a prior U(0.01; 0.04)
p_pop_change_newpop : the ratio of the changed population size to the present population size. It is drawn from a prior U(0.3; 0.5).

***********************************************
IMPORTANT: theta AND rho values in the command line refer to a locus of length 100 (as it is specified by the nsites).
----------
This number is rescaled for each of the loci (since their length may vary).
***********************************************

Summaries (averages and variances are calculated among loci).
s_average_segs  : average number of segregating sites.
s_variance_segs : variance of the number of segregating sites.
s_average_pi    : average pi.
s_variance_pi       : variance of pi.
s_average_w     : average Watterson's theta.
s_variance_w    : variance of Watterson's theta.
s_average_tajd  : average Tajima's D.
s_variance_tajd : variance of Tajima's D.
s_average_ZnS   : average ZnS.
s_variance_ZnS  : variance of ZnS.
s_average_FayWuH   : average Fay and Wu's H.
s_variance_FayWuH  : variance of Fay and Wu's H.
 
" > example4/readme.txt
./msABC 12 10 -t 10 -r -U 10 20 100 -eN -U 0.01 0.04 -U 0.3 0.5 --frag-begin --finp fragments.txt --frag-end --verbose -seeds 1 2 3 > example4/out.txt; 
mv log.txt example4; # just mv the header because it will be overwritten 
mv tempoutput_ms.txt example4/
mv seedms example4/

echo "example 5"
# multiple populations
echo -e "
#5

\t *** objective: Simulate 10 times, 6 (independent) loci from a sample of 12 lines.

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

" > example5/readme.txt
./msABC 12 10 -t 10 -r -U 10 20 1000 -I 2 7 5 0.5 -eN 0.01 0.3 --frag-begin --finp fragments2.txt --frag-end --verbose -seeds 1 2 3 > example5/out.txt; 
mv log.txt example5/
mv tempoutput_ms.txt example5/
mv seedms example5/


echo "example 6"
echo -e "
#6 
\t *** objective: Simulate 10 times, 6 (independent) loci. Sample size
is variable among loci

The length of the loci is variable as well.


The file fragments.txt describes the configuration of the loci:


id	n	length	pop
f1	12	1000	0
f2	10	500	1
f2	9	500	2
f3	12	100	0
f4	12	2200	0
f5	12	1000	0	
f6	6	100	1
f6	17	100	2

NOTICE!!!! for the loci that the sample size is 12 we use pop=0. 
\t\t for the rest of the loci (where n != 12) we should specify the sample size
for each of the populations. The order of the information is important.

Polymorphism data are stored in tempoutput_ms.txt.
The parameter values with the summary statistics are stored in out.txt.

Each line represents a different replication. The first line is a header.

A logfile log.txt is created which contains the command line.




Parameters:
------------

p_theta : the value of parameter theta used in the simulations. It is constant.
p_rho   : the value of parameter rho. It is drawn from U(10; 20)
p_npop  : the number of sub-populations
p_isl_totmig  : the total migration rate   
p_em_time_pops_1_2    : time that the migration size changes. It is constant (=0.1)  
p_mig_1_2       : the new migration rate for sub-populations 1 and 2. It drawn from a prior U(0.1; 0.5)
p_pop_change_time       : the time that the total population size changes. It is drawn from a prior U(0.01; 0.02)
p_pop_change_newpop     : the ratio of the changed population to the population at the present



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

" > example6/readme.txt

./msABC 12 10 -t 10 -r -U 10 20 1000 -I 2 7 5 0.1 -em 0.1 1 2 -U 0.1 0.5 -eN -U 0.01 0.02 0.3 --frag-begin --finp fragments3.txt --frag-end  --verbose -seeds 1 2 3 > example6/out.txt
mv tempoutput_ms.txt example6/
mv log.txt example6/
mv seedms example6/


echo "example 7"
echo -e "
#6 
\t *** objective: Simulate 10 times, 2 (independent) loci. Sample size
is variable among loci. 

The samples originate from two populations. For the first locus 24 lines are from each population.
For the second locus 24 lines are derived from the first population and 23 lines from the second.

Data includes incomplete information, i.e. in the observed fasta alignment there are 'N's.

workflow
--------
1. use the script missing_ms_report.pl to obtain the information about the missing data. 
\t this script will transform in continuous form the missing positions.
\t this script should run for each of the loci:

\t missing_ms_report.pl in=locus1.fa > locus1.missing
\t missing_ms_report.pl in=locus2.fa > locus2.missing

Then, we create a list that contains the files that describe the missing information:
\t echo \"locus1.missing\" > missing.list
\t echo \"locus2.missing\" >> missing.list

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

" > example7/readme.txt

echo "# example07: run the script to obtain the missing information"
./missing_ms_report.pl in=locus1.fa > locus1.missing
./missing_ms_report.pl in=locus2.fa > locus2.missing

echo "# example07: create the list of files that contain the missing information"
echo -e "locus1.missing\nlocus2.missing" > missing.list

./msABC 48 10 -I 2 24 24 -U 0.5 0.6 -t 10 -r 10 1000 -en -U 0.01 0.02 2 -U 0.1 0.2 --missing missing.list --frag-begin --finp fragments4.txt --frag-end --verbose -seeds 1 2 3 >  example7/out.txt
mv tempoutput_ms.txt example7/
mv log.txt example7/
mv seedms example7/



echo "example 8"
echo -e "
#8
\t *** objective: Simulate 10 times, 2 (independent) loci. Sample size
is variable among loci. 
\t \t print out only the Tajima's D as a summary statistic (for each subpopulation)

Define the summary statistics
------------------------------
You may define in a file (here options.txt) that is denoted after the --options
flag, the summary statistics that you wish to print out. 
Notice that if private summary statistics for each subpopulation are printed,
then the summary statistic for the global sample is printed as well.


The samples originate from two populations. For the first locus 24 lines are from each population.
For the second locus 24 lines are derived from the first population and 23 lines from the second.

Data includes incomplete information, i.e. in the observed fasta alignment there are 'N's.

workflow
--------
1. use the script missing_ms_report.pl to obtain the information about the missing data. 
\t this script will transform in continuous form the missing positions.
\t this script should run for each of the loci:

\t missing_ms_report.pl in=locus1.fa > locus1.missing
\t missing_ms_report.pl in=locus2.fa > locus2.missing

Then, we create a list that contains the files that describe the missing information:
\t echo \"locus1.missing\" > missing.list
\t echo \"locus2.missing\" >> missing.list

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
s_average_tajd_1 : the average Tajima's D for sub-population 1
s_variance_tajd_1 : the variance of Tajima's D for sub-population 1
s_average_tajd_2 : the average of Tajima's D for sub-population 2
s_variance_tajd_2 : the variance of Tajima's D for sub-population 2
s_average_tajd : the average of Tajima's D for the total sample
s_variance_tajd : the variance of Tajima's D for the total sample


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

" > example8/readme.txt

echo "# example08: run the script to obtain the missing information"
./missing_ms_report.pl in=locus1.fa > locus1.missing
./missing_ms_report.pl in=locus2.fa > locus2.missing

echo "# example08: create the list of files that contain the missing information"
echo -e "locus1.missing\nlocus2.missing" > missing.list

./msABC 48 10 -I 2 24 24 -U 0.5 0.6 -t 10 -r 10 1000 -en -U 0.01 0.02 2 -U 0.1 0.2 --missing missing.list --frag-begin --finp fragments4.txt --frag-end --verbose -seeds 1 2 3 --options options.txt >  example8/out.txt
mv tempoutput_ms.txt example8/
mv log.txt example8/
mv seedms example8/
cp options.txt example8/
