-e 
#3
	 *** objective: Simulate 10 times 4 (independent) loci from a sample of 12 lines.

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
 

