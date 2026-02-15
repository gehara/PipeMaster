-e 
#2
	 *** objective: Simulate 10 times one fragment from a sample of 12 lines,
	    theta is drawn from a prior U(10; 20). 

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
 

