#First create an "observed data file". Here, I create that using ms,
# but I pretend that the data is observed. THe data is stored in
# obs.out
# NOTICE: In order to run the following command, ms should be in the PATH

ms 10 5 -t 10 > obs.out

#Now use the msABC but instead of simulating data, 
#calculate summary statistics from the obs.out.
#Notice that the command line is the same like the normal msABC
#There is no need to give theta, rho, effective population size
#However, still you need to give the number of sequences (10)
# and the number of datasets (5)

../msABC 10 5 --obs obs.out

