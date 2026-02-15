# Assume the observed file seq.FA in fasta format

#first run the perl script in the msABC directory, fas2ms.pl
#The seq.FA contains 10 sequences. 9 of them are the ingroup, and the
#last, called "out" is the outgroup. The outgroup is used by fas2ms.pl
#to polarize the data. 
#The file seq.ms will be created which is the ms-like transformation
#of the seq.FA file. Some other files will be created, but they are
#not important for this application

#IMPORTANT: seq.ms contains 9 sequences, as many as the ingroup

../../fas2ms.pl fas=seq.FA outgroup=out > seq.ms 2>fas2ms.LOG

#Next run the msABC on the seq.ms
#

../../msABC 9 1 --obs seq.ms
