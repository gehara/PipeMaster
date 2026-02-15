/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>

#ifdef R_PACKAGE_BUILD
#include <R.h>
#include <Rinternals.h>
#endif

         double
ran1()
{
        double drand48();
        return( drand48() );
}


#ifdef R_PACKAGE_BUILD

/* In R package mode: manage seed in memory, no file I/O */
static unsigned short msABC_seed[3] = {3579, 27011, 59243};

void seedit_r(unsigned short *seedv) {
    msABC_seed[0] = seedv[0];
    msABC_seed[1] = seedv[1];
    msABC_seed[2] = seedv[2];
    seed48(msABC_seed);
}

void get_seed_r(unsigned short *seedv) {
    unsigned short *pseed = seed48(msABC_seed);
    seedv[0] = pseed[0];
    seedv[1] = pseed[1];
    seedv[2] = pseed[2];
}

void seedit( char *flag )
{
    unsigned short *seed48(), *pseed;

    if( flag[0] == 's' ) {
        seed48( msABC_seed );
    }
    else {
        /* "end" mode: just save current state, no file write */
        pseed = seed48(msABC_seed);
        msABC_seed[0] = pseed[0];
        msABC_seed[1] = pseed[1];
        msABC_seed[2] = pseed[2];
    }
}

int
commandlineseed( char **seeds)
{
    unsigned short seedv[3], *seed48();

    seedv[0] = atoi( seeds[0] );
    seedv[1] = atoi( seeds[1] );
    seedv[2] = atoi( seeds[2] );

    msABC_seed[0] = seedv[0];
    msABC_seed[1] = seedv[1];
    msABC_seed[2] = seedv[2];

    seed48(seedv);
    return(3);
}

#else

/* Original standalone mode with file I/O */
	void seedit( char *flag )
{
	FILE *pfseed;
	unsigned short seedv[3], seedv2[3], *pseed ;
	int i;

	if( flag[0] == 's' ) {
	   pfseed = fopen("seedms","r");
	   if( pfseed == NULL ) {
           seedv[0] = 3579 ; seedv[1] = 27011; seedv[2] = 59243;
	   }
	   else {
	       seedv2[0] = 3579; seedv2[1] = 27011; seedv2[2] = 59243;
           for(i=0;i<3;i++){
		       if(  fscanf(pfseed," %hd",seedv+i) < 1 )
		            seedv[i] = seedv2[i] ;
		   }
	       fclose( pfseed);
	   }
	   seed48( seedv );
	}
	else {
	     pfseed = fopen("seedms","w");
         pseed = seed48(seedv);
         fprintf(pfseed,"%d %d %d\n",pseed[0], pseed[1],pseed[2]);
	}
}

	int
commandlineseed( char **seeds)
{
  FILE *pfseed;
  pfseed = fopen("seedms","w");
	unsigned short seedv[3], *seed48();

	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
	fprintf(pfseed, "%d %d %d\n", seedv[0], seedv[1], seedv[2] );

	seed48(seedv);
	return(3);
}

#endif
