/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>

#ifdef R_PACKAGE_BUILD
#include <R.h>
#include <Rinternals.h>
#endif

/* mingw-w64 (Rtools) does not ship drand48/seed48. Provide a POSIX-compatible
 * 48-bit LCG fallback on Windows. Linux/macOS continue to use libc's versions
 * via <stdlib.h>. */
#if defined(_WIN32) || defined(__MINGW32__) || defined(__MINGW64__)

#include <stdint.h>

static uint64_t pm_drand48_state = (0x1234ABCDULL << 16) | 0x330EULL;

static double drand48(void) {
    pm_drand48_state = (0x5DEECE66DULL * pm_drand48_state + 0xBULL)
                       & ((1ULL << 48) - 1);
    return (double)pm_drand48_state / (double)(1ULL << 48);
}

static unsigned short *seed48(unsigned short xseed[3]) {
    static unsigned short prev[3];
    prev[0] = (unsigned short)(pm_drand48_state & 0xFFFFULL);
    prev[1] = (unsigned short)((pm_drand48_state >> 16) & 0xFFFFULL);
    prev[2] = (unsigned short)((pm_drand48_state >> 32) & 0xFFFFULL);
    pm_drand48_state = ((uint64_t)xseed[2] << 32) |
                       ((uint64_t)xseed[1] << 16) |
                       (uint64_t)xseed[0];
    return prev;
}

#endif /* _WIN32 / MinGW fallback */

double ran1(void)
{
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
    unsigned short *pseed;

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
    unsigned short seedv[3];

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
	unsigned short seedv[3];

	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
	fprintf(pfseed, "%d %d %d\n", seedv[0], seedv[1], seedv[2] );

	seed48(seedv);
	return(3);
}

#endif
