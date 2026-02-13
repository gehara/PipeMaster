/***** msABC.c     ************************************************

A modification of Hudson's ms program to facilitate multi-locus ABC analysis.

For more information about ms see the web page of Dick Hudson:
http://home.uchicago.edu/~rhudson1/source/mksamples.html

***************************************************************************/





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "ms.h"
#include "data_sumstat.h"
#include <time.h>

#ifdef R_PACKAGE_BUILD
#include <setjmp.h>
#include "msABC_capture.h"
#include "sfs_sim.h"
jmp_buf msABC_jmpbuf;
static int msABC_jmpbuf_set = 0;
#define EXIT_MSABC(code) do { if (msABC_jmpbuf_set) longjmp(msABC_jmpbuf, (code) ? (code) : 1); else exit(code); } while(0)
#else
#define EXIT_MSABC(code) exit(code)
#endif


     
#define SITESINC 10 
#define NDIST 5
/* maximun number of parameters that the prior distribution may have*/
#define MAXDISTPARS 5
/* maximum number of summary statistics */
#define MAXSUMSTAT 20
#define MAXFRPARS 10
#define MAXPOPS 1000
/*this should be the same in data_sumstat.c*/
#define MISSING '2'
/* the minimum number of sequences that some statistics (e.g. Tajima's D) require */
#define MINSEQ 2
unsigned maxsites = SITESINC ;
/* flag variable referenced by streec.c */
int flag;
/* a variable that stores the greatest sample size among loci */
int largestnsam;
int derived;
/*the name of the fragment file*/
char fragfile[100];
/*name of the options file */
char optionsfile[100];
/* name of the list of the missing data */
char missingfile[100];
/* name of the file with the observed data */
char obfile[100];
char poldatafile[100];
/*proportion of meaningful fragments for calculation of pairwise summary statistics
  if a summary statistics is undefined in more than this proportion then it will be reported a nan
  as average and variance
*/
double propfragments;


int* missingdata;
/* a flag that gets the value 1 if we use fragments */
int fragmode;
/* stores the number of loci we want to simulate */
int nofragments;
/* index of the fragment i.e. 1,2 .. it is used in the missing data */
int fragindex, missing_nof_fragments;
/* number of parameters for each locus */
int fcolumns;
/* the weight for each locus has been set to 1. Alternatively 
   it may be used a weight equal to the length of the fragment,
   because larger fragments have higher rec rate, i.e. less variance in sum. stats
*/
int weightmode;
/* denotes what measure of Fst we need to calculate (hsm, hbk, slatkin) */
char fst_type[10];
int tip;
int misdata;

/* loci parameters */
char **fname;
int *fsize;
double *frecbp;
double *fmubp;
int *flength;
int *fpop;
/* mu/rec overrides for batch mode (set by sfs_sim.c / msABC_wrapper.c) */
double *frag_mu_override = NULL;
int frag_mu_override_len = 0;
double *frag_rec_override = NULL;
int frag_rec_override_len = 0;
/****************/

double *timeevents; // store the events (-e)
/* parallel restrictions in the order of the events */
int** conditional_timeevents; 
int event; // an index for the total number of the events
int nof_resampling;


/* for private sstats */
int* samples_b; /* stores the beginning of the sample for each subpop */
int* samples_e; /* stores the end of the sample for each subpop */


int frpars[MAXFRPARS];
double effectiveN;

typedef struct {
  // this just denotes ch sum. stats should be used.
  int segs;
  int thetaw;
  int  thetapi;
  int tajimasD;
  int ZnS;
  int Fst;
  int shared;
  int private;
  int fixed_dif;
  int fst_pops;
  int* Fst_pops ;
  int fwh;
  int dvstat;
  int thomson_est;
  int thomson_var;
  /* private for each sub-pop summary statistics */
  int prisegs;
  int prithetaw;
  int prithetapi;
  int pritajimasD;
  int priZnS;
  int prifwh;
  int pridvstat;
  int prithomson_est;
  int prithomson_var;
  
  int pristats; // just a flag for all the private statistics
  

  int howmany;
} sumstats;

int filter;
int resample;
double bm; // back mutation probability

int sstats_array[MAXSUMSTAT];
double sstats_weights[MAXSUMSTAT];
double sstats_denominator_weights[MAXSUMSTAT];

int* allelic_map = NULL;
int* missing_map = NULL;
int** npop_allelic_map = NULL;
int** npop_missing_map = NULL;
// node structure of the coalescent tree
struct node{
	int abv;
	int ndes;
	float time;
	};

// segment structure of the coalescent tree
struct segl {
	int beg;
	struct node *ptree;
	int next;
	};

double *posit ;
double segfac ;
int count, ntbs, nseeds ;
struct params pars ;


/*added by Jeff. Ross-Ibara Aug 2007 */                                                                         
char ** rematrix(int nsam, int largestnsam, int howlong, char ** m);
//char ** reinitialize(int nsam, int largestnsam, int howlong, char ** m);
int ** rematrix_int(int nsam, int largestnsam, int howlong, int ** m);

/** the distribution flags, for the moment the Uniform is implemented
    It's trivial to add some more
*/
char distros[NDIST] = {'U', 'N', 'G', 'L', 'R'};	

/* stores the distributions used here.
   e.g. {U}niform, {N}ormal
*/
char *param_distr;

/* stores the parameters of the distribution. e.g. for the
   uniform it stores the min and max valu
*/
double** param_distr_values;
double* global_values;

// a flag: if 1 then the duration mode is used i.e. the time of each event is defined
// relatively to the time of the previous event.
int durationmode;
int printall;

// to imitate resulsts from NGS
int outgroup;
int phasemode;

/**
   options for the Fst
*/

/**
   The weights of the population. 
   This is usually proportional to the size of the populations
*/
double* weights; 

/**********************************************************************/

int initialize_global(){
  int i,j;
  derived = 0;

  propfragments = 0;
  nofragments = 0;
  misdata = 0;
  durationmode = 0;
  printall = 0;
  for(i=0; i<MAXFRPARS; i++)
    frpars[i] = 0;

  weightmode = 0;
  largestnsam = 0;
  strcpy(optionsfile,"");
  strcpy(missingfile, "");
  missingdata = NULL;
  missing_nof_fragments=0;

  bm = 0.;
  strcpy(fst_type, "hbk");

  /*
     the weights for the summary statistics for each fragment.
     It's set up to 1. Alternatively one can use the length of the fragment as weight
     for the variance and mean of the variance statistics.
     The length may be useful as a weight because large fragments contain less variance for the summary statistics
  */
  for(i=0; i<MAXSUMSTAT; i++)
    sstats_weights[i] = 1.;
  for(i=0; i<MAXSUMSTAT; i++)
    sstats_denominator_weights[i] = 0.;


  samples_b = NULL;
  samples_e = NULL;

#ifdef R_PACKAGE_BUILD
  /* Reset additional globals for safe repeated calls */
  count = 0;
  ntbs = 0;
  nseeds = 0;
  effectiveN = -1.;
  maxsites = SITESINC;
  segfac = 0.0;
  event = 0;
  nof_resampling = 0;
  filter = 0;
  resample = 0;
  fragmode = 0;
  fragindex = 0;
  fcolumns = 0;
  tip = 0;
  outgroup = 0;
  phasemode = 0;
  strcpy(fragfile, "");
  strcpy(obfile, "");
  strcpy(poldatafile, "");
  posit = NULL;
  param_distr = NULL;
  param_distr_values = NULL;
  global_values = NULL;
  weights = NULL;
  allelic_map = NULL;
  missing_map = NULL;
  npop_allelic_map = NULL;
  npop_missing_map = NULL;
  timeevents = NULL;
  conditional_timeevents = NULL;
  fname = NULL;
  fsize = NULL;
  frecbp = NULL;
  fmubp = NULL;
  flength = NULL;
  fpop = NULL;
  memset(&pars, 0, sizeof(struct params));
#endif

  return 1;

}



int modifyList(int outgroup, int phasemode, char** list, int segsites, int seqs)
{
  if(outgroup == 1)
    return 0;
  if(phasemode != 2) 
    return 0;

  int i=0, j=0;
  
  
  if(outgroup == 0 && phasemode == 2)
    {
      assert(seqs % 2 == 0);
      int * nones = calloc(segsites, sizeof(int));

      for(i=0; i<segsites; ++i)
	for(j=0; j<seqs; ++j)
	  if(list[j][i] == '1')
	    nones[i]+=1;

      for(i=0; i<segsites; ++i)
	{
	  

	  if(nones[i] < seqs/2)
	    {
	      for(j=0; j<seqs; ++j)
		{
		  if(list[j][i] == '1')
		    list[j][i] = '0';
		  else if(list[j][i] == '0')
		    list[j][i] = '1';
		}
	    }
	
	  /* shuffle them within pairs */
	  
	  for(j=0; j<seqs; j+=2)
	    {
	      
	      double ran = (double)rand();

	      double r = ran/(double)((RAND_MAX)+1.);
	     
	      fprintf(stderr, "R is %e\n", r);
 
	      if(r < 0.5)
		{
		  char temp = list[j][i];
		  list[j][i] = list[j+1][i];
		  list[j+1][i] = temp;
		}
	    }
	  
	}
      
    }
  return 1;
}

 
  
int msABC_main(int argc, char** argv){
  
  effectiveN = -1.;
  //set the timer
  //time_t t1,t2, t_a, t_b, t3, t4;
  //double timedif=0, timedif1=0.;
  //time(&t1);
  //FILE* timefile = fopen("time_measurements.txt", "w"); 
  
  initialize_global();
  
  // indices
  int i,j, npopi, npopj, ii;
  /* consider the summary statistic that are marked with 1 */
  // the last number denotes the number of summary statistics
  sumstats sstats = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, NULL, 1, 1, 1, 1,// global and pairwise
		     1, 1, 1, 1, 1, 1, 1, 1, 1,// private
		     1, 14}; // general flags
  

  /** here denote which summary statistics you want to use
     I can modif that to set up by the user
  */
  for(i=0; i<sstats.howmany; i++)
    sstats_array[i] = 1;

  
  param_distr = (char*)malloc( argc*sizeof(char));
  global_values = (double*)malloc( argc*sizeof(double));
  param_distr_values = (double**)malloc( argc*sizeof(double*));


  /*
    -------------------- LOCAL VARIABLES -----------------------------
  */

  /* Missing data */
  int** misseq;
  double **misstart, **misend;

  int zpop, zpop1, zpop2, z;
  FILE* fragin;
  FILE* optionsin=NULL; // the stream with the options
  FILE* missingin=NULL;
  FILE* obin=NULL; /* the file with the observed data */

  // the second letter in the -e. options, i.e. -eN, the N is stored in secarg
  char secarg;
  // initial values of parameters. The ones that are given in the command line
  double init_theta=0., init_rho=0.; //, output_theta=0., output_rho=0.; $ these are not needed anymore
  int init_nsites = 0, init_nsam=0;
  //  int total_nsites = 0; // $ this is not needed any more
  // the number of total replications, number of segregating sites
  int  k, howmany, segsites ;
  // the polymorphic table in list
  char **cmatrix(int, int);
  char **list, **tbsparamstrs ;
  // output files. the logfile stores the information about the command line and the random seeds.
  FILE *pf, *logfile, *tempdatafile=NULL ;
  double probss, tmrca, ttot ;
  void seedit( const char * ) ;
  void getpars( int argc, char *argv[], int *howmany )  ;
  int check_timeevents( int argc, double* timeevents);
  void getfragpars(FILE* f, int r, int c);
  void getfragdimensions(FILE* f, int* r, int* c);
  void get_timeevents(int argc, char* argv[]);
  
  void print_eventlist();
  int check_conditional_timeevents(int nof_resampling, int** conditional_timeevents, double* timeevents);
  double varstat(double sumsq, double ave, double den);
  double skewstat(double sum3, double sum2, double ave, double den);
  double kurtstat(double sum4, double sum3, double sum2, double ave, double den);
  int getmissing(FILE* missingin, int*** misseq, double*** misstart, double*** misend, int** missingdata);
  int getobservations(FILE* obin, char*** data);

  int put_missing(char** list, double* posit, int nsam, int segsites, int** misseq, double** misstart, double** misend, int missingdata, int fragindex);
  int getoptions(FILE* optionsin, 
		 int* segs,
		 int* thetaw,
		 int* thetapi,
		 int* tajimasD,
		 int* ZnS,
		 int* Fst,
		 int* shared,
		 int* private,
		 int* fixed_dif,
		 int* fst_pops,
		 int* Fst_pops,
		 int* fwh,
		 int* dvstat,
		 int* thomson_est,
		 int* thomson_var,
		 // private
		 int* prisegs,
		 int* prithetaw,
		 
		 
		 int*prithetapi,
		 int* pritajimasD,
		 int* priZnS,
		 int* prifwh,
		 int* pridvstat,
		 int* prithomson_est,
		 int* prithomson_var,
		 
		 int* pristats);
  
  double Fst(char fst[10], double p1, double p2, double p3, double p4);
  int gensam( char **list, double *probss, double *ptmrca, double *pttot, int largestnsam, int old_maxnsam) ;
  int modifyList(int outgroup, int phasemode, char** list, int segsites, int seqs);
  
  int printdata = 0;
  int sumstat = 1;
  int distribution_arguments( char distr); 
  double total_stats_weights[MAXSUMSTAT];
  int largesample = 100, largestinsegsites=100; // for initial memory allocation



  /*
    memory for the streec
  */
  int old_maxnsam = 0;

  /* the sample size in the previous round;
     It is  needed for the calculation of denominators in the summary statistics
  */
  int old_nsam = 0;
  
  /*
    For optimizing a bit the ZnS calculations
  */
  int* maxder;


  /**
     Fst calculations
  */

  

  double piT=0., 
    piS=0.,
    piB=0.,
    piD=0.,
    ppiT=0.,
    ppiS=0.,
    ppiB=0.,
    ppiD=0.;
  
  
  double fst=0.;
  double pairwise_fst = 0.;


  /*
    multiple populations 
  */
  int f_total_sample_size = 0; 
  int calculations_done = 0;
  int pairwise_calculations_done = 0;
  int ipop, jpop, pop1, pop2, pop, c;
  int fst_pops;
  int* Fst_pops;



  /* for the denominators of the statistics */
  double hn, sqhn, bn;

  /* for the denominators of the private statistics */
  double*prihn, *prisqhn, *pribn;

  
  
  
  
  
  
     
  
  /*
    ---------------------------- LOCAL VARIABLES -------------------------------
  */
  
  timeevents = malloc(argc*sizeof(double));
  conditional_timeevents = malloc( (unsigned)argc * sizeof(int*)); // 20091015
  for(i=0; i<argc; i++)
    conditional_timeevents[i] = malloc( (unsigned)argc * sizeof(int));
  
  /* put default values for each parameter 
     definitely there are less parameters than argc, 
     so allocating space for argc should be more than enough 
  */
  for(i=0; i<argc; i++)
    param_distr[i] = '-';
  
  /* initialize the vector of parameter values 
     I consider that a distribution cannot have more than 5 parameters 
  */
  for(i=0; i<argc; i++){
    param_distr_values[i] = (double*)malloc( MAXDISTPARS*sizeof(double));
    int ii = 0;
    for(ii=0; ii<5; ii++)
      param_distr_values[i][ii] = 0.0;
  }
  

  ntbs = 0 ;   /* these next few lines are for reading in parameters from a file (for each sample) */
  tbsparamstrs = (char **)malloc( argc*sizeof(char *) ) ;
  
#ifdef R_PACKAGE_BUILD
  /* Suppress log.txt file writes in R package mode */
  logfile = NULL;
#else
  /* open the logfile and write there the command line */
  if( (logfile = fopen("log.txt", "w") ) == NULL){
    fprintf(stderr, "\n\nERROR: cannot open the log.txt for output.\n");
    fprintf(stderr, "Program will exit now.\n");
    EXIT_MSABC(1);
  }
  for( i=0; i<argc; i++) fprintf(logfile, "%s ",argv[i]);
  fclose(logfile);
#endif
  
  if(tempdatafile != NULL)
    for( i=0; i<argc; i++) fprintf(tempdatafile,"%s ",argv[i]);
  
  for( i =0; i<argc; i++) tbsparamstrs[i] = (char *)malloc(30*sizeof(char) ) ;
  for( i = 1; i<argc ; i++)
    if( strcmp( argv[i],"tbs") == 0 )  argv[i] = tbsparamstrs[ ntbs++] ;
  
  
  
  if( ntbs > 0 )  for( k=0; k<ntbs; k++)  scanf(" %s", tbsparamstrs[k] );

  
 
  /**
     get the parameters 
     in case of resampling, check the order of the events
  */
  getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */
  
  /* initial sample size that is given from the command line */
  init_nsam = pars.cp.nsam;

  while( (resample == 1 && count && !check_timeevents(event, timeevents) ) ||
	 (
	  resample == 2 && count && 
	  (check_conditional_timeevents(nof_resampling, conditional_timeevents, timeevents) != -1)
	  )
	 ){
    get_timeevents(argc, argv);
  }
  

  if( !pars.commandlineseedflag ) seedit( "s");
#ifdef R_PACKAGE_BUILD
  pf = msABC_get_output_stream();
#else
  pf = stdout ;
#endif
  
  
  /* I guess that here the space for the polymorphic table is allocated */
  if( pars.mp.segsitesin ==  0 ) {
    list = cmatrix(pars.cp.nsam, maxsites+1);
    posit = (double *)malloc( (unsigned)( maxsites*sizeof( double)) ) ;
  }
  else {
    list = cmatrix(pars.cp.nsam, pars.mp.segsitesin+1 ) ;
    posit = (double *)malloc( (unsigned)( pars.mp.segsitesin*sizeof( double)) ) ;
    if( pars.mp.theta > 0.0 ){
      segfac = 1.0 ;
      for(  i= pars.mp.segsitesin; i > 1; i--) segfac *= i ;
    }
  }
  largestnsam = pars.cp.nsam;
  
  
  int* pritruesegsites;
  if(sumstat && sstats.pristats == 1 && pars.cp.npop > 1){
    pritruesegsites = (int*)malloc(pars.cp.npop * sizeof(int));
  }
  
  fst_pops = sstats.fst_pops;
  Fst_pops = (int*)malloc(pars.cp.npop*sizeof(int));

  /*denominators of private statistics */
  if(sumstat && pars.cp.npop > 1 && sstats.pristats == 1){
    prihn = (double*)malloc(pars.cp.npop * sizeof(double));
    prisqhn = (double*)malloc(pars.cp.npop * sizeof(double));
    pribn = (double*)malloc(pars.cp.npop * sizeof(double));
  }

  /* private summary statistics */
  samples_b = (int*)malloc(pars.cp.npop*sizeof(int));
  samples_e = (int*)malloc(pars.cp.npop*sizeof(int));


  double** sstats_popden_weights;
  if( (sstats.pristats==1) && (pars.cp.npop>1) ){
    /* for the denominator of private averages */
    sstats_popden_weights = malloc(MAXSUMSTAT*sizeof(double*));
    for(i=0; i<MAXSUMSTAT; i++)
      sstats_popden_weights[i] = malloc(pars.cp.npop * sizeof(double));
  }
  /*denominators for private averages */

 
  double ***sstats_popdenpair_weights; /* pp20100422 */
      /* for the denominator of pairwise statistics */

  /* inlude  pp20100422 */
  if(pars.cp.npop > 1){
    sstats_popdenpair_weights = malloc(MAXSUMSTAT*sizeof(double**));
    for(i=0; i<MAXSUMSTAT; i++){
      sstats_popdenpair_weights[i] = malloc(pars.cp.npop * sizeof(double*));
      for(j=0; j<pars.cp.npop; j++)
	sstats_popdenpair_weights[i][j] = malloc(pars.cp.npop * sizeof(double));
    }
  }
  

  if(strcmp(missingfile, "") && (missingin=fopen(missingfile, "r")) == NULL){
    fprintf(stderr, "\nERROR opening the missing info file\n");
    EXIT_MSABC(1);
  }
  if(missingin != NULL){
    // initial memory allocation; later realloc
    
    getmissing( missingin, &misseq, &misstart, &misend, &missingdata);
  }
  
  /* open the observed data */
  if(strcmp(obfile, "") && (obin = fopen(obfile, "r")) == NULL){
    fprintf(stderr, "\nERROR opening the observed files\n");
    EXIT_MSABC(1);
  }


  if( strcmp(poldatafile, "") && (tempdatafile = fopen(poldatafile, "w") ) == NULL){
    fprintf(stderr, "\nERROR opening the output file\n");
    EXIT_MSABC(1);
  
  }
  


		
	

  if(strcmp(optionsfile, "") && (optionsin=fopen( optionsfile, "r")) == NULL){
    fprintf(stderr, "\nERROR opening the options file\n");
    EXIT_MSABC(1);
  }
  if(optionsin != NULL){
    getoptions( optionsin, 
		&sstats.segs, 
		&sstats.thetaw,
		&sstats.thetapi,
		&sstats.tajimasD,
		&sstats.ZnS,
		&sstats.Fst,
		&sstats.shared,
		&sstats.private,
		&sstats.fixed_dif,
		&sstats.fst_pops,
		Fst_pops,
		&sstats.fwh,
		&sstats.dvstat,
		&sstats.thomson_est,
		&sstats.thomson_var,
		
		&sstats.prisegs,
		&sstats.prithetaw,
		&sstats.prithetapi,
		&sstats.pritajimasD,
		&sstats.priZnS,
		&sstats.prifwh,
		&sstats.pridvstat,
		&sstats.prithomson_est,
		&sstats.prithomson_var,
		
		&sstats.pristats
		);
    sstats_array[0] = sstats.segs;
    sstats_array[1] = sstats.thetaw;
    sstats_array[2] = sstats.thetapi;
    sstats_array[3] = sstats.tajimasD;
    sstats_array[4] = sstats.ZnS;
    sstats_array[5] = sstats.Fst;
    sstats_array[6] = sstats.shared;
    sstats_array[7] = sstats.private;
    sstats_array[8] = sstats.fixed_dif;
    sstats_array[9] = (sstats.fst_pops > 0) ? 1:0;
    sstats_array[10] = sstats.fwh;
    sstats_array[11] = sstats.dvstat;
    sstats_array[12] = sstats.thomson_est;
    sstats_array[13] = sstats.thomson_var;
    
    fst_pops = sstats.fst_pops;
    
  }

  /**
     Fst calculations, private, shared and fixed_dif matrix
     npop x npop
  */
  

  
  /* initialize Fst_pops */
  
  
  if( sstats.fst_pops == 1){
    fst_pops = pars.cp.npop;
    for(i=0; i<fst_pops; i++)
      Fst_pops[i] = i;
  }
  else if(sstats.fst_pops > 0){
    for(i=0; i<fst_pops; i++)
      Fst_pops[i] = i;
  }
  else{
    fst_pops = 0;
  }

  
  //printf("fstpops: %d\n", fst_pops);

  double** shared = (double**)malloc(pars.cp.npop * sizeof(double*));
  double** private = (double**)malloc(pars.cp.npop * sizeof(double*));
  double** fixed_dif = (double**)malloc(pars.cp.npop * sizeof(double*));
  

  for(i=0; i<pars.cp.npop; ++i){
    shared[i] = (double*)malloc(pars.cp.npop * sizeof(double));
    private[i] = (double*)malloc(pars.cp.npop * sizeof(double));
    fixed_dif[i] = (double*)malloc(pars.cp.npop * sizeof(double));
    for( j=0; j<pars.cp.npop; j++){
      shared[i][j] = 0.;
      private[i][j] = 0.;
      fixed_dif[i][j] = 0.;
    }
  }

  //save the initial configuration
  for(i=0; i<pars.cp.npop; ++i){
    //printf("initconfig: %d, config: %d\n", pars.cp.initconfig[i], pars.cp.config[i]);
    (pars.cp.initconfig)[i] = (pars.cp.config)[i];
  }
  
  /*
    the weights for the Fst calculation equal to the initial population size proportions
  */
  weights = (double*)malloc(pars.cp.npop * sizeof(double));

  /* p20100422 find the number of populations from which we have a sample */
  int sampledpops = pars.cp.npop;
  for(i=0; i<pars.cp.npop; i++){
    if( (pars.cp.config)[i] < 1)
      sampledpops--;
  }

  for(i=0; i<pars.cp.npop; i++){
    /* p20100422 if for a population we have no sample, make its weight 0 */
    if( (pars.cp.config)[i] > 0)
      weights[i] = 1./(double)sampledpops; // alternatively it could be pars.cp.size[i];
    else weights[i] = 0.;
  }
  

  double** ave_shared=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** ave_private=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** ave_fixed_dif=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** ave_pairwise_Fst = (double**)malloc(pars.cp.npop * sizeof(double*));
  double** var_shared=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** var_private=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** var_fixed_dif=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** var_pairwise_Fst = (double**)malloc(pars.cp.npop * sizeof(double*));
  double** skew_shared=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** skew_private=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** skew_fixed_dif=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** skew_pairwise_Fst = (double**)malloc(pars.cp.npop * sizeof(double*));
  double** kurt_shared=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** kurt_private=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** kurt_fixed_dif=(double**)malloc(pars.cp.npop * sizeof(double*));
  double** kurt_pairwise_Fst = (double**)malloc(pars.cp.npop * sizeof(double*));
  for(npopi=0; npopi<pars.cp.npop; ++npopi){
      ave_shared[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      ave_private[npopi]= (double*)malloc(pars.cp.npop * sizeof(double));
      ave_fixed_dif[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      ave_pairwise_Fst[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      var_shared[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      var_private[npopi]= (double*)malloc(pars.cp.npop * sizeof(double));
      var_fixed_dif[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      var_pairwise_Fst[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      skew_shared[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      skew_private[npopi]= (double*)malloc(pars.cp.npop * sizeof(double));
      skew_fixed_dif[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      skew_pairwise_Fst[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      kurt_shared[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      kurt_private[npopi]= (double*)malloc(pars.cp.npop * sizeof(double));
      kurt_fixed_dif[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
      kurt_pairwise_Fst[npopi] = (double*)malloc(pars.cp.npop * sizeof(double));
  }

  /* allocate memory for the summary statistics */
  
  int** ders=NULL;
  int*** npop_ders=NULL;
  if(sumstat){
    allelic_map = (int*)malloc(largestinsegsites*sizeof(int));
    missing_map = (int*)malloc(largestinsegsites*sizeof(int));
    ders = (int**)malloc(largesample*sizeof(int*)); //largesample is just an arbitrary number
    for(j=0; j<largesample; j++)
      ders[j] = (int*)malloc(largestinsegsites*sizeof(int));
  }

  /* if multiple populations then calculate the allelic map for all of them */
  if(sumstat && pars.cp.npop > 1 && sstats.pristats == 1){
    npop_allelic_map = malloc(pars.cp.npop*sizeof(int*));
    npop_missing_map = malloc(pars.cp.npop*sizeof(int*));
    //printf("largestinsegsites: %d -1\n", largestinsegsites);
    //EXIT_MSABC(5);
    for(j=0; j<pars.cp.npop; j++){
      npop_allelic_map[j] = malloc(largestinsegsites*sizeof(int));
      npop_missing_map[j] = malloc(largestinsegsites*sizeof(int));
    }

    npop_ders = (int***)malloc(pars.cp.npop*sizeof(int**));
    for(ipop=0; ipop<pars.cp.npop; ipop++){
      npop_ders[ipop]=(int**)malloc(largesample*sizeof(int*));
      for(j=0; j<largesample; j++)
	npop_ders[ipop][j] = (int*)malloc(largestinsegsites*sizeof(int));
    }
  }



  double *prisegs=NULL,
    *prithetapi=NULL, 
    *prithetaw = NULL, 
    *pritajd = NULL, 
    *prithetah = NULL, 
    *prifwh = NULL, 
    *prizns = NULL,
    *pridvk = NULL,
    *pridvh = NULL,
    *prithomson_est = NULL,
    *prithomson_var = NULL;
  
  
  
  //* *for private statistics *\/ */
/*   if(sumstat && sstats.pristats == 1 && pars.cp.npop > 1){ */
/*     prisegs = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     prithetapi = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     prithetaw = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     pritajd = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     prithetah = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     prifwh = (double*)malloc(pars.cp.npop * sizeof(double)); */
/*     prizns = (double*)malloc(pars.cp.npop * sizeof(double)); */
    
    
/*   } */
  
  
 


  double* priave_segs;
  double* priave_thetapi;
  double* priave_thetaw;
  double* priave_fwh;
  double* priave_tajimasD;
  double* priave_ZnS;
  double* priave_dvk;
  double* priave_dvh;
  double* priave_thomson_est;
  double* priave_thomson_var;

  double* privar_segs;
  double* privar_thetapi;
  double* privar_thetaw;
  double* privar_fwh;
  double* privar_tajimasD;
  double* privar_ZnS;
  double* privar_dvk;
  double* privar_dvh;
  double* privar_thomson_est;
  double* privar_thomson_var;

  double* priskew_segs;
  double* priskew_thetapi;
  double* priskew_thetaw;
  double* priskew_fwh;
  double* priskew_tajimasD;
  double* priskew_ZnS;
  double* priskew_dvk;
  double* priskew_dvh;
  double* priskew_thomson_est;
  double* priskew_thomson_var;

  double* prikurt_segs;
  double* prikurt_thetapi;
  double* prikurt_thetaw;
  double* prikurt_fwh;
  double* prikurt_tajimasD;
  double* prikurt_ZnS;
  double* prikurt_dvk;
  double* prikurt_dvh;
  double* prikurt_thomson_est;
  double* prikurt_thomson_var;

  /*for private statistics */
  if(sumstat && sstats.pristats == 1 && pars.cp.npop > 1){
    prisegs = (double*)malloc(pars.cp.npop * sizeof(double));
    prithetapi = (double*)malloc(pars.cp.npop * sizeof(double));
    prithetaw = (double*)malloc(pars.cp.npop * sizeof(double));
    pritajd = (double*)malloc(pars.cp.npop * sizeof(double));
    prithetah = (double*)malloc(pars.cp.npop * sizeof(double));
    prifwh = (double*)malloc(pars.cp.npop * sizeof(double));
    prizns = (double*)malloc(pars.cp.npop * sizeof(double));
    pridvk = (double*)malloc(pars.cp.npop * sizeof(double));
    pridvh = (double*)malloc(pars.cp.npop * sizeof(double));
    prithomson_est = (double*)malloc(pars.cp.npop * sizeof(double));
    prithomson_var = (double*)malloc(pars.cp.npop * sizeof(double));
    
    priave_segs = (double*)malloc(pars.cp.npop * sizeof(double));
    priave_thetapi= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_thetaw= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_fwh= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_tajimasD= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_ZnS= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_dvk= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_dvh= (double*)malloc(pars.cp.npop * sizeof(double));
    priave_thomson_est = (double*)malloc(pars.cp.npop * sizeof(double));
    priave_thomson_var = (double*)malloc(pars.cp.npop * sizeof(double));
    
    privar_segs= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_thetapi= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_thetaw= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_tajimasD= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_fwh= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_ZnS= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_dvk= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_dvh= (double*)malloc(pars.cp.npop * sizeof(double));
    privar_thomson_est = (double*)malloc(pars.cp.npop * sizeof(double));
    privar_thomson_var = (double*)malloc(pars.cp.npop * sizeof(double));

    priskew_segs= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_thetapi= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_thetaw= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_tajimasD= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_fwh= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_ZnS= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_dvk= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_dvh= (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_thomson_est = (double*)malloc(pars.cp.npop * sizeof(double));
    priskew_thomson_var = (double*)malloc(pars.cp.npop * sizeof(double));

    prikurt_segs= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_thetapi= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_thetaw= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_tajimasD= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_fwh= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_ZnS= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_dvk= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_dvh= (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_thomson_est = (double*)malloc(pars.cp.npop * sizeof(double));
    prikurt_thomson_var = (double*)malloc(pars.cp.npop * sizeof(double));

  }


  
  
  
  /* repeat for every replication ( the second option in ms )
     some people translate it as independent loci
     this is the main loop
  */
  
  count=0;
  while( howmany-count++ ){
    /* fprintf(stderr, "count: %d, argv[8]: %s\n", count, argv[8]); */
    
    
    calculations_done = 0;
    pairwise_calculations_done = 0;
    //printf("count: %d\n", count); // $p
    // reset the total sample size
    f_total_sample_size = 0;

    // timer
    /* int interval=100; */
/*     if( (count % interval) == 0){ */
/*       time(&t2); */
/*       timedif = difftime(t2, t1); */
/*       fprintf(stderr, "%d\t%e\n", count, timedif); */
/*     } */
    //    printf("count: %d\n", count); // $p
    /* if tbs */
    if( (ntbs > 0) && (count >1 ) ){
      for( k=0; k<ntbs; k++){ 
	if( scanf(" %s", tbsparamstrs[k]) == EOF ){
#ifndef R_PACKAGE_BUILD
	  if( !pars.commandlineseedflag ) seedit( "end" );
#endif
	  EXIT_MSABC(0);
	}
      }
      
      
      /**
	 get the parameters 
	 in case of resampling, check the order of the events
      */
      getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */
      
      while( (resample == 1 && count && !check_timeevents(event, timeevents) ) ||
	     (
	      resample == 2 && count && 
	      (check_conditional_timeevents(nof_resampling, conditional_timeevents, timeevents) != -1)
	      )
	     ){
	get_timeevents(argc, argv);
	
      }
      
      
    }
      /* /\** */
/* 	 get the parameters  */
/* 	 in case of resampling, check the order of the events */
/*       *\/ */
/*       do{ */
/* 	getpars( argc, argv, &howmany) ;   /\* results are stored in global variable, pars *\/ */
/*       } */
/*       while(resample == 1 && count && !check_timeevents(event, timeevents)); */
      
/*     } */
    
    /* it should re-read the parameters in every replication even in the first one
       for the first one: it should produce random numbers for the variables according to the random number seeds
       for the next: it should reproduce a random number for the parameter
    */
    if( ntbs <= 0){
      getpars( argc, argv, &howmany) ;   /* results are stored in global variable, pars */
      
      
      /* for(i = 0; i<event; i++){ */
/* 	fprintf(stderr, "%e\t", timeevents[i]); */
/*       } */
/*       fprintf(stderr, "\n"); */
      
      //int wwww=0;
      while( (resample == 1 && count && !check_timeevents(event, timeevents) ) ||
	     (
	      resample == 2 && count && 
	      (check_conditional_timeevents(nof_resampling, conditional_timeevents, timeevents) != -1)
	      )
	     ){
	//fprintf(stderr, "%d\n", wwww);
	get_timeevents(argc, argv);
	//wwww++;
	//printf("%d", 3);
      }
    }
    
    /* ******************************
       
       FRAGMENT MODE ... you want to simulate multiple loci
       --frag-begin flag has found and a file with the fragments will follow 
    */
    if(fragmode == 1){
      
      
      
      double ave_thetaw = 0.;
      double ave_thetapi = 0.;
      double ave_tajimasD = 0.;
      double ave_fwh = 0.;
      double ave_segs = 0.;
      double ave_ZnS = 0.;
      double ave_fst = 0.;
      double ave_dvk = 0.;
      double ave_dvh = 0.;
      double ave_thomson_est = 0.;
      double ave_thomson_var = 0.;
      
      double var_thetaw = 0.,
	var_thetapi = 0.,
	var_tajimasD = 0.,
	var_segs = 0.,
	var_ZnS = 0.,
	var_fst = 0.,
	var_fwh = 0.,
	var_dvk = 0.,
	var_dvh = 0.,
	var_thomson_est = 0.,
	var_thomson_var = 0.;

      double skew_thetaw=0., skew_thetapi=0., skew_tajimasD=0., skew_segs=0.,
	skew_ZnS=0., skew_fst=0., skew_fwh=0., skew_dvk=0., skew_dvh=0.,
	skew_thomson_est=0., skew_thomson_var=0.;
      double kurt_thetaw=0., kurt_thetapi=0., kurt_tajimasD=0., kurt_segs=0.,
	kurt_ZnS=0., kurt_fst=0., kurt_fwh=0., kurt_dvk=0., kurt_dvh=0.,
	kurt_thomson_est=0., kurt_thomson_var=0.;


      for(ipop=0; sumstat && sstats.pristats==1 && pars.cp.npop > 1 && ipop < pars.cp.npop; ipop++){
	priave_segs[ipop] = priave_thetapi[ipop] = priave_thetaw[ipop]=priave_fwh[ipop] 
	  = priave_tajimasD[ipop] = priave_ZnS[ipop] = priave_dvk[ipop] = priave_dvh[ipop] =0.;
        privar_segs[ipop] = privar_thetapi[ipop] = privar_thetaw[ipop] =0.;
	privar_tajimasD[ipop] = privar_fwh[ipop] = privar_ZnS[ipop] =0.;
	privar_dvk[ipop] = privar_dvh[ipop] = priave_thomson_est[ipop] = privar_thomson_est[ipop] = priave_thomson_var[ipop] = privar_thomson_var[ipop] = 0.;
	priskew_segs[ipop] = priskew_thetapi[ipop] = priskew_thetaw[ipop] = priskew_fwh[ipop] = priskew_tajimasD[ipop] = priskew_ZnS[ipop] = 0.;
	priskew_dvk[ipop] = priskew_dvh[ipop] = priskew_thomson_est[ipop] = priskew_thomson_var[ipop] = 0.;
	prikurt_segs[ipop] = prikurt_thetapi[ipop] = prikurt_thetaw[ipop] = prikurt_fwh[ipop] = prikurt_tajimasD[ipop] = prikurt_ZnS[ipop] = 0.;
	prikurt_dvk[ipop] = prikurt_dvh[ipop] = prikurt_thomson_est[ipop] = prikurt_thomson_var[ipop] = 0.;
      }
      
      
      for(npopi=0; npopi<pars.cp.npop; ++npopi){
	for(npopj=0; npopj<pars.cp.npop; ++npopj){
	  ave_shared[npopi][npopj] = 0.;
	  ave_private[npopi][npopj] = 0.;
	  ave_fixed_dif[npopi][npopj]= 0.;
	  var_shared[npopi][npopj] = 0.;
	  var_private[npopi][npopj] = 0.;
	  var_fixed_dif[npopi][npopj]= 0.;
	  ave_pairwise_Fst[npopi][npopj]=0.;
	  var_pairwise_Fst[npopi][npopj] = 0.;
	  skew_shared[npopi][npopj] = 0.;
	  skew_private[npopi][npopj] = 0.;
	  skew_fixed_dif[npopi][npopj] = 0.;
	  skew_pairwise_Fst[npopi][npopj] = 0.;
	  kurt_shared[npopi][npopj] = 0.;
	  kurt_private[npopi][npopj] = 0.;
	  kurt_fixed_dif[npopi][npopj] = 0.;
	  kurt_pairwise_Fst[npopi][npopj] = 0.;
	}
      }
           /* if the first replication then allocate memory etc */
      if(count <=1){

	/*
	  get the dimensions of the data
	*/
	
	if((fragin=fopen( fragfile ,"r"))==NULL){
	  fprintf(stderr, "\nERROR OPENNING THE FRAGMENT FILE\n");
	  EXIT_MSABC(1);
	}

	getfragdimensions(fragin, &nofragments, &fcolumns);
	
	fclose(fragin);

	/******** Dimensions OK **************************************/

	
	
	//allocate memory for the arrays

	
	
	/* memory for the name of the fragments */
	fname = (char**)malloc( nofragments * sizeof(char* ) );
	for( i=0; i<nofragments; i++)
	  fname[i] = (char*)malloc( 100 * sizeof(char) );
	
	/* for the size of the sample for each fragment */
	fsize = (int*)malloc( nofragments*sizeof(int) );
	frecbp = (double*)malloc( nofragments*sizeof(double));
	fmubp = (double*)malloc( nofragments*sizeof(double));
	flength = (int*)malloc( nofragments*sizeof(int));
	fpop = (int*)malloc( nofragments*sizeof(int));

	
	for(i=0; i<nofragments; i++)
	  fpop[i] = 1;

	if(pars.cp.npop > 1)
	  for(i=0; i<nofragments; i++)
	    fpop[i] = 0;
	
	/* open and read the file */
	if((fragin=fopen( fragfile ,"r"))==NULL){
	  fprintf(stderr, "\nERROR OPENNING THE FRAGMENT FILE\n");
	  EXIT_MSABC(1);
	}
	
	/* Store the information from a fragment file
	   input file, number of fragments, number of columns
	*/
	getfragpars( fragin, nofragments, fcolumns );

	/* Override mu/rec values if set by batch caller */
	if (frag_mu_override != NULL) {
	    int ov;
	    for (ov = 0; ov < nofragments && ov < frag_mu_override_len; ov++)
	        fmubp[ov] = frag_mu_override[ov];
	}
	if (frag_rec_override != NULL) {
	    int ov;
	    for (ov = 0; ov < nofragments && ov < frag_rec_override_len; ov++)
	        frecbp[ov] = frag_rec_override[ov];
	}

	for(ii=0; ii<MAXSUMSTAT; ii++)
	  total_stats_weights[ii] = 0;
	  
	
	/*
	  different modes of weights may be defined here
	*/
	fragindex=0;
	for(i=0; i<nofragments; i++){
	  /* check if fpop is always smaller or equal to pars.cp.npop */
	  if(fpop[i] > pars.cp.npop){
	    fprintf(stderr, "------------- ERROR -----------\n");
	    fprintf(stderr, "\nfpop at line %i is greater than the number of populations. This is not allowed.\n", i+1);
	    fprintf(stderr, "Please check that your command line, especially argument -I, and the locus file are correct\n\n");
	    EXIT_MSABC(1);
	  }
	    
	    
	  // skip if information is about the specific populations
	  if( fpop[i] && (fpop[i] != pars.cp.npop) )
	    continue;
	  
	  /* this is the real index of the fragment */
	  


	  /* calculate average (weighted) values for all the 
	     summary statistics.
	     For the moment I have set the weights to 1, so it's like a normal average
	  */
	  for(ii=0; ii<MAXSUMSTAT; ii++){
	    /* set the weights */
	    /* before it was only for the first dataset, which is not correct 
	       for the sstats_denominator_weights, It should be reseted every time
	    */
	    
	    if(weightmode == 1) // this is not used now
	      sstats_weights[ii] = (double)flength[i];
	  
	    /* and the denominator */
	    total_stats_weights[ii] += sstats_weights[ii];
	  }
	}
	
      } // if(count <=1)
      
      //      printf("count: %d\n", count); // $p
      /* initial values for theta and rho
	 they are needed in order to put the 
	 proper values of theta and rho when 
	 simulating each fragment
      */
      init_theta = pars.mp.theta;
      init_rho = pars.cp.r;
      init_nsites = pars.cp.nsites;
      
      //printf("count: %d\n", count); // (P
      // reset the initial denominator for the calculation of averages and variances
      for(ii=0; ii<MAXSUMSTAT; ii++){     
	sstats_denominator_weights[ii] = total_stats_weights[ii];
	for(ipop=0; ipop<pars.cp.npop && sstats.pristats==1 && pars.cp.npop>1; ipop++){
	  sstats_popden_weights[ii][ipop] = total_stats_weights[ii];
	  //printf("ipop: %d, pars.cp.npop: %d\n", ipop, pars.cp.npop);
	}

	/* pairwise statistics pp20100422 */
	if(pars.cp.npop > 1){
	  for(ipop =0; ipop<(pars.cp.npop); ipop++){
	    for(jpop = 0; jpop < pars.cp.npop; jpop++){
	      sstats_popdenpair_weights[ii][ipop][jpop] = total_stats_weights[ii];
	    }
	  }
	}
      }
      
      
      /* simulate each fragment  MAIN LOOP FOR THE FRAGMENTS *********************************************/
      /**
	 keep in mind that one line now does not mean one simulations.
	 It is possible that one line describes the 
	 features of fragmetns in different populations.
      */
      for(i=0; i<nofragments; i++){
	//fprintf(stderr, "fragment: %d\n", i);
	calculations_done = 0;
	pairwise_calculations_done = 0;
	/*
	  if specific information is given for the fragment in populations
	  then we should proceed differently
	*/
	if(fpop[i] != 0){
	  // the total sample size for this fragment 
	  f_total_sample_size += fsize[i];
	  (pars.cp.config)[ fpop[i] - 1 ] = fsize[i];
	  if(fpop[i] == pars.cp.npop){ // if information for the last population is read
	    pars.cp.nsam = f_total_sample_size;
	    f_total_sample_size = 0; // reset the total sample size info
	  }
	  else // go to the next round of the loop without doing the rest of the calculations
	    continue;
	}
	
	else{ // if information for the total sample is given
	  if(init_nsam != fsize[i]){
	    fprintf(stderr, "sample size in locus file does not match with sample size in command line. ERROR\n");
	    fprintf(stderr, "check the line %i in the locus file\n", i);
	    fprintf(stderr, "fsize: %d \t fpop: %d\n", fsize[i], fpop[i]);
	    EXIT_MSABC(1);
	  }
	  pars.cp.nsam = fsize[i];
	  for(npopi=0; npopi<pars.cp.npop; npopi++)
	    (pars.cp.config)[npopi] = (pars.cp.initconfig)[npopi];
	  
	} 	
    	
	/* Added by pavlos Jan 11 2009.  Re-allocates memory to list if new sample size is bigger than old 
	   when a larger sample size then 
	   reallocate the memory
	*/
	if( pars.cp.nsam > largestnsam && (largestnsam > 0) ){
	  
	  if( pars.mp.segsitesin ==	0 ){
	    list = rematrix(pars.cp.nsam, largestnsam, maxsites+1, list ); 
	  }
	  else{ 
	    list = rematrix(pars.cp.nsam, largestnsam, pars.mp.segsitesin+1, list );
	  }
	  largestnsam=pars.cp.nsam; 
	}
	
	
	/* recompute theta and rho for each fragment */
	for(j=2; j<=3; j++){ // j=2 is theta and j=3 is rho
	  /* 
	     In the case theta and rho are NOT
	     given in the file, then recalculate
	     from the command line, taking into account
	     the length of the fragment
	  */
	  
	  if((frpars[j] == 0)){
	
	    switch(j){
	    case 2:
	      /*
		divide the given theta with the number of nsites
		This is given in the command line in the rec parameter. 
	      */
	      pars.mp.theta = init_theta/(double)init_nsites * (double)flength[i];
	      break;
	    case 3:	      
	      pars.cp.r = init_rho*((double)flength[i]/(double)init_nsites);
	      break;
	    }
	  }
	  
	  /*
	    If the new theta and rho are given in 
	    the file, then consider these values
	    as the real ones
	  */
	  else if(frpars[j] != 0 ){
	    
	    switch(j){
	      
	    case 2:
	      if(effectiveN < 0){
		fprintf(stderr, " effective N should be specified in the fragment mode\n" );
		EXIT_MSABC(-1);
	      }
	      /*
		in the input file there is 
		theta per site. Maybe this should
		change in order to have theta
		per locus.
	      */
	      pars.mp.theta = 4. * fmubp[i] * flength[i] * effectiveN; // !! notice fthetabp is theta per site
	      break;
	    case 3:
	      if(effectiveN < 0){
		fprintf(stderr, " effective N should be specified in the fragment mode\n" );
		EXIT_MSABC(-1);
	      }
	      /*
		same explanations for the theta ( see above )
	      */
	      pars.cp.r = 4. * frecbp[i] * flength[i] * effectiveN;
	      break;
	    }
	  }
	}
	
	/*
	  use the fragments length as nsites
	*/
	pars.cp.nsites = flength[i];
	

	if(obin != NULL){
	  /*
	    !!!!!!!!!!!!!!!!!!!!   READ THE OBSERVED DATA  !!!!!!!!!!!!!!!!!!!!!!!!! 
	  */
	  segsites = getobservations(obin, &list);
	  
	  modifyList(outgroup, phasemode, list, segsites, pars.cp.nsam);
	  
	  if(segsites < 0){
	    fprintf(stderr, "segsites should be positive integer. it is %d\n", segsites);
	    EXIT_MSABC(1);
	  }
	}
	else{
	  /* 
	   !!!!!!!!!!!!!!!!!!!! SIMULATIONS FOR EACH LOCUS !!!!!!!!!!!!!!!!!!!!!!!!!! 
	  */
	  segsites = gensam( list, &probss, &tmrca, &ttot, largestnsam, old_maxnsam);
	  
	  modifyList(outgroup, phasemode, list, segsites, pars.cp.nsam);
	}

#ifdef R_PACKAGE_BUILD
	/* SFS accumulation mode in fragment loop: compute SFS and skip summary stats */
	if (sfs_mode_active) {
	  if (!sfs_expected_mode) {
	    sfs_accumulate(list, pars.cp.nsam, segsites, pars.cp.config, pars.cp.npop);
	  }
	  continue;  /* skip summary stats and text output for this fragment */
	}
#endif

	//add missing data if applicable
	if( (missingin != NULL) && (missingdata[fragindex++] > 0) ){
	  if( ! put_missing(list, posit, largestnsam, segsites, misseq,  misstart, misend, missingdata[fragindex-1], fragindex-1)){
	    fprintf(stderr, "Couldn't embed missing data\n");
	  }
	}
	
	if( old_maxnsam < pars.cp.nsam)
	  old_maxnsam = pars.cp.nsam;
		
	/* print the polymorphic table */
	/* if we are working on the observed data don't print the data */
	if( (obin == NULL) && tempdatafile!=NULL ){ 

	  int k;
	  
	  fprintf(tempdatafile,"\n//");
	  fprintf(tempdatafile, "\t%d\t%d\t%e\t%e\t%d", count, i, pars.mp.theta, pars.cp.r, flength[i]);
	  fprintf(tempdatafile, "\n");
	  
	  
	  if( pars.mp.timeflag ) fprintf(pf,"*time:\t%lf\t%lf\n",tmrca, ttot ) ;
	  if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
	    if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) 
	      fprintf(tempdatafile,"prob: %g\n", probss ) ;
	    
	    fprintf(tempdatafile,"segsites: %d\n",segsites);
	    
	    if( segsites > 0 )	fprintf(tempdatafile,"positions: ");
	    for( k=0; k<segsites; k++)
	      /* fprintf(pf,"%6.4lf ",posit[i] ); */
	      fprintf(tempdatafile,"%e ",posit[k] );
	    /* fprintf(pf,"\n"); */
	    fprintf(tempdatafile,"\n");
	    
	    for(k=0;k<pars.cp.nsam; k++) 
	      /* fprintf(pf,"%s\n", list[i] );  */
	      fprintf(tempdatafile, "%s\n", list[k] );
	  }
	}
	

	/* 
	   print the summary statistics 
	*/
	if(sumstat){ // This is the default that is used here
	  //time(&t3);
	  int insegsites = (segsites > 0) ? segsites:1;
	  if(insegsites > largestinsegsites || pars.cp.nsam > largesample){
	    
	    if(insegsites > largestinsegsites){
	      /* array that stores the class of SNP e.g. 1 for singleton, 2 for doubleton etc */
	      allelic_map = (int*)realloc(allelic_map, (unsigned)(insegsites*sizeof(int)));
	      missing_map = (int*)realloc(missing_map, (unsigned)(insegsites*sizeof(int)));
	      
	      for( j = 0; sstats.pristats==1 && j<pars.cp.npop && pars.cp.npop > 1; j++){
		/* the same as the allelic map, but now for each population */
		npop_allelic_map[j] = (int*)realloc(npop_allelic_map[j], (unsigned)(insegsites*sizeof(int)));
		npop_missing_map[j] = (int*)realloc(npop_missing_map[j], (unsigned)(insegsites*sizeof(int)));
	      }
	      largestinsegsites = insegsites;
	    }
	    
	    if( pars.cp.nsam <= largesample){
	      //printf("realloc ders[i]\n");
	      for( j=0; j<largesample; j++) {
		if(   ! (  ders[j] = (int *) realloc( ders[j], (unsigned) insegsites*sizeof(int) )  )   ) {
		  perror("realloc error in rematrix. ders");
		}
	      }
	      for(ipop=0; sstats.pristats ==1 && ipop <pars.cp.npop && pars.cp.npop > 1; ipop++){
		for( j=0; j<largesample; j++) {
		  if(   ! (  npop_ders[ipop][j] = 
			     (int *) realloc( npop_ders[ipop][j], (unsigned) insegsites*sizeof(int) )  )   ) {
		    perror("realloc error in rematrix. npop_ders");
		  }
		}
	      }
	      
	    }
	    
	    if( (pars.cp.nsam > largesample)){
	      //printf("realloc ders\n");
	      ders = rematrix_int(pars.cp.nsam, largesample, largestinsegsites+1, ders );
	      
	      for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++)
		npop_ders[ipop] = rematrix_int(pars.cp.nsam, largesample, largestinsegsites+1, npop_ders[ipop]);
	      
	      largesample = pars.cp.nsam;
	    }
	  }  
	  	  
	  int segs;
	  double thetapi=0., thetaw = 0., tajd = 0., thetah = 0., fwh = 0., zns=0., dvk=0., dvh=0., thomson_est=0., thomson_var = 0.;
	  for(ipop=0; ipop<pars.cp.npop && sumstat && pars.cp.npop > 1 && sstats.pristats == 1 ;  ipop++){
	    prisegs[ipop] = prithetapi[ipop] = prithetaw[ipop] = pritajd[ipop] 
	      = prithetah[ipop] = prifwh[ipop] = prizns[ipop] = pridvk[ipop] = pridvh[ipop] = 0.;
	  }
	  
	  
	  /* 
	     for the denominators of the statistics
	  */
	  
	  if(pars.cp.nsam != old_nsam){
	    /* denominators of Tajima's D, H etc */
	    denominators(pars.cp.nsam, &hn, &sqhn, &bn);
	  }
	  old_nsam = pars.cp.nsam;
	  
	  
	  for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++)
	    denominators( (pars.cp.config)[ipop], &prihn[ipop], &prisqhn[ipop], &pribn[ipop]);
	 

	  //fprintf(stderr, "pars.cp.nsam: %d\n", pars.cp.nsam); //20100712

	  set_almap(allelic_map, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);
	  set_missing(missing_map, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);
	  //printf("ok\n"); //.;
	  set_ders(ders, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);


	  samples_b[0] = 0;
	  samples_e[0] = (pars.cp.config)[0];
	  for(ipop=1; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++){
	    samples_b[ipop] = samples_b[ipop-1] + (pars.cp.config)[ipop-1];
	    samples_e[ipop] = samples_b[ipop] + (pars.cp.config)[ipop];
	  }
	  /* for(ipop=0; ipop<pars.cp.npop && pars.cp.npop > 1 && sstats.pristats == 1; ipop++){ */
/* 	    printf("pop: %d, b:%d, e:%d, config:%d\n", ipop+1, samples_b[ipop], samples_e[ipop], pars.cp.config[ipop]); */
/* 	  } */
	  
	  for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++){
	    set_almap(npop_allelic_map[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
	    set_missing(npop_missing_map[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
	    set_ders(npop_ders[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
	  }
	  
	  int truesegsites = seg_sites(allelic_map, pars.cp.nsam, segsites);
	  
	  
	  for(ipop=0;  sstats.pristats==1 && pars.cp.npop> 1 && ipop < pars.cp.npop; ipop++){
	    pritruesegsites[ipop] = seg_sites(npop_allelic_map[ipop], 
					      (pars.cp.config)[ipop],
					      segsites);
	  }
	  
	  /* CALCULATE SUMMARY STATISTICS PER LOCUS*/
	  int ii;
	  for(ii=0; ii<MAXSUMSTAT; ii++){
	    



	    //printf("-------------------------------, stat: %d, fragment %d\n", ii, i);
	    /* this is a flag that denotes if a certain summary statistic is used */
	    if(sstats_array[ii] == 1){
	      //printf("sstats %d\n", ii);
	      switch(ii){
	      case 0: // segsites
		if(sstats.pristats == 1 && pars.cp.npop > 1 ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    if( (pars.cp.config)[ipop] < MINSEQ ){
		      pritruesegsites[ipop] = 0;
		    }
		    prisegs[ipop] = pritruesegsites[ipop];
		    priave_segs[ipop] += sstats_weights[ii] * prisegs[ipop];
		    privar_segs[ipop] += prisegs[ipop] * prisegs[ipop];
		    priskew_segs[ipop] += prisegs[ipop] * prisegs[ipop] * prisegs[ipop];
		    prikurt_segs[ipop] += prisegs[ipop] * prisegs[ipop] * prisegs[ipop] * prisegs[ipop];
		    //printf("segsites for pop %d is %e\n", ipop, prisegs[ipop]);
		  }
		}
		segs = truesegsites;
		ave_segs += sstats_weights[ii] * (double)segs;
		var_segs += (double)segs*(double)segs;
		skew_segs += (double)segs*(double)segs*(double)segs;
		kurt_segs += (double)segs*(double)segs*(double)segs*(double)segs;

		break;
	      case 1: // theta pi
		if(sstats.pristats == 1 && pars.cp.npop > 1 ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    if( (pars.cp.config)[ipop] >= MINSEQ){
		      theta_pi(&prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]);
		      priave_thetapi[ipop] += sstats_weights[ii] * prithetapi[ipop];
		      privar_thetapi[ipop] += prithetapi[ipop] * prithetapi[ipop];
		      priskew_thetapi[ipop] += prithetapi[ipop] * prithetapi[ipop] * prithetapi[ipop];
		      prikurt_thetapi[ipop] += prithetapi[ipop] * prithetapi[ipop] * prithetapi[ipop] * prithetapi[ipop];
		    }
		    else
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		  }
		}

		if( theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map)){
		  ave_thetapi += sstats_weights[ii] * thetapi;
		  var_thetapi += thetapi*thetapi;
		  skew_thetapi += thetapi*thetapi*thetapi;
		  kurt_thetapi += thetapi*thetapi*thetapi*thetapi;
		}
		break;
	      case 2: // theta w
		
		if(sstats.pristats == 1 && pars.cp.npop > 1){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    theta_w(&prithetaw[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, prihn[ipop], npop_missing_map[ipop]);
		      priave_thetaw[ipop] += sstats_weights[ii] * prithetaw[ipop];
		      privar_thetaw[ipop] += prithetaw[ipop] * prithetaw[ipop];
		      priskew_thetaw[ipop] += prithetaw[ipop] * prithetaw[ipop] * prithetaw[ipop];
		      prikurt_thetaw[ipop] += prithetaw[ipop] * prithetaw[ipop] * prithetaw[ipop] * prithetaw[ipop];

		  }
		}

		if( theta_w(&thetaw, allelic_map, pars.cp.nsam, segsites, hn, missing_map)){
		  //thetaw /= (double)flength[i];
		  ave_thetaw += sstats_weights[ii] * thetaw;
		  var_thetaw += thetaw*thetaw;
		  skew_thetaw += thetaw*thetaw*thetaw;
		  kurt_thetaw += thetaw*thetaw*thetaw*thetaw;
		}

		break;
	      case 3: /* Tajima's D*/
		
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.pritajimasD ){
		  
		  for(ipop=0; ipop<pars.cp.npop; ipop++){

		    //if thetaw and thetapi are known
		    if( (prithetaw[ipop] && prithetapi[ipop] && pritruesegsites[ipop] ) &&
			tajD(&pritajd[ipop], pritruesegsites[ipop], (pars.cp.config)[ipop], 
			     prithetaw[ipop], prithetapi[ipop], prihn[ipop], prisqhn[ipop]) )
		      {
			
			priave_tajimasD[ipop] += sstats_weights[ii]*pritajd[ipop];
			privar_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop];
			priskew_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop]*pritajd[ipop];
			prikurt_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop]*pritajd[ipop]*pritajd[ipop];
			//printf("count %d, tajimas D for pop %d is %e\n", count, ipop, pritajd[ipop]);
		      }
		    else if(pritruesegsites[ipop] && (!sstats.prithetaw || !sstats.prithetapi)){

		      if( theta_w( &prithetaw[ipop], npop_allelic_map[ipop],
				   (pars.cp.config)[ipop], segsites, prihn[ipop], npop_missing_map[ipop]) &&
			  theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
			  tajD(&pritajd[ipop], pritruesegsites[ipop], (pars.cp.config)[ipop],
			       prithetaw[ipop], prithetapi[ipop], prihn[ipop], prisqhn[ipop]) )
			{

			  priave_tajimasD[ipop] += sstats_weights[ii]*pritajd[ipop];
			  privar_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop];
			  priskew_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop]*pritajd[ipop];
			  prikurt_tajimasD[ipop] += pritajd[ipop]*pritajd[ipop]*pritajd[ipop]*pritajd[ipop];
			}
		    }
		    else if(pritruesegsites[ipop] == 0) // Tajimas' D isn't defined
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		    
		  }
		}
		
		if( (thetaw && thetapi && truesegsites) &&
		    tajD(&tajd, truesegsites, pars.cp.nsam, thetaw, thetapi, hn, sqhn))
		  {
		    ave_tajimasD += sstats_weights[ii] * tajd;
		    var_tajimasD += tajd*tajd;
		    skew_tajimasD += tajd*tajd*tajd;
		    kurt_tajimasD += tajd*tajd*tajd*tajd;
		  }
		else if(truesegsites && (!thetaw || !thetapi)){
		  if(theta_w(&thetaw, allelic_map, pars.cp.nsam, segsites, hn, missing_map) &&
		     theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) &&
		     tajD(&tajd, truesegsites, pars.cp.nsam, thetaw, thetapi, hn, sqhn)
		     ){
		    ave_tajimasD += sstats_weights[ii] * tajd;
		    var_tajimasD += tajd*tajd;
		    skew_tajimasD += tajd*tajd*tajd;
		    kurt_tajimasD += tajd*tajd*tajd*tajd;
		  }
		}
		else if(truesegsites == 0){
		  sstats_denominator_weights[ii] -= sstats_weights[ii];
		}

		break;
		/*
		  calculate ZnS
		*/
	      case 4:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.priZnS){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    
		    if(pritruesegsites[ipop] > 1 && 
		       ZnS(&prizns[ipop], list, (pars.cp.config)[ipop], segsites, // segsites is the dimension of almap
			   filter, npop_allelic_map[ipop], npop_ders[ipop], npop_missing_map[ipop]) )
		      {
			priave_ZnS[ipop] += sstats_weights[ii] * prizns[ipop];
			privar_ZnS[ipop] += prizns[ipop] * prizns[ipop];
			priskew_ZnS[ipop] += prizns[ipop] * prizns[ipop] * prizns[ipop];
			prikurt_ZnS[ipop] += prizns[ipop] * prizns[ipop] * prizns[ipop] * prizns[ipop];
		      }
		    else if(pritruesegsites[ipop] <= 1)
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		  }
		}

		if(truesegsites > 1 && ZnS(&zns, list, pars.cp.nsam, truesegsites, filter, allelic_map, ders, missing_map)){
		  ave_ZnS += sstats_weights[ii] * zns;
		  var_ZnS += zns*zns;
		  skew_ZnS += zns*zns*zns;
		  kurt_ZnS += zns*zns*zns*zns;
		}
		else
		  sstats_denominator_weights[ii] -= sstats_weights[ii];

		break;
		
		/* 
		   in the case of multiple populations 
		   calculate the Fst
		*/
	      case 5:
		
		if(truesegsites > 0 &&
		   pars.cp.npop > 1 ){
		  
		  calculations(weights, 
			       list,
			       pars.cp.config,
			       pars.cp.npop,
			       truesegsites,
			       pars.cp.nsam,
			       
			       &piT,
			       &piS,
			       &piB,
			       &piD,
			       
			       shared,
			       private,
			       fixed_dif,
			       derived);
		   calculations_done=1;
		   fst = Fst(fst_type, piT, piS, piD, piB);
		   ave_fst += sstats_weights[ii] * fst;
		   var_fst += fst * fst;
		   skew_fst += fst * fst * fst;
		   kurt_fst += fst * fst * fst * fst;
		   //printf("Fst is %e\n", fst_hbk);
		 }
		 else
		   sstats_denominator_weights[ii] -= sstats_weights[ii];
		 break;
	       case 6: //shared polymorphism between populations
		 if(pars.cp.npop > 1 ){
		  
		  if(calculations_done != 1){
		    calculations(weights, 
				 list,
				 pars.cp.config,
				 pars.cp.npop,
				 truesegsites,
				 pars.cp.nsam,
				 
				 &piT,
				 &piS,
				 &piB,
				 &piD,
				 
				 shared,
				 private,
				 fixed_dif,
				 derived);
		    calculations_done = 1;
		  }
		  
		  for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      /* pp20100422 if there is no sample take care of the denominator */
		      if( (pars.cp.config)[npopi] > MINSEQ && (pars.cp.config)[npopj] > MINSEQ &&
			  !isnan(shared[npopi][npopj]) ){
			ave_shared[npopi][npopj] += shared[npopi][npopj];
			var_shared[npopi][npopj] += shared[npopi][npopj]*shared[npopi][npopj];
			skew_shared[npopi][npopj] += shared[npopi][npopj]*shared[npopi][npopj]*shared[npopi][npopj];
			kurt_shared[npopi][npopj] += shared[npopi][npopj]*shared[npopi][npopj]*shared[npopi][npopj]*shared[npopi][npopj];
		      }
		      else
			sstats_popdenpair_weights[ii][npopi][npopj] -= sstats_weights[ii];
		    }
		  }
		 }
		 else{
		   sstats_denominator_weights[ii] -= sstats_weights[ii];
		 }
		 break;

	      case 7: // Proportion of private polymorphisms between two populations
		
		if(pars.cp.npop > 1 && sstats.fixed_dif ){
		  if(calculations_done != 1){
		    calculations(weights, 
				 list,
				 pars.cp.config,
				 pars.cp.npop,
				 truesegsites,
				 pars.cp.nsam,
				 
				 &piT,
				 &piS,
				 &piB,
				 &piD,
				 
				 shared,
				 private,
				 fixed_dif,
				 derived);
		    calculations_done = 1;
		  }
		  
		  for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		    /* pp20100425 if there is no sample take care of the denominator */
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      if( (pars.cp.config)[npopi] > MINSEQ && (pars.cp.config)[npopj] > MINSEQ
			  && !isnan(private[npopi][npopj]) ){
			ave_private[npopi][npopj] += private[npopi][npopj];
			var_private[npopi][npopj] += private[npopi][npopj]*private[npopi][npopj];
			skew_private[npopi][npopj] += private[npopi][npopj]*private[npopi][npopj]*private[npopi][npopj];
			kurt_private[npopi][npopj] += private[npopi][npopj]*private[npopi][npopj]*private[npopi][npopj]*private[npopi][npopj];
		      }
		      else
			sstats_popdenpair_weights[ii][npopi][npopj] -= sstats_weights[ii];
		    }
		  }
		}
		else{
		  sstats_denominator_weights[ii] -= sstats_weights[ii];
		}
		break;
		
	      case 8: // Proportion of fixed polymorphisms between two populations
		
		if(pars.cp.npop > 1){
		  if(calculations_done != 1){
		    calculations(weights, 
				 list,
				 pars.cp.config,
				 pars.cp.npop,
				 truesegsites,
				 pars.cp.nsam,
				 
				 &piT,
				 &piS,
				 &piB,
				 &piD,
				 
				 shared,
				 private,
				 fixed_dif,
				 derived);
		    calculations_done = 1;
		  }
		  
		  for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		    /* pp20100425 if there is no sample take care of the denominator */
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      if( (pars.cp.config)[npopi] > MINSEQ&& (pars.cp.config)[npopj] > MINSEQ
			  && !isnan(fixed_dif[npopi][npopj])){
			ave_fixed_dif[npopi][npopj] += fixed_dif[npopi][npopj];
			var_fixed_dif[npopi][npopj] += fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj];
			skew_fixed_dif[npopi][npopj] += fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj];
			kurt_fixed_dif[npopi][npopj] += fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj]*fixed_dif[npopi][npopj];
		      }
		      else
			sstats_popdenpair_weights[ii][npopi][npopj] -= sstats_weights[ii];
		    }
		  }
		}
		else{
		  sstats_denominator_weights[ii] -= sstats_weights[ii];
		}
		
		break;
	      case 9: // Pairwise Fst
		
		//printf("fst_pops: %d\n", fst_pops);
		if(pars.cp.npop > 1 && pairwise_calculations_done != 1){
		  for(ipop=0; ipop<fst_pops-1; ipop++){
		    /* pp20100425 if there is no sample take care of the denominator */
		    for(jpop=ipop+1; jpop<fst_pops; jpop++){
		      if( (pars.cp.config)[ipop] > MINSEQ && (pars.cp.config)[jpop] > MINSEQ ){
			pop1=Fst_pops[ ipop ];
			pop2=Fst_pops[ jpop ];
			
			pairwiseFstcalculations(
						pop1, pop2,
						weights,
						list,
						pars.cp.config,
						pars.cp.npop,
						truesegsites,
						pars.cp.nsam,
						
						&ppiT,
						&ppiS,
						&ppiB,
						&ppiD);
			
			pairwise_fst = Fst(fst_type, ppiT, ppiS, ppiD, ppiB);
			if( !isnan(pairwise_fst) ){
			  ave_pairwise_Fst[pop1][pop2] += sstats_weights[ii]*pairwise_fst;
			  var_pairwise_Fst[pop1][pop2] += pairwise_fst * pairwise_fst;
			  skew_pairwise_Fst[pop1][pop2] += pairwise_fst * pairwise_fst * pairwise_fst;
			  kurt_pairwise_Fst[pop1][pop2] += pairwise_fst * pairwise_fst * pairwise_fst * pairwise_fst;
			}
			else{
			
			  sstats_popdenpair_weights[ii][pop1][pop2] -= sstats_weights[ii];
			}
		      }
		      else
			sstats_popdenpair_weights[ii][pop1][pop2] -= sstats_weights[ii];
		    }
		    
		    pairwise_calculations_done = 1;
		  }
		}
		else{
		  sstats_denominator_weights[ii] -= sstats_weights[ii];
		}
		
		break;		  
		
	      case 10: /* Fay Wu H */
		
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prifwh){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){

		    // if we know the thetapi and thetah
		    if( (pritruesegsites[ipop] && prithetapi[ipop] && prithetah[ipop] ) &&
			htest( &prifwh[ipop], (pars.cp.config)[ipop], 
				pritruesegsites[ipop], // the number of true segsites is needed 
			       prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
		      {
			priave_fwh[ipop] += sstats_weights[ii] * prifwh[ipop];
			privar_fwh[ipop] += prifwh[ipop]*prifwh[ipop];
			priskew_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			prikurt_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
		      }
		    else if(pritruesegsites[ipop] && prithetapi[ipop] && (!prithetah[ipop]) ){
		      if(thetaH( &prithetah[ipop], list,
				 pars.cp.nsam, // here the pars.cp.nsam is the dimension of the list
				 segsites, // this is the other dimension of the list
				 bm,
				 npop_allelic_map[ipop], samples_b[ipop], samples_e[ipop], npop_missing_map[ipop])  &&
			 htest( &prifwh[ipop], (pars.cp.config)[ipop],
				pritruesegsites[ipop], // the number of true segsites is needed
				prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
			{
			  priave_fwh[ipop] += sstats_weights[ii] * prifwh[ipop];
			  privar_fwh[ipop] += prifwh[ipop]*prifwh[ipop];
			  priskew_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			  prikurt_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			}
		    }
		    else if(pritruesegsites[ipop] && prithetapi[ipop] && (!prithetah[ipop]) ){
		      if(theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
			 htest( &prifwh[ipop], (pars.cp.config)[ipop],
				pritruesegsites[ipop], // the number of true segsites is needed
				prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
			{
			  priave_fwh[ipop] += sstats_weights[ii] * prifwh[ipop];
			  privar_fwh[ipop] += prifwh[ipop]*prifwh[ipop];
			  priskew_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			  prikurt_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			}
		    }
		    else if(pritruesegsites[ipop] && (!prithetapi[ipop]) && (!prithetah[ipop]) ){
		      if( theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
			  thetaH( &prithetah[ipop], list,
				  pars.cp.nsam, // here the pars.cp.nsam is the dimension of the list
				  segsites, // this is the other dimension of the list
				  bm,
				  npop_allelic_map[ipop], samples_b[ipop], samples_e[ipop], npop_missing_map[ipop] ) &&
			  htest( &prifwh[ipop], (pars.cp.config)[ipop],
				 pritruesegsites[ipop], // the number of true segsites is needed
				 prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
			{
			  priave_fwh[ipop] += sstats_weights[ii] * prifwh[ipop];
			  privar_fwh[ipop] += prifwh[ipop]*prifwh[ipop];
			  priskew_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			  prikurt_fwh[ipop] += prifwh[ipop]*prifwh[ipop]*prifwh[ipop]*prifwh[ipop];
			}
		    }
		    else if(pritruesegsites[ipop] <= 0)
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		    //printf("fwh:%e, priave_fwh[%d]:%e\n", prifwh[ipop], ipop, priave_fwh[ipop]);
		  }
		}
			    

		
		if( ( truesegsites && thetapi && thetah) &&
		    htest( &fwh, pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		  ave_fwh += sstats_weights[ii] * fwh;
		  var_fwh += fwh*fwh;
		  skew_fwh += fwh*fwh*fwh;
		  kurt_fwh += fwh*fwh*fwh*fwh;
		}
		else if( truesegsites && thetapi && (!thetah) ){
		  if(thetaH( &thetah, list, pars.cp.nsam, truesegsites, bm, allelic_map, 0, pars.cp.nsam, missing_map) &&
		     htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		    ave_fwh += sstats_weights[ii] * fwh;
		    var_fwh += fwh*fwh;
		    skew_fwh += fwh*fwh*fwh;
		    kurt_fwh += fwh*fwh*fwh*fwh;
		  }
		}
		else if( truesegsites && (!thetapi) && (thetah) ){
		  if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) &&
		     htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		    ave_fwh += sstats_weights[ii] * fwh;
		    var_fwh += fwh*fwh;
		    skew_fwh += fwh*fwh*fwh;
		    kurt_fwh += fwh*fwh*fwh*fwh;
		  }
		}
		else if( truesegsites && (!thetapi) && (!thetah) ){
		  if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) &&
		     thetaH(&thetah, list, pars.cp.nsam, truesegsites, bm, allelic_map, 0, pars.cp.nsam, missing_map) &&
		     htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		    ave_fwh += sstats_weights[ii] * fwh;
		    var_fwh += fwh*fwh;
		    skew_fwh += fwh*fwh*fwh;
		    kurt_fwh += fwh*fwh*fwh*fwh;
		  }
		}
		else if(truesegsites == 0){
		  sstats_denominator_weights[ii] -= sstats_weights[ii];
		}
		break;
		     
		  
	      case 11: /* Depaulis Veuille H, and K */
		
		if(sstats.pristats == 1 && pars.cp.npop > 1 && (sstats.pridvstat)  ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    if((pars.cp.config)[ipop] > MINSEQ){
		      dvstat(list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop], &pridvk[ipop], &pridvh[ipop]);
		      
		      priave_dvk[ipop] += sstats_weights[ii] * pridvk[ipop];
		      priave_dvh[ipop] += sstats_weights[ii] * pridvh[ipop];

		      privar_dvk[ipop] += pridvk[ipop] * pridvk[ipop];
		      privar_dvh[ipop] += pridvh[ipop] * pridvh[ipop];
		      priskew_dvk[ipop] += pridvk[ipop] * pridvk[ipop] * pridvk[ipop];
		      priskew_dvh[ipop] += pridvh[ipop] * pridvh[ipop] * pridvh[ipop];
		      prikurt_dvk[ipop] += pridvk[ipop] * pridvk[ipop] * pridvk[ipop] * pridvk[ipop];
		      prikurt_dvh[ipop] += pridvh[ipop] * pridvh[ipop] * pridvh[ipop] * pridvh[ipop];
		    }
		    else{
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		    }
		  }
		  
		}
		if( dvstat(list, pars.cp.nsam, segsites, 0, pars.cp.nsam, &dvk, &dvh)){
		  ave_dvk += sstats_weights[ii] * dvk;
		  ave_dvh += sstats_weights[ii] * dvh;
		  var_dvk += dvk*dvk;
		  var_dvh += dvh*dvh;
		  skew_dvk += dvk*dvk*dvk;
		  skew_dvh += dvh*dvh*dvh;
		  kurt_dvk += dvk*dvk*dvk*dvk;
		  kurt_dvh += dvh*dvh*dvh*dvh;
		}

		break;
		
	      case 12:
		if(sstats.pristats == 1 && pars.cp.npop > 1 ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    if( 
		       (pars.cp.config)[ipop] >= MINSEQ &&
		       thomsonEst(&prithomson_est[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop])
			)
		      
		      {
			priave_thomson_est[ipop] += sstats_weights[ii] * prithomson_est[ipop];
			privar_thomson_est[ipop] += prithomson_est[ipop] * prithomson_est[ipop];
			priskew_thomson_est[ipop] += prithomson_est[ipop] * prithomson_est[ipop] * prithomson_est[ipop];
			prikurt_thomson_est[ipop] += prithomson_est[ipop] * prithomson_est[ipop] * prithomson_est[ipop] * prithomson_est[ipop];
		      }
		    
		    else
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		  }
		}
		
		if( thomsonEst(&thomson_est, allelic_map, pars.cp.nsam, segsites, missing_map)){
		  ave_thomson_est += sstats_weights[ii] * thomson_est;
		  var_thomson_est += thomson_est*thomson_est;
		  skew_thomson_est += thomson_est*thomson_est*thomson_est;
		  kurt_thomson_est += thomson_est*thomson_est*thomson_est*thomson_est;
		}
		break;

		case 13:
		if(sstats.pristats == 1 && pars.cp.npop > 1 ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    if( (pars.cp.config)[ipop] >= MINSEQ){
		      thomsonVar(&prithomson_var[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites);
		      priave_thomson_var[ipop] += sstats_weights[ii] * prithomson_var[ipop];
		      privar_thomson_var[ipop] += prithomson_var[ipop] * prithomson_var[ipop];
		      priskew_thomson_var[ipop] += prithomson_var[ipop] * prithomson_var[ipop] * prithomson_var[ipop];
		      prikurt_thomson_var[ipop] += prithomson_var[ipop] * prithomson_var[ipop] * prithomson_var[ipop] * prithomson_var[ipop];
		    }
		    else
		      sstats_popden_weights[ii][ipop] -= sstats_weights[ii];
		  }
		}
	      
		if( thomsonVar(&thomson_var, allelic_map, pars.cp.nsam, segsites)){
		  ave_thomson_var += sstats_weights[ii] * thomson_var;
		  var_thomson_var += thomson_var*thomson_var;
		  skew_thomson_var += thomson_var*thomson_var*thomson_var;
		  kurt_thomson_var += thomson_var*thomson_var*thomson_var*thomson_var;
		}
		break;
	
		
		
	      default:
		fprintf(stderr, "Statistic is not implemented yet, case: %d!\n", ii);
		break;
	      }
	      
	    }
	    
	  }
	}
      }
      
       
      /*
	print the average and the variance of the summary statistics
      */
      if(sumstat){
	if(count == 1){
	  for(i=3; i<argc; i++){
	    //fprintf(stderr, "argv[%d]: %s\n", i, argv[i]); 
	    if(argv[i][0] != '-') continue;
	    switch(argv[i][1]){
	    case '-':
	      if(!strcmp(argv[i], "--N"))
		if(printall || param_distr[i] != '-')
		  fprintf(pf, "p_Ne\t");
	      break;
	    case 't':
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_theta\t");
	      break;
	    case 'n':
	      if(argv[i][2] != '\0') break;
	      i+=1; /* i+1 because the first one is the name of the population */
	      if(printall || param_distr[i] != '-') 
		fprintf(pf, "p_init_size_pop_%d\t", atoi(argv[i]) );
	      break;
	    case 'g':
	      if(argv[i][2] != '\0') break;
	      i+=1; /* i+1 because the first one is the name of the population */
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_g_pop_%d\t", atoi(argv[i]) );
	      break;
	    case 'r':
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_rho\t");
	      break;
	    case 's':
	      if(printall || param_distr[i] != '-'){ 
		if( argv[i][2] != 'e' ){
		  fprintf(pf, "p_insegs\t");
		}
		else{
		  i+=nseeds;
		}
	      }
	      break;
	    case 'I':
	      c=0;
	      if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_npop\t");
	      }
	      i += pars.cp.npop+1 ; // skip the number of populations (done already) and the configuration
	      if( (printall && (i+1 < argc) && argv[i+1][0] != '-') ||
		  ( (i+1 < argc) && (c = dist_argcheck(i+1,argc,argv)) ) )
		{
		  i++;
		  fprintf(pf, "p_isl_totmig\t");
		}
	      c = 0;
	      break;
	    case 'm':
	      c=0;
	      
	      if( argv[i][2] == 'a' ){
		i++;
		
		for(pop = 0; pop < pars.cp.npop; pop++){
		  for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		    c = 0;
		   
		    if( (c=  dist_argcheck(i, argc, argv) ) || printall ){
			fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
			i+=c;
			if(c) continue;
		    }
		    i++; 
		  }
		}
		i--; // go one back because the index will be increased in the initial loop
	      }
	      else{
		c=0;
		i++;
		pop = atoi(argv[i++])-1;
		pop2 = atoi(argv[i++])-1;
		if( (c = dist_argcheck(i, argc, argv) ) || printall ){
		    fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
		    if(c) i+=c-1;
		}
	      }
	      break;
	    case 'G':
	      if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_G\t");
	      }
	      break;
	    case 'e':
	      secarg = argv[i][2];
	      switch(secarg){
	      case 'M':
		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "p_globalMigRate_change_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_change_globalMigRate\t");
		}
		break;
	      case 'N':
		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "p_pop_change_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "p_pop_change_newpop\t");
	      }
	      break;
	      case 'n':
		if(argv[i][3] != '\0') break;
		z = i; 
		z += distribution_arguments( param_distr[i] ) + 2;
		zpop = atoi(argv[z]);
		
		//fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_en_time_pop_%d\t", zpop);
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		i++;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_en_size_pop_%d\t", atoi(argv[i]) );
		}
		break;
	      case 'g':
		if(argv[i][3] != '\0') break;
		z = i; 
		z += distribution_arguments( param_distr[i] ) + 2;
		zpop = atoi(argv[z]);
		
		//fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eg_time_pop_%d\t", zpop);
		}
		
		i += distribution_arguments( param_distr[i] ) + 1;
		i++;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eg_pop_%d\t", atoi(argv[i]) );
		}
		break;
	      case 'G':
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eG_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eG_rate\t");
		}
		break;
	      case 'j':
		c = distribution_arguments( param_distr[i] );
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_ej_time_pop_%d_%d\t", atoi(argv[i+2+c]), atoi(argv[i+3+c]) );
		}
		c=0;
		break;
	      case 's':
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_es_time\t");
		}
		break; 
	      case 'm':
		
		c = 0;
		if( argv[i][3] == 'a'){ // -ema
		
		  if(printall || param_distr[i] != '-'){
		      fprintf(pf, "p_ema_time\t");
		  }
		  i += distribution_arguments( param_distr[i] ) + 1;
		  i+=2; /* the time */
		  for(pop = 0; pop < pars.cp.npop; pop++){
		    for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		      c = 0;
		      if( (c=  dist_argcheck(i, argc, argv) ) || printall){
			fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
			i+=c;
			if(c) continue;
		      }
		      i++;
		    }
		  }
		  i--; // go one back because the index will be increased in the initial loop
		}
		else{
		  z = i; 
		  z += distribution_arguments( param_distr[i] ) + 2;
		  zpop1 = atoi(argv[z]);
		  zpop2 = atoi(argv[z+1]);
		//printf("here\n"); 
		  
		  if(printall || param_distr[i] != '-'){
		    fprintf(pf, "p_em_time_pops_%d_%d\t", zpop1, zpop2);
		  }
		  i += distribution_arguments( param_distr[i] ) + 1;
		  c=0;
		  pop = atoi(argv[++i])-1;
		  pop2 = atoi(argv[++i])-1;
		  
		  if(printall || (c = dist_argcheck(i, argc, argv) ) ){
		      fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
		      if(c) i+=c-1;
		  }
		}
		break;
	      }
	    }
	  }

	  /* print headers for sumstat */
	  for(i=0; i<MAXSUMSTAT; i++){
	    if(sstats_array[i] == 1){
	      switch(i){
	      case 0: // segsites

		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prisegs){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_segs_%d\t", ipop+1);
		    fprintf(pf, "s_var_segs_%d\t", ipop+1);
		    fprintf(pf, "s_skew_segs_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_segs_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_segs\t");
		fprintf(pf, "s_var_segs\t");
		fprintf(pf, "s_skew_segs\t");
		fprintf(pf, "s_kurt_segs\t");
		break;

	      case 1: // theta pi
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetapi){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_pi_%d\t", ipop+1);
		    fprintf(pf, "s_var_pi_%d\t", ipop+1);
		    fprintf(pf, "s_skew_pi_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_pi_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_pi\t");
		fprintf(pf, "s_var_pi\t");
		fprintf(pf, "s_skew_pi\t");
		fprintf(pf, "s_kurt_pi\t");
		break;

	      case 2: // theta w
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetaw){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_w_%d\t", ipop+1);
		    fprintf(pf, "s_var_w_%d\t", ipop+1);
		    fprintf(pf, "s_skew_w_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_w_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_w\t");
		fprintf(pf, "s_var_w\t");
		fprintf(pf, "s_skew_w\t");
		fprintf(pf, "s_kurt_w\t");
		break;

	      case 3: // tajima's D
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.pritajimasD){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_tajd_%d\t", ipop+1);
		    fprintf(pf, "s_var_tajd_%d\t", ipop+1);
		    fprintf(pf, "s_skew_tajd_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_tajd_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_tajd\t");
		fprintf(pf, "s_var_tajd\t");
		fprintf(pf, "s_skew_tajd\t");
		fprintf(pf, "s_kurt_tajd\t");
		break;

	      case 4:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.priZnS){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_ZnS_%d\t", ipop+1);
		    fprintf(pf, "s_var_ZnS_%d\t", ipop+1);
		    fprintf(pf, "s_skew_ZnS_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_ZnS_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_ZnS\t");
		fprintf(pf, "s_var_ZnS\t");
		fprintf(pf, "s_skew_ZnS\t");
		fprintf(pf, "s_kurt_ZnS\t");
		break;

	      case 5:
		if(pars.cp.npop > 1){
		  fprintf(pf, "s_mean_Fst\t");
		  fprintf(pf, "s_var_Fst\t");
		  fprintf(pf, "s_skew_Fst\t");
		  fprintf(pf, "s_kurt_Fst\t");
		}
		break;

	      case 6:
		for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		  //if((pars.cp.config)[npopi] < MINSEQ) continue;
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		    //if((pars.cp.config)[npopj] < MINSEQ) continue;
		    fprintf(pf, "s_mean_shared_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_var_shared_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_skew_shared_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_kurt_shared_%d_%d\t", npopi+1, npopj+1);
		  }
		}
		break;
	      case 7:
		for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		  //if((pars.cp.config)[npopi] < MINSEQ) continue;
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		    //if((pars.cp.config)[npopj] < MINSEQ) continue;
		    fprintf(pf, "s_mean_private_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_var_private_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_skew_private_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_kurt_private_%d_%d\t", npopi+1, npopj+1);
		  }
		}
		break;

	      case 8:
		for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		  //if((pars.cp.config)[npopi] < MINSEQ) continue;
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		    //if((pars.cp.config)[npopj] < MINSEQ) continue;
		    fprintf(pf, "s_mean_fixed_dif_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_var_fixed_dif_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_skew_fixed_dif_%d_%d\t", npopi+1, npopj+1);
		    fprintf(pf, "s_kurt_fixed_dif_%d_%d\t", npopi+1, npopj+1);
		  }
		}
		break;

	      case 9:
		for(ipop=0; ipop < fst_pops-1; ipop++){
		  //if((pars.cp.config)[ipop] < MINSEQ) continue;
		  for(jpop = ipop+1; jpop < fst_pops; jpop++){
		    //if((pars.cp.config)[jpop] < MINSEQ) continue;
		    pop1=Fst_pops[ ipop ];
		    pop2=Fst_pops[ jpop ];
		    fprintf(pf, "s_mean_pairwise_fst_%d_%d\t", ipop+1, jpop+1);
		    fprintf(pf, "s_var_pairwise_fst_%d_%d\t", ipop+1, jpop+1);
		    fprintf(pf, "s_skew_pairwise_fst_%d_%d\t", ipop+1, jpop+1);
		    fprintf(pf, "s_kurt_pairwise_fst_%d_%d\t", ipop+1, jpop+1);
		  }
		}
		break;

	      case 10:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prifwh){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_fwh_%d\t", ipop+1);
		    fprintf(pf, "s_var_fwh_%d\t", ipop+1);
		    fprintf(pf, "s_skew_fwh_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_fwh_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_FayWuH\t");
		fprintf(pf, "s_var_FayWuH\t");
		fprintf(pf, "s_skew_FayWuH\t");
		fprintf(pf, "s_kurt_FayWuH\t");
		break;

	      case 11: // dvk, dvh
		if(sstats.pristats == 1 && pars.cp.npop > 1 && (sstats.pridvstat)){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_dvk_%d\t", ipop+1);
		    fprintf(pf, "s_var_dvk_%d\t", ipop+1);
		    fprintf(pf, "s_skew_dvk_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_dvk_%d\t", ipop+1);
		    fprintf(pf, "s_mean_dvh_%d\t", ipop+1);
		    fprintf(pf, "s_var_dvh_%d\t", ipop+1);
		    fprintf(pf, "s_skew_dvh_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_dvh_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_dvk\t");
		fprintf(pf, "s_var_dvk\t");
		fprintf(pf, "s_skew_dvk\t");
		fprintf(pf, "s_kurt_dvk\t");
		fprintf(pf, "s_mean_dvh\t");
		fprintf(pf, "s_var_dvh\t");
		fprintf(pf, "s_skew_dvh\t");
		fprintf(pf, "s_kurt_dvh\t");

		break;

	      case 12: // thomson_est
		if(sstats.pristats == 1 && pars.cp.npop > 1){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_thomson_est_%d\t", ipop+1);
		    fprintf(pf, "s_var_thomson_est_%d\t", ipop+1);
		    fprintf(pf, "s_skew_thomson_est_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_thomson_est_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_thomson_est\t");
		fprintf(pf, "s_var_thomson_est\t");
		fprintf(pf, "s_skew_thomson_est\t");
		fprintf(pf, "s_kurt_thomson_est\t");
		break;


	      case 13: // thomson_var
		  if(sstats.pristats == 1 && pars.cp.npop > 1 ){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_mean_thomson_var_%d\t", ipop+1);
		    fprintf(pf, "s_var_thomson_var_%d\t", ipop+1);
		    fprintf(pf, "s_skew_thomson_var_%d\t", ipop+1);
		    fprintf(pf, "s_kurt_thomson_var_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_mean_thomson_var\t");
		fprintf(pf, "s_var_thomson_var\t");
		fprintf(pf, "s_skew_thomson_var\t");
		fprintf(pf, "s_kurt_thomson_var\t");
		break;
		
		
		
	      default:
		fprintf(stderr, "Statistic is not implemented yet,  case: %d!\n", i);
		EXIT_MSABC(1);
	      }
	    }
	  }
	  
	  fprintf(pf, "\n");
	}/* if count == 1 */
	


	for(i=3; i<argc; i++){
	  
	  if(argv[i][0] != '-') continue;
	  switch(argv[i][1]){
	    /* pavlos options */
	  case '-':
	    if(!strcmp(argv[i], "--N")){
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	    }
	    break;
	  case 't':
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 'n':
	    if(argv[i][2] != '\0') break;
	    i++;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 'g':
	    if(argv[i][2] != '\0') break;
	    i++;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 'r':
	    //fprintf(stderr, "%s\n", argv[i]);
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 's':
	    if(printall || param_distr[i] != '-'){ 
	      if( argv[i][2] != 'e' ){
		fprintf(pf, "%d\t",  pars.mp.segsitesin);
	      }
	      else
		i+=nseeds;
	    }
	    break;
	  case 'I':
	    c=0;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%d\t", pars.cp.npop);
	    }
	    i += pars.cp.npop+1 ; // skip the number of populations (done already) and the configuration
	    if( (printall && (i+1 < argc) && argv[i+1][0] != '-') ||
		( (i+1 < argc) && (c = dist_argcheck(i+1,argc,argv)) ) )
	      {
		i++;
		
		fprintf(pf, "%e\t", pars.cp.mig_mat[0][1] * (pars.cp.npop - 1));
	      }
	    c = 0;
	    break;
	  case 'm':
	    c=0;
	    
	    if( argv[i][2] == 'a' ){
	      i++;
	      for(pop = 0; pop < pars.cp.npop; pop++){
		for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		  c = 0;
		  
		  if(  (c=  dist_argcheck(i, argc, argv) ) || printall){
		    fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		    i+=c;
		    if(c) continue;
		  }
		  i++;
		}
	      }
	      i--; // go one back because the index will be increased in the initial loop
	    }
	    else{
	      c=0;
	      i++;
	      pop = atoi(argv[i++])-1;
	      pop2 = atoi(argv[i++])-1;
	      if((c = dist_argcheck(i, argc, argv) ) || printall){
		fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		
		if(c) i+=c-1;
	      }
	    }
	    break;
		
	  case 'G':
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", pars.cp.alphag[0]);
	    }
	    break;
	  case 'e':
	    
	    secarg = argv[i][2];
	    
	    switch(secarg){
	    case 'M': 
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      i += distribution_arguments( param_distr[i] ) + 1;
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break;
	    case 'N':
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      i += distribution_arguments( param_distr[i] ) + 1;
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break;
	    case 'n':
	      if(argv[i][3] != '\0') break;
	      z = i; 
	      z += distribution_arguments( param_distr[i] ) + 2;
	      zpop = atoi(argv[z]);

	      //fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      
	      i += distribution_arguments( param_distr[i] ) + 1;
	      i++;
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break;
	    case 'g':
	      if(argv[i][3] != '\0') break;
	      z = i; 
	      z += distribution_arguments( param_distr[i] ) + 2;
	      zpop = atoi(argv[z]);
	      
	      //fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      
	      i += distribution_arguments( param_distr[i] ) + 1;
	      i++;
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break;

	    case 'G':
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      i += distribution_arguments( param_distr[i] ) + 1;
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break;
	    case 'j':
	      c = distribution_arguments( param_distr[i] );
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      c=0;
	      break;
	      
	    case 's':
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break; 
	    case 'm':
	      
	      c = 0;
	      if( argv[i][3] == 'a'){ // -ema
		
		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "%e\t", global_values[i]);
		}
		
		i += distribution_arguments( param_distr[i] ) + 1;
		i+=2; /* the number of populations */
		
		for(pop = 0; pop < pars.cp.npop; pop++){
		  for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		    
		    c = 0;
		  
		    if((c=  dist_argcheck(i, argc, argv) ) || printall){
		      fprintf(pf, "%e\t", global_values[i-1]);
		      fprintf(stderr, "%d, %d\n", c, i);
		      //fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		      i+=c;
		      if(c) continue;
		    }
		    i++;
		  }
		}
		i--; // go one back because the index will be increased in the initial loop
	      }
	      else{
		z = i; 
		z += distribution_arguments( param_distr[i] ) + 2;
		zpop1 = atoi(argv[z]);
		zpop2 = atoi(argv[z+1]);
		//printf("here\n"); 

		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "%e\t", global_values[i]);
		}
		//fprintf(stderr, "*argv%s\n", argv[i]);
		i += distribution_arguments( param_distr[i] ) + 1;
		//fprintf(stderr, "*argv%s\n", argv[i]);

		c=0;
		
		pop = atoi(argv[++i])-1;
		pop2 = atoi(argv[++i])-1;
		//printf("*pop:%d, pop2:%d\n", pop, pop2);
		
		if((c = dist_argcheck(i, argc, argv) ) || printall){
		  fprintf(pf, "%e\t", global_values[i]);
		  /* fprintf(stderr, "\n\n%d\t%e\t%e\n\n", i, global_values[i], pars.cp.mig_mat[pop][pop2]); 
		     EXIT_MSABC(1);*/
		  if(c) i+=c-1;
		}
	      }
	      break;
	      
	    }
	  }
	}
	
	
	for(i=0; i<MAXSUMSTAT; i++){
	  if(sstats_array[i] == 1){
	    switch(i){
	    case 0: // segsites
	      
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prisegs){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		    priave_segs[ipop] /= sstats_popden_weights[i][ipop];
		    { double sk = skewstat(priskew_segs[ipop], privar_segs[ipop], priave_segs[ipop], sstats_popden_weights[i][ipop]);
		      double ku = kurtstat(prikurt_segs[ipop], priskew_segs[ipop], privar_segs[ipop], priave_segs[ipop], sstats_popden_weights[i][ipop]);
		      privar_segs[ipop] = varstat(privar_segs[ipop], priave_segs[ipop], sstats_popden_weights[i][ipop]);
		      fprintf(pf, "%e\t", priave_segs[ipop]);
		      fprintf(pf, "%e\t", privar_segs[ipop]);
		      fprintf(pf, "%e\t", sk);
		      fprintf(pf, "%e\t", ku);
		    }
		    //}
		}
	      }
	      ave_segs /= sstats_denominator_weights[i];

	      /* calculate the sample variance */
	      { double sk_segs = skewstat(skew_segs, var_segs, ave_segs, sstats_denominator_weights[i]);
	        double ku_segs = kurtstat(kurt_segs, skew_segs, var_segs, ave_segs, sstats_denominator_weights[i]);
	        var_segs = (var_segs/sstats_denominator_weights[i]
			    - ave_segs*ave_segs)/(sstats_denominator_weights[i] - 1)*sstats_denominator_weights[i] ;
	        fprintf(pf, "%e\t", ave_segs);
	        fprintf(pf, "%e\t", var_segs);
	        fprintf(pf, "%e\t", sk_segs);
	        fprintf(pf, "%e\t", ku_segs);
	      }
	      break;
	    case 1: // theta pi
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetapi){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		    priave_thetapi[ipop] /= sstats_popden_weights[i][ipop];
		    { double sk = skewstat(priskew_thetapi[ipop], privar_thetapi[ipop], priave_thetapi[ipop], sstats_popden_weights[i][ipop]);
		      double ku = kurtstat(prikurt_thetapi[ipop], priskew_thetapi[ipop], privar_thetapi[ipop], priave_thetapi[ipop], sstats_popden_weights[i][ipop]);
		      privar_thetapi[ipop] = varstat(privar_thetapi[ipop], priave_thetapi[ipop], sstats_popden_weights[i][ipop]);
		      fprintf(pf, "%e\t", priave_thetapi[ipop]);
		      fprintf(pf, "%e\t", privar_thetapi[ipop]);
		      fprintf(pf, "%e\t", sk);
		      fprintf(pf, "%e\t", ku);
		    }
		    //}
		}
	      }

	      ave_thetapi /= sstats_denominator_weights[i];
	      { double sk_thetapi = skewstat(skew_thetapi, var_thetapi, ave_thetapi, sstats_denominator_weights[i]);
	        double ku_thetapi = kurtstat(kurt_thetapi, skew_thetapi, var_thetapi, ave_thetapi, sstats_denominator_weights[i]);
	        var_thetapi = (var_thetapi/sstats_denominator_weights[i]
			       - ave_thetapi*ave_thetapi)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_thetapi);
	        fprintf(pf, "%e\t", var_thetapi);
	        fprintf(pf, "%e\t", sk_thetapi);
	        fprintf(pf, "%e\t", ku_thetapi);
	      }
	      break;

	    case 2: // theta w
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetaw){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		    priave_thetaw[ipop] /= sstats_popden_weights[i][ipop];
		    { double sk = skewstat(priskew_thetaw[ipop], privar_thetaw[ipop], priave_thetaw[ipop], sstats_popden_weights[i][ipop]);
		      double ku = kurtstat(prikurt_thetaw[ipop], priskew_thetaw[ipop], privar_thetaw[ipop], priave_thetaw[ipop], sstats_popden_weights[i][ipop]);
		      privar_thetaw[ipop] = varstat(privar_thetaw[ipop], priave_thetaw[ipop], sstats_popden_weights[i][ipop]);
		      fprintf(pf, "%e\t", priave_thetaw[ipop]);
		      fprintf(pf, "%e\t", privar_thetaw[ipop]);
		      fprintf(pf, "%e\t", sk);
		      fprintf(pf, "%e\t", ku);
		    }
		    //}
		}
	      }

	      ave_thetaw /= sstats_denominator_weights[i];
	      { double sk_thetaw = skewstat(skew_thetaw, var_thetaw, ave_thetaw, sstats_denominator_weights[i]);
	        double ku_thetaw = kurtstat(kurt_thetaw, skew_thetaw, var_thetaw, ave_thetaw, sstats_denominator_weights[i]);
	        var_thetaw = (var_thetaw/sstats_denominator_weights[i]
			      - ave_thetaw*ave_thetaw)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_thetaw);
	        fprintf(pf, "%e\t", var_thetaw);
	        fprintf(pf, "%e\t", sk_thetaw);
	        fprintf(pf, "%e\t", ku_thetaw);
	      }
	      break;
	      
	    case 3: // tajima's D
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.pritajimasD){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		    priave_tajimasD[ipop] /= sstats_popden_weights[i][ipop];
		    { double sk = skewstat(priskew_tajimasD[ipop], privar_tajimasD[ipop], priave_tajimasD[ipop], sstats_popden_weights[i][ipop]);
		      double ku = kurtstat(prikurt_tajimasD[ipop], priskew_tajimasD[ipop], privar_tajimasD[ipop], priave_tajimasD[ipop], sstats_popden_weights[i][ipop]);
		      privar_tajimasD[ipop] = varstat(privar_tajimasD[ipop], priave_tajimasD[ipop], sstats_popden_weights[i][ipop]);
		      fprintf(pf, "%e\t", priave_tajimasD[ipop]);
		      fprintf(pf, "%e\t", privar_tajimasD[ipop]);
		      fprintf(pf, "%e\t", sk);
		      fprintf(pf, "%e\t", ku);
		    }
		    //}
		}
	      }

	      ave_tajimasD /= sstats_denominator_weights[i];
	      { double sk_tajD = skewstat(skew_tajimasD, var_tajimasD, ave_tajimasD, sstats_denominator_weights[i]);
	        double ku_tajD = kurtstat(kurt_tajimasD, skew_tajimasD, var_tajimasD, ave_tajimasD, sstats_denominator_weights[i]);
	        var_tajimasD = (var_tajimasD/sstats_denominator_weights[i] - ave_tajimasD*ave_tajimasD)/(sstats_denominator_weights[i] - 1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_tajimasD);
	        fprintf(pf, "%e\t", var_tajimasD);
	        fprintf(pf, "%e\t", sk_tajD);
	        fprintf(pf, "%e\t", ku_tajD);
	      }
	      break;
	      
	    case 4: 
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.priZnS){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		    priave_ZnS[ipop] /= sstats_popden_weights[i][ipop];
		    { double sk = skewstat(priskew_ZnS[ipop], privar_ZnS[ipop], priave_ZnS[ipop], sstats_popden_weights[i][ipop]);
		      double ku = kurtstat(prikurt_ZnS[ipop], priskew_ZnS[ipop], privar_ZnS[ipop], priave_ZnS[ipop], sstats_popden_weights[i][ipop]);
		      privar_ZnS[ipop] = varstat(privar_ZnS[ipop], priave_ZnS[ipop], sstats_popden_weights[i][ipop]);
		      fprintf(pf, "%e\t", priave_ZnS[ipop]);
		      fprintf(pf, "%e\t", privar_ZnS[ipop]);
		      fprintf(pf, "%e\t", sk);
		      fprintf(pf, "%e\t", ku);
		    }
		    //}
		}
	      }

	      ave_ZnS /= sstats_denominator_weights[i];
	      { double sk_ZnS = skewstat(skew_ZnS, var_ZnS, ave_ZnS, sstats_denominator_weights[i]);
	        double ku_ZnS = kurtstat(kurt_ZnS, skew_ZnS, var_ZnS, ave_ZnS, sstats_denominator_weights[i]);
	        var_ZnS = (var_ZnS/sstats_denominator_weights[i] - ave_ZnS*ave_ZnS)/(sstats_denominator_weights[i] - 1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_ZnS);
	        fprintf(pf, "%e\t", var_ZnS);
	        fprintf(pf, "%e\t", sk_ZnS);
	        fprintf(pf, "%e\t", ku_ZnS);
	      }
	      break;
	    
	    case 5:
	      if(pars.cp.npop > 1){

		ave_fst /= sstats_denominator_weights[i];
		{ double sk_fst = skewstat(skew_fst, var_fst, ave_fst, sstats_denominator_weights[i]);
		  double ku_fst = kurtstat(kurt_fst, skew_fst, var_fst, ave_fst, sstats_denominator_weights[i]);
		  var_fst  = (var_fst/sstats_denominator_weights[i] - ave_fst*ave_fst) /
		    (sstats_denominator_weights[i] - 1.)* sstats_denominator_weights[i];
		  fprintf(pf, "%e\t", ave_fst);
		  fprintf(pf, "%e\t", var_fst);
		  fprintf(pf, "%e\t", sk_fst);
		  fprintf(pf, "%e\t", ku_fst);
		}
	      }
	      break;

	    case 6:
	      

	      for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		//if((pars.cp.config)[npopi] < MINSEQ) continue;
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  //if( (pars.cp.config)[npopj] < MINSEQ) continue;


		  /* divide by the number of fragments where the statistic is defined and is not nan */

		  if(sstats_popdenpair_weights[i][npopi][npopj] < propfragments * total_stats_weights[i]){
		    ave_shared[npopi][npopj] = 0.;
		    sstats_popdenpair_weights[i][npopi][npopj] = 0.;
		  }

		  ave_shared[npopi][npopj] /= sstats_popdenpair_weights[i][npopi][npopj];
		  { double sk_tmp = skewstat(skew_shared[npopi][npopj], var_shared[npopi][npopj], ave_shared[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    double ku_tmp = kurtstat(kurt_shared[npopi][npopj], skew_shared[npopi][npopj], var_shared[npopi][npopj], ave_shared[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    skew_shared[npopi][npopj] = sk_tmp;
		    kurt_shared[npopi][npopj] = ku_tmp;
		  }
		  var_shared[npopi][npopj] = (var_shared[npopi][npopj]/sstats_popdenpair_weights[i][npopi][npopj]
					      - ave_shared[npopi][npopj]*ave_shared[npopi][npopj])
		    /(sstats_popdenpair_weights[i][npopi][npopj] - 1.)*sstats_popdenpair_weights[i][npopi][npopj];

		}
	      }
	      for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  fprintf(pf, "%e\t", ave_shared[npopi][npopj]);
		  fprintf(pf, "%e\t", var_shared[npopi][npopj]);
		  fprintf(pf, "%e\t", skew_shared[npopi][npopj]);
		  fprintf(pf, "%e\t", kurt_shared[npopi][npopj]);
		}
	      }
	      break;  

	    case 7:
	      

	      for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		//if((pars.cp.config)[npopi] < MINSEQ) continue;
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  //if((pars.cp.config)[npopj] < MINSEQ) continue;
		  if(sstats_popdenpair_weights[i][npopi][npopj] < propfragments * total_stats_weights[i]){
		    ave_private[npopi][npopj] = 0.;
		    sstats_popdenpair_weights[i][npopi][npopj] = 0.;
		  }

		  ave_private[npopi][npopj] /= sstats_popdenpair_weights[i][npopi][npopj];
		  { double sk_tmp = skewstat(skew_private[npopi][npopj], var_private[npopi][npopj], ave_private[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    double ku_tmp = kurtstat(kurt_private[npopi][npopj], skew_private[npopi][npopj], var_private[npopi][npopj], ave_private[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    skew_private[npopi][npopj] = sk_tmp;
		    kurt_private[npopi][npopj] = ku_tmp;
		  }
		  var_private[npopi][npopj] = (var_private[npopi][npopj]/sstats_popdenpair_weights[i][npopi][npopj]
					      - ave_private[npopi][npopj]*ave_private[npopi][npopj])
		    /(sstats_popdenpair_weights[i][npopi][npopj] - 1.)*sstats_popdenpair_weights[i][npopi][npopj];

		}
	      }
	      for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  fprintf(pf, "%e\t", ave_private[npopi][npopj]);
		  fprintf(pf, "%e\t", var_private[npopi][npopj]);
		  fprintf(pf, "%e\t", skew_private[npopi][npopj]);
		  fprintf(pf, "%e\t", kurt_private[npopi][npopj]);
		}
	      }
	      break;    

	    case 8:
	      for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		//if((pars.cp.config)[npopi] < MINSEQ) continue;
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  //if((pars.cp.config)[npopj] < MINSEQ) continue;
		  if(sstats_popdenpair_weights[i][npopi][npopj] < propfragments * total_stats_weights[i]){
		    ave_fixed_dif[npopi][npopj] = 0.;
		    sstats_popdenpair_weights[i][npopi][npopj] = 0.;
		  }
		  ave_fixed_dif[npopi][npopj] /=sstats_popdenpair_weights[i][npopi][npopj];
		  { double sk_tmp = skewstat(skew_fixed_dif[npopi][npopj], var_fixed_dif[npopi][npopj], ave_fixed_dif[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    double ku_tmp = kurtstat(kurt_fixed_dif[npopi][npopj], skew_fixed_dif[npopi][npopj], var_fixed_dif[npopi][npopj], ave_fixed_dif[npopi][npopj], sstats_popdenpair_weights[i][npopi][npopj]);
		    skew_fixed_dif[npopi][npopj] = sk_tmp;
		    kurt_fixed_dif[npopi][npopj] = ku_tmp;
		  }
		  var_fixed_dif[npopi][npopj] = (var_fixed_dif[npopi][npopj]/sstats_popdenpair_weights[i][npopi][npopj]
						 - ave_fixed_dif[npopi][npopj]*ave_fixed_dif[npopi][npopj])
		    /(sstats_popdenpair_weights[i][npopi][npopj]-1.)*sstats_popdenpair_weights[i][npopi][npopj];


		}
	      }
	       for(npopi=0; npopi<pars.cp.npop-1; npopi++){
		for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		  fprintf(pf, "%e\t", ave_fixed_dif[npopi][npopj]);
		  fprintf(pf, "%e\t", var_fixed_dif[npopi][npopj]);
		  fprintf(pf, "%e\t", skew_fixed_dif[npopi][npopj]);
		  fprintf(pf, "%e\t", kurt_fixed_dif[npopi][npopj]);
		}
	      }


	      break; 
	      
	    case 9:
	      
	      for(ipop=0; ipop < fst_pops-1; ipop++){
		//if((pars.cp.config)[ipop] < MINSEQ) continue;
		for(jpop = ipop+1; jpop<fst_pops; jpop++){
		  //if((pars.cp.config)[jpop] < MINSEQ) continue;
		  pop1=Fst_pops[ ipop ];
		  pop2=Fst_pops[ jpop ];
		  if(sstats_popdenpair_weights[i][pop1][pop2] < propfragments * total_stats_weights[i]){
		    ave_pairwise_Fst[pop2][pop1] = 0.;
		    sstats_popdenpair_weights[i][pop1][pop2] = 0.;
		  }
		  ave_pairwise_Fst[pop1][pop2] /=sstats_popdenpair_weights[i][pop1][pop2];
		  { double sk_tmp = skewstat(skew_pairwise_Fst[pop1][pop2], var_pairwise_Fst[pop1][pop2], ave_pairwise_Fst[pop1][pop2], sstats_popdenpair_weights[i][pop1][pop2]);
		    double ku_tmp = kurtstat(kurt_pairwise_Fst[pop1][pop2], skew_pairwise_Fst[pop1][pop2], var_pairwise_Fst[pop1][pop2], ave_pairwise_Fst[pop1][pop2], sstats_popdenpair_weights[i][pop1][pop2]);
		    skew_pairwise_Fst[pop1][pop2] = sk_tmp;
		    kurt_pairwise_Fst[pop1][pop2] = ku_tmp;
		  }
		  var_pairwise_Fst[pop1][pop2] = (
						     var_pairwise_Fst[pop1][pop2]/
						     sstats_popdenpair_weights[i][pop1][pop2] -
						     ave_pairwise_Fst[pop1][pop2]*
						     ave_pairwise_Fst[pop1][pop2]
						     )/(sstats_popdenpair_weights[i][pop1][pop2] - 1.)*
		    sstats_popdenpair_weights[i][pop1][pop2];
		}
	      }

	      for(ipop=0; ipop < fst_pops-1; ipop++){
		//if((pars.cp.config)[ipop] < MINSEQ) continue;
		for(jpop = ipop+1; jpop<fst_pops; jpop++){
		  //if((pars.cp.config)[jpop] < MINSEQ) continue;
		  pop1=Fst_pops[ ipop ];
		  pop2=Fst_pops[ jpop ];
		  fprintf(pf, "%e\t", ave_pairwise_Fst[ pop1 ][ pop2 ]);
		  fprintf(pf, "%e\t", var_pairwise_Fst[ pop1 ][ pop2 ]);
		  fprintf(pf, "%e\t", skew_pairwise_Fst[ pop1 ][ pop2 ]);
		  fprintf(pf, "%e\t", kurt_pairwise_Fst[ pop1 ][ pop2 ]);
		}
	      }
	      break;
	    case 10:
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prifwh){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] < MINSEQ) continue;
		  priave_fwh[ipop] /= sstats_popden_weights[i][ipop];
		  //fprintf(pf, "pop:%d,den:%e,varfwh:%e\t", ipop, sstats_popden_weights[i][ipop], privar_fwh[ipop]);
		  { double sk = skewstat(priskew_fwh[ipop], privar_fwh[ipop], priave_fwh[ipop], sstats_popden_weights[i][ipop]);
		    double ku = kurtstat(prikurt_fwh[ipop], priskew_fwh[ipop], privar_fwh[ipop], priave_fwh[ipop], sstats_popden_weights[i][ipop]);
		    privar_fwh[ipop] = varstat(privar_fwh[ipop], priave_fwh[ipop], sstats_popden_weights[i][ipop]);
		    fprintf(pf, "%e\t", priave_fwh[ipop]);
		    fprintf(pf, "%e\t", privar_fwh[ipop]);
		    fprintf(pf, "%e\t", sk);
		    fprintf(pf, "%e\t", ku);
		  }
		}
	      }


	      ave_fwh /= sstats_denominator_weights[i];
	      { double sk_fwh = skewstat(skew_fwh, var_fwh, ave_fwh, sstats_denominator_weights[i]);
	        double ku_fwh = kurtstat(kurt_fwh, skew_fwh, var_fwh, ave_fwh, sstats_denominator_weights[i]);
	        var_fwh = varstat(var_fwh, ave_fwh, sstats_denominator_weights[i]);
	        fprintf(pf, "%e\t", ave_fwh);
	        fprintf(pf, "%e\t", var_fwh);
	        fprintf(pf, "%e\t", sk_fwh);
	        fprintf(pf, "%e\t", ku_fwh);
	      }
	      break;
	      
	    case 11: // dvk, dvh
	      if(sstats.pristats == 1 && pars.cp.npop > 1 && (sstats.pridvstat)){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] < MINSEQ) continue;
		  priave_dvk[ipop] /= sstats_popden_weights[i][ipop];
		  { double sk_k = skewstat(priskew_dvk[ipop], privar_dvk[ipop], priave_dvk[ipop], sstats_popden_weights[i][ipop]);
		    double ku_k = kurtstat(prikurt_dvk[ipop], priskew_dvk[ipop], privar_dvk[ipop], priave_dvk[ipop], sstats_popden_weights[i][ipop]);
		    privar_dvk[ipop] = varstat(privar_dvk[ipop], priave_dvk[ipop], sstats_popden_weights[i][ipop]);
		    priave_dvh[ipop] /= sstats_popden_weights[i][ipop];
		    double sk_h = skewstat(priskew_dvh[ipop], privar_dvh[ipop], priave_dvh[ipop], sstats_popden_weights[i][ipop]);
		    double ku_h = kurtstat(prikurt_dvh[ipop], priskew_dvh[ipop], privar_dvh[ipop], priave_dvh[ipop], sstats_popden_weights[i][ipop]);
		    privar_dvh[ipop] = varstat(privar_dvh[ipop], priave_dvh[ipop], sstats_popden_weights[i][ipop]);

		    fprintf(pf, "%e\t", priave_dvk[ipop]);
		    fprintf(pf, "%e\t", privar_dvk[ipop]);
		    fprintf(pf, "%e\t", sk_k);
		    fprintf(pf, "%e\t", ku_k);

		    fprintf(pf, "%e\t", priave_dvh[ipop]);
		    fprintf(pf, "%e\t", privar_dvh[ipop]);
		    fprintf(pf, "%e\t", sk_h);
		    fprintf(pf, "%e\t", ku_h);
		  }
		}
	      }


	      ave_dvk /= sstats_denominator_weights[i];
	      { double sk_dvk = skewstat(skew_dvk, var_dvk, ave_dvk, sstats_denominator_weights[i]);
	        double ku_dvk = kurtstat(kurt_dvk, skew_dvk, var_dvk, ave_dvk, sstats_denominator_weights[i]);
	        var_dvk = (var_dvk/sstats_denominator_weights[i]
			   - ave_dvk*ave_dvk)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        ave_dvh /= sstats_denominator_weights[i];
	        double sk_dvh = skewstat(skew_dvh, var_dvh, ave_dvh, sstats_denominator_weights[i]);
	        double ku_dvh = kurtstat(kurt_dvh, skew_dvh, var_dvh, ave_dvh, sstats_denominator_weights[i]);
	        var_dvh = (var_dvh/sstats_denominator_weights[i]
			      - ave_dvh*ave_dvh)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_dvk);
	        fprintf(pf, "%e\t", var_dvk);
	        fprintf(pf, "%e\t", sk_dvk);
	        fprintf(pf, "%e\t", ku_dvk);
	        fprintf(pf, "%e\t", ave_dvh);
	        fprintf(pf, "%e\t", var_dvh);
	        fprintf(pf, "%e\t", sk_dvh);
	        fprintf(pf, "%e\t", ku_dvh);
	      }
	      break;
	      
	    case 12: // thomson_est
	      if(sstats.pristats == 1 && pars.cp.npop > 1){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		  priave_thomson_est[ipop] /= sstats_popden_weights[i][ipop];
		  { double sk = skewstat(priskew_thomson_est[ipop], privar_thomson_est[ipop], priave_thomson_est[ipop], sstats_popden_weights[i][ipop]);
		    double ku = kurtstat(prikurt_thomson_est[ipop], priskew_thomson_est[ipop], privar_thomson_est[ipop], priave_thomson_est[ipop], sstats_popden_weights[i][ipop]);
		    privar_thomson_est[ipop] = varstat(privar_thomson_est[ipop], priave_thomson_est[ipop], sstats_popden_weights[i][ipop]);
		    fprintf(pf, "%e\t", priave_thomson_est[ipop]);
		    fprintf(pf, "%e\t", privar_thomson_est[ipop]);
		    fprintf(pf, "%e\t", sk);
		    fprintf(pf, "%e\t", ku);
		  }
		  //}
		}
	      }

	      ave_thomson_est /= sstats_denominator_weights[i];
	      { double sk_te = skewstat(skew_thomson_est, var_thomson_est, ave_thomson_est, sstats_denominator_weights[i]);
	        double ku_te = kurtstat(kurt_thomson_est, skew_thomson_est, var_thomson_est, ave_thomson_est, sstats_denominator_weights[i]);
	        var_thomson_est = (var_thomson_est/sstats_denominator_weights[i] - ave_thomson_est*ave_thomson_est)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_thomson_est);
	        fprintf(pf, "%e\t", var_thomson_est);
	        fprintf(pf, "%e\t", sk_te);
	        fprintf(pf, "%e\t", ku_te);
	      }
	      break;
	      
	    case 13: // thomson_var
	      if(sstats.pristats == 1 && pars.cp.npop > 1){
		for(ipop=0; ipop<pars.cp.npop; ipop++){
		  //if((pars.cp.config)[ipop] >= MINSEQ){
		  priave_thomson_var[ipop] /= sstats_popden_weights[i][ipop];
		  { double sk = skewstat(priskew_thomson_var[ipop], privar_thomson_var[ipop], priave_thomson_var[ipop], sstats_popden_weights[i][ipop]);
		    double ku = kurtstat(prikurt_thomson_var[ipop], priskew_thomson_var[ipop], privar_thomson_var[ipop], priave_thomson_var[ipop], sstats_popden_weights[i][ipop]);
		    privar_thomson_var[ipop] = varstat(privar_thomson_var[ipop], priave_thomson_var[ipop], sstats_popden_weights[i][ipop]);
		    fprintf(pf, "%e\t", priave_thomson_var[ipop]);
		    fprintf(pf, "%e\t", privar_thomson_var[ipop]);
		    fprintf(pf, "%e\t", sk);
		    fprintf(pf, "%e\t", ku);
		  }
		  //}
		}
	      }

	      ave_thomson_var /= sstats_denominator_weights[i];
	      { double sk_tv = skewstat(skew_thomson_var, var_thomson_var, ave_thomson_var, sstats_denominator_weights[i]);
	        double ku_tv = kurtstat(kurt_thomson_var, skew_thomson_var, var_thomson_var, ave_thomson_var, sstats_denominator_weights[i]);
	        var_thomson_var = (var_thomson_var/sstats_denominator_weights[i] - ave_thomson_var*ave_thomson_var)/(sstats_denominator_weights[i] -1)*sstats_denominator_weights[i];
	        fprintf(pf, "%e\t", ave_thomson_var);
	        fprintf(pf, "%e\t", var_thomson_var);
	        fprintf(pf, "%e\t", sk_tv);
	        fprintf(pf, "%e\t", ku_tv);
	      }
	      break;
	    
  
	      
	    default:
	      fprintf(stderr, "Statistic is not implemented yet,  case: %d!\n", i);
	      break;
	    }
	  }
	  
	}
	fprintf(pf, "\n");
/* 	if(count == 1){ */
/* 	  fprintf(pf, "\n"); */
	  
/* 	} */
      }

      // free the space for the fragments
      
      if(count == howmany){
	
	for(i=0; i<nofragments; i++)
	  free(fname[i]);
	free(fname);
	
	free(fsize);
	free(frecbp);
	free(fmubp);
	free(flength);
	free(fpop);
      }
      
    }

            
    if(fragmode == 1) continue;

   // added by Jeff Ross-Ibara August 2007
    if( pars.cp.nsam>largestnsam ){
      if( pars.mp.segsitesin ==	0 ){
	list = rematrix(pars.cp.nsam, largestnsam, maxsites+1, list ); 
      }
      else{ 
	list = rematrix( pars.cp.nsam, largestnsam, pars.mp.segsitesin+1, list );
      }
      largestnsam=pars.cp.nsam; 
    }
    
    
    if(obin != NULL){
      /* !!!!!!!!!!!!!!!!!!!! READ OBSERVED DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
      segsites = getobservations(obin, &list);
      modifyList(outgroup, phasemode, list, segsites, pars.cp.nsam);
      
      if(segsites < 0){
	fprintf(stderr, "segsites should be positive. it is %d\n", segsites);
	EXIT_MSABC(1);
      }
    }
    else{
      /* !!!!!!!!!!!!!!!!!!!!! SIMULATE DATA !!!!!!!!!!!!!!!!!!!!!!!!!!*/
      segsites = gensam( list, &probss, &tmrca, &ttot, largestnsam, old_maxnsam) ;
      modifyList(outgroup, phasemode, list, segsites, pars.cp.nsam);
    }

#ifdef R_PACKAGE_BUILD
    /* SFS accumulation mode: compute SFS from haplotypes and skip text output */
    if (sfs_mode_active) {
      sfs_accumulate(list, pars.cp.nsam, segsites, pars.cp.config, pars.cp.npop);
      continue;  /* skip summary stats and text output for this replicate */
    }
#endif

    if( (missingin !=NULL) && (missingdata[0] > 0) ){
      /* we have only one fragment i.e. the last parameter should be 0 */
      if( ! put_missing(list, posit, largestnsam, segsites, misseq,  misstart, misend, missingdata[0], 0)){ 
	fprintf(stderr, "Couldn't embed missing data\n");
      }
    }


    if(old_maxnsam < pars.cp.nsam)
      old_maxnsam = pars.cp.nsam;
    if( (obin == NULL) && tempdatafile!=NULL ){  
      
      fprintf(tempdatafile, "\n//");
      fprintf(tempdatafile, "\t%d\t%e\t%e", count, pars.mp.theta, pars.cp.r);
      if( ntbs >0 ){
	for(k=0; k< ntbs; k++) 
	  printf("\t%s", tbsparamstrs[k] ) ;

      }
      
      fprintf(tempdatafile, "\n");
      
      if( pars.mp.timeflag ) 
	//fprintf(pf,"*time:\t%lf\t%lf\n",tmrca, ttot ) ;
	fprintf(tempdatafile,"*time:\t%lf\t%lf\n",tmrca, ttot ) ;
      
      if( (segsites > 0 ) || ( pars.mp.theta > 0.0 ) ) {
	if( (pars.mp.segsitesin > 0 ) && ( pars.mp.theta > 0.0 )) 
	  /* fprintf(pf,"prob: %g\n", probss ) ;
	   */
	  fprintf(tempdatafile,"prob: %g\n", probss ) ;
	//fprintf(pf,"segsites: %d\n",segsites);
	fprintf(tempdatafile,"segsites: %d\n",segsites);
	if( segsites > 0 )	
	  //fprintf(pf,"positions: ");
	  fprintf(tempdatafile,"positions: ");
	for( i=0; i<segsites; i++)
	  /* fprintf(pf,"%6.4lf ",posit[i] ); */
	  fprintf(tempdatafile,"%e ",posit[i] );
	/* fprintf(pf,"\n"); */
	fprintf(tempdatafile,"\n");
	
	  for(i=0;i<pars.cp.nsam; i++) { /* fprintf(pf,"%s\n", list[i] ); */ 
	    fprintf(tempdatafile,"%s\n", list[i] );  
	  }
      }
    }
    
    
    if(sumstat){
      if(count == 1){
 	//headerfile = fopen("header.txt", "w");
	denominators(pars.cp.nsam, &hn, &sqhn, &bn);
      }

      /* print the parameters and the summary statistics */
      if(count == 1){
	for(i=3; i<argc; i++){
	  
	  //fprintf(stderr, "argv[9]: %s\n", argv[9]);
	  //fprintf(stderr, "*****argv[8]: %s\t argv: %s\n", argv[8], argv[i]);
	  
	  if(argv[i][0] != '-') continue;
	  switch(argv[i][1]){
	    /* case '-': */
	    /* 	      if(!strcmp(argv[i], "--N")) */
/* 		if(printall || param_distr[i] != '-') */
/* 		  fprintf(pf, "p_Ne\t"); */
/* 	      break; */
	    case 't':
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_theta\t");
	      break;
	    case 'n':
	      if(argv[i][2] != '\0') break;
	      i++;
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_init_size_pop_%d\t", atoi(argv[i]) );
	      break;
	    case 'g':
	      if(argv[i][2] != '\0') break;
	      i++;
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_g_pop_%d\t", atoi(argv[i]) );
	      break;
	    case 'r':
	      if(printall || param_distr[i] != '-')
		fprintf(pf, "p_rho\t");
	      break;
	    case 's':
	      if(printall || param_distr[i] != '-'){ 
		if( strcmp(argv[i], "-seeds") ){
		  fprintf(pf, "p_insegs\t");
		}
		else 
		  i+=nseeds;
	      }
	      break;
	    case 'I':
	      c=0;
	      if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_npop\t");
	      }
	      i += pars.cp.npop+1 ; // skip the number of populations (done already) and the configuration
	      if( (printall && (i+1 < argc) && argv[i+1][0] != '-') ||
		  ( (i+1 < argc) && (c = dist_argcheck(i+1,argc,argv)) ) )
		{
		  i++;
		  fprintf(pf, "p_isl_totmig\t");
		}
	      c = 0;
	      break;
	    case 'm':
	      c=0;
	      
	      if( argv[i][2] == 'a' ){
		i++;
		for(pop = 0; pop < pars.cp.npop; pop++){
		  for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		    c = 0;
		    
		    if( (c=  dist_argcheck(i, argc, argv) ) || printall){
			fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
			i+=c;
			if(c) continue;
		    }
		    i++; 
		  }
		}
		i--; // go one back because the index will be increased in the initial loop
	      }
	      else{
		c=0;
		i++;
		pop = atoi(argv[i++])-1;
		pop2 = atoi(argv[i++])-1;
		if((c = dist_argcheck(i, argc, argv) ) || printall){
		    fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
		    if(c) i+=c-1;
		}
	      }
	      break;
	    case 'G':
	      if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_G\t");
	      }
	      break;
	    case 'e':
	      secarg = argv[i][2];
	      switch(secarg){
	      case 'M':
		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "p_globalMutRate_change_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		if(printall || param_distr[i] != '-'){ 
		  fprintf(pf, "p_change_globalMutRate\t");
		}
		break;
	      case 'N':
		if(printall || param_distr[i] != '-'){
		  fprintf(pf, "p_pop_change_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "p_pop_change_newpop\t");
	      }
	      break;
	      case 'n':
		if(argv[i][3] != '\0') break;
		z = i; 
		z += distribution_arguments( param_distr[i] ) + 2;
		zpop = atoi(argv[z]);
		
		//fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_en_time_pop_%d\t", zpop);
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		i++;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_en_size_pop_%d\t", atoi(argv[i]) );
		}
		break;
	      case 'g':
		if(argv[i][3] != '\0') break;
		z = i; 
		z += distribution_arguments( param_distr[i] ) + 2;
		zpop = atoi(argv[z]);
		
		//fprintf(stderr, "z: %d, %s, %d\n", zpop, argv[z], distribution_arguments( param_distr[i] ) + 1);
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eg_time_pop_%d\t", zpop);
		}
		
		i += distribution_arguments( param_distr[i] ) + 1;
		i++;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eg_pop_%d\t", atoi(argv[i]) );
		}
		break;
	      case 'G':
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eG_time\t");
		}
		i += distribution_arguments( param_distr[i] ) + 1;
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_eG_rate\t");
		}
		break;
	      case 'j':
		c = distribution_arguments( param_distr[i] );
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_ej_time_pop_%d_%d\t", atoi(argv[i+2+c]), atoi(argv[i+3+c]) );
		}
		c=0;
		break;
	      case 's':
		if(printall || param_distr[i] != '-'){ 
		    fprintf(pf, "p_es_time\t");
		}
		break; 
	      case 'm':
		
		c = 0;
		if( argv[i][3] == 'a'){ // -ema
		
		  if(printall || param_distr[i] != '-'){
		      fprintf(pf, "p_ema_time\t");
		  }
		  i += distribution_arguments( param_distr[i] ) + 1;
		  i+=2; /* the time */
		  for(pop = 0; pop < pars.cp.npop; pop++){
		    for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		      c = 0;
		      if( (c=  dist_argcheck(i, argc, argv) ) || printall){
			fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
			i+=c;
			if(c) continue;
		      }
		    }
		  }
		  i--; // go one back because the index will be increased in the initial loop
		}
		else{
		  z = i; 
		  z += distribution_arguments( param_distr[i] ) + 2;
		  zpop1 = atoi(argv[z]);
		  zpop2 = atoi(argv[z+1]);
		//printf("here\n"); 
		  
		  if(printall || param_distr[i] != '-'){
		    fprintf(pf, "p_em_time_pops_%d_%d\t", zpop1, zpop2);
		  }
		  i += distribution_arguments( param_distr[i] ) + 1;
		  c=0;
		  pop = atoi(argv[++i])-1;
		  pop2 = atoi(argv[++i])-1;
		  
		  if((c = dist_argcheck(i, argc, argv) ) || printall){
		      fprintf(pf, "p_mig_%d_%d\t", pop+1, pop2+1);
		      if(c) i+=c-1;
		  }
		}
		break;
	      }
	    }
	    
	    //fprintf(stderr, "argv[9]: %s\n", argv[9]);
	    //fprintf(stderr, "*argv[8]: %s\t argv[10]: %s\n", argv[8], argv[10]);
	  }
	  
	  /* print the summary statistics headers */
	  
	  for(i=0; i<MAXSUMSTAT; i++){
	    //fprintf(stderr, "argv[9]: %s\n", argv[9]);
	    //fprintf(stderr, "*argv[8]: %s\n", argv[8]);
	    //printf("sstat %i\n", i);
	    if(sstats_array[i] == 1){
	      switch(i){
	      case 0: // segsites
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prisegs){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_segs_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_segs\t");
		break;

	      case 1:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetapi){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_pi_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_theta_pi\t");
		break;


	      case 2: // theta w
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetaw){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_theta_w_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_theta_w\t");
		break;

	      case 3:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.pritajimasD){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_tajimasD_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_tajimasD\t");
		break;

	      case 4:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.priZnS){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_ZnS_%d\t", ipop+1);
		  }
		}

		fprintf(pf,  "s_ZnS\t");
		break;

	      case 5:
		if(pars.cp.npop > 1)
		  fprintf(pf, "s_fst\t");
		break;

	      case 6:
		if(pars.cp.npop > 1){
		  for(npopi=0; npopi < pars.cp.npop - 1; npopi++){
		    //if((pars.cp.config)[npopi] < MINSEQ) continue;
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      //if((pars.cp.config)[npopj] < MINSEQ) continue;
		      fprintf(pf, "s_perc_shared_%d_%d\t", npopi+1, npopj+1);
		    }
		  }
		}
		break;

	      case 7:
		if(pars.cp.npop > 1){
		  for(npopi=0; npopi < pars.cp.npop - 1; npopi++){
		    //if((pars.cp.config)[npopi] < MINSEQ) continue;
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      //if((pars.cp.config)[npopj] < MINSEQ) continue;
		      fprintf(pf, "s_perc_private_%d_%d\t", npopi+1, npopj+1);
		    }
		  }
		}
		break;

	      case 8:
		if(pars.cp.npop > 1){
		  for(npopi=0; npopi < pars.cp.npop - 1; npopi++){
		    //if((pars.cp.config)[npopi] < MINSEQ) continue;
		    for(npopj=npopi+1; npopj<pars.cp.npop; npopj++){
		      //if((pars.cp.config)[npopj] < MINSEQ) continue;
		      fprintf(pf, "s_perc_fixed_dif_%d_%d\t", npopi+1, npopj+1);
		    }
		  }
		}
		break;

	      case 9:
		if(pars.cp.npop > 1){
		  for(ipop=0; ipop<fst_pops-1; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    for(jpop=ipop+1; jpop<fst_pops; jpop++){
		      //if((pars.cp.config)[jpop] < MINSEQ) continue;
		      pop1=Fst_pops[ ipop ];
		      pop2=Fst_pops[ jpop ];
		      fprintf(pf, "s_pairwise_fst_%d_%d\t", pop1+1, pop2+1);
		    }
		  }
		}
		break;

	      case 10:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prifwh){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_FayWuH_%d\t", ipop+1);
		  }
		}
		fprintf(pf, "s_FayWuH\t");
		break;

	      case 11:
		if(sstats.pristats == 1 && pars.cp.npop > 1 && (sstats.pridvstat)){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_dvk_%d\t", ipop+1);
		    fprintf(pf, "s_dvh_%d\t", ipop+1);

		  }
		}
		fprintf(pf, "s_dvk\t");
		fprintf(pf, "s_dvh\t");
		break;

	      case 12:
		if(sstats.pristats == 1 && pars.cp.npop > 1){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_thomson_est_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_thomson_est\t");
		break;

	      case 13:
		if(sstats.pristats == 1 && pars.cp.npop > 1){
		  for(ipop=0; ipop<pars.cp.npop; ipop++){
		    //if((pars.cp.config)[ipop] < MINSEQ) continue;
		    fprintf(pf, "s_thomson_var_%d\t", ipop+1);
		  }
		}
		fprintf(pf,  "s_thomson_var\t");
		break;
		

	      default:
		fprintf(stderr, "Statistic is not implemented yet!\n");
		EXIT_MSABC(1);
	      }
	    }
	  }
	  fprintf(pf, "\n");
      } /* if count == 1 */
      
       
      
      for(i=3; i<argc; i++){
	/* fprintf(stderr, "argv[8]: %s\n", argv[8]); */
/* 	fprintf(stderr, "-argv[9]: %s\n", argv[9]); */
/* 	fprintf(stderr, "-*argv[8]: %s\t argv: %s\n", argv[8], argv[i]); */
	if(argv[i][0] != '-') continue;
	switch(argv[i][1]){
	  
	case 't':
	  if(printall || param_distr[i] != '-'){ 
	    fprintf(pf, "%e\t", global_values[i]);
	  }
	  
	  
	  break;
	case 'n':
	  
	  if(argv[i][2] != '\0') break;
	  i++;
	  //fprintf(stderr, "**%c\t%c\n", argv[i][2], param_distr[i]);
	  if(printall || param_distr[i] != '-'){
	    fprintf(pf, "%e\t", global_values[i]);
	  }
	  break; 
	case 'g':
	  
	  if(argv[i][2] != '\0') break;
	  i++;
	  if(printall || param_distr[i] != '-'){ 
	    fprintf(pf, "%e\t", global_values[i]);
	  }
	  break; 
	case 'r':
	  //fprintf(stderr, "%s\n", argv[i]);
	  if(printall || param_distr[i] != '-'){
	    fprintf(pf, "%e\t", global_values[i]);
	  }
	  break;
	case 's':
	  if(printall || param_distr[i] != '-'){ 
	    if( argv[i][2] != 'e' ){
	      fprintf(pf, "%d\t",  pars.mp.segsitesin);
	    }
	    else{
	      i+=nseeds;
	    }
	  }
	  break;
	case 'I':
	  c=0;
	  if(printall || param_distr[i] != '-'){ 
	    fprintf(pf, "%d\t", pars.cp.npop);
	  }
	  /* 20091026 */
	  i += pars.cp.npop+1 ; // skip the number of populations (done already) and the configuration /
	  /* 20091026 */
	  if( (printall && (i+1 < argc) && argv[i+1][0] != '-') ||
	      /* 20091026 */
	      ( (i+1 < argc) && (c = dist_argcheck(i+1,argc,argv)) ) )
	    {
	      /* 20091026 */
	      i++;
	      fprintf(pf, "%e\t", pars.cp.mig_mat[0][1] * (pars.cp.npop - 1));
	      }
	  c = 0;
	  break;
	  
	
	case 'G':
	  if(printall || param_distr[i] != '-'){ 
	    fprintf(pf, "%e\t", pars.cp.alphag[0]);
	  }
	  break;
	  
	case 'm':
	  c=0;
	  
	  if( argv[i][2] == 'a' ){
	    i++;
	    for(pop = 0; pop < pars.cp.npop; pop++){
	      for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		c = 0;
		
		if( (c=  dist_argcheck(i, argc, argv) ) || printall){
		  fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		  i+=c;
		  if(c) continue;
		}
		i++;
	      }
	    }
	    i--; // go one back because the index will be increased in the initial loop
	  }
	  else{
	    c=0;
	    i++;
	    pop = atoi(argv[i++])-1;
	    pop2 = atoi(argv[i++])-1;
	    if((c = dist_argcheck(i, argc, argv) ) || printall){
	      fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
	      
	      if(c) i+=c-1;
	    }
	    
	  }
	  break;
	  
	  
	case 'e':
	  
	  secarg = argv[i][2];
	  
	  switch(secarg){
	  case 'N':
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }

	    i += distribution_arguments( param_distr[i] ) + 1;
	    

	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }

	    i += distribution_arguments( param_distr[i] ) + 1;
	    

	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 'n':
	    if(argv[i][3] != '\0') break;
	    z = i; 
	    z += distribution_arguments( param_distr[i] ) + 2;
	    zpop = atoi(argv[z]);
	    if(printall || param_distr[i] != '-'){
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    i += distribution_arguments( param_distr[i] ) + 1;
	    i++;
	    if(printall || param_distr[i] != '-'){
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;

	  case 'g':
	    if(argv[i][3] != '\0') break;
	    z = i; 
	    z += distribution_arguments( param_distr[i] ) + 2;
	    zpop = atoi(argv[z]);
	    if(printall || param_distr[i] != '-'){
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    i += distribution_arguments( param_distr[i] ) + 1;
	    i++;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	    
	  case 'G':
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    i += distribution_arguments( param_distr[i] ) + 1;
	    if(printall || param_distr[i] != '-'){ 
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    break;
	  case 'j':
	    c = distribution_arguments( param_distr[i] );
	    if(printall || param_distr[i] != '-'){;
	      fprintf(pf, "%e\t", global_values[i]);
	    }
	    c=0;
	    break;
	    
	    case 's':
	      if(printall || param_distr[i] != '-'){ 
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      break; 
	      
	  case 'm':
	    
	    c = 0;
	    if( argv[i][3] == 'a'){ // -ema
	      
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      
	      i += distribution_arguments( param_distr[i] ) + 1;
	      i+=2;
	      /* the time */
	      
	      for(pop = 0; pop < pars.cp.npop; pop++){
		for(pop2 = 0; pop2 < pars.cp.npop; pop2++){
		  c = 0;
		  
		  if( (c=  dist_argcheck(i, argc, argv) ) || printall){
		    //fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		    fprintf(pf, "%e\t", global_values[i-1]);
		    i+=c;
		    if(c) continue;
		  }
		  i++;
		}
	      }
	      i--; // go one back because the index will be increased in the initial loop
	    }
	    else{
	      
	      z = i; 
	      z += distribution_arguments( param_distr[i] ) + 2;
	      zpop1 = atoi(argv[z]);
	      zpop2 = atoi(argv[z+1]);
	      //printf("here\n"); 
	      
	      if(printall || param_distr[i] != '-'){
		fprintf(pf, "%e\t", global_values[i]);
	      }
	      //fprintf(stderr, "*argv%s\n", argv[i]);
	      i += distribution_arguments( param_distr[i] ) + 1;
	      //fprintf(stderr, "*argv%s\n", argv[i]);
	      
	      c=0;
	      
	      pop = atoi(argv[++i])-1;
	      pop2 = atoi(argv[++i])-1;
	      //printf("*pop:%d, pop2:%d\n", pop, pop2);
	      
	      if((c = dist_argcheck(i, argc, argv) ) || printall){
		fprintf(pf, "%e\t", global_values[i]);
		//fprintf(pf, "%e\t", pars.cp.mig_mat[pop][pop2]);
		if(c) i+=c-1;
	      }
	      
	    }
	   
	    break;
	  }
	}
	
      }


      
      int insegsites = (segsites > 0) ? segsites:1;
      if(insegsites > largestinsegsites || pars.cp.nsam > largesample){
	
	if(insegsites > largestinsegsites){
	  allelic_map = (int*)realloc(allelic_map, (unsigned)(insegsites*sizeof(int)));
	  missing_map = (int*)realloc(missing_map, (unsigned)(insegsites*sizeof(int)));
	  
	  for( j = 0; sstats.pristats==1 && j<pars.cp.npop && pars.cp.npop > 1; j++){
	    npop_allelic_map[j] = (int*)realloc(npop_allelic_map[j], (unsigned)(insegsites*sizeof(int)));
	    npop_missing_map[j] = (int*)realloc(npop_missing_map[j], (unsigned)(insegsites*sizeof(int)));
	  }
	  
	  largestinsegsites = insegsites;
	}
	
	if( pars.cp.nsam <= largesample){
	  //printf("realloc ders[i]\n");
	  for( j=0; j<largesample; j++) {
	    if(   ! (  ders[j] = (int *) realloc( ders[j], (unsigned) insegsites*sizeof(int) )  )   ) {
	      perror("realloc error in rematrix. ders");
	    }
	  }
	  for(ipop=0; sstats.pristats ==1 && ipop <pars.cp.npop && pars.cp.npop > 1; ipop++){
	    for( j=0; j<largesample; j++) {
	      if(   ! (  npop_ders[ipop][j] = 
			 (int *) realloc( npop_ders[ipop][j], (unsigned) insegsites*sizeof(int) )  )   ) {
		perror("realloc error in rematrix. npop_ders");
	      }
	    }
	  }
	  
	}
	
	if( (pars.cp.nsam > largesample)){
	  //printf("realloc ders\n");
	  ders = rematrix_int(pars.cp.nsam, largesample, largestinsegsites+1, ders );
	  
	  for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++)
	    npop_ders[ipop] = rematrix_int(pars.cp.nsam, largesample, largestinsegsites+1, npop_ders[ipop]);
	  
	  largesample = pars.cp.nsam;
	}
      }  
      
      int segs;
      double thetapi=0., thetaw = 0., tajd = 0., thetah = 0., fwh = 0., zns=0., dvh=0., dvk=0., thomson_est = 0., thomson_var = 0.;
      for(ipop=0; ipop<pars.cp.npop && sumstat && pars.cp.npop > 1 && sstats.pristats == 1 ;  ipop++){
	prisegs[ipop] = prithetapi[ipop] = prithetaw[ipop] = pritajd[ipop] 
	  = prithetah[ipop] = prifwh[ipop] = prizns[ipop] = pridvk[ipop] = pridvh[ipop] = 0.;
      }
      
      
      /* 
	 for the denominators of the statistics
      */
      
      if(pars.cp.nsam != old_nsam){
	//printf("calculation of denominators\n");
	denominators(pars.cp.nsam, &hn, &sqhn, &bn);
      }
      old_nsam = pars.cp.nsam;
      
     


      for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++)
	    denominators( (pars.cp.config)[ipop], &prihn[ipop], &prisqhn[ipop], &pribn[ipop]);
      
      samples_b[0] = 0;
      samples_e[0] = (pars.cp.config)[0];
      
      for(ipop=1; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop > 1; ipop++){
	samples_b[ipop] = samples_b[ipop-1] + (pars.cp.config)[ipop-1];
	samples_e[ipop] = samples_b[ipop] + (pars.cp.config)[ipop];
      }
      

      set_ders(ders, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);
      set_almap(allelic_map, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);
      

      set_missing(missing_map, list, pars.cp.nsam, segsites, 0, pars.cp.nsam);

      

      for(ipop=0; sstats.pristats==1 && ipop<pars.cp.npop && pars.cp.npop> 1; ipop++){
	set_almap(npop_allelic_map[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
	set_missing(npop_missing_map[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
	set_ders(npop_ders[ipop], list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop]);
      }

      int truesegsites = seg_sites(allelic_map, pars.cp.nsam, segsites);
      
      


      for(ipop=0;  sstats.pristats==1 && pars.cp.npop> 1 && ipop < pars.cp.npop; ipop++){
	pritruesegsites[ipop] = seg_sites(npop_allelic_map[ipop], 
					  (pars.cp.config)[ipop],
					  segsites);
      }


      for(i=0; i<MAXSUMSTAT; i++){
	//printf("sstat %i\n", i);
	if(sstats_array[i] == 1){
	  switch(i){
	  case 0: // segsites
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prisegs){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		prisegs[ipop] = pritruesegsites[ipop];
		fprintf(pf, "%e\t", prisegs[ipop]);
	      }
	    }
	    segs = truesegsites;
	    fprintf(pf, "%i\t", segs);
	    break;

	  case 1: // theta pi
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetapi){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		theta_pi(&prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]);
		fprintf(pf, "%e\t", prithetapi[ipop]);
	      }
	    }
	    
	    if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map)){
	      //thetapi /= pars.cp.nsites;
	      fprintf (pf, "%e\t", thetapi);
	      //thetapi *=  pars.cp.nsites;//add by SL 0104

	    }
	    else
	      fprintf(pf, "nan\t");
	    break;

	  case 2: // theta w
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prithetaw){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		theta_w(&prithetaw[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, prihn[ipop], npop_missing_map[ipop]);
		fprintf(pf, "%e\t", prithetaw[ipop]);
	      }
	    }
	    if(theta_w(&thetaw, allelic_map, pars.cp.nsam, segsites, hn, missing_map)){
	      //thetaw /= pars.cp.nsites;
	      fprintf(pf, "%e\t", thetaw);
	      //thetaw *= pars.cp.nsites;//add by SL 0104

	    }
	    else
	      fprintf(pf, "nan\t");
	    break;

	  case 3:
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.pritajimasD){
	      
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		if( (prithetaw[ipop] && prithetapi[ipop] && pritruesegsites[ipop] ) &&
		    tajD(&pritajd[ipop], pritruesegsites[ipop], (pars.cp.config)[ipop], 
			 prithetaw[ipop], prithetapi[ipop], prihn[ipop], prisqhn[ipop]) )
		  {
		    /* printf("\n*pop %d, thetaw: %e\n", ipop, prithetaw[ipop]); */
/* 		    printf("pop %d, thetapi: %e\n", ipop, prithetapi[ipop]); */
/* 		    printf("pop %d, tajima's D: %e\n", ipop, pritajd[ipop]); */
		    
		    fprintf(pf, "%e\t", pritajd[ipop]);
		  }
	       
		else if(pritruesegsites[ipop] && (!prithetaw[ipop] || !prithetapi[ipop])){
		  if( theta_w( &prithetaw[ipop], npop_allelic_map[ipop], 
			       (pars.cp.config)[ipop], segsites, prihn[ipop], npop_missing_map[ipop]) &&
		      theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
		      tajD(&pritajd[ipop], pritruesegsites[ipop], (pars.cp.config)[ipop], 
			   prithetaw[ipop], prithetapi[ipop], prihn[ipop], prisqhn[ipop]) )
		    {
		      /* printf("pop %d, thetaw: %e", ipop, prithetaw[ipop]); */
/* 		      printf("pop %d, thetapi: %e", ipop, prithetapi[ipop]); */
/* 		      printf("pop %d, tajima's D: %e", ipop, pritajd[ipop]); */
		      fprintf(pf, "%e\t", pritajd[ipop]);
		    }
		}
		else if(pritruesegsites[ipop] == 0) // Tajimas' D isn't defined
		  fprintf(pf, "nan\t");
	      }
	    }
	    if( (thetaw && thetapi && truesegsites) && 
		tajD(&tajd, truesegsites, pars.cp.nsam, thetaw, thetapi, hn, sqhn)
		)
	      fprintf(pf, "%e\t", tajd);
	    else if(truesegsites && (!thetaw || !thetapi)){
	      if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) && 
		 theta_w(&thetaw, allelic_map, pars.cp.nsam, segsites, hn, missing_map) &&
		 tajD(&tajd, truesegsites, pars.cp.nsam, thetaw, thetapi, hn, sqhn)
		 )
		fprintf(pf, "%e\t", tajd);
	    }
	    else
	      fprintf(pf, "nan\t");
	    
	    break;
	    
	  case 4:
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.priZnS){
	        
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		
		if(pritruesegsites[ipop] > 1 && 
		   ZnS(&prizns[ipop], list, (pars.cp.config)[ipop], segsites, // segsites is the dimension of almap
		       filter, npop_allelic_map[ipop], npop_ders[ipop], npop_missing_map[ipop]) )
		  {
		    fprintf(pf, "%e\t", prizns[ipop]);
		  }
		else{
		  fprintf(pf, "nan\t");
		}
	      }
	    }
	    
	    
	    if(ZnS(&zns, list, pars.cp.nsam, truesegsites, filter, allelic_map, ders, missing_map)){
	      fprintf(pf, "%e\t", zns);
	    }
	    else
	      fprintf(pf, "nan\t");
	    break;
	    /* 
		   in the case of multiple populations 
		   calculate the Fst
		*/
	  case 5:
	    if(pars.cp.npop > 1){
	      
	      if(truesegsites > 0 ){
		if(calculations_done != 1){
		  calculations(weights, 
			       list,
			       pars.cp.config,
			       pars.cp.npop,
			       truesegsites,
			       pars.cp.nsam,
			       
			       &piT,
			       &piS,
			       &piB,
			       &piD,
			       
			       shared,
			       private,
			       fixed_dif,
			       derived);
		  calculations_done = 1;
		}
		fst = Fst(fst_type, piT, piS, piD, piB);
		fprintf(pf, "%e\t", fst);
	      }
	      else
		fprintf(pf, "nan\t");
	    }
	    break;
	    
	  case 6:
	    if(pars.cp.npop > 1){
	      
	      if(pars.cp.npop > 1){
		if(calculations_done != 1){
		  calculations(weights, 
			       list,
			       pars.cp.config,
			       pars.cp.npop,
			       truesegsites,
			       pars.cp.nsam,
			       
			       &piT,
			       &piS,
			       &piB,
			       &piD,
			       
			       shared,
			       private,
			       fixed_dif,
			       derived);
		  calculations_done = 1;
		}
		for(npopi=0; npopi < pars.cp.npop - 1; npopi++)
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++)
		    fprintf(pf, "%e\t", shared[npopi][npopj]);
	      }
	      else
		fprintf(pf, "nan\t");
	    }
	    break;

	  case 7:
	    if(pars.cp.npop > 1){
	     
	      
	      if(pars.cp.npop > 1){
		if(calculations_done != 1){
		  calculations(weights, 
			       list,
			       pars.cp.config,
			       pars.cp.npop,
			       truesegsites,
			       pars.cp.nsam,
			       
			       &piT,
			       &piS,
			       &piB,
			       &piD,
			       
			       shared,
			       private,
			       fixed_dif,
			       derived);
		  calculations_done = 1;
		}
		for(npopi=0; npopi < pars.cp.npop - 1; npopi++)
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++)
		    fprintf(pf, "%e\t", private[npopi][npopj]);
	      }
	      else
		fprintf(pf, "nan\t");
	    }
	    break;
	    
	  
	  case 8:
	    if(pars.cp.npop > 1){
	      
	      
	      if(pars.cp.npop > 1){
		if(calculations_done != 1){
		  calculations(weights, 
			       list,
			       pars.cp.config,
			       pars.cp.npop,
			       truesegsites,
			       pars.cp.nsam,
			       
			       &piT,
			       &piS,
			       &piB,
			       &piD,
			       
			       shared,
			       private,
			       fixed_dif,
			       derived);
		  calculations_done = 1;
		  
		}
		for(npopi=0; npopi < pars.cp.npop - 1; npopi++)
		  for(npopj=npopi+1; npopj<pars.cp.npop; npopj++)
		    fprintf(pf, "%e\t", fixed_dif[npopi][npopj]);
	      }
	      else
		fprintf(pf, "nan\t");
	    }
	    break;
		 
	  case 9:
	    if(pars.cp.npop > 1){
	      
	      
	      
	      if(pars.cp.npop > 1){
		if(pairwise_calculations_done != 1){
		  for(ipop=0; ipop<fst_pops-1; ipop++){
		    for(jpop=ipop+1; jpop<fst_pops; jpop++){
		      pop1=Fst_pops[ ipop ];
		      pop2=Fst_pops[ jpop ];
		      pairwiseFstcalculations(
					      pop1, pop2,
					      weights,
					      list,
					      pars.cp.config,
					      pars.cp.npop,
					      truesegsites,
					      pars.cp.nsam,
					      
					      &ppiT,
					      &ppiS,
					      &ppiB,
					      &ppiD);
		      
		      pairwise_fst = Fst(fst_type,ppiT, ppiS, ppiD, ppiB);
		      fprintf(pf, "%e\t", pairwise_fst);
		    }
		  }
		  
		  pairwise_calculations_done = 1;
		}
	      }
	      else{
		fprintf(pf, "nan\t");
	      }
	    }
	  
	    break;
	    
	  case 10:
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && sstats.prifwh){
	      
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		
		// if we know the thetapi and thetah
		if( (pritruesegsites[ipop] && prithetapi[ipop] && prithetah[ipop] ) &&
		    htest( &prifwh[ipop], (pars.cp.config)[ipop], 
			   pritruesegsites[ipop], // the number of true segsites is needed 
			   prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
		  {
		    fprintf(pf, "%e\t", prifwh[ipop]);
		  }
		else if(pritruesegsites[ipop] && prithetapi[ipop] && (!prithetah[ipop]) ){
		  
		  if(thetaH( &prithetah[ipop], list, 
			     pars.cp.nsam, // here the pars.cp.nsam is the dimension of the list
			     segsites, // this is the other dimension of the list
			     bm,
			     npop_allelic_map[ipop], samples_b[ipop], samples_e[ipop], npop_missing_map[ipop])  &&
		     htest( &prifwh[ipop], (pars.cp.config)[ipop], 
			    pritruesegsites[ipop], // the number of true segsites is needed 
			    prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
		    {
		      fprintf(pf, "%e\t", prifwh[ipop]);
		    }
		}
		else if(pritruesegsites[ipop] && prithetapi[ipop] && (!prithetah[ipop]) ){
		  if(theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
		     htest( &prifwh[ipop], (pars.cp.config)[ipop], 
			    pritruesegsites[ipop], // the number of true segsites is needed 
			    prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
		    {
		      fprintf(pf, "%e\t", prifwh[ipop]);
		    }
		}
		else if(pritruesegsites[ipop] && (!prithetapi[ipop]) && (!prithetah[ipop]) ){
		  if( theta_pi( &prithetapi[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]) &&
		      thetaH( &prithetah[ipop], list, 
			      pars.cp.nsam, // here the pars.cp.nsam is the dimension of the list
			      segsites, // this is the other dimension of the list
			      bm,
			      npop_allelic_map[ipop], samples_b[ipop], samples_e[ipop], npop_missing_map[ipop] ) &&
		      htest( &prifwh[ipop], (pars.cp.config)[ipop], 
			     pritruesegsites[ipop], // the number of true segsites is needed 
			     prithetapi[ipop], prithetah[ipop], prihn[ipop], prisqhn[ipop], pribn[ipop] ) )
		    {
		      fprintf(pf, "%e\t", prifwh[ipop]);
		    }
		}
		else if(pritruesegsites[ipop] == 0)
		  fprintf(pf, "nan\t");
		//printf("fwh:%e, priave_fwh[%d]:%e\n", prifwh[ipop], ipop, priave_fwh[ipop]);
	      }
	    }

		    
	    
	    if( ( truesegsites && thetapi && thetah) && 
		htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
	      fprintf(pf, "%e\t", fwh);
	      
	    }
	    else if( truesegsites && thetapi && (!thetah) ){
	      if(thetaH(&thetah, list, pars.cp.nsam, truesegsites, bm, allelic_map, 0, pars.cp.nsam, missing_map) &&
		 htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		fprintf(pf, "%e\t", fwh);
	      }
	    }
	    else if( truesegsites && (!thetapi) && (thetah) ){
	      if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) &&
		 htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		fprintf(pf, "%e\t", fwh);
	      }
	    }
	    else if( truesegsites && (!thetapi) && (!thetah) ){
	      if(theta_pi(&thetapi, allelic_map, pars.cp.nsam, segsites, missing_map) &&
		 thetaH(&thetah, list, pars.cp.nsam, truesegsites, bm, allelic_map, 0, pars.cp.nsam, missing_map) &&
		 htest( &fwh,  pars.cp.nsam, truesegsites, thetapi, thetah, hn, sqhn, bn) ){
		fprintf(pf, "%e\t", fwh);
	      }
	    }
	    else if(truesegsites == 0){
	      fprintf(pf, "nan\t");
	    }
	    break;
	    
	  case 11: // dvk, dvh
	    if(sstats.pristats == 1 && pars.cp.npop > 1 && (sstats.pridvstat)){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		dvstat(list, pars.cp.nsam, segsites, samples_b[ipop], samples_e[ipop], &pridvk[ipop], &pridvh[ipop]);
		fprintf(pf, "%e\t", pridvk[ipop]);
		fprintf(pf, "%e\t", pridvh[ipop]);
	      }
	    }
	    if(dvstat(list, pars.cp.nsam, segsites, 0, pars.cp.nsam, &dvk, &dvh)){
	      //thetaw /= pars.cp.nsites;
	      fprintf(pf, "%e\t", dvk);
	      fprintf(pf, "%e\t", dvh);
	      
	      //thetaw *= pars.cp.nsites;//add by SL 0104
	      
	    }
	    else
	      fprintf(pf, "nan\t");
	    break;
	    

	  case 12: // thomson_est
	    if(sstats.pristats == 1 && pars.cp.npop > 1){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		thomsonEst(&prithomson_est[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites, npop_missing_map[ipop]);
		fprintf(pf, "%e\t", prithomson_est[ipop]);
	      }
	    }
	    
	    if(thomsonEst(&thomson_est, allelic_map, pars.cp.nsam, segsites, missing_map)){
	      //thetapi /= pars.cp.nsites;
	      fprintf (pf, "%e\t", thomson_est);
	      //thetapi *=  pars.cp.nsites;//add by SL 0104
	      
	    }
	    else
	      fprintf(pf, "nan\t");
	    break;

	    
	  case 13: // thomson_var
	    if(sstats.pristats == 1 && pars.cp.npop > 1){
	      for(ipop=0; ipop<pars.cp.npop; ipop++){
		thomsonVar(&prithomson_est[ipop], npop_allelic_map[ipop], (pars.cp.config)[ipop], segsites);
		fprintf(pf, "%e\t", prithomson_var[ipop]);
	      }
	    }
	    
	    if(thomsonVar(&thomson_var, allelic_map, pars.cp.nsam, segsites)){
	      fprintf (pf, "%e\t", thomson_var);
	    }
	    else
	      fprintf(pf, "nan\t");
	    break;

	    
	  default:
	    fprintf(stderr, "Statistic is not implemented yet!\n");
	    break;
	  }
	  
	}
	
      }
      /* fprintf(stderr, "argv[8]: %s\n", argv[8]); */
      fprintf(pf, "\n");
      /* if(count == 1) */
/* 	fclose(headerfile); */
    }
    
  
  }


  /*********************** FREE MEMORY ****************************/
  
  if(pars.cp.npop > 1){
    for(i=0; i<MAXSUMSTAT; i++){
      for(j=0; j<pars.cp.npop; j++)
	free(sstats_popdenpair_weights[i][j]);
      free(sstats_popdenpair_weights[i]);
    }
    free(sstats_popdenpair_weights);
  }


  if(missingin != NULL){
    /* missing data */
    for(i=0; i<missing_nof_fragments; i++){
      free(misseq[i]);
      free(misstart[i]);
      free(misend[i]);
    }
    free(misseq);
    free(misstart);
    free(misend);
    free(missingdata);
    /****************/
  }
  


  if(sumstat && pars.cp.npop > 1 && sstats.pristats == 1){

    for(i=0; i<MAXSUMSTAT; i++)
      free(sstats_popden_weights[i]);
    free(sstats_popden_weights);

    for(i=0; i<pars.cp.npop; i++){
      free(npop_allelic_map[i]);
      free(npop_missing_map[i]);
    }

    free(npop_allelic_map);
    free(npop_missing_map);

    for(ipop=0; ipop<pars.cp.npop; ipop++){
      for(i=0; i<largesample; i++)
	free(npop_ders[ipop][i]);
    }
    for(ipop=0; ipop<pars.cp.npop; ipop++)
      free(npop_ders[ipop]);
    
    free(npop_ders);
    
    free(prisegs);
    free(prithetapi);
    free(prithetaw);
    free(pritajd);
    free(prithetah);
    free(prifwh);
    free(prizns);
    free(prithomson_est);
    free(prithomson_var);
    
    free(priave_segs);
    free(priave_thetapi);
    free(priave_thetaw);
    free( priave_fwh);
    free( priave_tajimasD);
    free( priave_ZnS);
    free( privar_segs);
    free( privar_thetapi);
    free( privar_thetaw);
    free( privar_fwh);
    free( privar_tajimasD);
    free( privar_ZnS);
    free(priave_dvk);
    free(priave_dvh);
    free(privar_dvk);
    free(privar_dvh);
    free(privar_thomson_est);
    free(privar_thomson_var);

    free(priskew_segs);
    free(priskew_thetapi);
    free(priskew_thetaw);
    free(priskew_fwh);
    free(priskew_tajimasD);
    free(priskew_ZnS);
    free(priskew_dvk);
    free(priskew_dvh);
    free(priskew_thomson_est);
    free(priskew_thomson_var);

    free(prikurt_segs);
    free(prikurt_thetapi);
    free(prikurt_thetaw);
    free(prikurt_fwh);
    free(prikurt_tajimasD);
    free(prikurt_ZnS);
    free(prikurt_dvk);
    free(prikurt_dvh);
    free(prikurt_thomson_est);
    free(prikurt_thomson_var);

    free(pritruesegsites);
    free(pribn);
    free(prisqhn);
    free(prihn);
    free(pridvk);
    free(pridvh);
    
  }
  
  free(samples_e);
  free(samples_b);
  
  free(param_distr);
  free(weights);
  free(global_values);

  

  free(Fst_pops);
  for(i=0; i<pars.cp.npop; i++){
    free(shared[i]);
    free(private[i]);
    free(fixed_dif[i]);
    free(ave_shared[i]);
    free(ave_private[i]);
    free(ave_pairwise_Fst[i]);
    free(var_shared[i]);
    free(var_private[i]);
    free(var_fixed_dif[i]);
    free(var_pairwise_Fst[i]);
    free(ave_fixed_dif[i]);
    free(skew_shared[i]);
    free(skew_private[i]);
    free(skew_fixed_dif[i]);
    free(skew_pairwise_Fst[i]);
    free(kurt_shared[i]);
    free(kurt_private[i]);
    free(kurt_fixed_dif[i]);
    free(kurt_pairwise_Fst[i]);
  }
  free(shared);
  free(private);
  free(fixed_dif);
  free(ave_shared);
  free(ave_private);
  free(ave_pairwise_Fst);
  free(var_shared);
  free(var_private);
  free(var_fixed_dif);
  free(var_pairwise_Fst);
  free(ave_fixed_dif);
  free(skew_shared);
  free(skew_private);
  free(skew_fixed_dif);
  free(skew_pairwise_Fst);
  free(kurt_shared);
  free(kurt_private);
  free(kurt_fixed_dif);
  free(kurt_pairwise_Fst);
  
  for(i=0; i<argc; i++){
    free(tbsparamstrs[i]);
    free(param_distr_values[i]);
  }
  free(tbsparamstrs);
  free(param_distr_values); param_distr_values = NULL;

  for(i=0; i<largestnsam; i++)
      free(list[i]);
  free(list);

  if(sumstat){
    for(i=0; i<largesample; i++)
      free(ders[i]);
    free(ders);
    free(allelic_map); allelic_map = NULL;
    free(missing_map); missing_map = NULL;
  }

  free(posit); posit = NULL;

  free(timeevents); timeevents = NULL;
  
#ifndef R_PACKAGE_BUILD
  if( !pars.commandlineseedflag ) seedit( "end" );
#endif
  if(tempdatafile != NULL)
    fclose(tempdatafile);

  for(i=0; i<argc; i++)
    free(conditional_timeevents[i]);
  free(conditional_timeevents); conditional_timeevents = NULL;

  //fclose(timefile);

}

#ifdef R_PACKAGE_BUILD
void msABC_set_jmpbuf_active(int active) {
    msABC_jmpbuf_set = active;
}
#endif




#ifdef R_PACKAGE_BUILD
static void sfs_from_tree(struct node *ptree, int nsam, double theta_seg,
                           int *config, int npop);
#endif

int gensam( char **list, double *pprobss, double *ptmrca, double *pttot, int largestnsam, int old_maxnsam)
{

  //fprintf(stderr, "pars.mp.theta: %e\n", pars.mp.theta);
  int nsegs, h, i, k, j, seg, ns, start, end, len, segsit ;
  struct segl *seglst, *segtre_mig(struct c_params *p, int *nsegs) ; /* used to be: [MAXSEG];  */
  double nsinv,  tseg, tt, ttime(struct node *, int nsam), ttimemf(struct node *, int nsam, int mfreq) ;
  double *pk;
  int *ss;
  int segsitesin,nsites;
  double theta, es ;
  int nsam, mfreq ;
  void prtree( struct node *ptree, int nsam);
  int make_gametes(int nsam, int mfreq,  struct node *ptree, double tt, int newsites, int ns, char **list );
  void ndes_setup( struct node *, int nsam );
  
  //  fprintf(stderr, "gensam starts...\n"); //$$p
  //fprintf(stderr, "#P pars.cp.nsams %d\n", pars.cp.nsam); //$$p
  
  nsites = pars.cp.nsites ;
  nsinv = 1./nsites;
 
  seglst = segtre_mig(&(pars.cp),  &nsegs) ; //******************
 
  nsam = pars.cp.nsam;
  segsitesin = pars.mp.segsitesin ;
  theta = pars.mp.theta ;
  
  

  mfreq = pars.mp.mfreq ;
  //  fprintf(stderr, "#P pars.mp.theta: %e\t%e\n", pars.mp.theta, theta); //$$p
  //fprintf(stderr, "#P pars.cp.nsites: %d\n", pars.cp.nsites);//$$p
  
  

/*   fprintf(stderr, "nsites: %d\nnsiv: %e\nnsam: %d\nsegsitesin: %d\ntheta: %e\nmfreq: %d\n\n", */
/*                      nsites, nsinv,  nsam, segsitesin, pars.mp.theta, mfreq); //$$p */
  
  if( pars.mp.treeflag ) {
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
      if( (pars.cp.r > 0.0 ) || (pars.cp.f > 0.0) ){
	end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
	start = seglst[seg].beg ;
	len = end - start + 1 ;
	fprintf(stdout,"[%d]", len);
      }
      prtree( seglst[seg].ptree, nsam ) ;
      if( (segsitesin == 0) && ( theta == 0.0 ) && ( pars.mp.timeflag == 0 ) ) 
	free(seglst[seg].ptree) ;
    }
  }
  
  if( pars.mp.timeflag ) {
    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      if( (nsegs==1) || ( ( start <= nsites/2) && ( end >= nsites/2 ) ) )
	*ptmrca = (seglst[seg].ptree + 2*nsam-2) -> time ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) tt += ttime(seglst[seg].ptree,nsam)*tseg ;
      else tt += ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
      if( (segsitesin == 0) && ( theta == 0.0 )  ) 
	free(seglst[seg].ptree) ;
    }
    *pttot = tt ;
  }	
  
  if( (segsitesin == 0) && ( theta > 0.0)   ) {
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) {
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len*(theta/nsites) ;
#ifdef R_PACKAGE_BUILD
      if (sfs_expected_mode && sfs_mode_active) {
        ndes_setup(seglst[seg].ptree, nsam);
        sfs_from_tree(seglst[seg].ptree, nsam, tseg,
                      pars.cp.config, pars.cp.npop);
        free(seglst[seg].ptree);
        continue;
      }
#endif
      if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
      if( mfreq == 1) tt = ttime(seglst[seg].ptree, nsam);
      else tt = ttimemf(seglst[seg].ptree, nsam, mfreq );
      segsit = poisso( tseg*tt );
      if( (segsit + ns) >= maxsites ) {
	maxsites = segsit + ns + SITESINC ;
	posit = (double *)realloc(posit, maxsites*sizeof(double) ) ;
	/*modified by JRI: needs to be largestnsam to cover all values in list*/
	biggerlist(largestnsam, list) ;
      }
      make_gametes(nsam,mfreq,seglst[seg].ptree,tt, segsit, ns, list );
      free(seglst[seg].ptree) ;
      locate(segsit,start*nsinv, len*nsinv,posit+ns);
      ns += segsit;
    }
  }
  else if( segsitesin > 0 ) {
    
    pk = (double *)malloc((unsigned)(nsegs*sizeof(double)));
    ss = (int *)malloc((unsigned)(nsegs*sizeof(int)));
    if( (pk==NULL) || (ss==NULL) ) perror("malloc error. gensam.2");
    
    
    tt = 0.0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      if( mfreq > 1 ) ndes_setup( seglst[seg].ptree, nsam );
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites ;
      if( mfreq == 1 ) pk[k] = ttime(seglst[seg].ptree,nsam)*tseg ;
      else pk[k] = ttimemf(seglst[seg].ptree,nsam, mfreq)*tseg ;
      tt += pk[k] ;
    }
    if( theta > 0.0 ) { 
      es = theta * tt ;
      *pprobss = exp( -es )*pow( es, (double) segsitesin) / segfac ;
    }
    if( tt > 0.0 ) {
      for (k=0;k<nsegs;k++) pk[k] /= tt ;
      mnmial(segsitesin,nsegs,pk,ss);
    }
    else
      for( k=0; k<nsegs; k++) ss[k] = 0 ;
    ns = 0 ;
    for( seg=0, k=0; k<nsegs; seg=seglst[seg].next, k++) { 
      end = ( k<nsegs-1 ? seglst[seglst[seg].next].beg -1 : nsites-1 );
      start = seglst[seg].beg ;
      len = end - start + 1 ;
      tseg = len/(double)nsites;
      make_gametes(nsam,mfreq,seglst[seg].ptree,tt*pk[k]/tseg, ss[k], ns, list);
      
      free(seglst[seg].ptree) ;
      locate(ss[k],start*nsinv, len*nsinv,posit+ns);   
      ns += ss[k] ;
    }
    free(pk);
    free(ss);
    
  }
  for(i=0;i<nsam;i++) list[i][ns] = '\0' ;
  //fprintf(stderr, "#P segsites from sims: %d\n", ns);
  return( ns ) ;
}

void 
ndes_setup(struct node *ptree, int nsam )
{
	int i ;

	for( i=0; i<nsam; i++) (ptree+i)->ndes = 1 ;
	for( i = nsam; i< 2*nsam -1; i++) (ptree+i)->ndes = 0 ;
	for( i= 0; i< 2*nsam -2 ; i++)  (ptree+((ptree+i)->abv))->ndes += (ptree+i)->ndes ;

}

#ifdef R_PACKAGE_BUILD
/*
 * ndes_per_pop_setup - compute per-population descendant counts for each node.
 * pop_ndes is a (2*nsam-1) x npop array (row-major): pop_ndes[node*npop + pop].
 * Tips 0..config[0]-1 belong to pop 0, config[0]..config[0]+config[1]-1 to pop 1, etc.
 */
static void
ndes_per_pop_setup(struct node *ptree, int nsam, int *config, int npop, int *pop_ndes)
{
	int i, p, tip_start;

	/* Zero the array */
	memset(pop_ndes, 0, (2*nsam - 1) * npop * sizeof(int));

	/* Assign tips to populations */
	tip_start = 0;
	for (p = 0; p < npop; p++) {
		for (i = tip_start; i < tip_start + config[p]; i++) {
			pop_ndes[i * npop + p] = 1;
		}
		tip_start += config[p];
	}

	/* Propagate upward */
	for (i = 0; i < 2*nsam - 2; i++) {
		int parent = (ptree+i)->abv;
		for (p = 0; p < npop; p++) {
			pop_ndes[parent * npop + p] += pop_ndes[i * npop + p];
		}
	}
}

/*
 * sfs_from_tree - compute expected SFS directly from tree branch lengths.
 * For each branch i: E[SFS[ndes_i]] += theta_seg * b_i
 * For one.snp mode: SFS[k] += L_k / T_total (probability a random mutation is in class k)
 */
static void
sfs_from_tree(struct node *ptree, int nsam, double theta_seg,
              int *config, int npop)
{
	int i, p;

	if (npop == 1) {
		/* ---- 1D SFS ---- */
		if (sfs_one_snp) {
			/* Accumulate branch lengths per frequency class, then normalize */
			double *temp = (double *)calloc(nsam - 1, sizeof(double));
			if (temp == NULL) return;
			double t_total = 0.0;
			for (i = 0; i < 2*nsam - 2; i++) {
				double b_i = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
				int nd = (ptree+i)->ndes;
				t_total += b_i;
				if (nd >= 1 && nd <= nsam - 1) {
					temp[nd - 1] += b_i;
				}
			}
			if (t_total > 0.0) {
				for (i = 0; i < nsam - 1; i++) {
					sfs_accumulator[i] += temp[i] / t_total;
				}
			}
			free(temp);
		} else {
			for (i = 0; i < 2*nsam - 2; i++) {
				double b_i = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
				int nd = (ptree+i)->ndes;
				if (nd >= 1 && nd <= nsam - 1) {
					sfs_accumulator[nd - 1] += theta_seg * b_i;
				}
			}
		}
	} else {
		/* ---- Joint SFS (multi-population) ---- */
		int *pop_ndes = (int *)malloc((2*nsam - 1) * npop * sizeof(int));
		if (pop_ndes == NULL) return;

		ndes_per_pop_setup(ptree, nsam, config, npop, pop_ndes);

		/* Compute strides for flat indexing (pop 0 varies fastest) */
		int strides[npop];
		strides[0] = 1;
		for (p = 1; p < npop; p++) {
			strides[p] = strides[p - 1] * (config[p - 1] + 1);
		}

		if (sfs_one_snp) {
			/* Accumulate branch lengths per joint class, then normalize */
			double *temp = (double *)calloc(sfs_accum_len, sizeof(double));
			if (temp == NULL) { free(pop_ndes); return; }
			double t_total = 0.0;
			for (i = 0; i < 2*nsam - 2; i++) {
				double b_i = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;
				t_total += b_i;

				/* Compute flat index from per-pop descendants */
				int flat_idx = 0;
				int valid = 1;
				for (p = 0; p < npop; p++) {
					int nd_p = pop_ndes[i * npop + p];
					if (nd_p < 0 || nd_p > config[p]) { valid = 0; break; }
					flat_idx += nd_p * strides[p];
				}
				if (valid && flat_idx >= 0 && flat_idx < sfs_accum_len) {
					temp[flat_idx] += b_i;
				}
			}
			if (t_total > 0.0) {
				for (i = 0; i < sfs_accum_len; i++) {
					sfs_accumulator[i] += temp[i] / t_total;
				}
			}
			free(temp);
		} else {
			for (i = 0; i < 2*nsam - 2; i++) {
				double b_i = (ptree + (ptree+i)->abv)->time - (ptree+i)->time;

				int flat_idx = 0;
				int valid = 1;
				for (p = 0; p < npop; p++) {
					int nd_p = pop_ndes[i * npop + p];
					if (nd_p < 0 || nd_p > config[p]) { valid = 0; break; }
					flat_idx += nd_p * strides[p];
				}
				if (valid && flat_idx >= 0 && flat_idx < sfs_accum_len) {
					sfs_accumulator[flat_idx] += theta_seg * b_i;
				}
			}
		}

		free(pop_ndes);
	}
}
#endif /* R_PACKAGE_BUILD */

/*modified by JRI*/
int biggerlist(int largestnsam, char ** list ){
	int i;
	for( i=0; i<largestnsam; i++){	
		list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
		if( list[i] == NULL ) perror( "realloc error. bigger");
	}
}

		   


/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned) nsam*sizeof( char* ) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) len*sizeof( char ) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}


/* Added by JRI Aug 2007.  Reallocates space for gametes (character strings) */
char ** rematrix(int nsam, int largestnsam, int howlong, char ** m){   
  int i, j;
  
  if( ! (m = (char **) realloc(m,  ((unsigned) nsam*sizeof( char* )) ) ) ) {
    perror("realloc error in rematrix");
  }
  
  for( i=0; i<largestnsam; i++) {
    if(   ! (  m[i] = (char *) realloc( m[i], (unsigned) howlong*sizeof(char) )  )   ) {
      perror("realloc error in rematrix. 2");
    }
  }
  for( j=largestnsam; j<nsam; j++) {
    if( ! ( m[j] = (char *) malloc( ((unsigned) howlong*sizeof( char )) ))){
      perror("malloc error in rematrix. 3");
    }
  }
  return( m );
}

char** reinitialize(int nsam, int previousnsam, int howlong, char** m){
  int i,j;
  for(i=0; i<previousnsam; i++){
    if(m[i] != NULL)
      free(m[i]);
  }
  if(m != NULL) free(m);
  
  m=(char**)malloc((unsigned)nsam*sizeof(char*));
  for(i=0; i<nsam; i++)
    m[i] = (char*)malloc((unsigned)howlong*sizeof(char));
  return(m);
}
  
    


int ** rematrix_int(int nsam, int largestnsam, int howlong, int ** m){   
	int i, j;
	if( ! (m = (int**) realloc(m,  ((unsigned) nsam*sizeof( int* )) ) ) ) {
		perror("realloc error in rematrix");
	}
	for( i=0; i<largestnsam; i++) {
		if(   ! (  m[i] = (int *) realloc( m[i], (unsigned) howlong*sizeof(int) )  )   ) {
			perror("realloc error in rematrix. 2");
		}
	}
 	for( j=largestnsam; j<nsam; j++) {
 		if( ! ( m[j] = (int*) malloc( ((unsigned) howlong*sizeof( int )) ))){
 			perror("malloc error in rematrix. 3");
 		}
	}
	return( m );
}


double ** rematrix_double(int newy, int oldy, int x, double ** m){   
	int i, j;
	if( ! (m = realloc(m,  ((unsigned)newy*sizeof( double* )) ) ) ) {
	  perror("realloc error in rematrix");
	}
	for( i=0; i<oldy; i++) {
	  if(   ! (  m[i] = realloc( m[i], (unsigned)x*sizeof(double) )  )   ) {
	    perror("realloc error in rematrix. 2");
	  }
	}
 	for( j=oldy; j<newy; j++) {
	  if( ! ( m[j] = malloc( ((unsigned)x*sizeof( double )) ))){
	    perror("malloc error in rematrix. 3");
	  }
	}
	return( m );
}



	int
locate(n,beg,len,ptr)
	int n;
	double beg, len, *ptr;
{
	int i;

	ordran(n,ptr);
	for(i=0; i<n; i++)
		ptr[i] = beg + ptr[i]*len ;

}

int NSEEDS = 3 ;

  void
getpars(int argc, char *argv[], int *phowmany )
{

  int durationmode_events = 0;
  int durationmode_counter = 0;
  int arg, carg, i, j, sum , pop , argstart, npop , npop2, pop2 ;
  int recgiven=0; // in the persite ms model the number of recombinant fragments should be given, even if the recombination is 0
  double migr=0.;
  double mij, psize, palpha ;
  void addtoelist( struct devent *pt, struct devent *elist ); 
  void argcheck( int arg, int argc, char ** ) ;
  int commandlineseed( char ** ) ;
  void free_eventlist( struct devent *pt, int npop );
  void checkarguments(char**, char*, int, int);
  
  double get_param( int arg, char distr);
  int dist_argcheck(int arg, int  argc, char** argv);
  int distribution_arguments( char distr );
  
  struct devent *ptemp , *pt ;
  FILE *pf ;
  char ch3 ;
  
  // for storing the order of the events;
  event=0;

  /** countis a global variable */
  /** if count == 0 then this is before the first replication */
  if( count == 0 ) {
    
    if( argc < 4 ){ 
      
      for(arg=0; arg < argc; arg++){
	if(!strcmp(argv[arg], "-ms")){
	  usage1();
	  EXIT_MSABC(1);
	}
	else if(!strcmp(argv[arg], "-h")){
	  usage2();
	  EXIT_MSABC(1);
	}
	else if(!strcmp(argv[arg], "-p")){
	  usage3();
	  EXIT_MSABC(1);
	}
	else if(!strcmp(argv[arg], "-all")){
	  usage1();
	  usage2();
	  usage3();
	  EXIT_MSABC(1);
	}
      }
      fprintf(stderr,"Too few command line arguments\n"); usage();
      usage();
    }

    
    /** read the sample size and store it */
    pars.cp.nsam = atoi( argv[1] );
    //largestnsam = pars.cp.nsam;
    if( pars.cp.nsam <= 0 ) { fprintf(stderr,"First argument error. nsam <= 0. \n"); usage();}
    
    /** read the number of replications and store it */
    *phowmany = atoi( argv[2] );
    if( *phowmany  <= 0 ) { fprintf(stderr,"Second argument error. howmany <= 0. \n"); usage();}
    
    
    pars.commandlineseedflag = 0 ;
    pars.cp.r = pars.mp.theta =  pars.cp.f = 0.0 ;
    pars.cp.track_len = 0. ;
    pars.cp.npop = npop = 1 ;
    
    /** alocate space for the migration matrix 
	just denote that it will be a 2d array. real allocations hasn't happened yet */
    pars.cp.mig_mat = (double **)malloc( (unsigned) sizeof( double *) );
    pars.cp.mig_mat[0] = (double *)malloc( (unsigned)sizeof(double ));
    
    
    pars.cp.mig_mat[0][0] =  0.0 ;
    pars.mp.segsitesin = 0 ;
    pars.mp.treeflag = 0 ;
    pars.mp.timeflag = 0 ;
    
    /** just print out only polymorphic sites */
    /** maybe this could be user defined... maybe... */
    pars.mp.mfreq = 1 ; // TODO
    
    pars.cp.config = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    (pars.cp.config)[0] = pars.cp.nsam ;
    pars.cp.initconfig = (int *) malloc( (unsigned)(( pars.cp.npop +1 ) *sizeof( int)) );
    (pars.cp.initconfig)[0] = pars.cp.nsam ;
    pars.cp.size= (double *) malloc( (unsigned)( pars.cp.npop *sizeof( double )) );
    (pars.cp.size)[0] = 1.0  ;
    pars.cp.alphag = (double *) malloc( (unsigned)(( pars.cp.npop ) *sizeof( double )) );
    (pars.cp.alphag)[0] = 0.0  ;
    pars.cp.nsites = 2 ;
  }
  else{
    /* next two lines added by JRI Aug. 2007.  Sets pop size to new tbs popsize */
    pars.cp.nsam = atoi( argv[1] ); 	
    (pars.cp.config)[0] = pars.cp.nsam ; 
    npop = pars.cp.npop ;
    free_eventlist( pars.cp.deventlist, npop ); // 20091015
  }

  //free_eventlist( pars.cp.deventlist, npop ); // 20091015

  pars.cp.deventlist = NULL ;
  
  arg = 3 ;
  double current_time = 0.;

  // initialize fragmode (i.e. no fragment mode)
  fragmode = 0;
  fragindex=0;

  fcolumns = 0;
  
  filter = 0;
  resample = 0;
  
  int cr = 0; // conditional resampling

  while( arg < argc ){
    /* fprintf(stderr, "argv[9]: %s \t argv[arg] : %s\t\targv[8]: %s, %c, %c\n", argv[9], argv[arg], argv[8], argv[8][2], argv[8][4]); */
    
    int c;
    if(!count)
      checkarguments(argv, argv[arg], arg, argc);

    if(!strcmp(argv[arg], "--outgroup"))
      {
	outgroup = atoi(argv[++arg]);
	++arg; 
	continue;
      }

    if(!strcmp(argv[arg], "--phasemode"))
      {
	phasemode = atoi(argv[++arg]);
	++arg;
	continue;
      }

    if(!strcmp(argv[arg], "--DER")){
      derived = atoi(argv[++arg]);
      assert(derived == 0 || derived ==1 || derived == 2);
      ++arg;
      continue;
    }
    
    if(!strcmp(argv[arg], "--obs")){
      strcpy(obfile, argv[++arg]);
      ++arg;
      continue;
    }

    if(!strcmp(argv[arg], "--poldata")){
      strcpy(poldatafile, argv[++arg]);
      ++arg;
      continue;
    }

    if(!strcmp(argv[arg], "--pf")){
      propfragments = atof(argv[++arg]);
      ++arg;
      continue;
    }

    if(!strcmp(argv[arg], "--missing")){
      strcpy(missingfile, argv[++arg]);
      ++arg;
      continue;
    }

    if(!strcmp(argv[arg], "--fst")){		
      strcpy(fst_type, argv[++arg]);
      ++arg;
      continue;
    }
    
    if(!strcmp(argv[arg], "--resample")){
      c=0;
      
      if( (arg + 2 < argc) && (argv[arg+1][0] != '-') ){
	
	c = atoi( argv[++arg] ) + 1; // +1 because the next number will denote how many events will be there
	conditional_timeevents[cr][0] = c - 1;
	for( j=1; j<c; j++){ // it starts from 1 becasue pos 1 is allocated by the number of elements; c = elements + 1
	  assert(argv[++arg][0] != '-');
	  conditional_timeevents[cr][j] = atoi(argv[arg]);
	}
	cr++;
	resample = 2;
	nof_resampling = cr;
      }
      else
	resample = 1;

      arg++;
      //fprintf(stderr, "c: resample: %d\t%s\n", resample, argv[arg]);
      
      continue;
    }

    if(!strcmp(argv[arg], "--bm")){
      bm = atof(argv[++arg]);
      ++arg;
      continue;
    }
    if(!strcmp(argv[arg], "--options")){ 
      strcpy( optionsfile, argv[++arg]);
      ++arg;
      continue;
    }
    if(!strcmp(argv[arg], "--verbose")){
      printall = 1;
      ++arg;
      continue;
    }
    
    if(!strcmp(argv[arg], "--filter")){ filter = 1; ++arg; continue; }
    
    //fprintf(stderr, "arg: %d\targv[arg]: %s\t argc: %d\n", arg, argv[arg], argc); //#p
    if(!strcmp(argv[arg], "--dur-mode")){ 
      resample = 0; 
      durationmode = 1; 
      if( (arg + 1< argc) && argv[arg+1][0] != '-' )
	{
	  durationmode_events = atoi(argv[++arg]);
	  durationmode_counter = 0;
	}
      
      ++arg;
      
      continue;
    }
    
    //fragments

    if(!strcmp(argv[arg], "--frag-begin")){ 
      //fprintf(stderr, "fragments... command line\n");
      fragmode = 1; 
      while(strcmp( argv[++arg], "--frag-end")){
	//fprintf(stderr, "i: %d, arg: %s\n", arg, argv[arg]);
	if(!strcmp( argv[arg], "--finp")){
	  strcpy(fragfile, argv[++arg]); 
	  continue;
	}
	/* number of rows (fragments) in the file, 
	   without the header
	*/
	if(!strcmp( argv[arg], "--nf")){
	  nofragments = atoi(argv[++arg]); 
	  continue;
	}
	/* number of columns in the file
	 */
	if(!strcmp( argv[arg], "--np")){
	  fcolumns = atoi(argv[++arg]);
	  continue;
	}
	
	/* what kind of weights to use */
	if(!strcmp( argv[arg], "--w")){
	  weightmode = atoi(argv[++arg]); 
	  continue;
	}
	
	if(!strcmp( argv[arg], "--N")){
	  c = 0;
	  arg++;
	  //printf("-count %d, param-distr: %c\n", count, param_distr[arg-1]);
	  if( !count && (c = dist_argcheck(arg, argc, argv)) ){
	    
	  }
	  else if( count && param_distr[arg - 1] != '-'){
	    //fprintf(stderr, "now I check theta.3 - distr: %c\n", argv[arg][1]); 
	    effectiveN = get_param(arg - 1, argv[arg][1]);
	    global_values[arg-1] = effectiveN;
	    arg = arg + distribution_arguments( param_distr[arg - 1] ) + 1;
	    //fprintf(stderr, "-effective: %e\n", effectiveN);
	    
	  }
	  else{
	    argcheck( arg, argc, argv);
	    effectiveN = atof(  argv[arg] );
	    global_values[arg-1] = effectiveN;
	    arg++;
	  }
      	  //fprintf(stderr, "c: %d\n", c); //#p
	  //fprintf(stderr, "effective: %e\n", effectiveN);
	  arg += c-1;
	  continue;
	}
	//fprintf(stderr, "i: %d, arg: %s\n", arg, argv[arg]);
	fprintf(stderr, "argument %s is not supported. Check the manual -- FRAG HELP\n", argv[arg]);
	EXIT_MSABC(-1);
      }
      //fprintf(stderr, "*i: %d, arg: %s, effectiveN: %e\n", arg, argv[arg], effectiveN);
      if(!strcmp(argv[arg], "--frag-end")){ arg++; continue; }
      
    }
    if(!strcmp(argv[arg], "--printall")){ 
      arg++; 
      printall = 1; 
      continue; 
    }
   
    

    

    if( argv[arg][0] != '-' ) { fprintf(stderr," argument should be -%s ?\n", argv[arg]); usage();}
    switch ( argv[arg][1] ){
    case 'f' :
      if( ntbs > 0 ) { fprintf(stderr," can't use tbs args and -f option.\n"); EXIT_MSABC(1); }
      arg++;
      argcheck( arg, argc, argv);
      pf = fopen( argv[arg], "r" ) ;
      if( pf == NULL ) {fprintf(stderr," no parameter file %s\n", argv[arg] ); EXIT_MSABC(0);}
      arg++;
      argc++ ;
      argv = (char **)malloc(  (unsigned)(argc+1)*sizeof( char *) ) ;
      argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      argstart = arg ;
      while( fscanf(pf," %s", argv[arg]) != EOF ) {
	arg++;
	argc++;
	argv = (char **)realloc( argv, (unsigned)argc*sizeof( char*) ) ;
	argv[arg] = (char *)malloc( (unsigned)(20*sizeof( char )) ) ;
      }
      fclose(pf);
      
      argc--;
      arg = argstart ;
      break;
    case 'r' : 
      recgiven = 1;
      arg++;
      c = 0;
      //fprintf(stderr, "parsing r\n"); //#p
      if( !count && (c = dist_argcheck(arg, argc, argv)) ){
	//fprintf(stderr, "parsing r.1\n"); //#p
	/* argv[arg][1] should denote the distributio i.e. U (for uniform), G (for gamma) etc */
	//fprintf(stderr, "now I check r.2 - distr: %c\n", argv[arg][1]); //#p
	//pars.mp.theta = get_param(arg-1, argv[arg][1]);
	
      }
      else if( count && param_distr[arg - 1] != '-'){
	//fprintf(stderr, "parsing r.2\n"); //#p
	pars.cp.r = get_param(arg - 1, argv[arg][1]);
	global_values[arg-1] = pars.cp.r;
	arg = arg + distribution_arguments( param_distr[arg - 1] ) + 1;
	  
      }
      else{
	//fprintf(stderr, "parsing r.3\n"); //#p
    	argcheck( arg, argc, argv);
	pars.cp.r = atof(  argv[arg++] );
	global_values[arg-2] = pars.cp.r;
      }
      
      //fprintf(stderr, "c: %d\n", c); //#p
      arg += c;
      //fprintf(stderr, "arg: %d\targv[arg]: %s\t argc: %d\n", arg, argv[arg], argc); //#p
      //if(count)
	//fprintf(stderr, "r:%e\n", pars.cp.r); //#p

     
     
      
       if( arg >= argc || argv[arg][0] == '-' || atoi(argv[arg]) <= 1){
	fprintf(stderr,"with -r option must specify both rec_rate and nsites > 1 (argv[arg]: %s)\n", argv[arg]);
	usage();
      }
      
      pars.cp.nsites = atoi( argv[arg++]);
      //while(++arg < argc && argv[arg][0] != '-');
      break;		
    case 'c' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.cp.f = atof(  argv[arg++] );
      argcheck( arg, argc, argv);
      pars.cp.track_len = atof( argv[arg++]);
      if( pars.cp.track_len <1. ){
	fprintf(stderr,"with -c option must specify both f and track_len>0\n");
	usage();
      }
      break;		
    case 't' : 
      arg++;
      //fprintf(stderr, "now I check theta\n"); //#p
      c = 0;
      if( !count && (c = dist_argcheck(arg, argc, argv)) ){
	/* argv[arg][1] should denote the distributio i.e. U (for uniform), G (for gamma) etc */
	//fprintf(stderr, "now I check theta.2 - distr: %c\n", argv[arg][1]); //#p
	//pars.mp.theta = get_param(arg-1, argv[arg][1]);
	
      }
      else if( count && param_distr[arg - 1] != '-'){
	//fprintf(stderr, "now I check theta.3 - distr: %c\n", argv[arg][1]); 
	pars.mp.theta = get_param(arg - 1, argv[arg][1]);
	global_values[arg-1] = pars.mp.theta;
	arg = arg + distribution_arguments( param_distr[arg - 1] ) + 1;
      }
      else{
    	argcheck( arg, argc, argv);
	pars.mp.theta = atof(  argv[arg++] );
	global_values[arg-2] = pars.mp.theta;
      }
      
      //fprintf(stderr, "c: %d\n", c); //#p
      arg += c;
      //fprintf(stderr, "arg: %d\targv[arg]: %s\t argc: %d\n", arg, argv[arg], argc); //#p
      if(count){}
      //fprintf(stderr, "theta:%e\n", pars.mp.theta); //#p
      break;
    case 's' : 
      arg++;
      /* fprintf(stderr, "** argv[arg-1]: %s, argv[arg-1][2]: %c, argv[8]: %s\n", argv[arg-1], argv[arg-1][2], argv[8]); */
      argcheck( arg, argc, argv);
      if( argv[arg-1][2] == 'e' ){  /* command line seeds */
	pars.commandlineseedflag = 1 ;
	if( count == 0 ) nseeds = commandlineseed(argv+arg );
	arg += nseeds ;
      }
      else {
	pars.mp.segsitesin = atoi(  argv[arg++] );
      }
      break;
    case 'F' : 
      arg++;
      argcheck( arg, argc, argv);
      pars.mp.mfreq = atoi(  argv[arg++] );
      if( (pars.mp.mfreq < 2 ) || (pars.mp.mfreq > pars.cp.nsam/2 ) ){
	fprintf(stderr," mfreq must be >= 2 and <= nsam/2.\n");
	usage();
      }
      break;
    case 'T' : 
      pars.mp.treeflag = 1 ;
      arg++;
      break;
    case 'L' : 
      pars.mp.timeflag = 1 ;
      arg++;
      break;
    case 'I' : 
      c = 0;
      arg++;
      if( count == 0 ) {
	argcheck( arg, argc, argv);
	pars.cp.npop = atoi( argv[arg]);
	pars.cp.config = (int *) realloc( pars.cp.config, (unsigned)( pars.cp.npop*sizeof( int)));
	pars.cp.initconfig = (int *) realloc( pars.cp.initconfig, (unsigned)( pars.cp.npop*sizeof( int)));
	npop = pars.cp.npop ;
      }
      arg++;
      for( i=0; i< pars.cp.npop; i++) {
	argcheck( arg, argc, argv);
	pars.cp.config[i] = atoi( argv[arg++]);
	pars.cp.initconfig[i] = pars.cp.config[i];
      }
      if( count == 0 ){
	pars.cp.mig_mat = 
	  (double **)realloc(pars.cp.mig_mat, (unsigned)(pars.cp.npop*sizeof(double *) )) ;
	pars.cp.mig_mat[0] = 
	  (double *)realloc(pars.cp.mig_mat[0], (unsigned)( pars.cp.npop*sizeof(double)));
	for(i=1; i<pars.cp.npop; i++) 
	  pars.cp.mig_mat[i] = (double *)malloc( (unsigned)( pars.cp.npop*sizeof(double)));
	
	pars.cp.size = (double *)realloc( pars.cp.size, (unsigned)( pars.cp.npop*sizeof( double )));

	pars.cp.alphag =  (double *) realloc( pars.cp.alphag, (unsigned)( pars.cp.npop*sizeof( double )));
	for( i=1; i< pars.cp.npop ; i++) {
	  (pars.cp.size)[i] = (pars.cp.size)[0]  ;
	  (pars.cp.alphag)[i] = (pars.cp.alphag)[0] ;
	}
      }
      
      if(!count && (c = dist_argcheck(arg, argc, argv) ) ){
	//printf("here!\n");
      }
      else if( count && param_distr[arg - 1] != '-'){
	migr = get_param(arg - 1, argv[arg][1]);
	global_values[arg-1] = migr;
	arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
      }
      
      else if(arg<argc && argv[arg][0] != '-'){
	argcheck( arg, argc, argv);
	
	migr = atof(  argv[arg++] );
	//printf("\nmigration is %e\nnext argument is %d", migr, arg);
	global_values[arg-2] = migr;
      } 
      // default value of migr is 0;
      for( i=0; i<pars.cp.npop; i++) 
	for( j=0; j<pars.cp.npop; j++) 
	  pars.cp.mig_mat[i][j] = migr/(pars.cp.npop-1) ;

      for( i=0; i< pars.cp.npop; i++) pars.cp.mig_mat[i][i] = migr ;
      arg += c;
      break;
    case 'm' :
      c = 0;
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
      if( argv[arg][2] == 'a' ) {
	arg++;
	for( pop = 0; pop <npop; pop++)
	  for( pop2 = 0; pop2 <npop; pop2++){
	  
	    //fprintf(stderr, "arg: %d, argv %s\n", arg, argv[arg]);

	    /* if(pop == pop2){ */
/* 	      pars.cp.mig_mat[pop][pop2]= 0.; */
/* 	      arg++; */
/* 	      continue; */
/* 	    } */

	    c = 0;
	    //printf("%s, c:%d, pop:%d, pop2:%d, count:%d\n", argv[arg], c, pop, pop2, count);
	    //if( !count && (arg <argc) && ( argv[arg][0] != '-' ) ) {
	    if(!count && (c = dist_argcheck(arg, argc, argv) ) ){
	      //printf("%s, c:%d, pop:%d, pop2:%d\n", argv[arg], c, pop, pop2);
	      arg+=c;  // 
	      //fprintf(stderr, "WARNING: Ln3607 test 20090820 it was arg+=c; before. I think the arg+=c+1 is correct\n");
	    }
	    else if( count && param_distr[arg - 1] != '-'){
	      //printf("now I check migrate for pop %d, pop %d - dist: %c\n", pop, pop2, argv[arg][1]);
	      pars.cp.mig_mat[pop][pop2]  = get_param(arg - 1, argv[arg][1]);
	      global_values[arg-1] = pars.cp.mig_mat[pop][pop2];
	      arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
	    }
	    else if( arg < argc && argv[arg][0] != '-'){
	      argcheck( arg, argc, argv);
	      pars.cp.mig_mat[pop][pop2] = atof(  argv[arg++] );
	      global_values[arg-2] = pars.cp.mig_mat[pop][pop2];
	    }
	    else if( arg < argc && argv[arg][0] == '-'){
	      warning("-ma option", arg, argv[arg] );
	    }
	  }
	    
	
	for( pop = 0; pop < npop; pop++) {
	  pars.cp.mig_mat[pop][pop] = 0.0 ;
	  for( pop2 = 0; pop2 < npop; pop2++){
	    
	    if( pop2 != pop ) pars.cp.mig_mat[pop][pop] += pars.cp.mig_mat[pop][pop2] ;
	    //fprintf(stderr, "mig %d %d: %e\n", pop, pop2, pars.cp.mig_mat[pop][pop2]);
	  }
	}
	//fprintf(stderr, "WARNING: Ln3632 test 20090820 rm arg += c\n");
	//arg += c;
      }
      else {
	arg++;
	argcheck( arg, argc, argv);
	i = atoi( argv[arg++] ) -1;
	argcheck( arg, argc, argv);
	j = atoi( argv[arg++] ) -1;
	if(!count && (c = dist_argcheck(arg, argc, argv) ) ){}
	else if( count && param_distr[arg - 1] != '-'){
	  mij  = get_param(arg - 1, argv[arg][1]);
	  /*printf("mij: %e\n", mij);*/
	  pars.cp.mig_mat[i][i] += mij - pars.cp.mig_mat[i][j];
	  pars.cp.mig_mat[i][j] = mij;
	  global_values[arg-1] = mij;
	  arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
	}
	else if( arg < argc && argv[arg][0] != '-'){
	  argcheck( arg, argc, argv);
	  mij = atof( argv[arg++] );
	  //printf("%f\n", mij);
	  //EXIT_MSABC(1);
	  pars.cp.mig_mat[i][i] += mij - pars.cp.mig_mat[i][j];
	  pars.cp.mig_mat[i][j] = mij;
	  global_values[arg-2] = mij;
	}
	arg+=c;
      }
      
      break;
    case 'n' :
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
      arg++;
      c = 0;
      argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      if(!count && (c = dist_argcheck(arg, argc, argv)) ){}
      else if(count && param_distr[arg - 1] != '-'){
	psize = get_param(arg -1, argv[arg][1]);
	pars.cp.size[pop] = psize;
	global_values[arg-1] = psize;
	arg = arg + distribution_arguments( param_distr[arg-1]) +1;
      }
      else{
	argcheck( arg, argc, argv);
	psize = atof( argv[arg++] );
	pars.cp.size[pop] = psize ;
	global_values[arg-2] = psize;
      }
      arg += c;
      break;
      
    case 'g' :
      c = 0;
      if( npop < 2 ) { fprintf(stderr,"Must use -I option first.\n"); usage();}
      arg++;
      argcheck( arg, argc, argv);
      pop = atoi( argv[arg++] ) -1;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
      if( ! count && (c = dist_argcheck(arg, argc, argv)) ) {}
      else if(count && param_distr[arg - 1] != '-'){
	palpha = get_param(arg -1, argv[arg][1]);
	pars.cp.alphag[pop] = palpha;
	global_values[arg-1] = palpha;
	arg = arg + distribution_arguments(param_distr[arg - 1] ) + 1;
      }
      else{
	palpha = atof( argv[arg++] );
	pars.cp.alphag[pop] = palpha ;
	global_values[arg-2] = palpha;
      }
      arg += c;
      break;
    case 'G' :
      arg++;
      c=0;
      if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -G.\n"); usage(); }
      if( !count && (c = dist_argcheck(arg, argc, argv)) ){}
      else if(count && param_distr[arg - 1] != '-'){
	palpha = get_param(arg - 1, argv[arg][1]);
	for( i=0; i<pars.cp.npop; i++) 
	  pars.cp.alphag[i] = palpha ;
	global_values[arg-1] = palpha;
	arg = arg + distribution_arguments( param_distr[arg - 1] ) + 1;
      }
      else{
    	palpha = atof( argv[arg++] );
	for( i=0; i<pars.cp.npop; i++) 
	  pars.cp.alphag[i] = palpha ;
	global_values[arg-2] = palpha;
      }
      arg+=c;
      break;
    case 'e' :
      /* set the current time. This increases every time that an appropriate events takes place */
      
      
      pt = (struct devent *)malloc( sizeof( struct devent) ) ;
      pt->detype = argv[arg][2] ;
      
      
      ch3 = argv[arg][3] ;
      arg++;
      c = 0;
      //fprintf(stderr, "0\n");
      /* if duration mode */
      //fprintf(stderr, "arg: %d, duration mode: %i, duration_counter: %d, duration_events: %d\n", arg, durationmode, durationmode_counter, durationmode_events);
      
      if( (durationmode && durationmode_events == 0) ||
	  (durationmode && (durationmode_events > durationmode_counter++))
	 && (pt->detype == 'N' || pt->detype == 'G' || 
			  pt->detype == 'n' || pt->detype == 'g' ||
			  pt->detype == 'j' || pt->detype == 'M' ||
			  pt->detype == 'm'
			  )
	  ){
	//fprintf(stderr, "DUR-MODE\n");
	/**
	   here the time parameter does not refer to the beginning of the event
	   It refers to how long the event has lasted. That means that when the event begun depends on
	   when the previous events finished. Thus, the user must be careful to put the events in the right
	   backwards order 
	*/
	c=0;
	double duration = 0.;
	//fprintf(stderr, "BEFORE current-time: %e for the event %c, duration is %e\n", current_time,pt->detype, duration);
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	  carg = arg-1;
	}
	else if( count && param_distr[arg-1] !='-'){
	  duration = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = current_time + duration;
	  carg = arg ;
	  arg = arg + distribution_arguments( param_distr[arg -1] ) + 1; 
	  
	}
	else{
	  argcheck( arg, argc, argv);
	  duration = atof( argv[arg++] ) ;
	  global_values[arg-2] = current_time + duration;
	  carg = arg;
	}
	arg += c;
	
	current_time += duration;
	pt->time = current_time;
	//fprintf(stderr, "time: %e for the event %c, duration is %e\n", pt->time, pt->detype, duration);
	
      }
      else if( (durationmode == 0) || (durationmode  == 1 && durationmode_counter > durationmode_events ) ){
	durationmode = 0;
	//fprintf(stderr, "2\n");
	/* read the time */
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	  pt->time = 0.;
	  carg = arg;
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->time = current_time = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->time;
	  carg = arg;
	  //fprintf(stderr, "carg:%d\n", carg);
	  arg = arg + distribution_arguments( param_distr[arg -1] ) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  pt->time = current_time = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->time;
	  carg = arg;
	}
	arg += c;
      }
      
      //      printf("count %d, time: %e\n", count, pt->time); //$$$
      
      timeevents[event]=pt->time;
      event++;
      //argcheck( arg, argc, argv);
      //pt->time = atof( argv[arg++] ) ;
      pt->nextde = NULL ;
      //fprintf(stderr, "time: pt->time: %e\n", pt->time);
      /* sort the events */
      if( pars.cp.deventlist == NULL ){
	//fprintf(stderr, "OK...arg %d\n", arg);
	pars.cp.deventlist = pt ;
      }
      else if (  pt->time < pars.cp.deventlist->time ) { 
	//fprintf(stderr, "arg: %d\n", arg);
	ptemp = pars.cp.deventlist ;
	pars.cp.deventlist = pt ;
	pt->nextde = ptemp ;	
      }	
      else{
	//fprintf(stderr, "arg: %d\n", arg);
	addtoelist( pt, pars.cp.deventlist ) ;
      }
      // arg++ has been done a few lines above
      switch( pt->detype ) {
      case 'N' :
	
	//fprintf(stderr, "carg:%d\n", carg);
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	break;
      case 'G' :
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	break;
      case 'M' :
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	break;
      case 'n' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	//argcheck( arg, argc, argv);
	
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	//pt->paramv = atof( argv[arg++] ) ;
	
	break;
      case 'g' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
	
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	//pt->paramv = atof( argv[arg++] ) ;
	break;
      case 's' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;

	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	//argcheck( arg, argc, argv);
	//pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'm' :
	c = 0;
	if( ch3 == 'a' ) {
	  pt->detype = 'a' ;
	  argcheck( arg, argc, argv);
	  npop2 = atoi( argv[arg++] ) ;
	  pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
	  for( pop =0; pop <npop2; pop++){
	    (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
	    for( i=0; i<npop2; i++){
	      c = 0;
	      if( i == pop ) arg++;
	      else {
		
		if(!count && (c = dist_argcheck(arg, argc, argv) ) ){ 
		  arg+=c; 
		}
		else if( count && param_distr[arg - 1] != '-'){
		  (pt->mat)[pop][i] = get_param(arg - 1, argv[arg][1]);
		  global_values[arg-1] = (pt->mat)[pop][i];
		  arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
		}
		else if(arg < argc && argv[arg][0] != '-'){
		  
		  argcheck( arg, argc, argv); 
		  (pt->mat)[pop][i] = atof( argv[arg++] ) ;
		  global_values[arg - 2] = (pt->mat)[pop][i];
		}
		else if(arg<argc && argv[arg][0] == '-'){
		  warning("-ema", arg, argv[arg]);
		}
	      }
	    }
	  }
	  
	  for( pop = 0; pop < npop2; pop++) {
	    (pt->mat)[pop][pop] = 0.0 ;
	    for( pop2 = 0; pop2 < npop2; pop2++){
	      if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
	    }
	  }	
	}
	
	else {
	  argcheck( arg, argc, argv);
	  pt->popi = atoi( argv[arg++] ) -1 ;
	  argcheck( arg, argc, argv);
	  pt->popj = atoi( argv[arg++] ) -1 ;
	  
	  if(!count && (c = dist_argcheck(arg, argc, argv) ) ){
	  }
	  else if( count && param_distr[arg - 1] != '-'){
	    pt->paramv = get_param(arg - 1, argv[arg][1]);
	    global_values[arg-1] = pt->paramv;
	    arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
	  }
	  else if( arg < argc && argv[arg][0] != '-'){
	    argcheck( arg, argc, argv);
	    pt->paramv = atof( argv[arg++] ) ;
	    global_values[arg-2] = pt->paramv;
	    /* fprintf(stderr, "%d\t%e\n", arg-2, global_values[arg-2]); */
	    /* EXIT_MSABC(1); */
	  }
	  arg+=c;
	}
	
	break;
      case 'j' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	argcheck( arg, argc, argv);
	pt->popj = atoi( argv[arg++] ) -1 ;
	break;
      default: fprintf(stderr,"e event\n");  usage();
      }
      break;
    default: 
      fprintf(stderr," option %s (at %d) is not supported\n", argv[arg], arg); 
      EXIT_MSABC(1);
    }
  }
  
  /* if there is no observed data */
  if(!strcmp(obfile, "") && count > 0 && recgiven == 0 && effectiveN < 0 && fragmode){
    fprintf(stderr, "Please provide the recombination rate and the number of fragments that can recombine!\n\n");
    fprintf(stderr, "-- This is needed in order to rescale theta (and rho if != 0) for each of the locus.\n");
    fprintf(stderr, "-- If recombination rate = 0 and theta (in -t) refers to a locus of length 100, then please write -r 0 100\n");
    fprintf(stderr, "-- Consult the manual for more details or contact the authors (pavlidis@bio.lmu.de)\n");
    EXIT_MSABC(-1);
  }
  
  /* the next check should be done only when count > 0 since the values for the parameters are set now at every step */
  if( count && (pars.mp.theta == 0.0 && effectiveN < 0) && ( pars.mp.segsitesin == 0 ) && ( pars.mp.treeflag == 0 ) && (pars.mp.timeflag == 0) && !strcmp(obfile, "") ) {
    fprintf(stderr," either -s or -t or -T option must be used. \n");
    usage();
    EXIT_MSABC(1);
  }
  sum = 0 ;
  for( i=0; i< pars.cp.npop; i++) sum += (pars.cp.config)[i] ;
  if( sum != pars.cp.nsam ) {
    fprintf(stderr," sum sample sizes != nsam\n");
    usage();
    EXIT_MSABC(1);
  }
}

double Fst(char fst[10], double p1, double p2, double p3, double p4){
  if(!strcmp(fst, "hbk"))
    return Fst_HBK(p2,p1);
  else if(!strcmp(fst, "hsm"))
    return Fst_HSM(p3,p2);
  else if(!strcmp(fst, "slatkin"))
    return Fst_Slatkin(p3, p2);
}


int getobservations(FILE* obin, char*** data){
  int i=0, j=0;
  int limx = 10000, limy = 1000;
  char line[1000000]; 
  int segsites;
  int first = 0;
    
  /* read the next dataset */
  while( fscanf(obin, "%s", line) != EOF){
    //printf("%s\n", line);
    if(!strcmp(line, "segsites:")){
      first = 1;
      /* read the segsites */
      if(fscanf(obin, "%d", &segsites) == EOF){
	fprintf(stderr, "cannot read the number of segregating sites\n");
	EXIT_MSABC(1);
      }
      else{
	//fprintf(stderr, "%d, %d, %d\n", pars.cp.nsam, largestnsam, segsites+1);
	(*data) = rematrix(largestnsam, pars.cp.nsam, segsites+1, (*data) ); 
	// read the rest of the line
	fgets(line, sizeof(line), obin);
      }
      
      
      /* skip the next line */
      if(fgets(line, sizeof(line), obin) == NULL){
	fprintf(stderr, "cannot read the next line (after segsites)\n");
	EXIT_MSABC(1);
      }
      
      continue;
      
    } 
    
    
    if(first && !strcmp(line, "//")){
      //fprintf(stderr, "%d\n", segsites);
      return segsites;
    }
    
    
    if( first && (line[0] == '1' || line[0] == '0' || line[0] == MISSING)){
      assert(j < 1000000 && i < pars.cp.nsam);
      j=0;
      
      /* save the data */
      while(line[j] != '\0' && line[j] != '\n'){
	
	(*data)[i][j] = line[j];
	j++;
      }
      
      i++;
      
    }
    //printf("%s\n", line);
  }
  
  return segsites;
}



int getmissing(FILE* missingin, int*** misseq, double*** misstart, double*** misend, int** missingdata){
  
  int i = 0, j = 0, jj=0;
  int limx=1000, limy=100;
  int m1;
  double m2,m3;
  char filename[100];
  FILE* fragmissing=NULL;
  
  *misseq = malloc((unsigned)limy * sizeof(int*));
  *misstart = malloc((unsigned)limy * sizeof(double*));
  *misend = malloc((unsigned)limy * sizeof(double*));
  *missingdata = malloc((unsigned)limy * sizeof(int));
  
  for(i=0; i<limy; i++){
    (*misseq)[i] = malloc((unsigned)limx * sizeof(int));
    (*misstart)[i] = malloc((unsigned)limx * sizeof(double));
    (*misend)[i] = malloc((unsigned)limx * sizeof(double));
  }
  

  //missingdata=0;
  while(fscanf(missingin, "%s", filename)!=EOF){
    
    //fprintf(stderr, "file is: %s\n", filename);
    if( (fragmissing = fopen( filename, "r")) == NULL){
      fprintf(stderr, "\nERROR opening file %s\n", filename);
      EXIT_MSABC(1);
    }

    
  
    if(limy < j+1){
      //fprintf(stderr, "rematrix y direction\n");
      *misseq = rematrix_int(limy+100, limy, limx, *misseq);
      *misstart = rematrix_double(limy+100, limy, limx, *misstart);
      *misend = rematrix_double(limy+100, limy, limx, *misend);
      if( !(*missingdata = realloc( *missingdata, (unsigned)(limy+100)*sizeof(int))))
	perror("realloc error in missingdata\n");
      limy += 100;
    }

    i=0;
    (*missingdata)[j]=0;
    
    while(fscanf(fragmissing, "%d%lf%lf", &m1, &m2, &m3) != EOF){
      (*misseq)[j][i] = m1;
      (*misstart)[j][i] = m2;
      (*misend)[j][i] = m3;
      //fprintf(stderr, "mising data: %d\t%e\t%e\n", m1, m2, m3);
      ++i;
      (*missingdata)[j]+=1;

      if(limx < i+1){
	//fprintf(stderr, "realloc in missing data arrays -- x direction\n");
	//realloc memory for the arrays
	limx += 1000;
	for(jj=0; jj<limy; jj++){
	  if( !( (*misseq)[jj] = realloc((*misseq)[jj], (unsigned)limx * sizeof(int)) ) )
	    perror("realloc error in misseq\n");
	  if( !( (*misstart)[jj] = realloc((*misstart)[jj], (unsigned)limx * sizeof(double)) ) )
	    perror("realloc error in misstart\n");
	  if( !( (*misend)[jj] = realloc((*misend)[jj], (unsigned)limx * sizeof(double)) ) )
	    perror("realloc error in misend\n");
	}
      }
    }
    //fprintf(stderr, "fragment %d contains %d missing\n", j, (*missingdata)[j]);
    
    /* for(i=0; i<90; i++){ */
/*       fprintf(stderr, "** misseq[%d][%d]: %d\n", j, i, (*misseq)[j][i]); */
/*     } */
    j++;
  }
  missing_nof_fragments = limy;
  return 1;
}

int put_missing(char** list, double* posit, int nsam, int segsites, int** misseq, double** misstart, double** misend, int missingdata, int fragment){
  int i=0, j=0, s=0;

  //fprintf(stderr, "fragment missing: %d\n", fragment);
  
  /* for(i=0; i<9000; i++){ */
/*     fprintf(stderr, "misseq[%d]: %d\n", i, misseq[i]); */
/*   } */
  
//  fprintf(stderr, "missingdata: %d\tfragment %d\n", missingdata, fragment);

  for(i=0; i<missingdata; i++){ 
    //fprintf(stderr, "****** i: %d\tmissing: %d\tmisseq: %d\n", i, missingdata, misseq[fragment][i]);
    if( (misseq[fragment][i] > s)){
      j=0;      
      s = misseq[fragment][i];
    }
    while( (j < segsites) && (posit[j] < misend[fragment][i])){
      if(posit[j] > misstart[fragment][i])
	list[s][j] = MISSING;
      j++;
    }
  }
  return 1;
}
  
	  
int getoptions(FILE* optionsin, 
	       int* segs,
	       int* thetaw,
	       int* thetapi,
	       int* tajimasD,
	       int* ZnS,
	       int* Fst,
	       int* shared,
	       int* private,
	       int* fixed_dif,
	       int* fst_pops,
	       int* Fst_pops,
	       int* fwh,
	       int* dvstat,
	       int* thomson_est,
	       int* thomson_var,
	       // private
	       int* prisegs,
	       int* prithetaw,
	       int*prithetapi,
	       int* pritajimasD,
	       int* priZnS,
	       int* prifwh,
	       int* pridvstat,
	       int* prithomson_est,
	       int* prithomson_var,

	       int* pristats
	       

	       ){

  char opt[100];
  char line[1000];
  int n=0;
  int d=0, p=0, i;

  *segs = *thetaw = *thetapi = *tajimasD = *ZnS = *Fst = *shared = *private = *fixed_dif = *fst_pops 
    = *fwh = *prisegs = * prithetaw = *prithetapi = *priZnS = *prifwh = *pridvstat = *pristats = 0;
  //Fst_pops = NULL;

  
  
  while(fscanf(optionsin, "%s", opt) != EOF){// read until the end of file
    
    //printf("opt: %s\n", opt);
    if(!strcmp(opt, "tajd")){ 
      fscanf(optionsin, "%d", &d);
      if(d == 1) *tajimasD = 1;
      continue;
    }
    else if(!strcmp(opt, "segs")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *segs = 1;
    }
    else if(!strcmp(opt, "thetaw")){
      fscanf(optionsin, "%d", &d);
      if(d == 1) *thetaw = 1;
      }
    else if(!strcmp(opt, "thetapi")){
      fscanf(optionsin, "%d", &d);
      if(d == 1) *thetapi = 1;
    }
    else if(!strcmp(opt, "zns")){
      fscanf(optionsin, "%d", &d);
      if(d == 1) *ZnS = 1;
    }
    else if(!strcmp(opt, "shared")){
      fscanf(optionsin, "%d", &d);
      if(d == 1) *shared = 1;
    }
    else if(!strcmp(opt, "private")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *private = 1;
      }
      else if(!strcmp(opt, "fixed")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *fixed_dif = 1;
      }
      else if(!strcmp(opt, "fst")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *Fst = 1;
      }
      else if(!strcmp(opt, "pairfst")){
	fscanf(optionsin, "%d", &d);
	*fst_pops = d;
	//printf("d:%d\n", d);
	if(d > 1)
	  for(i=0; i<d; i++){
	    fscanf(optionsin, "%d", &p);
	    //printf("p:%d\n", p);
	    assert(p>0);
	    Fst_pops[i] = p-1;
	    //printf("p:%d\n", p);
	  }
      }
      else if(!strcmp(opt, "fwh")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *fwh = 1;
      }
      else if(!strcmp(opt, "dvstat")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *dvstat = 1;
      }
      else if(!strcmp(opt, "pri-stats")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats = *prisegs = *prithetaw = *prithetapi = *pritajimasD = *priZnS = *prifwh = 
		     *segs = *thetaw = *thetapi = *tajimasD = *ZnS = *fwh = 1;
      }
      else if(!strcmp(opt, "pri-segs")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats = *prisegs =  *segs = 1;
      }
      else if(!strcmp(opt, "pri-thetaw")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*prithetaw = *thetaw = 1;
      }
      else if(!strcmp(opt, "pri-thetapi")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*prithetapi = *thetapi = 1;
      }
      else if(!strcmp(opt, "pri-tajd")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*pritajimasD = *tajimasD = 1;
	
      }
      else if(!strcmp(opt, "pri-ZnS")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*priZnS = *ZnS = 1;
      }
      else if(!strcmp(opt, "pri-fwh")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*prifwh = *fwh = 1;
      }
      else if(!strcmp(opt, "pri-dvstat")){
	fscanf(optionsin, "%d", &d);
	if(d == 1) *pristats =*pridvstat = *dvstat = 1;
      }
    
  }
  
  return 1;
}


void getfragdimensions(FILE* ifile, int* rows, int* columns){
  int r=0, c=0, i=0;
  int maxcol=100;
  char head[200];
  char line[10000];
  
  for(i=0; i<maxcol; i++){
    
    if(fscanf(ifile, "%[\n]", head) || fscanf(ifile, "%*[ \t]%[\n]", head)){
      c = i;
      break;
    }
    if(!fscanf(ifile, "%s", head)){
      fprintf(stderr, "Could not read properly the heading in the fragment file\n");
      EXIT_MSABC(1);
    }
    
  }

  while(fgets(line, 9999, ifile)!=NULL){
    r++;
  }

  *rows = r;
  *columns = c;
  /*fprintf(stderr, "r: %d\t c: %d\n", r, c);*/
  
}
  
 
  
  
  


void getfragpars(FILE* ifile, int rows, int columns){
  //  printf("ok1\n"); // $$p
  int c=0;
  char head[200];
  char* name;
  name = (char*)malloc(100*sizeof(char));
  int k;
  double d;
  int j,i,ii;
  
  for(j=0; j<columns; j++){

    
    
    if(!fscanf(ifile, "%s", head)){
      fprintf(stderr, "Could not read properly the heading in the fragment file\n");
      EXIT_MSABC(1);
    }

      
    if(!strcmp(head, "id")) {
      c++;
      frpars[0] = j+1;
      continue;
    }
    if(!strcmp(head, "n")){
      c++;
      frpars[1] = j+1;
      continue;
    }
    if(!strcmp(head, "mu")){
      c++;
      frpars[2] = j+1;
      continue;
    }
    if(!strcmp(head, "rec")){
      c++;
      frpars[3] = j+1;
      continue;
    }
    if(!strcmp(head, "length")){
      c++;
      frpars[4] = j+1;
      continue;
    }
    if(!strcmp(head, "pop")){
      c++;
      frpars[5] = j+1;
      continue;
    }
    fprintf(stderr, "heading %s is not allowed\n", head);
    EXIT_MSABC(1);
  }

  if(!c){
    fprintf(stderr, "You should provide the header in the file with the loci information\n");
    fprintf(stderr, "allowed headings are: id, n, mu, rec, length, pop\n");
    EXIT_MSABC(1);
  }
  
      
     
  for(i=0; i<rows; i++){
    for(j=0; j<columns; j++){
      
      for(ii=0; ii<MAXFRPARS; ii++){
	//printf("col:%d, par:%d\n", j, ii);
	if(frpars[ii] == j+1){

	  switch(ii){
	  case 0:
	    fscanf(ifile, "%s", name); 
	    strcpy(fname[i], name);
	    //	    printf("fname function: %s\n", fname[i]); //$$p
	    break;
	  case 1:
	    fscanf(ifile, "%d", &k);
	    fsize[i] = k;
	    //	    printf("fsize function: %d\n", fsize[i]); //$$p
	  
	    break;
	  case 2: 
	    
	    fscanf(ifile,  "%lf", &d);
	    fmubp[i] = d;
	    
	    
	    break;
	  case 3:
	    
	    fscanf(ifile, "%lf", &d);
	    frecbp[i] = d;
	    
	    break;
	  case 4:
	    
	    fscanf(ifile, "%d", &k);
	    flength[i] = k;
	    //	    printf("flength function: %d\n", flength[i]); //$$p
	    break;
	  case 5:
	    
	    fscanf(ifile, "%d", &k);
	    fpop[i] = k;
	    //	    printf("fpop function: %d\n", fpop[i]); //$$p
	    break;
	  }
	  
	  break;
	}
      }
    }
    //fprintf(stderr, "#P, fragpars, theta: %e\n", fthetabp[i]);
  }
  //  printf("ok2\n"); // $$p
  free(name);
}

	void
argcheck( int arg, int argc, char *argv[] )
{
  //fprintf(stderr, "arg: %d\targv[arg]: %s\n", arg, argv[arg]);
	if( (arg >= argc ) || ( argv[arg][0] == '-') ) {
	   fprintf(stderr,"not enough arguments after %s\n", argv[arg-1] ) ;
	   fprintf(stderr,"For usage type: ms<return>\n");
	   EXIT_MSABC(0);
	  }
}


int  distribution_arguments(char distr){
  switch(distr){
  case 'U':
    return 2;
  case 'N':
    return 2;
  case 'G':
    return 2;
  case 'L':
    return 2;
  case 'R':
    return 2;
 default: 
   return 0;
  }
}

int get_distro_parms(int arg, int argc, char** argv, char distr){
  int cnt=-1;
  while( (++cnt + arg) < argc ){
  
    //fprintf(stderr, "argv[arg+cnt][0]: %c\n", argv[arg+cnt][0]); //#p
    //fprintf(stderr, "cnt: %d\targc: %d\targ+cnt: %d\n", cnt, argc, arg+cnt); //#p
    //cnt++;
    if(cnt == 0 && arg > 0){
      /* arg-1 is the index of the flag i.e. -t */
      param_distr[arg-1] = distr;
      continue;
    }
    else if( arg+cnt >= argc || cnt > distribution_arguments(distr) ){
      
      //fprintf(stderr, "cnt: %d\n", cnt); //#p
      break;
    }

    assert(cnt < MAXDISTPARS);
    /* the first value of the paramter (cnt=1) should go to position 0 and so on */
    /* the first value of the parameter is the next value of the arg given that arg denotes the distribution */
    /* arg-1 is the index of the flag i.e. -t */
    param_distr_values[arg-1][cnt-1] = atof(argv[arg+cnt]);
    
    
    //assert(arg++ < argc);
  }
  return  cnt;
}
  

/**
   PP check if the value should be retrieved from a distribution
*/
int dist_argcheck( int arg, int argc, char* argv[] ){
  //fprintf(stderr, "arg: %d\targv[arg]: %s\n", arg, argv[arg]);
  int is_distribution = 0;
  char distr = '-';
  int c = 0;

  if( (arg < argc) ){
    if(argv[arg][0] == '-'){
      int i=0;
      for( i=0; i<NDIST; i++){
	//fprintf(stderr, "argv[arg]: %s\n", argv[arg]);
	if(argv[arg][1] == distros[i]){
	  distr = argv[arg][1];
	  
	  break;
	}
      }

      //fprintf(stderr, "distribution is: %c\n", distr);
      if(distr == '-'){ // probably then double option with or without the parameter value. Then just continue;
	//argcheck( arg, argc, argv);
	return 0; // probably meaningless
      }
    }
    
    if(distr != '-'){
      /* set the values for the distribution in the global arrays param_dist and param_distr_values */
      c = get_distro_parms(arg, argc, argv, distr);
      //fprintf(stderr, "c is: %d\n", c);
    }
  }
  //printf("c is %d\n", c); //$$$
  return c;
}// dist_argchec
      
      
     

int warning(const char* message, int i, char* a){
  fprintf(stderr, "Option %s\nShouldn't you specify the value for the parameter? I read '-' but value is expected at argument %d, which is %s\n", message, i, a);
  EXIT_MSABC(-1);
  return 1;
}
      
	  

int
usage()
{
  fprintf(stderr, "########################   Hudson's ms commands  ########################\n");
  fprintf(stderr,"usage: ms nsam howmany \n");
  fprintf(stderr,"  Options: \n"); 
  fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
  fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
  fprintf(stderr,"\t -T          (Output gene tree.)\n");
  fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
  fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
  fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
  fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
  fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
  fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
  fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
  fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
  fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
  fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
  fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
  fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
  fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
  fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
  fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
  fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
  fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
  fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
  fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
  fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
  fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
  fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
  fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
  fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
  fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

  fprintf(stderr, "\n\n\n#########################   msABC COMMANDS  ######################################\n\n\n");
  	
  fprintf(stderr, "Options\n");
  fprintf(stderr, "\t --missing <file> : file that contains a list of files with missing data\n");
  fprintf(stderr, "\t --fst < hbk|slatkin|hsm > : the type of Fst\n");
  fprintf(stderr, "\t --resample <# events under order control> <event_1> ... <event_n>\n");
  fprintf(stderr, "\t --bm <back mutation rate> : needed for calculation of Fay and Wu H (0.0)\n");
  fprintf(stderr, "\t --options <file> : file with summary statistics to print\n");
  fprintf(stderr, "\t --verbose : prints the polymorphic data and the constant parameters\n");
  fprintf(stderr, "\t --filter : it discards singletons from the calculation of ZnS\n");
  fprintf(stderr, "\t --dur-mode : enables the duration mode. i.e. all the times are relative to the previous event\n");
  fprintf(stderr, "\t --pf (new 2010-07-12): specify the proportion of fragments where a non-NaN value is required in the pairwise comparisons, e.g. pairwise fst, shared, private polymorphism\n"); 
  fprintf(stderr, "\t --frag-begin : initiates the 'section' of command line that refers to the loci (within locus 'section')\n");
  fprintf(stderr, "\t --finp <file> : the file that describes the loci properties (within locus 'section')\n");
  //fprintf(stderr, "\t --w  I should study this option better
  fprintf(stderr, "\t --N <float> : reads the effective population size (within locus 'section')\n");
  fprintf(stderr, "\t --frag-end : ends the locus 'section'\n");
  fprintf(stderr, "\t --printall : prints the polymorphic data. Similar to --verbose\n\n");
  fprintf(stderr, "\nHELP OPTIONS\n");
  fprintf(stderr, "\t msABC -ms : to get only the help that is associated with ms\n");
  fprintf(stderr, "\t msABC -h : to get only the help that is associated with msABC\n");
  fprintf(stderr, "\t msABC -p : to get only the help that is associated with random distributions\n");
  fprintf(stderr, "\t msABC -all : to get the previous help outputs\n\n\n");
  fprintf(stderr, "\t ------- For more questions contact: pavlidis@bio.lmu.de\n\n");
  
  
  
  EXIT_MSABC(1);
}


int usage1(){
  fprintf(stderr,"usage: ms nsam howmany \n");
  fprintf(stderr,"  Options: \n"); 
  fprintf(stderr,"\t -t theta   (this option and/or the next must be used. Theta = 4*N0*u )\n");
  fprintf(stderr,"\t -s segsites   ( fixed number of segregating sites)\n");
  fprintf(stderr,"\t -T          (Output gene tree.)\n");
  fprintf(stderr,"\t -F minfreq     Output only sites with freq of minor allele >= minfreq.\n");
  fprintf(stderr,"\t -r rho nsites     (rho here is 4Nc)\n");
  fprintf(stderr,"\t\t -c f track_len   (f = ratio of conversion rate to rec rate. tracklen is mean length.) \n");
  fprintf(stderr,"\t\t\t if rho = 0.,  f = 4*N0*g, with g the gene conversion rate.\n"); 
  fprintf(stderr,"\t -G alpha  ( N(t) = N0*exp(-alpha*t) .  alpha = -log(Np/Nr)/t\n");      
  fprintf(stderr,"\t -I npop n1 n2 ... [mig_rate] (all elements of mig matrix set to mig_rate/(npop-1) \n");    
  fprintf(stderr,"\t\t -m i j m_ij    (i,j-th element of mig matrix set to m_ij.)\n"); 
  fprintf(stderr,"\t\t -ma m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n"); 
  fprintf(stderr,"\t\t -n i size_i   (popi has size set to size_i*N0 \n");
  fprintf(stderr,"\t\t -g i alpha_i  (If used must appear after -M option.)\n"); 
  fprintf(stderr,"\t   The following options modify parameters at the time 't' specified as the first argument:\n");
  fprintf(stderr,"\t -eG t alpha  (Modify growth rate of all pop's.)\n");     
  fprintf(stderr,"\t -eg t i alpha_i  (Modify growth rate of pop i.) \n");    
  fprintf(stderr,"\t -eM t mig_rate   (Modify the mig matrix so all elements are mig_rate/(npop-1)\n"); 
  fprintf(stderr,"\t -em t i j m_ij    (i,j-th element of mig matrix set to m_ij at time t )\n"); 
  fprintf(stderr,"\t -ema t npop  m_11 m_12 m_13 m_21 m_22 m_23 ...(Assign values to elements of migration matrix.)\n");  
  fprintf(stderr,"\t -eN t size  (Modify pop sizes. New sizes = size*N0 ) \n");    
  fprintf(stderr,"\t -en t i size_i  (Modify pop size of pop i.  New size of popi = size_i*N0 .)\n");
  fprintf(stderr,"\t -es t i proportion  (Split: pop i -> pop-i + pop-npop, npop increases by 1.\n");    
  fprintf(stderr,"\t\t proportion is probability that each lineage stays in pop-i. (p, 1-p are admixt. proport.\n");
  fprintf(stderr,"\t\t Size of pop npop is set to N0 and alpha = 0.0 , size and alpha of pop i are unchanged.\n");
  fprintf(stderr,"\t -ej t i j   ( Join lineages in pop i and pop j into pop j\n");
  fprintf(stderr,"\t\t  size, alpha and M are unchanged.\n");  
  fprintf(stderr,"\t  -f filename     ( Read command line arguments from file filename.)\n");  
  fprintf(stderr," See msdoc.pdf for explanation of these parameters.\n");

}


int usage2(){
  fprintf(stderr, "Options\n");
  fprintf(stderr, "\t --missing <file> : file that contains a list of files with missing data\n");
  fprintf(stderr, "\t --fst < hbk|slatkin|hsm > : the type of Fst\n");
  fprintf(stderr, "\t --resample <# events under order control> <event_1> ... <event_n>\n");
  fprintf(stderr, "\t --bm <back mutation rate> : needed for calculation of Fay and Wu H (0.0)\n");
  fprintf(stderr, "\t --options <file> : file with summary statistics to print\n");
  fprintf(stderr, "\t --verbose : prints the polymorphic data and the constant parameters\n");
  fprintf(stderr, "\t --filter : it discards singletons from the calculation of ZnS\n");
  fprintf(stderr, "\t --dur-mode : enables the duration mode. i.e. all the times are relative to the previous event\n");
  fprintf(stderr, "\t --pf (new 2010-07-12): specify the proportion of fragments where a non-NaN value is required in the pairwise comparisons, e.g. pairwise fst, shared, private polymorphism\n"); 
  fprintf(stderr, "\t --frag-begin : initiates the 'section' of command line that refers to the loci (within locus 'section')\n");
  fprintf(stderr, "\t --finp <file> : the file that describes the loci properties (within locus 'section')\n");
  //fprintf(stderr, "\t --w  I should study this option better
  fprintf(stderr, "\t --N <float> : reads the effective population size (within locus 'section')\n");
  fprintf(stderr, "\t --frag-end : ends the locus 'section'\n");
  fprintf(stderr, "\t --printall : prints the polymorphic data. Similar to --verbose\n");
  
}

int usage3(){
  fprintf(stderr, "Parameters should be:\n");
  fprintf(stderr, "-U <min> <max> : U(min, max)\n");
  fprintf(stderr, "-N <mean> <sigma^2> : Normal(mean, sigma^2)\n");
  fprintf(stderr, "-G <shape> <scale> : Gamma(shape, scale)\n");
  fprintf(stderr, "-L <ln(mean)> <ln(sigma^2)> : Log-Normal(mean, sigma^2)\n");
  fprintf(stderr, "-R <min> <max>  . It will return a log-uniform value between min and max. For example if you provide min=10 and max=1000, then it will be the same probability to return a value in (10,100) or (100, 1000). The program gets the min and max, calulates its logs (logmin and logmax), then it gets a uniform number r, between logmin, logmax; then it returns the EXP(r).");
}


	void
addtoelist( struct devent *pt, struct devent *elist ) 
{
	struct devent *plast, *pevent, *ptemp  ;

	pevent = elist ;
	while(  (pevent != NULL ) && ( pevent->time <= pt->time ) )  {
		plast = pevent ;
		pevent = pevent->nextde ;
		}
	ptemp = plast->nextde ;
	plast->nextde = pt ;
	pt->nextde = ptemp ;
}

	void 
free_eventlist( struct devent *pt, int npop )
{
   struct devent *next ;
   int pop ;
   
   while( pt != NULL){
	  next = pt->nextde ;
	  if( pt->detype == 'a' ) {
	     for( pop = 0; pop < npop; pop++) free( (pt->mat)[pop] );
		 free( pt->mat );
	  }
	  free(pt);
	  pt = next ;
   }
}

	
/************ make_gametes.c  *******************************************
*
*
*****************************************************************************/

#define STATE1 '1'
#define STATE2 '0'

	int
make_gametes(int nsam, int mfreq, struct node *ptree, double tt, int newsites, int ns, char **list )
{
	int  tip, j,  node ;
        int pickb(int nsam, struct node *ptree, double tt), 
            pickbmf(int nsam, int mfreq, struct node *ptree, double tt) ;

	for(  j=ns; j< ns+newsites ;  j++ ) {
		if( mfreq == 1 ) node = pickb(  nsam, ptree, tt);
		else node = pickbmf(  nsam, mfreq, ptree, tt);
		for( tip=0; tip < nsam ; tip++) {
		   if( tdesn(ptree, tip, node) ) list[tip][j] = STATE1 ;
		   else list[tip][j] = STATE2 ;
		   }
		}
}


/***  ttime.c : Returns the total time in the tree, *ptree, with nsam tips. **/

	double
ttime( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i;

	t = (ptree + 2*nsam-2) -> time ;
	for( i=nsam; i< 2*nsam-1 ; i++)
		t += (ptree + i)-> time ;
	return(t);
}


	double
ttimemf( ptree, nsam, mfreq)
	struct node *ptree;
	int nsam, mfreq;
{
	double t;
	int i;

	t = 0. ;
	for( i=0;  i< 2*nsam-2  ; i++)
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		t += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	return(t);
}


	void
prtree( ptree, nsam)
	struct node *ptree;
	int nsam;
{
	double t;
	int i, *descl, *descr ;
	void parens( struct node *ptree, int *descl, int *descr, int noden );

	descl = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	descr = (int *)malloc( (unsigned)(2*nsam-1)*sizeof( int) );
	for( i=0; i<2*nsam-1; i++) descl[i] = descr[i] = -1 ;
	for( i = 0; i< 2*nsam-2; i++){
	  if( descl[ (ptree+i)->abv ] == -1 ) descl[(ptree+i)->abv] = i ;
	  else descr[ (ptree+i)->abv] = i ;
	 }
	parens( ptree, descl, descr, 2*nsam-2);
	free( descl ) ;
	free( descr ) ;
}

	void
parens( struct node *ptree, int *descl, int *descr,  int noden)
{
	double time ;

   if( descl[noden] == -1 ) {
	printf("%d:%5.3lf", noden+1, (ptree+ ((ptree+noden)->abv))->time );
	}
   else{
	printf("(");
	parens( ptree, descl,descr, descl[noden] ) ;
	printf(",");
	parens(ptree, descl, descr, descr[noden] ) ;
	if( (ptree+noden)->abv == 0 ) printf(");\n"); 
	else {
	  time = (ptree + (ptree+noden)->abv )->time - (ptree+noden)->time ;
	  printf("):%5.3lf", time );
	  }
        }
}

/***  pickb : returns a random branch from the tree. The probability of picking
              a particular branch is proportional to its duration. tt is total
	      time in tree.   ****/

	int
pickb(nsam, ptree, tt)
	int nsam;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
		if( y >= x ) return( i ) ;
		}
	return( i );
}

	int
pickbmf(nsam, mfreq, ptree, tt )
	int nsam, mfreq;
	struct node *ptree;
	double tt;
{
	double x, y, ran1();
	int i;

	x = ran1()*tt;
	for( i=0, y=0; i < 2*nsam-2 ; i++) {
	  if( ( (ptree+i)->ndes >= mfreq )  && ( (ptree+i)->ndes <= nsam-mfreq) )
		y += (ptree + (ptree+i)->abv )->time - (ptree+i)->time ;
	  if( y >= x ) return( i ) ;
	}
	return( i );
}

/****  tdesn : returns 1 if tip is a descendant of node in *ptree, otherwise 0. **/

	int
tdesn(ptree, tip, node )
	struct node *ptree;
	int tip, node;
{
	int k;

	for( k= tip ; k < node ; k = (ptree+k)->abv ) ;
	if( k==node ) return(1);
	else return(0);
}


/* pick2()  */

	int
pick2(n,i,j)
	int n, *i, *j;
{
	double ran1();

	*i = n * ran1() ;
	while( ( *j = n * ran1() ) == *i )
		;
	return(0) ;
}

/**** ordran.c  ***/

	int
ordran(n,pbuf)
	int n;
	double pbuf[];
{
	ranvec(n,pbuf);
	order(n,pbuf);
	return;
}


	int
mnmial(n,nclass,p,rv)
	int n, nclass, rv[];
	double p[];
{
	double ran1();
	double x, s;
	int i, j;

	for(i=0; i<nclass; i++) rv[i]=0;
	for(i=0; i<n ; i++) {
	   x = ran1();
	   j=0;
	   s = p[0];
	   while( (x > s) && ( j<(nclass-1) ) )  s += p[++j];
	   rv[j]++;
	   }
	return(j);
}

        int
order(n,pbuf)
        int n;
        double pbuf[];
{
        int gap, i, j;
        double temp;

        for( gap= n/2; gap>0; gap /= 2)
           for( i=gap; i<n; i++)
                for( j=i-gap; j>=0 && pbuf[j]>pbuf[j+gap]; j -=gap) {
                   temp = pbuf[j];
                   pbuf[j] = pbuf[j+gap];
                   pbuf[j+gap] = temp;
                   }
}


	int
ranvec(n,pbuf)
	int n;
	double pbuf[];
{
	int i;
	double ran1();

	for(i=0; i<n; i++)
		pbuf[i] = ran1();

	return;
}



	int
poisso(u)
	double u;
{
	double  cump, ru, ran1(), p, gasdev() ;
	int i=1;

	if( u > 30. ) return( (int)(0.5 + gasdev(u,u)) );
	ru = ran1();
	p = exp(-u);
	if( ru < p) return(0);
	cump = p;
	
	while( ru > ( cump += (p *= u/i ) ) )
		i++;
	return(i);
}


/* a slight modification of crecipes version */
/* returns a normally distributed deviate with zero mean and unit ariance, using ran1(idum) as the source*/
/* of uniform deviates */

double gasdev(m,v)
	double m, v;
{
	static int iset=0;
	static float gset;
	float fac,r,v1,v2;
	double ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran1()-1.0;
			v2=2.0*ran1()-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset= v1*fac;
		iset=1;
		return( m + sqrt(v)*v2*fac);
	} else {
		iset=0;
		return( m + sqrt(v)*gset ) ;
	}
}


/** return a uniform distributed value between [a, b) */
double uniform_ab(double a , double b){
  double ran1();
  double r = ran1();
  int i=0;
  
  /* for(i=0; i<10; i++){ */
    
  //fprintf(stderr, "a random number is: %e\n", r); 
/*   } */
  assert(a < b);
  return ( (b - a)*r + a );
}

double loguniform(double a, double b){
  double ran1();
  assert(a<b);
  assert(a > 0 && b > 0);
  double r = ran1();
  double ta = log(a);
  double tb = log(b);
  double tr = ((tb - ta)*r + ta);
  return exp(tr);
}


/**  Implements the Polar form of the Box-Muller
     Transformation
			 
     (c) Copyright 1994, Everett F. Carter Jr.
     Permission is granted by the author to use
     this software for any application provided this
     copyright notice is preserved.
*/
double normal_box_muller(double m, double s){
  double ran1();
  double x1, x2, y1, w;
  static double y2;
  static int use_last=0;


  if( use_last ){
    y1 = y2;
    use_last = 0;
  }
  else{
    do{
      x1 = 2. * ran1() - 1.;
      x2 = 2. * ran1() - 1.;
      w = x1*x1 + x2*x2;
    } while (w >= 1.);
    
    w = sqrt( ( -2 * log( w ))/ w);
    y1 = x1 * w;
    y2 = x2 * w;
    use_last = 1;
  }
    return (m + y1 * s);
  
}


double lognormal(double logmean, double logsd){
  return exp(normal_box_muller(logmean, logsd));
}


double gds(double parameter, double scale);   /* generator for gamma distribution */
double nsc(void){
  return normal_box_muller(0.0, 1.0);
}

double uniform(){
  double ran1();
  return ran1();
}


/*------------------------------------------------------------------*/

/***************************************************************
 * Gamma Distribution - Rejection algorithm gs combined with   *
 *                      Acceptance complement method gd        *
 ***************************************************************
 * 
 * REFERENCES: -   J.H. Ahrens, U. Dieter 1974:                *
 *                 Computer Methods For Sampling From Gamma,   *
 *                 Beta, Poisson and Binomial Distributions.   *
 *                 Computing 12, pp. 223 - 246.                *
 *             -   J.H. Ahrnes, U. Dieter 1982: Generating     *
 *                 Gamma By A Modified Rejection Technique.    *
 ***************************************************************/
double gds(parameter, scale)
double  parameter;
double scale;
{
  double a=parameter;                 /* a...scale parameter */
  assert(a > 0);

static double aa = -1.0, aaa = -1.0,
       b=0., c=0., d=0., e=0., r=0., s=0., si=0., ss=0., q0=0.,
       q1 = 0.0416666664, q2 =  0.0208333723, q3 = 0.0079849875,
       q4 = 0.0015746717, q5 = -0.0003349403, q6 = 0.0003340332,
       q7 = 0.0006053049, q8 = -0.0004701849, q9 = 0.0001710320,
       a1 = 0.333333333,  a2 = -0.249999949,  a3 = 0.199999867,
       a4= -0.166677482,  a5 =  0.142873973,  a6= -0.124385581,
       a7 = 0.110368310,  a8 = -0.112750886,  a9 = 0.104089866,
       e1 = 1.000000000,  e2 =  0.499999994,  e3 = 0.166666848,
       e4 = 0.041664508,  e5 =  0.008345522,  e6 = 0.001353826,
       e7 = 0.000247453;

double gds,p,q,t,sign_u,u,v,w,x;

if (a < 1.0)
  {               /* CASE A: Acceptance rejection algorithm gs */
   b = 1.0 + 0.36788794412 * a;                      /* Step 1 */

   for(;;)
      {p = b *  uniform();
       if (p <= 1.0)
         {gds = exp(log(p) / a);      /* Step 2. Case gds <= 1 */
          if (log( uniform()) <= -gds) return(gds);
         }
       else
         {gds = - log ((b - p) / a);   /* Step 3. Case gds > 1 */
          if (log( uniform()) <= ((a - 1.0) * log(gds)))
             return(scale * gds);
         }
      }
  }
else
  {if (a != aa)  /* CASE B: Acceptance complement algorithm gd */
     {aa = a;                          /* Step 1. Preparations */
      ss = a - 0.5;
      s = sqrt(ss);
      d = 5.656854249 - 12.0 * s;
     }

   t = nsc();                  /* Step 2. Normal deviate */
                               
   x = s + 0.5 * t; 
   gds = x * x;
   if (t >= 0.0)
      return(scale * gds);                     /* Immediate acceptance */

   u =  uniform();          /* Step 3. Uniform random number */
   if (d * u <= t * t * t) return(scale * gds);  /* Squeeze acceptance */

   if (a != aaa)
     {aaa = a;                  /* Step 4. Set-up for hat case */
      r = 1.0 / a;
      q0 = ((((((((q9 * r + q8) * r + q7) * r + q6) *  r + q5) *
           r + q4) * r + q3) * r + q2) * r + q1) * r;
      if (a > 3.686)
        {if (a > 13.022)
           {b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
           }
         else
           {b = 1.654 + 0.0076 * ss;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
           }
        }

      else
        {b = 0.463 + s - 0.178 * ss;
         si = 1.235;
         c = 0.195 / s - 0.079 + 0.016 * s;
        }
      }

    if (x > 0.0)                   /* Step 5. Calculation of q */
      {v = t / (s + s);                             /* Step 6. */
       if (fabs(v) > 0.25)
         {q = q0 - s * t + 0.25 * t * t + (ss+ss) * log(1.0 + v);
         }
       else
         {q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) *
              v + a7) * v + a6) *  v + a5) * v + a4) *
              v + a3) * v + a2) * v + a1) * v;
         }                      /* Step 7. Quotient acceptance */
       if (log(1.0 - u) <= q)
          return(scale* gds);
      }

    for(;;)
      {                /* Step 8. Double exponential deviate t */
       do
         {e = -log( uniform());
          u =  uniform();
          u = u + u - 1.0;
          sign_u = (u > 0)? 1.0 : -1.0;
          t = b + (e * si) * sign_u;
         }                           /* Step 9. Rejection of t */
       while (t <= -0.71874483771719);

    v = t / (s + s);                      /* Step 10. New q(t) */
    if (fabs(v) > 0.25)
      {q = q0 - s * t + 0.25 * t * t + (ss + ss) * log(1.0 + v);
      }
    else
      {q = q0 + 0.5 * t * t * ((((((((a9 * v + a8) * v + a7) *
           v + a6) * v + a5) * v + a4) * v + a3) * v + a2) *
           v + a1) * v;
      }
    if (q <= 0.0) continue;                        /* Step 11. */
    if (q > 0.5)
      {w = exp(q) - 1.0;
      }
    else
      {w = ((((((e7 * q + e6) * q + e5) * q + e4) * q + e3) *
           q + e2) * q + e1) * q;
      }                             /* Step 12. Hat acceptance */
    if ( c * u * sign_u <= w * exp(e - 0.5 * t * t))
      {x = s + 0.5 * t;
       return(scale*x*x);
      }
    }
  }
}



double get_param(int arg, char distr){
  switch(distr){
  case 'U':
    //fprintf(stderr, "low bound: %e \t upper bound: %e\n",  param_distr_values[arg][0], param_distr_values[arg][1]);
    return uniform_ab( param_distr_values[arg][0], param_distr_values[arg][1] );
    break;
  case 'N': 
    return normal_box_muller( param_distr_values[arg][0], param_distr_values[arg][1] );
    break;
  case 'G':
    return gds( param_distr_values[arg][0], param_distr_values[arg][1] );
    break;
  case 'L':
    return lognormal( param_distr_values[arg][0], param_distr_values[arg][1] );
    break;
  case 'R':
    return loguniform( param_distr_values[arg][0], param_distr_values[arg][1] );
    break;
  default:
    fprintf(stderr, "Error, there is no such distribution: %c\n",  distr);
    EXIT_MSABC(1);
    
  }
}


int check_timeevents(int events, double* timeevents){
  int i;
  for(i=1; i<events; i++){
    //printf("event %d, time: %e -- event: %d, time: %e\n", i-1, timeevents[i-1], i, timeevents[i]);
    if(timeevents[i] <= timeevents[i-1])
      return 0;
  }
  //printf("Events are fine\n");
  return 1;
}

/* if something is wrong in the order return the line number starting from zero.
   If everything is fine, then return -1
*/
int check_conditional_timeevents(int nof_resampling, int** conditional_timeevents, double* timeevents){
  int i,j, ev1, ev2, nofev;
  
  for(i=0; i< nof_resampling; i++){
    nofev = conditional_timeevents[i][0];
    assert(nofev > 1);
    for(j=2; j<=nofev; j++){
      ev1 = conditional_timeevents[i][j-1];
      ev2 = conditional_timeevents[i][j];
      if( timeevents[ev2-1] <= timeevents[ev1-1])
	return i;
    }
  }
  return -1;
}


double varstat( double sumsq, double ave, double den){
  return (sumsq/den - ave*ave)/(den - 1)* den;
}

double skewstat(double sum3, double sum2, double ave, double den) {
    double m2 = sum2/den - ave*ave;
    double m3 = sum3/den - 3.0*ave*(sum2/den) + 2.0*ave*ave*ave;
    if (m2 <= 0.0) return 0.0;
    return m3 / pow(m2, 1.5);
}

double kurtstat(double sum4, double sum3, double sum2, double ave, double den) {
    double m2 = sum2/den - ave*ave;
    double m4 = sum4/den - 4.0*ave*(sum3/den) + 6.0*ave*ave*(sum2/den) - 3.0*ave*ave*ave*ave;
    if (m2 <= 0.0) return 0.0;
    return m4 / (m2*m2);
}

 
void get_timeevents(int argc, char* argv[]){
  
  int arg, carg, i=0, j=0, sum , pop , argstart, npop , npop2, pop2 ;
  double migr=0.;
  double mij, psize, palpha ;
  void addtoelist( struct devent *pt, struct devent *elist ); 
  void argcheck( int arg, int argc, char ** ) ;
  void free_eventlist( struct devent *pt, int npop );
  
  double get_param( int arg, char distr);
  int dist_argcheck(int arg, int  argc, char** argv);
  int distribution_arguments( char distr );
  
  struct devent *ptemp , *pt ;
  FILE *pf ;
  char ch3 ;
  
  npop = pars.cp.npop;

  // free the eventlist anyway, since this function is not used the first time
  free_eventlist( pars.cp.deventlist, npop);

  pars.cp.deventlist = NULL;
  arg=3;
  while(arg < argc){
    int c;
    
    if( argv[arg][0]== '-' && argv[arg][1] == 'e'){
      /* set the current time. This increases every time that an appropriate events takes place */
      
      
      pt = (struct devent *)malloc( sizeof( struct devent) ) ;
      pt->detype = argv[arg][2] ;
      
      
      ch3 = argv[arg][3] ;
      arg++;
      c = 0;
      
      if(durationmode){
	fprintf(stderr, "durationmode==1 && resample is not compatible\n");
	EXIT_MSABC(1);
      }
      
      /* read the time */
      if( !count && (c = dist_argcheck(arg, argc, argv))){
	pt->time = 0.;
	carg = arg;
      }
      else if( count && param_distr[arg-1] !='-'){
	pt->time = get_param(arg-1, argv[arg][1]);
	global_values[arg-1] = pt->time;
	carg = arg;
	//fprintf(stderr, "carg:%d\n", carg);
	arg = arg + distribution_arguments( param_distr[arg -1] ) + 1;
      }
      else{
	argcheck( arg, argc, argv);
	pt->time = atof( argv[arg++] ) ;
	global_values[arg-2] = pt->time;
	carg = arg;
      }
      arg += c;
      
      
      //      printf("count %d, time: %e\n", count, pt->time); //$$$
      
      timeevents[i]=pt->time;
      i++;
      //argcheck( arg, argc, argv);
      //pt->time = atof( argv[arg++] ) ;
      pt->nextde = NULL ;
      
      /* sort the events */
      if( pars.cp.deventlist == NULL ) 
	pars.cp.deventlist = pt ;
      else if ( pt->time < pars.cp.deventlist->time ) { 
	ptemp = pars.cp.deventlist ;
	pars.cp.deventlist = pt ;
	pt->nextde = ptemp ;	
      }	
      else
	addtoelist( pt, pars.cp.deventlist ) ;
      // arg++ has been done a few lines above
      switch( pt->detype ) {
      case 'N' :
	
	//fprintf(stderr, "carg:%d\n", carg);
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	break;
      case 'G' :
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eG.\n"); usage(); }
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	break;
      case 'M' :
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	break;
      case 'n' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	//argcheck( arg, argc, argv);
	
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	//pt->paramv = atof( argv[arg++] ) ;
	
	break;
      case 'g' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	if( arg >= argc ) { fprintf(stderr,"Not enough arg's after -eg.\n"); usage(); }
	
	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  //argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	//pt->paramv = atof( argv[arg++] ) ;
	break;
      case 's' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;

	/* read the value */
	c = 0;
	if( !count && (c = dist_argcheck(arg, argc, argv))){
	}
	else if( count && param_distr[arg-1] !='-'){
	  pt->paramv = get_param(arg-1, argv[arg][1]);
	  global_values[arg-1] = pt->paramv;
	  arg = arg + distribution_arguments( param_distr[arg -1]) + 1;
	}
	else{
	  argcheck( arg, argc, argv);
	  pt->paramv = atof( argv[arg++] ) ;
	  global_values[arg-2] = pt->paramv;
	}
	arg += c;
	
	//argcheck( arg, argc, argv);
	//pt->paramv = atof( argv[arg++] ) ;
	break;
      case 'm' :
	c = 0;
	if( ch3 == 'a' ) {
	  pt->detype = 'a' ;
	  argcheck( arg, argc, argv);
	  npop2 = atoi( argv[arg++] ) ;
	  pt->mat = (double **)malloc( (unsigned)npop2*sizeof( double *) ) ;
	  for( pop =0; pop <npop2; pop++){
	    (pt->mat)[pop] = (double *)malloc( (unsigned)npop2*sizeof( double) );
	    for( i=0; i<npop2; i++){
	      c = 0;
	      if( i == pop ) arg++;
	      else {
		
		if(!count && (c = dist_argcheck(arg, argc, argv) ) ){ 
		  arg+=c; 
		}
		else if( count && param_distr[arg - 1] != '-'){
		  (pt->mat)[pop][i] = get_param(arg - 1, argv[arg][1]);
		  global_values[arg-1] = (pt->mat)[pop][i];
		  arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
		}
		else if(arg < argc && argv[arg][0] != '-'){
		  
		  argcheck( arg, argc, argv); 
		  (pt->mat)[pop][i] = atof( argv[arg++] ) ;
		  global_values[arg - 2] = (pt->mat)[pop][i];
		}
		else if(arg<argc && argv[arg][0] == '-'){
		  warning("-ema", arg, argv[arg]);
		}
	      }
	    }
	  }
	  
	  for( pop = 0; pop < npop2; pop++) {
	    (pt->mat)[pop][pop] = 0.0 ;
	    for( pop2 = 0; pop2 < npop2; pop2++){
	      if( pop2 != pop ) (pt->mat)[pop][pop] += (pt->mat)[pop][pop2] ;
	    }
	  }	
	}
	else {
	  argcheck( arg, argc, argv);
	  pt->popi = atoi( argv[arg++] ) -1 ;
	  argcheck( arg, argc, argv);
	  pt->popj = atoi( argv[arg++] ) -1 ;
	  
	  if(!count && (c = dist_argcheck(arg, argc, argv) ) ){
	  }
	  else if( count && param_distr[arg - 1] != '-'){
	    pt->paramv = get_param(arg - 1, argv[arg][1]);
	    global_values[arg-1] = pt->paramv;
	    arg = arg + distribution_arguments( param_distr[arg-1] ) + 1;
	  }
	  else if( arg < argc && argv[arg][0] != '-'){
	    argcheck( arg, argc, argv);
	    pt->paramv = atof( argv[arg++] ) ;
	    global_values[arg-2] = pt->paramv;
	  }
	  arg+=c;
	}
	
	break;
      case 'j' :
	argcheck( arg, argc, argv);
	pt->popi = atoi( argv[arg++] ) -1 ;
	argcheck( arg, argc, argv);
	pt->popj = atoi( argv[arg++] ) -1 ;
	break;
      default: fprintf(stderr,"e event\n");  usage();
      }
      //break;
      /*  default:  */
      /*       fprintf(stderr," option default\n");  usage() ; */
    }
    else
      arg++;
  }

}


void print_eventlist(){
  struct devent *visit = pars.cp.deventlist;
  
  
   while(visit != NULL) {
   printf("time: %f\ttype: %c\n", visit->time, visit->detype);
   visit = visit->nextde;
 }
 printf("\n");
}
    

void checkarguments(char** argv, char* arg, int i, int tot){
  if( strcmp(arg, "--missing") &&
      strcmp(arg, "--fst") &&
      strcmp(arg, "--obs") &&
      strcmp(arg, "--poldata") && 
      strcmp(arg, "--resample") &&
      strcmp(arg, "--bm") &&
      strcmp(arg, "--options") &&
      strcmp(arg, "--verbose") &&
      strcmp(arg, "--filter") &&
      strcmp(arg, "--dur-mode") &&
      strcmp(arg, "--frag-begin") &&
      strcmp(arg, "--pf") &&
      strcmp(arg, "--DER") && 
      strcmp(arg, "--outgroup") && 
      strcmp(arg, "--phasemode") && 
      /* the following lines should be commented out
	 because these flags belong to the section --frag-begin, --frag-end
      */

      /* strcmp(arg, "--finp") && */
/*       strcmp(arg, "--nf") && */
/*       strcmp(arg, "--np") && */
/*       strcmp(arg, "--w") && */
/*       strcmp(arg, "--N") && */
/*       strcmp(arg, "--frag-end") && */
      strcmp(arg, "--printall") &&
      strcmp(arg, "-f") &&
      strcmp(arg, "-r") &&
      strcmp(arg, "-c") &&
      strcmp(arg, "-t") &&
      strcmp(arg, "-s") &&
      strcmp(arg, "-seeds") &&
      strcmp(arg, "-F") &&
      strcmp(arg, "-L") &&
      strcmp(arg, "-I") &&
      strcmp(arg, "-m") &&
      strcmp(arg, "-n") &&
      strcmp(arg, "-g") &&
      strcmp(arg, "-G") &&
      strcmp(arg, "-eN") &&
      strcmp(arg, "-eG") &&
      strcmp(arg, "-eM") &&
      strcmp(arg, "-en") &&
      strcmp(arg, "-eg") &&
      strcmp(arg, "-es") && 
      strcmp(arg, "-ema") &&
      strcmp(arg, "-ma") &&
      strcmp(arg, "-em") &&
      strcmp(arg, "-ej")){
    
    fprintf(stderr, "\n\n!!!!!!!!!!!!!!!!!!!! argument %d: %s out of %d is invalid !!!!!!!!!!!!!!!!!!!!\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n", i, arg, tot);
    
    int k;
    for(k=0; k<i; ++k){
      fprintf(stderr, "%s ", argv[k]);
    }
    fprintf(stderr, "              * * *\n");
    for(k=i; k<tot; ++k){
      fprintf(stderr, "%s ", argv[k]);
    }
    fprintf(stderr, "* * * * * * * *\n");
  
    fprintf(stderr, "run msABC without parameters to get the possible options\n");
        EXIT_MSABC(1);
  }

}
	     
