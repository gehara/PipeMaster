#include "data_sumstat.h"
#include "assert.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define MISSING '2'

void set_ders(int** ders, char** list, int n, int segsites, int start, int end){
  assert(end <= n);
  int i,j,k;
  for(i=0; i<segsites; i++){
    
    k=0;
    for(j=start; j<end; j++){
      //printf("k:%d\n", k); //#p
      if(list[j][i] == '1'){
	ders[k][i] = j;
	k++;
      }
    }
  }
}
      
    


void set_maxder(int* maxder, char** list, int n, int segsites, int start, int end){
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }
  
  int i, j, max;
  for(i=0; i<segsites; i++){
    max=start;
    for(j=start; j<end; j++)
      if(list[j][i] == '1')  max=j;
    maxder[i] = max;
  }
}
    

/*set the allelic map */
void set_almap(int* almap, char** list, int n, int segsites, int start, int end){
  //printf("start: %d, end: %d, n: %d\n", start, end, n);
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }
  //printf("enn and n is %d and %d\n", end, n);
  
  int i,j;
  for(i=0; i<segsites; i++){
    int ones = 0;
    
    for( j=start; j<end; j++)
      ones += ( list[j][i] == '1' );
    
    
    almap[i] = ones;
    //fprintf(stderr, "i:%i\tsegsites:%i\talmap[]:%i\tstart:%d\tend:%d\n", i,segsites, almap[i], start, end);
  }
}


void set_missing(int* missing, char** list, int n, int segsites, int start, int end){
  if(end > n){
    fprintf(stderr, "end:%d, n:%d\n", end, n);
    assert(end <= n);
  }

  int i,j,mis=0;
  for(i=0; i<segsites; i++){
    mis=0;
    
    for( j=start; j<end; j++)
      mis += ( list[j][i] == MISSING );

    missing[i] = mis;
    //fprintf(stderr, "i:%i\tsegsites:%i\tmissing[]:%i\tstart:%d\tend:%d\n", i,segsites, missing[i], start, end);
  }
}


// number of segregating sites is based only on the '1's
int seg_sites(int* almap, int n, int segsites){
  int i;
  int s = 0;
  
  for(i=0; i<segsites; i++)
    if((almap[i] > 0) && (almap[i] < n))
      s++;
  return s;
}

/* this is the Thomson estimator as it is described in Hudson et al 2007. 
   Notice that I have adapted that to missing data. Thus, for every segregating sites I divide by
   the size of the site (excluding missing data).
*/
int thomsonEst(double* thomson, int* almap, int n, int segsites, int* missing){
  
  double te = 0.;
  int i,nsam;
  assert(n>1);
  
  for(i=0; i<segsites; ++i){
    nsam = n - missing[i];
    
    if(nsam <= 0) return 0;
    
    if(almap[i] == 0 || almap[i] == nsam) continue;

    te += ((double)almap[i]/(double)nsam);
  }
  *thomson = te;

  return 1;
}

/* this is the estimate of variance of Thomson estimator, 
   see Hudson 2007: The variance of coalescent times estimates
*/
int thomsonVar(double* thomson, int* almap, int n, int segsites){
  double tv = 0.;
  int i,ones;
  assert(n > 1);
  
  int* sfs = malloc((n-1)*sizeof(int));
  for(i=0; i<n-1; ++i)
    sfs[i] = 0;

  for(i=0; i<segsites; ++i){
    
    if(almap[i] == 0 ||almap[i] == n) continue;

    ones = almap[i];
    if(ones > 0 && ones < n)
      sfs[ones-1] ++;
  }
  
  for(i=1; i<n; ++i)
    tv += i*i*sfs[i-1];
  tv /= ((double)n*n);
  *thomson = tv;
  free(sfs);
  return 1;
}

/*
  calculates the theta pi of simulated data stored in list
*/
int theta_pi(double* theta, int* almap, int n, int segsites, int* missing){
  double pi = 0.;
  int i,j, ones,  nsam;
  
  assert(n > 1);
  double denom;

  //denom = n*(n - 1);
  for(i=0; i<segsites; i++){
    nsam = n - missing[i];
    
    denom = nsam*(nsam - 1.);
    
    if(denom <= 0) continue;
      
    double ssh = 0.; // sum of site homozygosity;
    ones = almap[i];
    
    /* in cases of non-polymorphic sites add 0 to pi and continue the loop */
    if( ones == nsam || ones == 0)
      continue;
    /* for the 1's class of alleles */
    ssh += ((double)ones * (double)(ones - 1))/denom;
    /* for the 0's class of alleles */
    ssh += ((double)(nsam - ones) * (double)( nsam -ones - 1.))/ denom;
    
    pi += (1.0 - ssh);
  }
  *theta = pi;
  return 1;
}

// these denominators get the maximum n for their calculation. 
int denominators(int n, double* hn, double* sqhn, double* bn){
  int i;
  *hn=0.;
  *sqhn = 0.;
  for(i=1; i<n; i++){
    *hn += 1./(double)i;
    *sqhn += 1./(double)(i*i);
  }
  *bn = *sqhn + 1./((double)n*n);
  return 1;
}


int theta_w(double* theta, int* almap, int n, int segsites, double denom, int* missing){
  double w = 0.;
  int i,j;
  int true_segsites = 0, nsam;
  
  for(i=0; i<segsites; i++){
    nsam = n - missing[i];
    if(nsam <= 0) continue;
    /* in cases of non-polymorphic sites add 0 to pi and continue the loop */
    if( almap[i] == nsam || almap[i] == 0)
      continue;
    
    /* increase the true segsites by one; */
    true_segsites++;
  }
  
  w = (true_segsites)/denom;
  *theta = w;
  //printf("*w:%e\n", w);
  return 1;
}

// this denominator works with the original sample size. i.e. does not know about missing data
double Dnominator(int n, int segsites, double hn, double sqhn){
  double epsilon = 1e-10;
  double a1 = 0., a2 = 0., b1 = 0., b2=0., c1=0., c2=0., e1=0., e2=0., denom=0.;
  int i;

  
  a1 = hn;
  a2 = sqhn;
  /* for(i=1; i<n; i++){ */
/*     a1 += 1./(double)i; */
/*     a2 += 1./pow( (double)i, 2.0 ); */
/*   } */
  b1 = ((double)n+1.)/(3.*((double)n-1.));
  b2 = (2.*((double)n*(double)n + (double)n + 3.))/(9.*(double)n*((double)n-1.));
  c1 = b1 - 1./a1;
  
  
  if(fabs(c1) < epsilon) c1 = 0.;
 
  c2 = b2 - (n+2)/(a1*n) + (a2/(a1*a1));
  if(fabs(c2) < epsilon) c2 = 0.;

  e1 = c1/a1;
  if(fabs(e1) < epsilon) e1 = 0.;

  e2 = c2/(a1*a1 + a2);
  if(fabs(e2) < epsilon) e2 = 0.;

    
  denom=(e1*segsites + e2*segsites*((double)segsites - 1.) );
  if(fabs(denom) < epsilon) denom = epsilon;

  if(denom <= 0){
    printf("sample size: %d\n", n);
    printf("segsites: %d\ne1: %e\ne2: %e\n", segsites, e1, e2);
    fprintf(stderr,
	    "a1: %e\na2: %e\nb1: %e\nb2: %e\nc1: %e\nc2: %e\ne1: %e\ne2: %e\n",
	    a1, a2, b1, b2, c1, c2, e1, e2);
  }
  assert(denom > 0);

  denom = sqrt(denom);
  return denom;
}

// missing data affect tajD trhough the numerator
int tajD(double* tajd, int segsites, int n, double thetaw, double thetap, double hn, double sqhn){
  if(n <= 0)
    return 0;
  
  double denom, td;
  denom = Dnominator(n, segsites, hn, sqhn);
  assert(denom > 0.);
  td = (thetap - thetaw)/denom;
  *tajd = td;
  return 1;
}



// missing data affect ZnS biasing it to higher values
int ZnS(double* ld, char** list, int n, int segsites, int filter, int* almap, int** ders, int* missing){
  
  double zns = 0.;
  int class = 0;
  int exc;
  
  //exc = (int*)malloc(segsites*sizeof(int));
  
  int i =0, j=0;
  int counter = 0;
  
  for(i=0; i<segsites-1; i++){
    //printf("i:%d\n", i);
    exc = ( (almap[i] == filter) || (almap[i] == n-filter) ) ? 1 : 0;
    
    if(exc) continue;
    for(j=i+1; j<segsites; j++){
      exc = ( (almap[j] == filter) || (almap[j] == n-filter) ) ? 1 : 0;

      if(!exc){
	//printf("*i,j: %d,%d - %d,%d\n", i,j,segsites -1, segsites);
	zns += r2(list, i, j, n, almap, ders, missing);
	//printf("i,j: %d,%d - %d,%d\n", i,j,segsites -1, segsites);
	counter ++;
      }
    }
  }
  
  
  
  if(counter == 0)
    return 0;
  
  zns /= counter;
  *ld = zns;
  //printf("zns %e\n", zns);
  //exit(-1);
  return 1;
}

int ZnA(double* ld, char** list, int n, int segsites, int filter){}

int FuLiD(){}

int  VarPi(){}


double r2(char** list, int x1, int x2, int n, int* almap, int** ders, int* missing){
  int i=0;
  int m1 = 0, der1=0;
  int m2 = 0, der2=0;
  int  m12 = 0;
  double  cor = 0.;
  
  double m1a, m2a, m12a;
  int use=x1;
  assert(n>0);
  
  double max_missing = (missing[x1] > missing[x2]) ? (double)missing[x1] : (double)missing[x2];
  double nsam = (double)n - max_missing;
   
  if(nsam <= 0) return 0.;
  
  //fprintf(stderr, "n: %d, x1: %d, x2: %d, nsam: %e, max_missing: %e\n", n, x1, x2, nsam, max_missing);

  m1 = almap[x1];
  m2 = almap[x2];
  // choose the min of the maxs
  int  min = m1;

  if(m2 < m1){
    min=m2;
    use = x2;
  }
    
  for(i=0; i<min; i++){ // notice <
          
    //
    m12 += ((list[ ders[i][use]  ][x1] == '1') && (list[ ders[i][use] ][x2] == '1'));
    /* if(m12 > 10 ){ */
/*       printf("m1:%d, m2:%d, m12:%d\n", m1, m2, m12); */
/*       printf("*: %d, %d, %d, %c, %c, %d, %d\n", i, ders[i][use], min, list[ ders[i][use]  ][x1], list[ ders[i][use] ][x2], x1, x2); */
/*     } */
  }
  
  if(max_missing > 0){
    m1 = m2 = 0;
    for(i=0; i<almap[x1]; i++)
      m1 += (list[ ders[i][x1] ][x1] == '1' && list[ ders[i][x1] ][x2] != MISSING);
    
    for(i=0; i<almap[x2]; i++)
      m2 += (list[ ders[i][x2] ][x2] == '1' && list[ ders[i][x2] ][x1] != MISSING);
    
  }
	
  
 
  //printf("max: %d\t%d\t%d\t%d\n", min, m1, m2, m12);
  if(m1 == 0 || m1 == nsam || m2 == nsam || m2 == 0){
    return 0.;
  }
  
  m1a = (double)m1/nsam;
  m2a = (double)m2/nsam;
  m12a = (double)m12/nsam;
  
  

  cor = (m12a - m1a*m2a)*(m12a - m1a*m2a)/(m1a*(1.-m1a)*m2a*(1.-m2a));
  //printf("cor: %f\n", cor);
  /* if(use == x1) */
/*     for(i=0; i<min; i++) */
/*       fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x1], list[ ders[i][use]  ][x2], ders[i][use], use); */
/*   else */
/*     for(i=0; i<min; i++) */
/*       fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x2], list[ ders[i][use]  ][x1], ders[i][use], use); */
 


  if(!(cor >=-0.0001 && cor <= 1.0001)){
    for(i=0; i<min; i++)
      fprintf(stderr, "%c,%c at %d from %d \n", list[ ders[i][use]  ][x1], list[ ders[i][use]  ][x2], ders[i][use], use);

    for(i=0; i<n; i++)
      fprintf(stderr, "%c,", list[i][x1]);
    fprintf(stderr, "\n");

    for(i=0; i<n; i++)
      fprintf(stderr, "%c,", list[i][x2]);
    fprintf(stderr, "\n");

    fprintf(stderr, "corr: %e\tm1: %d\tm2: %d\tm12: %d\n", cor, m1, m2, m12);
    exit(-1);
  }
  
  return cor;
}



    
/**
   Calculations of Fst
   It needs the weights for the populations because the total frequency is needed to be calculated
   The weights can be obtained by the initial command.   

   Code adapted from libsequence
*/
int calculations(double* weights, // larger pops have greater weight
		 char** list, // the polymorphic table
		 int* config, // the configuration of the sample
		 int npop,
		 int segsites,
		 int n,

		 double* piT,
		 double* piS,
		 double* piB,
		 double* piD,
		 
		 double** shared, // the fraction of share polymorphisms
		 double** private,
		 double** fixed_dif, 
		 int derived){ 

  

  //printf("Fst function\n");
  

  int i, 
    j, 
    k=0, 
    start=0, 
    end=0,
    ni, nj,
    segs_ij=0,
    start_i,
    end_i,
    start_j,
    end_j,
    zeros_i,
    zeros_j;

  double w_ii_sq =0.,
    denom=0.,
    pi=0.,
    pi_ij=0.,
    sum_wi_wj=0.,
    weighted_pi_ij = 0.,
    weighted_pi_ii = 0.;

  static int* samples_b = NULL;
  static int* samples_e = NULL;
  static int first=1;

  if(first ==1){
    samples_b = (int*)malloc(npop * sizeof(int));
    samples_e = (int*)malloc(npop * sizeof(int));
    first = 0;
  }
  

  
  for(i=0; i<npop; i++){
    for(j=0; j<npop; j++){
      shared[i][j] = 0.;
      private[i][j] = 0.;
      fixed_dif[i][j] = 0.;
    }
  }
  
  end=0;
  start=0;
  samples_b[0] = 0;
  samples_e[0] = config[0];
  for(i=1; i<npop; i++){
    samples_b[i] = samples_b[i-1]+config[i-1];
    samples_e[i] = samples_b[i]+config[i];
    //printf("config %d is %d\n", i, config[i]);
  }

  int* almap_i = (int*)malloc(segsites*sizeof(int));
  int* almap_j = (int*)malloc(segsites*sizeof(int));
  
  int* missing_i = (int*)malloc(segsites*sizeof(int));
  int* missing_j = (int*)malloc(segsites*sizeof(int));

  for(i=0; i<npop; i++){ // over pops
    /* printf("i %d, start %d, end %d\n", i, start, end); */

    
    start = samples_b[i];
    end = samples_e[i];
    /* printf("%d, %d, %d, %d, %e\n", i, segsites, start, end, weights[i]); */
    set_almap(almap_i, list, n, segsites, start, end);
    set_missing(missing_i, list, n, segsites, start, end);

    pi=0.;
    w_ii_sq += weights[i]*weights[i];
    for(j=0; j<segsites; j++){ // over pol sites
      
      // set the sample size
      ni=config[i] - missing_i[j];
      
      if(ni < 2) continue; /* pi for a population of a sample size 0 or 1 is 0 */
      
      double ssh=0.;
      denom = (double)ni * ((double)ni-1.);

      ssh+= (ni-almap_i[j]>0) ? (double)(ni-almap_i[j])*(double)(ni-almap_i[j] - 1)/denom : 0.;
      ssh+=(almap_i[j] > 0)?(double)almap_i[j] * (double)(almap_i[j] - 1)/denom : 0.;
      pi += (1. - ssh);
      
    }//segsites
    weighted_pi_ii += weights[i]*weights[i]*pi;
    //printf("*%d, %d, %e\n", i, segsites, weights[i]);
  }//pops
  
  //between population divergence
  for(i=0; i<npop-1; i++){
    start_i=samples_b[ i];
    end_i=samples_e[i];
    
    set_almap(almap_i, list, n, segsites, start_i, end_i);
    set_missing(missing_i, list, n, segsites, start_i, end_i);
    
    for (j = i+1; j<npop; j++){
      start_j=samples_b[j];
      end_j=samples_e[j];
      //printf("BETWEEN i %d, start %d, end %d\n", j, start_j, end_j);
      set_almap(almap_j, list, n, segsites, start_j, end_j);
      set_missing(missing_j, list, n, segsites, start_j, end_j);

      pi_ij=0.;
      sum_wi_wj += weights[i] * weights[j];
      segs_ij = 0;

      for(k=0; k<segsites; k++){ // over segsites
	//printf("\n***k=%d\n", k);
	ni=config[i] - missing_i[k];
	nj=config[j] - missing_j[k];
	zeros_i = ni - almap_i[k];
	zeros_j = nj - almap_j[k];
	
	if((ni <=0) || (nj<=0) )
	  continue;
	
	//fprintf(stderr, "ni: %d, nj: %d, zi: %d, zj: %d, oni: %d, onj: %d\n", ni, nj, zeros_i, zeros_j, almap_i[k], almap_j[k]);

	pi_ij += (double)zeros_i*(double)almap_j[k]/(double)ni/(double)nj;
	pi_ij += (double)almap_i[k]*(double)zeros_j/(double)ni/(double)nj;
	
	if( ( (zeros_j + zeros_i) < (ni + nj) ) &&
	    ( (zeros_j + zeros_i) > 0) ){
	  segs_ij++;
	}
	      
	//printf("site:%d\tzeros_i: %d, almap_i: %d,  zeros_j: %d, almap_j: %d\n", k, zeros_i, almap_i[k], zeros_j, almap_j[k]);
	if( (!zeros_i && zeros_j && almap_j[k]) ||
	    (!zeros_j && zeros_i && almap_i[k]) ||
	    (!almap_j[k] && almap_i[k] && zeros_i) ||
	    (!almap_i[k] && almap_j[k] && zeros_j) ){
	  private[i][j] = private[i][j] + 1.0;
	  //printf("private: %e\n", private[i][j]);
	}
	if( zeros_i && zeros_j && almap_i[k] && almap_j[k]){
	  shared[i][j] = shared[i][j] + 1.0;
	}
        
	if( ( 
	     (derived == 0) && 
	     ( 
	      (!zeros_i && !almap_j[k]) || 
	       (!zeros_j && !almap_i[k]) 
	       )
	      )
	    ||
	    ( 
	     (derived == 1) && 
	      (!zeros_i && !almap_j[k]) 
	      )
	    ||
	    ( 
	     (derived == 2) && 
	      (!zeros_j && !almap_i[k]) 
	      )
	    ){
	  fixed_dif[i][j] = fixed_dif[i][j] + 1.0;
	}
	
      }
      

      if(segs_ij > 0){
	shared[i][j] /= (double)segs_ij;
	fixed_dif[i][j] /= (double)segs_ij;
	private[i][j] /= (double)segs_ij;
      }
      else if(segs_ij == 0){
	shared[i][j] =	fixed_dif[i][j] = private[i][j]  = 0./0.;
      }
      weighted_pi_ij += weights[i]*weights[j]*pi_ij;
    }//pops_j
  }//pops_i

  *piT = weighted_pi_ii + 2.0*weighted_pi_ij;
  *piS = weighted_pi_ii / w_ii_sq;
  *piD = (*piT - *piS)/(2.0 * sum_wi_wj);
  *piB = weighted_pi_ij / sum_wi_wj;
  /* printf("\n"); */
/*   printf("piT: %e\n", *piT); */
/*   printf("piS: %e\n", *piS); */
/*   printf("piD: %e\n", *piD); */
/*   printf("piB: %e\n", *piB); */
/*   printf("sum_wi_wj: %e\n", sum_wi_wj); */
/*   printf("weighted_pi_ii: %e\n", weighted_pi_ii); */
/*   printf("weighted_pi_ij: %e\n", weighted_pi_ij); */


  free(almap_i);
  free(almap_j);
  free(missing_i);
  free(missing_j);

  return 1;
}
	

double Fst_HSM(double piD, 
	       double piS){
  return piD/(piS + piD);
}

double Fst_Slatkin(double piD,
		   double piS){
  return piD/(2.0*piS + piD);
}

double Fst_HBK(double piS,
	       double piT){
  //printf("piS:%e\tpiT:%e\n", piS, piT);
  return 1. - (piS/piT);
}

 


/**
   Calculations of pairwise Fst Values
   It needs the weights for the populations because the total frequency is needed to be calculated
   The weights can be obtained by the initial command.

   Code adapted from libsequence
*/
int pairwiseFstcalculations(int popi, int popj, // the two populations
			    double* weights, // larger pops have greater weight
			    char** list, // the polymorphic table
			    int* config, // the configuration of the sample
			    int npop,
			    int segsites,
			    int n,
			    
			    double* piT,
			    double* piS,
			    double* piB,
			    double* piD
			    ){ 

  

  //
  

  int i,
    j,
    k,
    start=0,
    end=0,
    ni, nj,
    segs_ij=0,
    start_i,
    end_i,
    start_j,
    end_j,
    zeros_i,
    zeros_j;

  double w_ii_sq =0.,
    denom=0.,
    pi=0.,
    pi_ij=0.,
    sum_wi_wj=0.,
    weighted_pi_ij = 0.,
    weighted_pi_ii = 0.;

  
  static int* samples_b = NULL;
  static int* samples_e = NULL;
  static int first = 1;
  
  int total_npop = npop; // save the total number of populations into npop;
  npop = 2;
  int popind[2]; // save the indices of populations
  popind[0] = popi;
  popind[1] = popj;

  end=0;
  start=0;

  if( first == 1){
    samples_b = (int*)malloc(total_npop * sizeof(int));
    samples_e = (int*)malloc(total_npop * sizeof(int));
    first = 0;
  }
  
  samples_b[0] = 0;
  samples_e[0] = samples_b[0] + config[0];
  for(i=1; i<total_npop; i++){
    samples_b[i] = samples_b[i-1] + config[i-1];
    samples_e[i] = samples_b[i] + config[i];
    //printf("config %d is %d\n", i, config[i]);
  }
  
  // printf("segsites: %d\n", segsites);

  int* almap_i = (int*)calloc(segsites,sizeof(int));
  int* almap_j = (int*)calloc(segsites,sizeof(int));

  int* missing_i = (int*)calloc(segsites, sizeof(int));
  int* missing_j = (int*)calloc(segsites, sizeof(int));
  
 

 
  //printf("Fst function\n");
  for(i=0; i<npop; i++){ // over pops
    

    start = samples_b[ popind[i] ];
    end = samples_e[ popind[i] ];
    /* printf("*i %d, start %d, end %d\n", popind[i], start, end); */
/*     printf("*%d, %d, %d, %d, %e\n", i, segsites, start, end, weights[i]); */
    set_almap(almap_i, list, n, segsites, start, end); // set the allelic map for the population
    set_missing(missing_i, list, n, segsites, start, end);
    pi=0.;
    w_ii_sq += ((double)total_npop/2.0 * weights[ popind[i] ] )* ((double)total_npop/2.0 *weights[ popind[i] ]);
    //w_ii_sq += .5 * .5;
    //fprintf(stderr, "weight %e, total: %d, weights: %e, %e\n", w_ii_sq, total_npop, weights[ popind[i] ], weights[ popind[i] ]);
    for(j=0; j<segsites; j++){ // over pol sites
      ni=config[ popind[i] ] - missing_i[j]; // configuration of the population
      if(ni < 2) continue;
      double ssh=0.;
      denom = (double)ni * ((double)ni-1.);


      ssh+= (ni-almap_i[j]>0) ?
	(double)(ni-almap_i[j])*(double)(ni-almap_i[j] - 1)/denom : 0.;
      ssh+=(almap_i[j] > 0)?(double)almap_i[j] * (double)(almap_i[j] - 1)/denom : 0.;
      pi += (1. - ssh);
      
    }//segsites
    weighted_pi_ii += weights[ popind[i] ]*weights[ popind[i] ]*pi;
    //printf("*%d, %d, %e\n", i, segsites, weights[i]);
  }//pops
  
  //between population divergence
  for(i=0; i<npop-1; i++){
    start_i=samples_b[ popind[i] ];
    end_i=samples_e[ popind[i] ];
    
    set_almap(almap_i, list, n, segsites, start_i, end_i);
    set_missing(missing_i, list, n, segsites, start_i, end_i);

    for (j = i+1; j<npop; j++){
      start_j=samples_b[ popind[j] ];
      end_j=samples_e[ popind[j] ];
      //printf("BETWEEN i %d, start %d, end %d\n", j, start_j, end_j);
      set_almap(almap_j, list, n, segsites, start_j, end_j);
      set_missing(missing_j, list, n, segsites, start_j, end_j);

      pi_ij=0.;
      sum_wi_wj += weights[ popind[i] ] * weights[popind[j] ];
      segs_ij = 0;

      for(k=0; k<segsites; k++){ // over segsites
	ni=config[ popind[i] ] - missing_i[k];
	nj=config[ popind[j] ] - missing_j[k];
	zeros_i = ni - almap_i[k];
	zeros_j = nj - almap_j[k];
	if((ni<=0) || (nj<=0) ) continue;

	pi_ij += (double)zeros_i*(double)almap_j[k]/(double)ni/(double)nj;
	pi_ij += (double)zeros_j*(double)almap_i[k]/(double)ni/(double)nj;
      }// sites
      weighted_pi_ij += weights[ popind[i] ] * weights[ popind[j] ] * pi_ij;
    }//pops_j
  }//pops_i

  *piT = weighted_pi_ii + 2.0*weighted_pi_ij;
  *piS = weighted_pi_ii / w_ii_sq;
  *piB = weighted_pi_ij / sum_wi_wj;
  *piD = (*piT - *piS)/(2.0 * sum_wi_wj);
  /* printf("*\n"); */
   /* printf("*piT: %e\n", *piT);  */
/*    printf("*piS: %e\n", *piS);  */
/*    printf("w_ii_sq: %e\n", w_ii_sq); */
/*   printf("*piD: %e\n", *piD); */
/*   printf("*piB: %e\n", *piB); */
/*   printf("*sum_wi_wj: %e\n", sum_wi_wj); */
/*   printf("*weighted_pi_ii: %e\n", weighted_pi_ii); */
/*   printf("*weighted_pi_ij: %e\n", weighted_pi_ij); */

  free(almap_i);
  free(almap_j);
  free(missing_i);
  free(missing_j);
  return 1;
}

	   
	
int thetaH(double *h, char** list, int n, int segs, double bm, int* almap, int start, int end, int* missing){
  int i, j, k;
  double f = 0.;
  int noc=0;

 
  //printf("******************\n");
        
  assert(start >=0 && end <= n);
  int nsam = end - start ; // start is zero and end the sample size


  if(bm <= 0.) // if no back mutation use the allelic map.
    for(i=0; i<segs; i++){
      nsam = end - start - missing[i];
      if(nsam < 2) continue;
      if(almap[i] == 0 || almap[i] == nsam) 
	continue;
      noc = almap[i];
      f += (2.*(double)noc*(double)noc)/((double)nsam*((double)nsam-1.));
      
      //fprintf(stderr, "noc: %d, nsam: %d, end : %d, start: %d\n", 
      //      noc, nsam, end, start);
      
      
    }
  /* if there is the possibility for back mutation then the allelic map
     is not exactly the one we observe.
     In this case the number of mis-inferred states is binomially distributed.
  */
  else if(bm > 0) {
    for(i=0; i<segs; i++){
      nsam = end - start - missing[i];
      noc = 0;
      for(j=start; j<end; j++){
	if(rand() > bm && (list[j][i] != MISSING) )
	  noc += (list[j][i] - '0');
	else if(list[j][i] != MISSING) 
	  noc += ('1' - list[j][i]);
      }
      if(noc == 0 || noc == nsam) continue;
      f += (2.*(double)noc*(double)noc)/((double)nsam*((double)nsam-1.));

      
    }
  }

  *h = f;
  return 1;
}


int hDenominator(double* hden, int n, int segs, double hn, double sqhn, double bn){
  double b1, b2, c1, c2, e1, e2, w;
  double den;
  b1 = (n+1.)/(3. * (n-1));
  b2 = (double)2*(n*n + n+3.)/(double)(9. * n * (n-1));
  c1 = b1 - 1./hn;
  c2 = b2 - (n + 2.)/(1. * n * hn) + sqhn/(hn*hn);
  e1 = c1/hn;
  e2 = c2/(hn*hn + sqhn);
  w = ((double) segs)/hn;
  
  *hden = sqrt((n - 2.)*w/(6*(n-1.)) + (18 * n*n*(n*3 + 2)*bn - (88*n*n*n + 9.*n*n - 13.*n+6)) * w*w/(9.*n*(n-1) * (n-1)));
  return 1;
}
  

int htest( double* h, int n, int segs, double thetaPi, double thetaH, double hn, double sqhn, double bn){
  double hden = 0;
  double fwh=0.;
    
  hDenominator(&hden, n, segs, hn, sqhn, bn);
  
  //fprintf(stderr, "thetapi: %e\tthetaH: %e\thden: %e\n", thetaPi, thetaH, hden);
  
  fwh = 0.5*(thetaPi - thetaH)/hden;
  
  *h = fwh;
  return 1;

}


int dvstat( char** list, int n, int segs, int start, int end, double* dvk, double* dvh){
  assert(end <= n);
  assert(start >= 0);
  int i,j, nsam=end-start, k=end-start;
  int* haplo = malloc((unsigned)n * sizeof(int));
  
  *dvk = k;
  *dvh = 1.0;

  for(i=0; i<n; i++)
    haplo[i] = 1;
  
  for(i=start; i<end-1; i++){
    if(!haplo[i]) continue;
    for(j=i+1; j<end; j++){
      if(!haplo[j]) continue;
      if(!mystrcmp(list[i], list[j], segs)){
	haplo[j] = 0;
	haplo[i]++;
	continue;
      } 
    }
  }
  
  
  for(i=start; i<end; i++){
    //printf("\n--i: %d, haplo: %d, nsam: %d, dvk: %e\n", i, haplo[i], nsam, *dvk);
    if(!haplo[i]){
      *dvk = *dvk - 1.;
      continue;
    }
    *dvh -= pow((double)haplo[i]/(double)nsam, 2);
    
  }
   
  free(haplo);
  return 1;
}

int mystrcmp(char* p1, char* p2, int length) 
{ 
  // A variable to hold the distance in ASCII charcters 
  // between the current two characters we're comparing. 
  int dist = 0, i=0;; 
  
  // Keep checking while the distance is 0, and we're 
  // not hitting NULL on one of the strings. 
  while (!dist && i<length){
    // Get the distance while incrementing the 
    // pointers to the next character. 
    dist = (*p2++) - (*p1++); 
    i++;
  }
  // Check the last distance and according to this 
  // return (1) if the first string is bigger, (-1) 
  // if the second string is bigger, or (0) if the 
  // strings are identical. 
  if (dist > 0) 
    return (-1); 
    else if (dist < 0) 
      return (1); 
  
  return (0); 
}
	
      
    
      
