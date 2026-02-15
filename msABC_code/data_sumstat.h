int segsites(int* almap, int n, int segsites, int start, int end);
int denominators(int n, double* hn, double* sqhn, double* bn);
int thetaH(double *h, char** list, int n, int segs, double bm, int* almap, int start, int end, int* missing);
int hDenominator(double* hden, int n, int segs, double hn, double sqhn, double bn);
int theta_pi(double* theta, int* almap, int n, int segsites, int* missing);
void set_maxder(int* maxder, char** list, int n, int segsites, int start, int end);
int theta_w(double* theta, int* almap, int n, int segsites, double denom, int* missing); 
void set_almap(int* almap, char** list, int n, int segsites, int start, int end);
void set_ders(int** ders, char** list, int n, int segsites, int start, int end);
void set_missing(int* missing, char** list, int n, int segsites, int start, int end);

double Dnominator(int n, int segsites, double hn, double sqhn);
int  tajD(double* tajd, int segsites, int n, double thetaw, double thetap, double hn, double sqhn);
int htest( double* h,  int n, int segs, double thetaPi, double thetaH, double hn, double sqhn, double bn);
int ZnS(double* ld, char** list, int n, int segsites, int filter, int* almap, int** ders, int* missing);
double r2(char** list, int x1, int x2, int n, int* almap, int** ders, int* missing);
int  ZnA(double* ld, char** list, int n, int segsites, int filter);
int  FuLiD();
int VarPi();

int calculations(double* weights,
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
		 int derived); // the fraction of private polymorphisms

double Fst_HSM(double piD, 
	       double piS);

double Fst_Slatkin(double piD,
		   double piS);

double Fst_HBK(double piS,
	       double piT);

int dvstat( char** list, int n, int segs, int start, int end, double* dvk, double* dvh);
int mystrcmp(char* p1, char* p2, int length);

/* this is the Thomson estimator as it is described in Hudson et al 2007. 
   Notice that I have adapted that to missing data. Thus, for every segregating sites I divide by
   the size of the site (excluding missing data).
*/
int thomsonEst(double* thomson, int* almap, int n, int segsites, int* missing);
int thomsonVar(double* thomson, int* almap, int n, int segsites);
