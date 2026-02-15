struct devent {
  double time;
  int popi;
  int popj;
  double paramv;
  double **mat ;
  char detype ;
  struct devent *nextde;
} ;
struct c_params {
  int npop;
  int nsam;
  int *config;
  int *initconfig;
  double **mig_mat;
  //double mig; // the global migration rate
  double r;
  /** hold the global value of r */
  double glob_r;
  int nsites;
  double f;
  double track_len;
  double *size;
  double *alphag;
  
  struct devent *deventlist ;
} ;
struct m_params {
  double theta;
  /** hold the global value of theta */
  double glob_theta;
  int segsitesin;
  
  int treeflag;
  int timeflag;
  int mfreq;
} ;
struct params { 
  struct c_params cp;
  struct m_params mp;
  int commandlineseedflag ;
};


