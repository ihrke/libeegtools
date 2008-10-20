/* helper.c
 *
 */	
#include "helper.h"

/* ---------------------------------------------------------------------------- 
   -- Helper functions                                                       -- 
   ---------------------------------------------------------------------------- */

/** deep copy of double ptrptr
 */
double** copy_double_ptrptr(const double **s, int N, int n){
  double **r;
  int i, j;

  r = (double**)malloc(N*sizeof(double*));
  for(i=0; i<N; i++){
    r[i] = (double*)malloc(n*sizeof(double));
	 for(j=0; j<n; j++){
		r[i][j] = s[i][j];
	 }
  }
  return r;
}

/** Deep-copy everything from source to dest. It is necessary, that
	 all number-of-element fields of dest are the same as those in
	 source, since no memory is reallocated. (comparison with
	 eegdata_cmp_settings() )
	 
	 \params dest,source
	 \return 0 success, failure else
 */
int      copy_similar_eegdata( EEGdata *dest, const EEGdata *source ){
  int i;
  if( eegdata_cmp_settings( dest, source ) ){
	 return -1;
  }
  memcpy( dest->markers, source->markers, (dest->nmarkers)*sizeof(unsigned long) );
  for( i=0; i<dest->nbchan; i++ ){
	 memcpy( dest->d[i], source->d[i], (dest->n)*sizeof( double ) );
  }
  return 0;
}

/** shallow copy of ModelData struct
 */
void    copy_modeldata(const ModelData *m1, ModelData *m2){
	m2->den_params=m1->den_params;
	m2->tw_params=m1->tw_params;
	m2->N=m1->N;
	m2->n=m1->n;
	m2->R=m1->R;
	m2->Ri=m1->Ri;
	m2->sRi=m1->sRi;
	m2->times=m1->times;
	m2->si=m1->si;
	m2->u=m1->u;
	m2->ui=m1->ui;
	m2->additional=m1->additional;
}

void free_warppath(WarpPath *p){
	free(p->upath);
	free(p->spath);
	free(p);
}

void print_modeldata(FILE *out, const ModelData *m){
	fprintf(out, "ModelData struct:\n");
	fprintf(out, "  den_params=\n");
	print_denoisingparameters(out, m->den_params);
	fprintf(out, "  tw_params=\n");
	print_timewarpparameters(out, m->tw_params);
	fprintf(out, "  N=%i\n", m->N);
	fprintf(out, "  n=%i\n", m->n);
	fprintf(out, "  R=%i\n", m->R);
	fprintf(out, "  Ri=%p\n", m->Ri);
	fprintf(out, "  sRi=%p\n", m->sRi);
	fprintf(out, "  times=%p\n", m->times);
	fprintf(out, "  si=%p\n", m->si);
	fprintf(out, "  u=%p\n", m->u);
	fprintf(out, "  ui=%p\n", m->ui);
	fprintf(out, "  additional=%p\n", m->additional);
}
void    print_denoisingparameters(FILE *out, const DenoisingParameters *p){
	fprintf(out, "   DenoisingParameters struct:\n"	);
	fprintf(out, "      L=%i\n", p->L);
	fprintf(out, "      sigextfct=%p\n", p->sigextfct);
	fprintf(out, "      cleanfct=%p\n", p->cleanfct);
	fprintf(out, "      eta=%p\n", p->eta);
}
void    print_timewarpparameters(FILE *out, const TimewarpParameters *p){
	fprintf(out, "   TimewarpParameters struct:\n"	);
	fprintf(out, "      theta1=%i\n", p->theta1);
	fprintf(out, "      theta1=%i\n", p->theta2);
}

		
void free_modeldata(ModelData *m){
	int i;
	if(!(m->den_params==NULL))	free(m->den_params);
	if(!(m->tw_params==NULL))	free(m->tw_params);
	if(!(m->additional==NULL))	free(m->additional);
	if(!(m->times==NULL)) free(m->times);
	if(!(m->sRi==NULL)) free(m->sRi);
	if(!(m->Ri==NULL)) free(m->Ri);
	if(!(m->u==NULL)) free(m->u);
	for(i=0; i<m->N; i++){
		if(m->si!=NULL)
			if(!(m->si[i]==NULL)) free(m->si[i]);
		if(m->ui!=NULL)
			if(!(m->ui[i]==NULL)) free(m->ui[i]);
	}
	if(m->si!=NULL) free(m->si);
	if(m->ui!=NULL) free(m->ui);
	free(m);
}

WarpPath* init_warppath(WarpPath *path, int J, int K){
  if(path==NULL){
	 path = (WarpPath*)malloc(sizeof(WarpPath));
  }
  path->J=J;
  path->K=K;
  path->upath = (int*)calloc(J+K, sizeof(int));
  path->spath = (int*)calloc(J+K, sizeof(int));

  return path;
}

EEGdata* init_eegdata(int nbchan, int nsamples, int nmarkers){
  int i;
  EEGdata *eeg;

  eeg = (EEGdata*) malloc( sizeof(EEGdata) );
  eeg->nbchan = nbchan;
  eeg->n = nsamples;
  eeg->d = (double**) malloc( nbchan * sizeof(double*) );
  for(i=0; i<nbchan; i++)
	 eeg->d[i] = (double*) malloc( nsamples*sizeof(double) );
  eeg->markers = (unsigned long*)malloc( nmarkers*sizeof(unsigned long) );
  eeg->nmarkers=nmarkers;

  return eeg;
}

EEGdata_trials* init_eegdata_trials(int nbtrials, int nmarkers_per_trial, int nbchan, int nbsamples){
  EEGdata_trials *eeg;
  int i;

  eeg = (EEGdata_trials*)malloc( sizeof( EEGdata_trials ) );
  eeg->nmarkers_per_trial = nmarkers_per_trial;
  eeg->ntrials = nbtrials;
  eeg->markers = (unsigned long **)malloc( nbtrials * sizeof( unsigned long* ) );
  for(i=0; i<nbtrials; i++){
	 eeg->markers[i] = (unsigned long*) malloc( nmarkers_per_trial * sizeof( unsigned long ) );
  }
  eeg->times = (double*)malloc( nbsamples*sizeof( double ) );
  eeg->data = (EEGdata**)malloc( nbtrials*sizeof( EEGdata* ) );
  for(i=0; i<nbtrials; i++)
	 eeg->data[i] = init_eegdata(nbchan, nbsamples, nmarkers_per_trial);
  
  return eeg;
}
void    print_eegdata(FILE *out, const EEGdata *eeg){
  fprintf(out, "EEGdata:\n"	);
	fprintf(out, "      nbchan  =%i\n", eeg->nbchan);
	fprintf(out, "      n       =%i\n", eeg->n);
	fprintf(out, "      nmarkers=%i\n", eeg->nmarkers);
}

void print_eegdata_trials(FILE *out, const EEGdata_trials *eeg){
  	fprintf(out, "EEGdata_trials:\n"	);
	fprintf(out, "      ntrials=%i\n", eeg->ntrials);
	fprintf(out, "      nmarkers_per_trial=%i\n", eeg->nmarkers_per_trial);
}

void free_eegdata_trials(EEGdata_trials *eeg){
  int i;
  for(i=0; i<eeg->ntrials; i++){
	 free_eegdata(eeg->data[i]);
	 free(eeg->markers[i]);
  }
  free(eeg->data);
  free(eeg->markers);
  free(eeg);
}

void free_eegdata(EEGdata *eeg){
  int i;
  for(i=0; i<eeg->nbchan; i++){
	 free(eeg->d[i]);
  }
  free(eeg->d);
  //  free(eeg);
}

int v_printf(int v, char *format, ...){
  extern int *verbosity;
  int tmp;
  if(!verbosity){
	 tmp=0;
	 verbosity=&tmp;
  }
  va_list arglist;

  va_start(arglist,format);

  if (v<=*verbosity)
    vfprintf(stderr,format,arglist);
  va_end(arglist);

  return 0;
}

size_t  ffread(void *ptr, size_t size, size_t nmemb, FILE *stream){
  int bread;
  /*  dprintf("size=%i, nmemb=%i\n", size, nmemb);*/
  bread=fread(ptr, size, nmemb, stream);
  /*  dprintf("feof=%i, ferror=%i\n", feof(stream), ferror(stream));
		dprintf("bread=%i, exp=%i\n", bread, size*nmemb);*/
  if(bread<nmemb){
	 errormsg(ERR_IO, 1);
  }
}


void    qsort_int_index( int *idx_idx, const int *idx, int n ){
  int i,j;
  int *tmp;
  tmp = (int*)malloc( n*sizeof( int ) );
  memcpy( tmp, idx, n*sizeof( int ) );

  qsort( tmp, n, sizeof( int ), compare_ints);
  for(i=0; i<n; i++){
	 for(j=0; j<n; j++){
		if( idx[i] == tmp[j] ){
		  idx_idx[i]=j;
		  tmp[j]=INT_MIN;
		  break;
		}
	 }
  }
  free(tmp);
}


/* well tested */
void swap_bytes(void *ptr, int nmemb){
  uint8_t tmp;
  uint8_t *nptr=(uint8_t*)ptr;
  int i;

  for(i=0; i<nmemb/2; i++){
	 tmp = nptr[i];
	 nptr[i]=nptr[nmemb-i-1];
	 nptr[nmemb-i-1]=tmp;
  }
}

int is_little_endian(){
  long l=1; 
  void *ptr=&l; 
  uint8_t t =*(uint8_t*)ptr;
  if(t==1) return 0;
  else if(t==0) return 1;
  else errormsg(ERR_ENDIAN, ERR_FATAL);
}

int compare_ints (const int *a, const int *b) {
  if (*a > *b)
	 return 1;
  else if (*a < *b)
	 return -1;
  else
	 return 0;
}

/** wrapper */
void wswap(void *ptr, int nmemb, int flag){
  if(flag)
	 swap_bytes(ptr, nmemb);
}
int vprint_vector(const char* name, double *v, int n){
  int i;
  fprintf(stderr, "%s = \n", name);
  for(i=0; i<n; i++)
    fprintf(stderr, " %2.2f ", v[i]);
  fprintf(stderr, "\n");
  return 0;
}


void   errormsg(int err_no, int fatal){
  switch(err_no){
  case ERR_IO:
    errprintf("IO Error\n");
    break;
  case ERR_GSL:
    errprintf("Error in the GSL-library\n");
    break;
  case ERR_PLOT:
    errprintf("Error in the Plot-library\n");
    break;
  case ERR_ENDIAN:
    errprintf( "Error in the Endianness\n");
    break;
  case ERR_PARSEMAT:
    errprintf("Error while parsing .mat-file! Continue at your own risk...\n");
    break;
  default:
    errprintf("Unknown Error number\n");
  }
  if(fatal){
    errprintf(0, "... Fatal\n");
    exit(err_no);
  } else errprintf("\n");
}

void      reset_warppath(WarpPath *P, int J, int K){
  memset( P->upath, 0, P->J*sizeof(int) );
  memset( P->spath, 0, P->K*sizeof(int) );
  P->J = J; P->K = K;
  memset( P->upath, 0, J*sizeof(int) );
  memset( P->spath, 0, K*sizeof(int) );
}


int       eegdata_cmp_settings( EEGdata *s1, EEGdata *s2 ){
  if( s1->nbchan != s2->nbchan || s1->n != s2->n || s1->nmarkers != s2->nmarkers ){
	 return 1;
  } else {
	 return 0;
  }
}

/** count the number of occurences of c in s
	 \param s - the string (\0 terminated)
	 \param c - the character
	 \return count( s==c )
*/
int strcount( const char *s, char c ){
  int i, count=0;
  for( i=0; i<strlen(s); i++ ){
	 if( s[i]==c ) count++;
  }
  return count;
}


/* ----------------------------------------------------------------------
	MATLAB
   ----------------------------------------------------------------------*/
#ifdef HAVE_MATLAB
Engine* ml_init(){
  Engine *m;
  if (!(m = engOpen(MATLAB_STARTUP_CMD))) {
    dprintf("\nCan't start MATLAB engine\n");
    return NULL;
  }
  dprintf("Matlab successfully started...\n");
  return m;
}

int ml_close(Engine *m){
  return engClose(m);
}


void ml_plot(Engine *matlab, const double *r, const double *v, int n, const char *color, int new){
  mxArray *a1, *a2;
  double *d1,*d2;
  char tmpstring[255];
 
  if(r!=NULL){
    a1 = mxCreateDoubleMatrix(1,n,mxREAL);
    d1 = mxGetPr(a1); 
    d1 = memcpy(d1, r, n*sizeof(double)); 
    if(engPutVariable(matlab, "times", a1))
       fprintf(stderr,"Error putting variable to MATLAB\n");
  }
  a2 = mxCreateDoubleMatrix(1,n,mxREAL);
  d2 = mxGetPr(a2); d2 = memcpy(d2, v, n*sizeof(double));

  if(engPutVariable(matlab, "sig", a2))
    fprintf(stderr,"Error putting variable to MATLAB\n");
  
  if(new) engEvalString(matlab, "figure;");
  if(r!=NULL)
    sprintf(tmpstring, "plot(times, sig, '%s');", color);
  else
    sprintf(tmpstring, "plot(sig, '%s');", color);

  dprintf("Evaluating: '%s'\n",tmpstring); 
  engEvalString(matlab, tmpstring);

  if(r!=NULL)  mxDestroyArray(a1);
  mxDestroyArray(a2);
  return;
} 

void ml_plot_path(Engine *matlab, const int *path, int K){
  int J, m, i, j, k;
  double *x, *y;
  mxArray *xp, *yp;
  char tmpstring[255];

  J = path[K-1];
  m = MAX(J, K);

  xp = mxCreateDoubleMatrix(1,m,mxREAL);
  x = mxGetPr(xp); 

  yp = mxCreateDoubleMatrix(1,m,mxREAL);
  y = mxGetPr(yp); 

  i=0; j=0; k=0; 
  for(i=0; i<m; i++){
    x[i] = (double)k;
    y[i] = (double)j;

    if(path[k]==j){  k++;  }
    else if(path[k]>j) j++;
  }

  if(engPutVariable(matlab, "xpath", xp) ||engPutVariable(matlab, "ypath", yp))
    fprintf(stderr,"Error putting variable to MATLAB\n");
  sprintf(tmpstring, "figure; plot(xpath, ypath); axis([0 %i 0 %i]);", K, J);
  engEvalString(matlab,tmpstring);
  
  mxDestroyArray(xp);  mxDestroyArray(yp);
  return;
}

void ml_wait(Engine *matlab){
  engEvalString(matlab, "waitforbuttonpress;");
}
#endif
