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

/** return number of occurences of c in f.
	 f is rewound to current position 
*/
int stream_count_char( FILE* f, char c ){
  long curpos;
  int numchar=0;
  curpos = ftell( f );
  while( (c=fgetc( f ))!=EOF ){
	 if( c=='\n' ) numchar++;
  }

  fseek( f, curpos, SEEK_SET );

  return numchar;
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
	 errprintf("dest and source not similar enough! Aborting!\n");
	 return -1;
  }
  memcpy( dest->markers, source->markers, (dest->nmarkers)*sizeof(unsigned long) );
  for( i=0; i<dest->nbchan; i++ ){
	 memcpy( dest->d[i], source->d[i], (dest->n)*sizeof( double ) );
  }
  return 0;
}


/** Deep-copy everything from source to a freshly allocated struct;
	 
	 \params source
	 \return dest
 */
EEGdata_trials*      clone_eegdata_trials( const EEGdata_trials *source ){
  int i;
  int nchan;
  int nsamples;
  EEGdata_trials *dest;
  
  if( source->ntrials<1 ){
	 errprintf(" empty eegdata_trials\n");
  }
  nchan = source->data[0]->nbchan;
  nsamples = source->data[0]->n;
  dest = init_eegdata_trials( source->ntrials, source->nmarkers_per_trial, nchan, nsamples, source->times);
  for( i=0; i<source->ntrials; i++ ){
	 copy_similar_eegdata( dest->data[i], source->data[i] );
  }

  return dest;
}

void free_warppath(WarpPath *p){
	free(p->upath);
	free(p->spath);
	free(p);
}

void print_channelinfo( FILE* out, const ChannelInfo *c ){
  fprintf( out, "Channel No. %i/%i - '%s' at (%.2f, %.2f, %.2f)\n", 
			  c->num, c->num_chans, c->label, c->x, c->y, c->z );
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

/** set all data->d[][] to zero
 */
void reset_eegdata( EEGdata* eeg ){
  int i;
  for( i=0; i<eeg->nbchan; i++ ){
	 memset( eeg->d[i], 0, eeg->n*sizeof(double) );
  }
}

/** initializes an eegdata-struct with values. 
	 \param nbchan - num channels
	 \param nsamples - num samples
	 \param nmarkers - num_markers
	 \return the pointer to the new eegdata-struct.
*/
EEGdata* init_eegdata(int nbchan, int nsamples, int nmarkers){
  int i;
  EEGdata *eeg;

  eeg = (EEGdata*) malloc( sizeof(EEGdata) );

  eeg->nbchan = nbchan;
  eeg->n = nsamples;
  eeg->d = (double**) malloc( nbchan * sizeof(double*) );
  for(i=0; i<nbchan; i++){
	 eeg->d[i] = (double*) malloc( nsamples*sizeof(double) );
  }
  eeg->markers = (unsigned long*)malloc( nmarkers*sizeof(unsigned long) );
  eeg->nmarkers=nmarkers;

  return eeg;
}

/** allocate new EEGdata_trials struct
	 \param times is NULL (no information) or a nbsamples long double array containing the times-entries
	              for the eeg-data
	 \return the struct
*/
EEGdata_trials* init_eegdata_trials(int nbtrials, int nmarkers_per_trial, int nbchan, int nbsamples, double *times){
  EEGdata_trials *eeg;
  int i;

  eeg = (EEGdata_trials*)malloc( sizeof( EEGdata_trials ) );
  eeg->nmarkers_per_trial = nmarkers_per_trial;
  eeg->ntrials = nbtrials;
  eeg->markers = (unsigned long **)malloc( nbtrials * sizeof( unsigned long* ) );
  for(i=0; i<nbtrials; i++){
	 eeg->markers[i] = (unsigned long*) malloc( nmarkers_per_trial * sizeof( unsigned long ) );
  }
  eeg->nsamples = nbsamples;
  eeg->times = (double*)malloc( nbsamples*sizeof( double ) );
  if( times!=NULL ){
	 memcpy( eeg->times, times, nbsamples*sizeof( double ) );
  }
  eeg->data = (EEGdata**)malloc( nbtrials*sizeof( EEGdata* ) );
  for(i=0; i<nbtrials; i++){
	 eeg->data[i] = init_eegdata(nbchan, nbsamples, nmarkers_per_trial);
  }
  
  /* compute sampling rate from sampling interval times[1]-times[0] */
  if( times!=NULL ){
	 eeg->sampling_rate = 1000.0/(times[1]-times[0]);
  } else {
	 eeg->sampling_rate = 0;
  }
  return eeg;
}
void    print_eegdata(FILE *out, const EEGdata *eeg){
  fprintf(out, "EEGdata:\n"	);
	fprintf(out, "      nbchan  =%i\n", eeg->nbchan);
	fprintf(out, "      n       =%i\n", eeg->n);
	fprintf(out, "      nmarkers=%i\n", eeg->nmarkers);
	if( eeg->nmarkers>0 ){
	  fprintf(out, "        markers[0]=%ld\n", eeg->markers[0]);
	  fprintf(out, "        markers[nm-1]=%ld\n",eeg->markers[eeg->nmarkers-1]);
	}
	if( eeg->n>0 ){
	  fprintf(out, "        d[0][0]=%f\n", eeg->d[0][0]);
	  fprintf(out, "        d[0][n-1]=%f\n",eeg->d[0][eeg->n-1]);
	}
}

void print_eegdata_trials(FILE *out, const EEGdata_trials *eeg){
  	fprintf(out, "EEGdata_trials:\n"	);
	fprintf(out, "      ntrials=%i\n", eeg->ntrials);
	fprintf(out, "      nmarkers_per_trial=%i\n", eeg->nmarkers_per_trial);
	if( eeg->nmarkers_per_trial>0 ){
	  fprintf(out, "        markers[0][0]=%ld\n", eeg->markers[0][0]);
	  fprintf(out, "        markers[0][nm-1]=%ld\n",eeg->markers[0][eeg->nmarkers_per_trial-1]);
	}
	fprintf(out, "      sampling_rate=%f\n", eeg->sampling_rate);
	if( eeg->times!=NULL ){
	  fprintf(out, "      times[0]=%f\n", eeg->times[0]);
	  fprintf(out, "      times[n-1]=%f\n", eeg->times[eeg->nsamples-1]);
	} else {
	  fprintf(out, "      times[0]=<NULL>\n");
	}
	if( eeg->data[0] && eeg->data[0]->d[0] && eeg->data[0]->d[0][0] ){
	  fprintf(out, "      data[0]->d[0][0]=%f\n", eeg->data[0]->d[0][0] );
	  fprintf(out, "      data[0]->d[0][n-1]=%f\n", eeg->data[0]->d[0][eeg->nsamples-1] );
	} else {
	  fprintf(out, "      data[0]->d[0][0]=<NULL>\n" );
	}
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
  free(eeg->markers);
  free(eeg);
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
  size_t bread;
  /*  dprintf("size=%i, nmemb=%i\n", size, nmemb);*/
  bread=fread(ptr, size, nmemb, stream);
  /*  dprintf("feof=%i, ferror=%i\n", feof(stream), ferror(stream));
		dprintf("bread=%i, exp=%i\n", bread, size*nmemb);*/
  if(bread<nmemb){
	 errormsg(ERR_IO, 1);
  }
  return bread;
}
size_t ffwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream){
  size_t bwrite;
  /*  dprintf("size=%i, nmemb=%i\n", size, nmemb);*/
  bwrite=fwrite(ptr, size, nmemb, stream);
  /*  dprintf("feof=%i, ferror=%i\n", feof(stream), ferror(stream));
		dprintf("bread=%i, exp=%i\n", bread, size*nmemb);*/
  if(bwrite<nmemb){
	 errormsg(ERR_IO, 1);
  }
  return bwrite;
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
  if(t==1){
	 return 0;
  } else if(t==0){ 
	 return 1;
  } else {
	 errormsg(ERR_ENDIAN, ERR_FATAL);
  }
  return 0;
}

int compare_ints (const void *a, const void *b) {
  int *i1, *i2;
  i1=(int*)a;
  i2=(int*)b;
  if (*i1 > *i2)
	 return 1;
  else if (*i1 < *i2)
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
    errprintf("... Fatal\n");
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


int       eegdata_cmp_settings( const EEGdata *s1, const EEGdata *s2 ){
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
