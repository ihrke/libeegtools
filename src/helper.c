/* helper.c
 *
 */	
#include "helper.h"

/* ---------------------------------------------------------------------------- 
   -- Helper functions                                                       -- 
   ---------------------------------------------------------------------------- */

/** returns a random integer in the range [from, from+1, ..., to]
 */
int randint( int from, int to ){
  return ((int)(drand48() * to))+from;
}


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
	

/** displays a progress-bar like this
	 \code
	 [ #####\                          ]
	 \endcode
	 \param flag one of PROGRESSBAR_INIT, PROGRESSBAR_CONTINUE_LONG, 
	             PROGRESSBAR_CONTINUE_SHORT
	 \param num meaning depends on flag
 */
void progressbar_rotating( int flag, int num ){
  int c, i;

  switch(flag){
  case PROGRESSBAR_INIT:
	 progress_status.max_progress = num;
	 progress_status.cur_progress = 0;
	 progress_status.prev_progress= 0;
	 fprintf( stdout, "[ " );
	 for( i=0; i<PROGRESSBAR_NUMCOLS; i++ ){
		fprintf( stdout, " " );
	 }
	 fprintf( stdout, " ]" );
	 for( i=0; i<PROGRESSBAR_NUMCOLS+2; i++ ){
		fprintf( stdout, "\b" );
	 }
	 break;
  case PROGRESSBAR_CONTINUE_LONG:
	 c = (num*PROGRESSBAR_NUMCOLS/progress_status.max_progress);
	 //printf("c=%i, cur_progress=%i, num=%i, p=%i\n", c, cur_progress, num, max_progress);
	 if( c>progress_status.cur_progress ){
		fprintf( stdout, "#" );
		progress_status.cur_progress++;
	 }
	 break;
  case PROGRESSBAR_CONTINUE_SHORT:
	 c = (progress_status.prev_progress++ % 4);
	 switch(c){
	 case 0: c = '/'; break;
	 case 1: c = '-'; break;
	 case 2: c = '\\'; break;
	 case 3: c = '|'; break;
	 }
	 fprintf( stdout, "%c", c);
	 fprintf( stdout, "\b" );
	 break;
  case PROGRESSBAR_FINISH:
	 fprintf( stdout, "\n");
	 break;
  } /* switch */
}
