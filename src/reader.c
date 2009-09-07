#include "reader.h"
#include "helper.h"
#include <matio.h>
#include "eeg.h"

#ifdef MATIO
double get_double_from_struct_field( matvar_t *eeg, const char *name, int struct_array_index );

/** This reader uses MatIO to parse EEGlab files and construct
	 a LibEEGTools EEG struct.
	 It was developed with EEGlab version 6.01b.
	 
	 Currently, the EEGlab set must be "epoched", that means that the data must be
	 three-dimensional with channels x n x trials. Continuous data is not yet supported
	 and the reader will fail.

	 \param eeglab .set file
	 \return the EEG struct
*/
EEG* read_eeglab_file( const char *file ){
  mat_t *mat;
  matvar_t *meeg; /* contains the struct 'EEG' from EEGlab */
  matvar_t *tmp, *tmp2, *event, *epoch;  
  int nfields; 
  int c,i,j;
  EEG *eeg;
  int nbchan, ntrials, n;

  dprintf("Reading file: '%s'\n", file);
  mat = Mat_Open( file, MAT_ACC_RDONLY);
  if( !mat ){
	 errprintf("Error opening MATLAB file '%s'\n", file );
	 return NULL;
  }
  meeg = Mat_VarRead( mat, "EEG" );
  if( meeg->class_type!=MAT_C_STRUCT ){
	 errprintf("EEG does not appear to be a struct\n" );
	 return NULL;
  }
  nfields = Mat_VarGetNumberOfFields( meeg );
#ifdef DEBUG
  dprintf( "There are %i fields in the EEG struct\n", nfields );
  for( i=1; i<=nfields; i++ ){ /* MATLAB 1-relative indices */
	 tmp = Mat_VarGetStructField( meeg, &i, BY_INDEX, 0 );
	 dprintf("Found field: '%s'\n", tmp->name);
  }
#endif

  /* dimensions */
  nbchan = (int)get_double_from_struct_field( meeg, "nbchan",0 );
  ntrials= (int)get_double_from_struct_field( meeg, "trials",0 );
  n      = (int)get_double_from_struct_field( meeg, "pnts",0 );
  eeg = eeg_init( nbchan, ntrials, n );

  /* filename  */
  eeg->filename=(char*)malloc( (strlen(file)+2)*sizeof(char) );
  strcpy( eeg->filename, file );

  /* comments */
  tmp = Mat_VarGetStructField( meeg, "comments", BY_NAME, 0 );
  eeg->comment=(char*)malloc( tmp->dims[0]*sizeof(char) );
  for( i=0; i<tmp->dims[0]+1; i++ ){
	 eeg->comment[i]='\0';
  }
  memcpy( eeg->comment, tmp->data, tmp->dims[0]*sizeof(char));

  /* sampling rate */
  eeg->sampling_rate = get_double_from_struct_field( meeg, "srate", 0);

  /* times */
  tmp = Mat_VarGetStructField( meeg, "times", BY_NAME, 0 );
  if( tmp->dims[1]!=n ){
	 errprintf("times-array should be of length n: %i != %i\n", tmp->dims[1], n );
	 eeg_free( eeg );
	 return NULL;
  } 
  if( tmp->data_size != sizeof(double) ){
	 errprintf("times is not double format, %i!=%li\n", tmp->data_size, sizeof(double));
	 eeg_free( eeg );
	 return NULL;	 
  }
  eeg->times = (double*) malloc( n*sizeof(double) );
  memcpy(eeg->times, tmp->data, n*sizeof(double) );

  /* channel info */
  eeg->chaninfo = (ChannelInfo*)malloc( nbchan*sizeof(ChannelInfo) );
  tmp = Mat_VarGetStructField( meeg, "chanlocs", BY_NAME, 0 );
  dprintf("chanlocs: %i,%i\n", tmp->dims[0], tmp->dims[1]);
  
  for( i=0; i<nbchan; i++ ){
	 tmp2 = Mat_VarGetStructField( tmp, "labels", BY_NAME, i );
	 eeg->chaninfo[i].num = i;
	 eeg->chaninfo[i].num_chans = nbchan;
	 strcpy(eeg->chaninfo[i].label, (char*)tmp2->data);
	 eeg->chaninfo[i].x = get_double_from_struct_field( tmp, "X", i);
	 eeg->chaninfo[i].y = get_double_from_struct_field( tmp, "Y", i);
	 eeg->chaninfo[i].z = get_double_from_struct_field( tmp, "Z", i);
  }
  
  
  /* data */
  tmp = Mat_VarGetStructField( meeg, "data", BY_NAME, 0 );
  if( tmp->dims[0]!=nbchan || tmp->dims[1]!=n || tmp->dims[2]!=ntrials ){
	 errprintf("(nbchan,ntrials,n)=(%i,%i,%i), should be (%i,%i,%i)\n",
				  tmp->dims[0], tmp->dims[2], tmp->dims[1], nbchan, ntrials, n );
	 eeg_free( eeg );
	 return NULL;
  }
  if( tmp->data_size != sizeof(float) ){
	 errprintf("data is not in float format, sizeof(data)=%i, sizeof(float)=%li\n", 
				  tmp->data_size, sizeof(float));
	 eeg_free( eeg );
	 return NULL;	 
  } 
  float x;
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<ntrials; i++ ){
		for( j=0; j<n; j++ ){
		  x=((float*)tmp->data)[ c*n*ntrials+i*n+j ];
		  eeg->data[c][i][j] = (double)x;
		}
	 }
  }
  
  /* markers */
  epoch = Mat_VarGetStructField( meeg, "epoch", BY_NAME, 0 );
  if( epoch->dims[0] == 0 || epoch->dims[1] < ntrials ){
	 errprintf("no epoch field, or wrong dimensions (%i,%i)\n", 
				  epoch->dims[0],epoch->dims[1]);
	 eeg_free( eeg );
	 return NULL;
  }
  eeg->nmarkers = (unsigned int*) malloc( ntrials*sizeof(unsigned int) );
  eeg->markers = (unsigned int**) malloc( ntrials*sizeof(unsigned int*) );
  eeg->marker_labels = (char***) malloc( ntrials*sizeof(char**) );

  for( i=0; i<ntrials; i++ ){
	 tmp  = Mat_VarGetStructField( epoch, "eventlatency", BY_NAME, i );
	 tmp2 = Mat_VarGetStructField( epoch, "eventtype", BY_NAME, i );
	 dprintf("%i, %i\n", i, tmp->dims[1]);
	 eeg->nmarkers[i] = tmp->dims[1];
	 eeg->markers[i] = (unsigned int*) malloc( eeg->nmarkers[i]*sizeof(unsigned int) );
	 eeg->marker_labels[i] = (char**) malloc( eeg->nmarkers[i]*sizeof(char*) );
	 for( j=0; j<eeg->nmarkers[i]; j++ ){
		/* label */
		event = Mat_VarGetCell( tmp2, j ); /* MATLAB index */
		eeg->marker_labels[i][j] = (char*)malloc( (strlen((char*)event->data)+1)*sizeof(char) );
		strcpy( eeg->marker_labels[i][j], (char*)event->data );
		/* latency */
		event = Mat_VarGetCell( tmp, j ); /* MATLAB index */
		eeg->markers[i][j] = closest_index( eeg->times, n, ((double*)event->data)[0] );
	 }
  }

  dprintf("Finished reading '%s'\n", file );

  return eeg;
}


double get_double_from_struct_field( matvar_t *eeg, const char *name, int struct_array_index ){
  matvar_t *tmp; 

  tmp = Mat_VarGetStructField( eeg, (char*)name, BY_NAME, struct_array_index );
  if( tmp->rank != 1 && tmp->dims[0]<1 ){
	 errprintf("field '%s' wrong, rank=%i,tmp->dims[0]=%i\n",name, tmp->rank,tmp->dims[0] );
	 return -1;
  }  
  dprintf( "found: %s=%f\n", name, (((double*)tmp->data)[0]) );
  return (((double*)tmp->data)[0]);
}

#endif /* MATIO */


/** reads from FILE* until '\n' and puts it into line.
	 '\0' is added to line.
	 \param f
	 \param line empty string (long enough), or NULL ->own memory
	 \return NULL if error, else ptr to line
 */
char* read_line( FILE *f, char *line ){
  long curpos;
  int length=0;
  char c;
  int i;

  if( line==NULL ){ /* find out how much memory is needed */
	 curpos = ftell( f );
	 while( (c=fgetc( f ))!='\n' ) 
		length++;
	 dprintf( "length=%i\n", length );
	 fseek( f, curpos, SEEK_SET );
	 line = (char*) malloc( (length+1)*sizeof(char) );
  }
  i = 0;
  while( (c=fgetc( f ))!='\n' ){
	 line[i++]=c;
  }
  line[i] = '\0';
  return line;
}

/** reads ChannelInfo struct from .ced file (EEGlab).
	 Format is: One line with labels (tab-sep), followed by lines of values (tab-sep).
	 This function reads:
	   - labels: "labels" "X", "Y", "Z"
		- values: char*  double double double
		\param fname filename
		\param chans ChannelInfo struct (large enough) or NULL (memory allocated in function)
		\return chaninfo read from file or NULL if error occured
*/
ChannelInfo* read_chaninfo_ced( const char *fname, ChannelInfo *chans ){
  int num_chans=0;
  FILE *f;
  char buf[MAX_LINE_LENGTH];
  int Label_idx=-1, 
	 X_idx=-1, Y_idx=-1, Z_idx=-1;
  char *strstart, *strend;
  int idx, i;

  if( (f=fopen( fname, "r" ))==NULL ){
	 errprintf("couldn't open file '%s'\n", fname );
	 return NULL;
  }
  
  read_line( f, buf );
  dprintf("labels='%s'\n", buf);

  /* labels present ?*/
  strstart=buf;
  idx = 0;
  while( (strend=strchr( strstart, '\t' ))!=NULL){
	 *strend='\0';
	 dprintf( "curlab='%s'\n", strstart );
	 if( !strcasecmp( strstart, "labels" ) ){
		Label_idx=idx;
	 } else if( !strcasecmp( strstart, "X" ) ){
		X_idx=idx;
	 } else if( !strcasecmp( strstart, "Y" ) ){
		Y_idx=idx;
	 } else if( !strcasecmp( strstart, "Z" ) ){
		Z_idx=idx;
	 }
	 strstart = strend+1;
	 idx++;
  }
  dprintf(" l=%i, X=%i, Y=%i, Z=%i\n", Label_idx, X_idx, Y_idx, Z_idx );
  if( X_idx<0 || Y_idx<0 || Z_idx<0 || Label_idx<0 ){
	 errprintf( " ERROR: did not find one of the labels: X=%i, Y=%i, Z=%i, Labels=%i\n", 
					X_idx, Y_idx, Z_idx, Label_idx );
	 return NULL;
  }

  num_chans = stream_count_char( f, '\n' );
  dprintf(" num_chans = %i\n", num_chans );

  if( chans==ALLOC_IN_FCT ){
	 chans = (ChannelInfo*) malloc( num_chans*sizeof(ChannelInfo) );
  } else { /* num_chans from chans-struct */
	 if( num_chans!=chans[0].num_chans ){
		errprintf( " ERROR: not enough memory in chans-struct\n");
		return NULL;
	 }
  }
  
  for( i=0; i<num_chans; i++ ){ /* for all channels */  
	 read_line( f, buf );
	 strstart=buf;
	 idx = 0;
	 chans[i].num_chans=num_chans;
	 chans[i].num = i+1;
	 while( (strend=strchr( strstart, '\t' ))!=NULL){
		*strend='\0';
		if( idx==Label_idx ){
		  strcpy( chans[i].label, strstart );
		} else if( idx==X_idx ){
		  chans[i].x = atof( strstart );
		} else if( idx==Y_idx ){
		  chans[i].y = atof( strstart );
		} else if( idx==Z_idx ){
		  chans[i].z = atof( strstart );
		}
		strstart = strend+1;
		idx++;
	 }
	 print_channelinfo( stderr, chans+i );
  }

  return chans;
}

/** read xdim x ydim matrix from ascii file with delimiter 'delim';
	 if NULL is passed as d, the functin allocates memory.
 */
double** read_double_matrix_ascii(const char *fname, int xdim, int ydim, double **d){
  FILE *f;
  int x,y;

  if((f=fopen(fname, "r"))==NULL){
	 errprintf("failed reading %s\n", fname );
	 return NULL;
  }

  if(d==NULL){
	 d = (double**) malloc( ydim*sizeof(double*) );
	 for( y=0; y<ydim; y++ )
		d[y] = (double*) malloc( xdim*sizeof(double) );
  }

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fscanf(f, " %lf ", &(d[y][x]));
	 }
  }
  fclose(f);
  return d;
}

/** read 1D vector frome file 
 */
double* read_double_vector_ascii( const char *fname, int N, double *v ){
  FILE *f;
  int i;

  if((f=fopen(fname, "r"))==NULL)
	 errormsg(ERR_IO, 1);

  if(v==NULL){
	 v = vector_init( NULL, N, 0.0 );
  }
  for( i=0; i<N; i++ ){
	 fscanf(f, " %lf ", &(v[i]));
  }
  fclose( f );
  return v;
}

/** Reads EEG-data from binary file. Format is Cxn doubles.
 * \param C - number of channels
 * \param n - number of samples per channel
 */
EEGdata* read_continuous_eeg_from_binfile(const char *file, int C, int n){
  int i;
  EEGdata *eeg;
  FILE *f;
  
  if((f = fopen(file, "rb"))==NULL) return NULL;

  eeg = (EEGdata*)malloc(sizeof(EEGdata));
  eeg->nbchan=C;
  eeg->n=n;
  eeg->d=(double**)malloc(C*sizeof(double*));
  for(i=0; i<C; i++){
	 eeg->d[i]=(double*)malloc(n*sizeof(double));
	 if(fread(eeg->d[i], sizeof(double), n, f)<n){
		dprintf("ERROR: read less bytes than requested\n");
	 }
  }
  return eeg;
}

/**  Reads EEG-data from binary file. See also \ref rawfileformat
	  The format is as follows (each segment
	  in the table is a 64 bits double):
	  \code
	  1) number of channels
	  2) number of trials
	  3) number of samples per trial segment
	  4) number of markers per trial
	  5) n-doubles giving the times-array (for each sample from 1...n this array
	     gives the corresponding time in ms)
 	  6) num_markers*num_trials-doubles giving markers in samples [1,...,n]
	  7) raw EEG-data in the format of:
	      - channels x trial x samples
			- i.e. first all trials of the first channel one after the other, than
			  the 2nd and so on
	  \endcode
\param file - name of file
*/
EEGdata_trials* read_eegtrials_from_raw(const char *file){
  EEGdata_trials* eeg;
  FILE *f;
  int i, j, c;
  double nbchan_d, ntrials_d, nsamples_d, nmarkers_d;
  unsigned int nbchan, ntrials, nsamples, nmarkers;

  if((f = fopen(file, "rb"))==NULL) errormsg(ERR_IO, 1);
  ffread( &nbchan_d,   sizeof(double), 1, f );
  ffread( &ntrials_d,  sizeof(double), 1, f );
  ffread( &nsamples_d, sizeof(double), 1, f );
  ffread( &nmarkers_d, sizeof(double), 1, f );

  nbchan   = (unsigned int) nbchan_d;
  ntrials  = (unsigned int) ntrials_d;
  nsamples = (unsigned int) nsamples_d;
  nmarkers = (unsigned int) nmarkers_d;

  dprintf("(nbchan, ntrials, nsamples, nmarkers) = (%i,%i,%i,%i)\n", nbchan, ntrials, nsamples, nmarkers);

  /* allocating all memory */
  eeg = init_eegdata_trials(ntrials, nmarkers, nbchan, nsamples, NULL);

  /* read times-array */
  ffread(eeg->times, sizeof(double), nsamples, f);
  eeg->sampling_rate=1000.0/(eeg->times[1]-eeg->times[0]);

  /* read markers */
  double *tmp;
  tmp = (double*) malloc( nmarkers * sizeof(double) );
  for( i=0; i<ntrials; i++ ){
	 ffread(tmp, sizeof(double), nmarkers, f);
	 /* dprintf("tmp[0]=%f, tmp[1]=%f\n", tmp[0], tmp[1]); */
	 for( j=0; j<nmarkers; j++){
		eeg->markers[i][j] = (unsigned long)tmp[j];
	 } 
	 memcpy( eeg->data[i]->markers, eeg->markers[i], nmarkers*sizeof(unsigned long) );
  }
  free(tmp);

  /* read data */
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<ntrials; i++ ){
		/*		dprintf("eeg->data[%i]->n = %i\n", i, eeg->data[i]->n);*/
		ffread( eeg->data[i]->d[c], sizeof( double ), nsamples, f );
	 }
  }

  return eeg;
}
