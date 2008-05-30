#include "reader.h"
#include "helper.h"


/** read xdim x ydim matrix from ascii file with delimiter 'delim';
	 if NULL is passed as d, the functin allocates memory.
 */
double** read_double_matrix_ascii(const char *fname, int xdim, int ydim, double **d){
  FILE *f;
  int x,y;

  if((f=fopen(fname, "r"))==NULL)
	 errormsg(ERR_IO, 1);

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

/**  Reads EEG-data from binary file. The format is as follows (each segment
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
  int i, c;
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

  eeg = init_eegdata_trials(ntrials, nmarkers, nbchan, nsamples);
  ffread(eeg->times, sizeof(double), nsamples, f);
  
  /* read markers */
  for( i=0; i<ntrials; i++ ){
	 ffread(eeg->markers[i], sizeof(double), nmarkers, f);
	 memcpy( eeg->data[i]->markers, eeg->markers[i], nmarkers*sizeof(double) );
  }

  /* read data */
  for( c=0; c<nbchan; c++ )
	 for( i=0; i<ntrials; i++ )
		ffread( eeg->data[i]->d[c], sizeof( double ), nsamples, f );

  return eeg;
}

/** parse Matfile Header. Assume file pointer is at beginning of file.
 */
void matfile_get_header(FILE *f, char *description, int *swapflag){ 
  uint16_t version; 
  uint16_t endian;
  uint8_t  *tmp;
  
  dprintf("islittle=%i\n", is_little_endian());
  /* header */
  ffread(description, sizeof(uint8_t), 124, f);
  dprintf("%s\n", description);
  ffread(&version, sizeof(uint16_t), 1, f);
  ffread(&endian, sizeof(uint16_t), 1, f);
  tmp = (uint8_t*)&endian;
  dprintf("Endianness: %c%c\n", tmp[0], tmp[1]);
  if(tmp[0]=='I' && tmp[1]=='M') *swapflag=1;
  else *swapflag=0;
  dprintf("version2=%#x\n", version);
}
/** parse Matfile Matrix-Header. Assume file pointer is at beginning
	 of a new miMATRIX-Element.
	 After the function, the following parameters are:
	 \param type - pointer to one integer according to mxTYPE 
	 \param ndims - number of dimensions
	 \param dims - ndims-long array holding the dimensions
	 \param namebuf - matlab's variable name for that matrix
 */
void matfile_matrix_header(FILE *f, uint8_t *type, int *ndims, uint32_t *dims, char *namebuf, int swapflag){
  uint8_t forget[16];
  uint32_t flags, first;
  uint8_t  *tmp8;
  uint16_t *tmp16;
  uint32_t nbytes;
  int namelen;
  int i;

  /* flag subelement, is never compressed */
  ffread(forget, sizeof(uint32_t), 2, f); /* header of flag-subelement */
  ffread(&flags, sizeof(uint32_t), 1, f);
  dprintf("flags = %#x\n", flags);  
  wswap(&flags, 4, swapflag);
  tmp8 = (uint8_t*)&flags;
  *type=tmp8[3];
  dprintf("class=%i\n", *type);
  ffread(forget, sizeof(uint32_t), 1, f); /* undefined */
  
  /* dimensions subelement, could be compressed */
  ffread(&first, sizeof(uint32_t), 1, f); /* dims-type */
  if(is_compressed_format(first, swapflag)){
	 dprintf("this is compressed!!! COULD GO WRONG!\n");
	 tmp16 = (uint16_t*)&first;
	 nbytes = (uint32_t)tmp16[1];
  } else {
	 ffread(&nbytes, sizeof(uint32_t), 1, f); 
  }
  *ndims=nbytes/4;
  dprintf("ndims=%i\n", *ndims);
  ffread(dims, sizeof(uint32_t), *ndims, f);
  for(i=0; i<*ndims; i++)
	 dprintf("ndims[%i]=%i\n", i+1, dims[i]);

  /* name subelement */
  namelen=0;
  ffread(&first, sizeof(uint32_t), 1, f); /* name-type */
  if(is_compressed_format(first, swapflag)){
	 tmp16  = (uint16_t*)&first;
	 nbytes = (uint32_t)tmp16[1];	 
	 namelen = 4;
  } else {
	 ffread(&nbytes, sizeof(uint32_t), 1, f); 
	 if(nbytes>0)
		namelen = nbytes+(nbytes/8+1)*8;
  }
  dprintf("nbytes=%i, namelen=%i\n", nbytes, namelen);
  ffread(namebuf, sizeof(uint8_t), namelen, f);
  namebuf[namelen]='\0';
}

/** assumes that f points to a struct-field in the matfile where 
	 exactly one int is stored (trials, nbchan...) in a miMATRIX
	 field.
*/
double matfile_get_single_value(FILE *f, int swapflag){
  uint32_t dtype, nbytes; /* these are for the miMATRIX */
  uint32_t dtypea, nbytesa; /* for miDOUBLE */
  uint16_t *tmp16;
  int ndims;
  uint32_t dims[20];
  uint8_t mtype;
  char namebuf[255];

  ffread( &dtype, sizeof( uint32_t ), 1, f );
  if( !is_compressed_format( dtype, swapflag ) ){
	 ffread(&nbytes, sizeof(uint32_t), 1, f);
  } else {
	 errormsg( ERR_PARSEMAT, ERR_NOFATAL );
  }

  if(dtype==miMATRIX){
	 dprintf("got matrix\n");
	 matfile_matrix_header(f, &mtype, &ndims, dims, namebuf, swapflag); 
	 if(ndims<2 || dims[0]!=1 || dims[1]!=1 || mtype!=mxDOUBLE_CLASS )
		errormsg( ERR_PARSEMAT, ERR_NOFATAL );
	 dprintf("namebuf=%s of (%i x %i), ndims=%i\n", namebuf, dims[0], dims[1], ndims);
	 ffread( &dtypea, sizeof( uint32_t ), 1, f);
	 if( is_compressed_format( dtypea, swapflag ) ){
		tmp16 = (uint16_t*) &dtype;
		nbytesa= (uint32_t) tmp16[ (swapflag?0:1) ];
	 } else {
		ffread( &nbytesa, sizeof( uint32_t ), 1, f);
	 }
	 dprintf("nbytesa=%i\n", nbytesa);
  } else {
	 errormsg(ERR_PARSEMAT, ERR_NOFATAL);
  }
}

/** parse .mat file pointed to by f and fill EEGdata_trials-struct eeg.
	 It is assumed that f points to the first byte after the name "EEG"
	 in the .mat file.
	 \b CAUTION: The matlab-reader is not yet functional! don't use!
*/
void parse_eeg_struct_from_matfile(FILE *f, EEGdata_trials *eeg, int swapflag){
  uint32_t thirtytwo;
  uint32_t dtype, nbytes, first;

  int nfields, i, j;
  char namebuf[255];
  enum fields { NBCHAN=0, TRIALS, PNTS, TIMES, DATA, EPOCH };
  int idx[6];   /* indices in mat-file for the fields:
						 - nbchan, trials, pnts, times,  data, epoch */
  int idx_idx[6]; /* index array for sorted idx-array */
  unsigned long byte_idx[6]; /* position in file stream, to be used with fseek() */
  int nbchan, trials, pnts;

  /* field name length subelement, always compressed,
	  should give 32 bytes for each name as restriction*/
  fseek(f, sizeof( uint32_t ), SEEK_CUR);
  ffread(&thirtytwo, sizeof( uint32_t ), 1, f);

  /* field names subelement */
  fseek(f, sizeof( uint32_t ), SEEK_CUR);
  ffread(&nbytes, sizeof( uint32_t ), 1, f);
  nfields = nbytes/thirtytwo;
  dprintf("nbytes=%i, nfields=%i\n", nbytes, nfields);

  for(i=1; i<=nfields; i++){
	 ffread(namebuf, sizeof( uint8_t ), thirtytwo, f);
	 v_printf(99, "read buf[%i]='%s'\n", i, namebuf);
	 if(strcmp(namebuf, "nbchan")==0){
		dprintf("reading nbchan\n");
		idx[NBCHAN] = i;
	 } else if( strcmp( namebuf, "trials" ) == 0 ){
		dprintf("reading trials\n");
		idx[TRIALS] = i;
	 } else if( strcmp( namebuf, "pnts" ) == 0 ){
		dprintf("reading pnts\n");
		idx[PNTS] = i;
	 } else if( strcmp( namebuf, "times" ) == 0 ){
		dprintf("reading times\n");
		idx[TIMES] = i;
	 } else if( strcmp( namebuf, "data" ) == 0 ){
		dprintf("reading data\n");
		idx[DATA] = i;
	 } else if( strcmp( namebuf, "epoch" ) == 0 ){
		dprintf("reading epoch\n");
		idx[EPOCH] = i;
	 } else {
		v_printf(99, "...ignoring\n");
	 }
  }

  qsort_int_index( idx_idx, idx, 6 );

  j=1;
  for(i=0; i<6; i++){ /* fill byte_idx array */
	 while( j < idx[idx_idx[i]] ){/* skip uninteresting fields */
		ffread(&first, sizeof( uint32_t ), 1, f);
		if( is_compressed_format( first, swapflag ) ){
		  fseek( f, sizeof( uint32_t ), SEEK_CUR );
		} else {
		  ffread( &nbytes, sizeof( uint32_t ), 1, f );
		  fseek( f, nbytes, SEEK_CUR );
		}
		j++;
	 }
	 
	 dprintf("j=%i\n",j);
	 byte_idx[idx_idx[i]] = ftell( f ); /* get byte position */

	 ffread( &first, sizeof( uint32_t ), 1, f ); /* and skip */
	 if( is_compressed_format( first, swapflag ) ){
		fseek( f, sizeof( uint32_t ), SEEK_CUR );
	 } else {
		ffread( &nbytes, sizeof( uint32_t ), 1, f );
		fseek( f, nbytes, SEEK_CUR );
	 }
	 j++;
  }	
  
  for(i=0; i<6; i++){
	 dprintf("idx[%i]=%i, byteidx=%i (%#x), idx_idx[%i]=%i\n", i, idx[i], byte_idx[i],byte_idx[i], i, idx_idx[i]);
  }

 /* get nbchans */ 
  fseek( f, byte_idx[NBCHAN], SEEK_SET );
  dprintf("at position=%#x\n", ftell( f ) );
  nbchan = (int) matfile_get_single_value( f, swapflag );

  /* get trials */
  fseek( f, byte_idx[TRIALS], SEEK_SET );
  dprintf("at position=%#x\n", ftell( f ) );
  trials = (int) matfile_get_single_value( f, swapflag );

 
}

/** reads segmented EEG-data from eeglab-set-file .
	 \b CAUTION: The matlab-reader is not yet functional! don't use!
 * \param file -- complete filepath to eeglab set-file
 */
EEGdata_trials* read_segmented_eeg_from_eeglabset(const char *file){
  FILE *f;
  EEGdata_trials *eeg;
  char description[124];
  uint32_t dtype, nbytes;
  uint8_t mtype;
  int ndims;
  uint32_t dims[10];
  char namebuf[255];

  int i;
  int swapf; /* flag; do byte swapping or not? */

  if((f = fopen(file, "rb"))==NULL) errormsg(ERR_IO, 1);
  eeg = (EEGdata_trials*)malloc(sizeof(EEGdata_trials));

  matfile_get_header(f, &description, &swapf);

  while(1){ /* find EEG-struct */
 	 ffread(&dtype, sizeof(uint32_t), 1, f);
	 ffread(&nbytes, sizeof(uint32_t), 1, f);

	 dprintf("compressed: %i\n", is_compressed_format(dtype, swapf));
	 dprintf("dt=%li\n", dtype);
	 dprintf("nb=%li\n", nbytes);

	 if(dtype==miMATRIX){
		dprintf("found matrix\n");
		matfile_matrix_header(f, &mtype, &ndims, dims, namebuf, swapf); 
		dprintf("type=%i, ndims=%i, name=%s\n", mtype, ndims, namebuf);
		if(mtype==mxSTRUCT_CLASS){
		  dprintf("it is a struct\n");		  
		  if(strcmp(namebuf, "EEG")==0){
			 dprintf("EEG-struct found\n");
			 parse_eeg_struct_from_matfile(f, eeg, swapf);
		  }
		}
	 }
	 fseek(f, nbytes, SEEK_CUR);
	 dprintf("fp=%li\n", ftell(f));
  }
  
  fclose(f);
  return eeg;
}

/** Determines whether Data element is in compressed format in 
 * the Matlab .mat-file.
 * \param first -  32 bit integer
 */
int is_compressed_format(uint32_t first, int swapflag){
  uint16_t *tmp;
  wswap(&first, 4, swapflag);
  tmp = (uint16_t*)&first;
  if(tmp[0]!=0){
	 dprintf("this is a compressed element\n"); 
	 return 1;
  } else 
	 return 0;
}
