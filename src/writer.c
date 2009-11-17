#include "writer.h"
#include "eeg.h"
#include <matio.h>

#ifdef MATIO
/** This write uses MatIO to create an EEGlab file.
	 It was developed with EEGlab version 6.01b.
	 
	 \param eeg the struct
	 \param file eeglab .set file
	 \return error code
*/
int write_eeglab_file( EEG* eeg, const char *file ){
  mat_t *mfile;
  matvar_t *meeg;
  if( !(mfile=Mat_Create( file, NULL )) ){
	 errprintf("Could not open '%s' for writing\n");
	 return -1;
  }

  int dims[2]={1,1}, rank=2;
  //meeg = Mat_VarCreate( "EEG", MAT_T_STRUCT, , rank, dims, NULL, 0 );
  
  Mat_VarWrite( mfile, meeg, 0 );
  return 0;
}
#endif //MATIO


/** writes header of raw-file.
	 \param f file pointer
	 \param nbchan number of channels
	 \param nbtrials number of trials
	 \param nsamples number of samples
	 \param nmarkers number of markers per trial
 */
void write_raw_header( FILE *f, int nbchan, int nbtrials, int nsamples, 
							  int nmarkers ){
  double nbchan_d, ntrials_d, nsamples_d, nmarkers_d;
  
  nbchan_d   = (double)nbchan;
  ntrials_d  = (double)nbtrials;
  nsamples_d = (double)nsamples;
  nmarkers_d = (double)nmarkers;

  ffwrite( &nbchan_d,   sizeof(double), 1, f );  
  ffwrite( &ntrials_d,  sizeof(double), 1, f );
  ffwrite( &nsamples_d, sizeof(double), 1, f );
  ffwrite( &nmarkers_d, sizeof(double), 1, f );
}

/** writes markers from EEGdata to RAW-file.
 */
void write_raw_markers( FILE *f, const EEGdata *eeg ){
  double *tmp;
  int i;

  tmp = (double*) malloc( eeg->nmarkers*sizeof(double) );
  for( i=0; i<eeg->nmarkers;  i++ ){ /* typecast */
	 tmp[i] = (double)eeg->markers[i];
  }
  ffwrite( tmp, sizeof(double), eeg->nmarkers, f );
  free( tmp );
}

/** write a file in RAW-file format from an EEGdata_trials struct.
 */
void write_eegtrials_to_raw( const EEGdata_trials *eeg, FILE *f ){
 
  int i, c, t;

  write_raw_header( f, eeg->data[0]->nbchan, eeg->ntrials, 
						  eeg->data[0]->n, eeg->nmarkers_per_trial );

  ffwrite( eeg->times, sizeof(double), eeg->data[0]->n, f );
  
  /* write markers */
  for( i=0; i<eeg->ntrials; i++ ){
	 write_raw_markers( f, eeg->data[i] );
  }
			  
  for( c=0; c<eeg->data[0]->nbchan; c++ ){
	 for( t=0; t<eeg->ntrials; t++ ){
		ffwrite( eeg->data[t]->d[c], sizeof(double), eeg->data[0]->n, f );
	 }
  }
}
/** write a file in RAW-file format from an EEGdata_trials struct.
	 calls  write_eegtrials_to_raw()
 */
void write_eegtrials_to_raw_file( const EEGdata_trials *eeg, const char *fname ){
  FILE *out;
  if((out=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_eegtrials_to_raw( eeg, out );

  fclose( out );
}

void write_double_matrix_ascii_file(const char *fname, const double **d, int xdim, int ydim){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);
  
  write_double_matrix_ascii(f, d, xdim, ydim);

  fclose(f);
}

void write_double_matrix_ascii(FILE *out, const double **d, int xdim, int ydim){
  int x, y;

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fprintf(out, "%f ", d[x][y]);
	 }
	 fprintf(out, "\n");
  }
}

		  
void write_double_vector_ascii(FILE *out, const double *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%f ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_double_vector_ascii_file(const char *fname, const double *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_double_vector_ascii(f, v, n);
  fclose(f);
}


		  
void write_int_vector_ascii(FILE *out, const int *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%i ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_int_vector_ascii_file(const char *fname, const int *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_int_vector_ascii(f, v, n);
  fclose(f);
}

/** ignore markers, currently. This is handeled by the raw-format 
 */
void write_eegdata_ascii(FILE *out, const EEGdata *eeg ){
  int c, i;
  for( c=0; c<eeg->nbchan; c++ ){
	 for( i=0; i<eeg->n; i++ ){
		fprintf( out, "%f ", eeg->d[c][i] );
	 }
	 fprintf( out, "\n" );
  }
}

void write_eegdata_ascii_file(const char *fname, const EEGdata *eeg ){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_eegdata_ascii(f, eeg );
  fclose(f);
}


int _eegtrials_to_double_getn( EEGdata_trials *eeg ){
  int n=0;
  
  n += 4; 							  /* fix fields in struct */
  n += eeg->nmarkers_per_trial*eeg->ntrials; /* markers */
  n += eeg->nsamples;								/* times */
  n += 1;												/* nbchan */
  n += (eeg->data[0]->nbchan)*(eeg->nsamples)*(eeg->ntrials); /* data */

  return n;
}

/** pack a complete EEGdata_trials struct (with all sub-structs) 
	 into a single double-array. Can be used, e.g. for MPISend().
	 \param eeg
	 \param n - output (number values in return vector)
	 \return vector-representation of eeg
*/
double* eegtrials_to_double( EEGdata_trials *eeg, int *n ){
  int i, j, c;
  double *v;
  int idx = 0;

  *n = _eegtrials_to_double_getn( eeg );
  v = (double*) malloc( (*n)*sizeof(double) );

  /* fix values */
  v[idx++] = (double) eeg->ntrials;
  v[idx++] = (double) eeg->nsamples;
  v[idx++] = (double) eeg->nmarkers_per_trial;
  v[idx++] = (double) eeg->data[0]->nbchan;
  v[idx++] = eeg->sampling_rate;

  /* times */
  for( i=0; i<eeg->nsamples; i++ ){
	 v[idx++] = eeg->times[i];
  }

  /* markers */
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=0; j<eeg->nmarkers_per_trial; j++ ){
		v[idx++] = (double) eeg->markers[i][j];
	 }
  }

  /* EEG-data */
  for( i=0; i<eeg->ntrials; i++ ){
	 for( c=0; c<eeg->data[0]->nbchan; c++ ){
		for( j=0; j<eeg->nsamples; j++ ){
		  v[idx++] = eeg->data[i]->d[c][j];
		}
	 }
  }

  return v;
}

/** reverse the process from eegtrials_to_double() (that is, 
	 construct a struct out of a stream of doubles).
	 Can be used, e.g. for MPIRecv().
*/
EEGdata_trials* eegtrials_from_double( double *stream, int n ){
  EEGdata_trials *eeg;
  int i,j,c;
  unsigned int nbtrials, nmarkers, nbchan, nbsamples;
  int idx;

  nbtrials = (unsigned int) stream[0];
  nbsamples= (unsigned int) stream[1];
  nmarkers = (unsigned int) stream[2];
  nbchan   = (unsigned int) stream[3];
  eeg = init_eegdata_trials( nbtrials, nmarkers, nbchan, nbsamples, stream+5 );
  
  /* copy markers */
  idx = 5+nbsamples;				  /* after times array */
  for( i=0; i<nbtrials; i++ ){
	 for( j=0; j<nmarkers; j++ ){
		eeg->markers[i][j] = (unsigned long)stream[idx++];
	 }
	 memcpy( eeg->data[i]->markers, eeg->markers[i], nmarkers*sizeof(unsigned long) );
  }

  /* copy data */
  for( i=0; i<nbtrials; i++ ){
	 for( c=0; c<nbchan; c++ ){
		for( j=0; j<nbsamples; j++ ){
		  eeg->data[i]->d[c][j]=stream[idx++]; 
		}
	 }
  }
  
  return eeg;
}

/** pack a complete EEGdata struct
	 into a single double-array. Can be used, e.g. for MPISend().
	 \param eeg
	 \param n - output (number values in return vector)
	 \return vector-representation of eeg
*/
double* eegdata_to_double( EEGdata *eeg, int *n ){
  double *v;
  int i,c;
  int idx=0;
  
  *n = (3+eeg->nmarkers+(eeg->nbchan*eeg->n));
  v = (double*) malloc( (*n)*sizeof(double) );

  v[idx++]=(double)eeg->nbchan;
  v[idx++]=(double)eeg->n;
  v[idx++]=(double)eeg->nmarkers;

  /* markers */
  for( i=0; i<eeg->nmarkers; i++ ){
	 v[idx++]=(double)eeg->markers[i];
  }

  /* data */
  for( c=0; c<eeg->nbchan; c++ ){
	 for( i=0; i<eeg->n; i++ ){
		v[idx++]=eeg->d[c][i];
	 }
  }

  return v;
}

/** reverse the process from eegdata_to_double() (that is, 
	 construct a struct out of a stream of doubles).
	 Can be used, e.g. for MPIRecv().
*/
EEGdata* eegdata_from_double( double *stream, int n ){
  EEGdata *eeg;
  int i,c;
  int nmarkers, nbchan, nbsamples;
  int idx;

  nbchan   = (int) stream[0];
  nbsamples= (int) stream[1];
  nmarkers = (int) stream[2];

  eeg = init_eegdata( nbchan, nbsamples, nmarkers );
  
  /* copy markers */
  idx=3;
  for( i=0; i<nmarkers; i++ ){
	 eeg->markers[i] = (unsigned long)stream[idx++];
  }

  /* copy data */
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<nbsamples; i++ ){
		eeg->d[c][i]=stream[idx++]; 
	 }
  }

  return eeg;
}
