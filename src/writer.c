#include "writer.h"



void write_eegtrials_to_raw( const EEGdata_trials *eeg, FILE *f ){
  double nbchan_d, ntrials_d, nsamples_d, nmarkers_d;
  int i,j;

  nbchan_d   = (double)eeg->data[0]->nbchan;
  ntrials_d  = (double)eeg->ntrials;
  nsamples_d = (double)eeg->data[0]->n;
  nmarkers_d = (double)eeg->nmarkers_per_trial;

  ffwrite( &nbchan_d,   sizeof(double), 1, f );  
  ffwrite( &ntrials_d,  sizeof(double), 1, f );
  ffwrite( &nsamples_d, sizeof(double), 1, f );
  ffwrite( &nmarkers_d, sizeof(double), 1, f );

  ffwrite( eeg->times, sizeof(double), eeg->data[0]->n, f );
  
  /* write markers */
  double *tmp;
  tmp = (double*) malloc( eeg->nmarkers_per_trial * sizeof(double) );
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=0; j<eeg->nmarkers_per_trial;  j++ ){ /* typecast */
		tmp[j] = (double)eeg->markers[i][j];
	 }
	 ffwrite( tmp, sizeof(double), eeg->nmarkers_per_trial, f );
  }
  free( tmp );
			  
  for( c=0; c<eeg->data[0]->nbchan; c++ ){
	 for( t=0; t<eeg->ntrials; t++ ){
		ffwrite( eeg->data[t]->d[c], sizeof(double), eeg->data[0]->n, f );
	 }
  }
}

void write_eegtrials_to_raw_file( const EEGdata_trials *eeg, const char *file ){
  FILE *out;
  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_eegtrials_to_raw( eeg, out );

  fclose( out );
}

void write_double_matrix_ascii_file(const char *fname, const double **d, int xdim, int ydim){
  FILE *f;
  int x, y;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_double_matrix_ascii(f, d, xdim, ydim);

  fclose(f);
}

void write_double_matrix_ascii(FILE *out, const double **d, int xdim, int ydim){
  int x, y;

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fprintf(out, "%f ", d[y][x]);
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
