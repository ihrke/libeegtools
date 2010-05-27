#include "writer.h"
#include "eeg.h"
#include "optarg.h"
#include "linalg.h"
#include <math.h>

#ifdef MATIO
#include <matio.h>
/** This write uses MatIO to create an EEGlab file.
	 It was developed with EEGlab version 6.01b.
	 
	 \param eeg the struct
	 \param file eeglab .set file
	 \return error code
*/
int write_eeglab_file( EEG* eeg, const char *file ){
  mat_t *mfile;
  matvar_t *meeg=NULL;
  if( !(mfile=Mat_Create( file, NULL )) ){
	 errprintf("Could not open '%s' for writing\n", file);
	 return -1;
  }

  //int dims[2]={1,1}, rank=2;
  //meeg = Mat_VarCreate( "EEG", MAT_T_STRUCT, , rank, dims, NULL, 0 );
  
  Mat_VarWrite( mfile, meeg, 0 );
  return 0;
}

/** \brief save Array-struct as MATLAB .mat file.

	 Currently, the function converts all arrays to
	 double-arrays. Should be easy to implement for arbitrary types if
	 needed.

	 \param a the Array
	 \param varname MATLAB-name of the array
	 \param file the filename
	 \param append if TRUE, the function tries to open() the file, else it is overwritten
	 \return error code
 */
int write_array_matlab( const Array *a, const char *varname, const char *file, bool append ){
  mat_t *mfile;
  matvar_t *marr=NULL;
  int i;
  
  if( append ){
	 if( !(mfile=Mat_Open( file, MAT_ACC_RDWR )) ){
		errprintf("Could not open '%s', creating new file\n", file);
	 } 
  }
  
  if( !mfile ){
	 if( !(mfile=Mat_Create( file, NULL )) ){
		errprintf("Could not open '%s' for writing\n", file);
		return -1;
	 }
  }

  int ndim = MAX(2, a->ndim);
  int *size = (int*)malloc( ndim*sizeof(int));
  if( a->ndim==1 ){
	 size[0]=1;
	 size[1]=a->ndim;
  } else {
	 memcpy( size, a->size, ndim*sizeof(int));
  }
  dprintf("Writing to file '%s' variable '%s', ndim=%i\n", file, varname, ndim);

  /* convert to column major for MATLAB */
  Array *b=array_convert_rowcolmajor( (Array*)a, TRUE );
  
  /* up-cast to DOUBLE - copy of b */
  Array *c=array_new( DOUBLE, b->ndim, b->size );
  for( i=0; i<array_NUMEL(b); i++ ){
	 array_dtype_to_double( array_INDEXMEM1(c,i), array_INDEXMEM1(b,i), b->dtype );
  }
  array_free( b );

  marr = Mat_VarCreate( varname, MAT_C_DOUBLE, MAT_T_DOUBLE, 
								ndim, size, c->data, MEM_CONSERVE /* Array remains owner */
								);
  
  dprintf("mfile=%p, marr=%p\n", mfile, marr );
  int succ=Mat_VarWrite( mfile, marr, 0 );
  dprintf("done writing with succ=%i\n", succ);

  Mat_Close( mfile );
  Mat_VarFree( marr );
  array_free( c );
  free( size );

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

/**
 write double matrix to file in ASCII format.
	 \param opts may contain
	 - "precision=int" number of significant digits (after comma); default=6;
*/
void write_double_dblpp_ascii_file(const char *fname, const double **d, int xdim, int ydim, OptArgList *opts){
  FILE *f;

  dprintf("opening '%s'\n", fname );
  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);
  
  write_double_dblpp_ascii(f, d, xdim, ydim, opts);

  fclose(f);
}

/** write double matrix to file in ASCII format.
	 \param opts may contain
	 - "precision=int" number of significant digits (after comma); default=6;
 */
void write_double_dblpp_ascii(FILE *out, const double **d, int xdim, int ydim, OptArgList *opts){
  int x, y;
  int precision;
  char tformat[20];
  double a;

  precision=6;
  if( optarglist_has_key( opts, "precision" ) ){
	 a = optarglist_scalar_by_key( opts, "precision" );
	 if( !isnan( a ) ) precision=(int)a;
  }
  
  sprintf( tformat, "%%.%if ", precision );
  dprintf("Printing in format: '%s'\n", tformat);

  for( y=0; y<ydim; y++ ){
	 for( x=0; x<xdim; x++ ){
		fprintf(out, tformat, d[x][y]);
	 }
	 fprintf(out, "\n");
  }
}

		  
void write_double_dblp_ascii(FILE *out, const double *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%f ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_double_dblp_ascii_file(const char *fname, const double *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_double_dblp_ascii(f, v, n);
  fclose(f);
}


		  
void write_int_dblp_ascii(FILE *out, const int *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%i ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_int_dblp_ascii_file(const char *fname, const int *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_int_dblp_ascii(f, v, n);
  fclose(f);
}




