/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke   *
 *   ihrke@nld.ds.mpg.de
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "io.h"
#include "linalg.h"

/* ------------- READER --------------------- */

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

/**\brief Read a matrix from blank-separated file.
	
	The format is something like this:
	\verbatim
-0.005344   0.696318   0.005792   589.041211   
-0.005348   0.691603   0.007035   589.074254   
-0.005352   0.686796   0.008324   589.096248   
-0.005516   0.680181   0.010164   589.118214   
...
	\endverbatim
	
	\param fname the name of the file
	\return the matrix
 */
Array* read_matrix_from_text( const char *fname ){ 
  char buf[MAX_LINE_LENGTH];
  int nrows=0;
  int ncols=0;
  FILE *f;

  if( (f=fopen( fname, "r" ))==NULL ){
	 errprintf("couldn't open file '%s'\n", fname );
	 return NULL;
  }
  
  nrows=stream_count_char( f, '\n' );
  rewind( f );
  read_line( f, buf );
  char c,prev='1';
  while( (c=fgetc(f))!='\n' ){
	 if( isblank(c) && !isblank(prev) ){
		ncols++;
	 }
	 prev=c; 
  }
  rewind( f );
  dprintf( "matrix is (%i,%i)\n", nrows, ncols );

  Array *a=array_new2( DOUBLE, 2, nrows, ncols );
  int i,j,k;
  for( i=0; i<nrows; i++ ){
	 for( k=0; k<ncols; k++ ){
		j=0;
		while( !isblank((c=fgetc(f))) ){
		  buf[j++]=c;
		} 
		buf[j]='\0';
		dprintf("%s\n", buf );
		mat_IDX( a, i, k )=atof( buf );
		while( isblank((c=fgetc(f))) ) ;
		fseek( f, -1, SEEK_CUR);
	 }
  }
  dprintf("Done\n");

  fclose( f );
  return a;
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
double** read_dblpp_ascii(const char *fname, int xdim, int ydim, double **d){
  FILE *f;
  int x,y;
  int flag;

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
		flag=fscanf(f, " %lf ", &(d[y][x]));
	 }
  }
  fclose(f);
  return d;
}

/** read 1D vector frome file 
 */
double* read_dblp_ascii( const char *fname, int N, double *v ){
  FILE *f;
  int i,flag;

  if((f=fopen(fname, "r"))==NULL)
	 errormsg(ERR_IO, 1);

  if(v==NULL){
	 v = dblp_init( NULL, N, 0.0 );
  }
  for( i=0; i<N; i++ ){
	 flag=fscanf(f, " %lf ", &(v[i]));
  }
  fclose( f );
  return v;
}



/* ------------- WRITER --------------------- */

/**
 write double matrix to file in ASCII format.
	 \param opts may contain
	 - "precision=int" number of significant digits (after comma); default=6;
*/
void write_dblpp_ascii_file(const char *fname, const double **d, int xdim, int ydim, OptArgList *opts){
  FILE *f;

  dprintf("opening '%s'\n", fname );
  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);
  
  write_dblpp_ascii(f, d, xdim, ydim, opts);

  fclose(f);
}

/** write double matrix to file in ASCII format.
	 \param opts may contain
	 - "precision=int" number of significant digits (after comma); default=6;
 */
void write_dblpp_ascii(FILE *out, const double **d, int xdim, int ydim, OptArgList *opts){
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

		  
void write_dblp_ascii(FILE *out, const double *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%f ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_dblp_ascii_file(const char *fname, const double *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_dblp_ascii(f, v, n);
  fclose(f);
}


		  
void write_intp_ascii(FILE *out, const int *v, int n){
  int x;
  for( x=0; x<n; x++ ){
	 fprintf(out, "%i ", v[x]);
  }
  fprintf(out, "\n");
}


		  
void write_intp_ascii_file(const char *fname, const int *v, int n){
  FILE *f;

  if((f=fopen(fname, "w"))==NULL)
	 errormsg(ERR_IO, 1);

  write_intp_ascii(f, v, n);
  fclose(f);
}
