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
double** read_double_dblpp_ascii(const char *fname, int xdim, int ydim, double **d){
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

/** Reads EEG-data from binary file. Format is Cxn doubles.
 * \param C - number of channels
 * \param n - number of samples per channel
 */
EEG* read_continuous_eeg_from_binfile(const char *file, int C, int n){
  int c;
  EEG *eeg;
  FILE *f;
  
  if((f = fopen(file, "rb"))==NULL) return NULL;

  eeg = eeg_init( C, 1, n );
  for(c=0; c<C; c++){
	 if(fread(eeg->data[c][0], sizeof(double), n, f)<n){
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
EEG* read_eeg_from_raw(const char *file){
  EEG* eeg;
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
  eeg = eeg_init( nbchan, ntrials, nsamples );
  eeg->times = (double*) malloc( nsamples*sizeof(double) );
  eeg->filename = (char*) malloc( (strlen( file )+1)*sizeof(char) );
  strcpy( eeg->filename, file );

  /* read times-array */
  ffread(eeg->times, sizeof(double), nsamples, f);
  eeg->sampling_rate=1000.0/(eeg->times[1]-eeg->times[0]);

  /* read markers */
  double *tmp;
  tmp = (double*) malloc( nmarkers * sizeof(double) );
  eeg->nmarkers = (unsigned int*) malloc( eeg->ntrials*sizeof(unsigned int) );
  eeg->markers  = (unsigned int**)malloc( eeg->ntrials*sizeof(unsigned int*) );

  for( i=0; i<ntrials; i++ ){
	 eeg->nmarkers[i] = nmarkers;
	 eeg->markers[i] = (unsigned int*) malloc( nmarkers*sizeof(unsigned int) );
	 ffread(tmp, sizeof(double), nmarkers, f);
	 /* dprintf("tmp[0]=%f, tmp[1]=%f\n", tmp[0], tmp[1]); */
	 for( j=0; j<nmarkers; j++){
		eeg->markers[i][j] = (unsigned int)tmp[j];
	 }
  }
  free(tmp);

  /* read data */
  for( c=0; c<nbchan; c++ ){
	 for( i=0; i<ntrials; i++ ){
		/*		dprintf("eeg->data[%i]->n = %i\n", i, eeg->data[i]->n);*/
		ffread( eeg->data[c][i], sizeof( double ), nsamples, f );
	 }
  }

  return eeg;
}




/* ------------- WRITER --------------------- */

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
