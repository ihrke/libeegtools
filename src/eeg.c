#include "eeg.h"
#include "helper.h"
#include "mathadd.h"

EEG* eeg_init            ( int nbchan, int ntrials, int nsamples ){
  EEG *eeg;
  int c,i,j;

  eeg = (EEG*)malloc( sizeof( EEG ) );
  eeg->nbchan = nbchan;
  eeg->ntrials= ntrials;
  eeg->n      = nsamples;
  eeg->data   = (double***) malloc( nbchan*sizeof( double** ));
  for( c=0; c<nbchan; c++ ){
	 eeg->data[c] = (double**) malloc( ntrials*sizeof(double**) );
	 for( i=0; i<ntrials; i++ ){
		eeg->data[c][i] = (double*) malloc( nsamples*sizeof(double) );
		for( j=0; j<nsamples; j++ ){
		  eeg->data[c][i][j] = 0.0;
		}
	 }
  }

  eeg->filename=NULL;

  eeg->comment=NULL;
  eeg->chaninfo=NULL;
  eeg->times=NULL;

  eeg->markers=NULL;
  eeg->nmarkers=NULL;
  eeg->marker_labels=NULL;

  return eeg;
}
/** allocate memory for markers.
	 \param nmarkers_per_trial
	 \param eeg
 */
EEG* eeg_init_markers    ( int nmarkers_per_trial, EEG *eeg ){
  int i;
  eeg->nmarkers= (unsigned int*) malloc( eeg->ntrials*sizeof(unsigned int) );
  eeg->markers = (unsigned int**)malloc( eeg->ntrials*sizeof(unsigned int*) );
  eeg->marker_labels=(char***)   malloc(eeg->ntrials*sizeof(char***) );
  for( i=0; i<eeg->ntrials; i++ ){
	 eeg->nmarkers[i] = nmarkers_per_trial;
	 eeg->markers[i] = (unsigned int*) malloc( nmarkers_per_trial*sizeof(unsigned int) );
	 eeg->marker_labels[i] = (char**)  malloc( nmarkers_per_trial*sizeof(char*) );
  }
  return eeg;
}

/** append text to the comment-field in the EEG-struct.
	 \param eeg the input struct
	 \param comment the text to append
 */
void eeg_append_comment( EEG *eeg, const char *comment ){
  if( !eeg->comment ){
	 eeg->comment = (char*) malloc( strlen(comment)*sizeof(char) );
	 strcpy( eeg->comment, comment );
  } else {
	 eeg->comment = (char*) realloc( eeg->comment, (strlen(eeg->comment)+strlen(comment))*sizeof(char) );
	 strcpy( eeg->comment+strlen(eeg->comment), comment );
  }
}

/** extract a list of channels from the EEG-struct eeg (data and channelinfo).
	 The resulting EEG struct is returned (either freshly allocated) or
	 overwritten in the struct.
	 \param eeg the input (and perhaps output) eeg
	 \param channels list of channels to extract
	 \param nchannels number of channels in channels
	 \param alloc TRUE: allocate new memory and return; FALSE: delete channels that have 
	               not been selected from the input dataset
 */
EEG* eeg_extract_channels( EEG* eeg, const int *channels, int nchannels, Boolean alloc ){
  EEG *outeeg;
  int c, i;
  int minus_nchannels;


  if( alloc ){
	 outeeg = eeg_clone( eeg, EEG_CLONE_ALL );
  } else {
	 outeeg = eeg;
  }

  minus_nchannels = nchannels;
  for( c=0; c<nchannels; c++ ){ /* check channel list */
	 if( channels[c]<0 || channels[c]>=outeeg->nbchan ){
		warnprintf("Skipping invalid channel '%i'\n", channels[c] );
		minus_nchannels--;
		continue;
	 }
  }

  for( c=0; c<outeeg->nbchan; c++ ){					 /* modify */
	 if( !isin_intarray( channels, nchannels, c ) ){ 
		if( outeeg->data ){
		  dprintf("free data for channel '%i'\n", c );
		  for( i=0; i<outeeg->ntrials; i++ )
			 safer_free( outeeg->data[c][i] );
		  safer_free( outeeg->data[c] );
		  outeeg->data[c] = NULL;
		}

		if( outeeg->chaninfo ){
		  dprintf("mark chaninfo for channel '%i' for deletion\n", c );
		  outeeg->chaninfo[c].num=-1;
		}
	 }
  }

  /* move empty fields such that there are no gaps in data and chaninfo */
  Boolean moved;
  for( c=0; c<outeeg->nbchan; c++ ){
	 dprintf("c=%i, data=%p, chan.num=%i\n",c, outeeg->data[c], outeeg->chaninfo[c].num);
	 moved=FALSE;
	 if( !outeeg->data[c] ){	  /* move all remaining data */	
		for( i=c+1; i<outeeg->nbchan; i++ ){
		  outeeg->data[i-1] = outeeg->data[i];
		}
		moved=TRUE;
	 }

	 if( outeeg->chaninfo[c].num==-1 ){ /* move remaining chaninfos */
		for( i=c+1; i<outeeg->nbchan; i++ ){
		  memcpy( &(outeeg->chaninfo[i-1]), &(outeeg->chaninfo[i]), sizeof(ChannelInfo) );
		}
		moved=TRUE;
	 }
	 if( moved ) {
		c--;
		outeeg->nbchan--;
	 }
  }

  outeeg->nbchan=minus_nchannels;

  return outeeg;
}

/** Check the eeg-struct for inconsitencies.
	 \warning This function can result in a segmentation fault.
	           The function tries to access all data that must be 
				  available due to the parameter settings within the struct.
				  If this it is not the case, the function will segfault.
				  This is useful to run before doing large calculation.
				  You will know that it will fail a lot earlier then :-)
	 \param eeg input eeg
	 \param flags one of EEG_CHECK_*
	 \return 0 if ok; something else if not															\
*/
int eeg_check( EEG *eeg, int flags ){

}

/** extract a list of trials from the EEG-struct eeg (data markers).
	 The resulting EEG struct is returned (either freshly allocated) or
	 overwritten in the struct.
	 \param eeg the input (and perhaps output) eeg
	 \param trials list of trials to extract
	 \param ntrials number of trials in trials-array
	 \param alloc TRUE: allocate new memory and return; FALSE: delete trials that have 
	               not been selected from the input dataset
 */
EEG* eeg_extract_trials  ( EEG* eeg, const int *trials,   int ntrials,   Boolean alloc ){
  EEG *outeeg;
  int c,i,j;
  int minus_ntrials;

  if( alloc ){
	 outeeg = eeg_clone( eeg, EEG_CLONE_ALL );
  } else {
	 outeeg = eeg;
  }

  minus_ntrials = ntrials;
  for( i=0; i<ntrials; i++ ){ /* check trial list */
	 if( trials[i]<0 || trials[i]>=outeeg->ntrials ){
		warnprintf("Skipping invalid Trial '%i'\n", trials[i] );
		minus_ntrials--;
		continue;
	 }
  }

  for( i=0; i<outeeg->ntrials; i++ ){
	 if( !isin_intarray( trials, ntrials, i) ){ /* remove */
		dprintf("Removing data for trial '%i'\n", i);
		for( c=0; c<outeeg->nbchan; c++ ){
		  safer_free( outeeg->data[c][i] );
		  outeeg->data[c][i]=NULL;
		}

		dprintf("Removing markers for trial '%i'\n", i);
		if( outeeg->markers ){
		  safer_free( outeeg->markers[i] );
		  outeeg->markers[i] = NULL;
		}

		if( outeeg->marker_labels ){
		  for( j=0; j<outeeg->nmarkers[i]; j++ ){
			 safer_free( outeeg->marker_labels[i][j] );
		  }
		  safer_free( outeeg->marker_labels[i] );
		  outeeg->marker_labels[i]=NULL;
		}
		if(outeeg->nmarkers )
		  outeeg->nmarkers[i] = 0;
	 }
  }

  /* shifting data/markers */
  Boolean moved;
  for( i=0; i<outeeg->ntrials; i++ ){
	 moved=FALSE;
	 for( c=0; c<outeeg->nbchan; c++ ){
		if( !outeeg->data[c][i] ){	  /* move all remaining data */	
		  for( j=i+1; j<outeeg->ntrials; j++ ){
			 outeeg->data[c][j-1] = outeeg->data[c][j];
		  }
		  moved=TRUE;
		}
	 }

	 if( outeeg->markers[i]==NULL ){
		for( j=i+1; j<outeeg->ntrials; j++ ){
		  outeeg->markers[j-1] = outeeg->markers[j];
		  outeeg->nmarkers[j-1] = outeeg->nmarkers[j];
		  if( outeeg->marker_labels ){
			 outeeg->marker_labels[j-1] = outeeg->marker_labels[j];
		  }
		}
		moved=TRUE;
	 }

	 if( moved ) {
		i--;
		outeeg->ntrials--;
	 }
  }

  outeeg->ntrials = minus_ntrials;
  return outeeg;
}

/** This function free's all memory in EEG-struct that it knows about,
	 including the eeg pointer itself.
 */
void eeg_free( EEG *eeg ){
  int c,i,j;
  safer_free( eeg->filename );
  safer_free( eeg->comment );
  safer_free( eeg->times );
  if( eeg->data ){
	 for( c=0; c<eeg->nbchan; c++ ){
		for( i=0; i<eeg->ntrials; i++ ){
		  safer_free( eeg->data[c][i] );
		}
		safer_free( eeg->data[c] );
	 }
	 safer_free( eeg->data );
  }

  if( eeg->markers ){
	 for( i=0; i<eeg->ntrials; i++ ){
		safer_free( eeg->markers[i] );
		if( eeg->marker_labels ){
		  for( j=0; j<eeg->nmarkers[i]; j++ ){
			 safer_free( eeg->marker_labels[i][j] );
		  }
		  safer_free( eeg->marker_labels[i] );
		}
	 }
	 safer_free( eeg->marker_labels );
	 safer_free( eeg->markers );
  }
  safer_free( eeg->nmarkers );
  safer_free( eeg->chaninfo );
  safer_free( eeg );
}

/** Pretty-print an EEG - struct.
	 \param out the output stream
	 \param eeg the struct
	 \param preview number of items to preview in data
*/
void eeg_print( FILE *out, const EEG *eeg, int preview ){
  int c,i,j;
  fprintf(out, 
			 "EEG:\n"
			 " filename      = '%s'\n"
			 " comment       = '%s'\n"
			 " nbchan        = %i\n"
			 " ntrials       = %i\n"
			 " n             = %i\n"
			 " sampling_rate = %f\n",
			 (eeg->filename)?(eeg->filename):"<NULL>", 
			 (eeg->comment)?(eeg->comment):"<NULL>",
			 eeg->nbchan, eeg->ntrials,
			 eeg->n, eeg->sampling_rate);  
  if( eeg->times!=NULL ){
	 fprintf(out, 
			 " times[0]      = %f\n", eeg->times[0]);
	 fprintf(out,
			 " times[n-1]    = %f\n", eeg->times[eeg->n-1]);
  } else {
	 fprintf(out, 
			 " times[0]      = <NULL>\n");
  }
  fprintf(out,   
			 " chaninfo      = \n");
  if( eeg->chaninfo ){
	 for( i=0; i<eeg->nbchan; i++ ){
		fprintf( out, "   ");
		print_channelinfo( out,  &(eeg->chaninfo[i]) );
	 }
  } else {
	 fprintf( out, "   <NULL>\n");
  }
  fprintf( out, " data = \n" );
  if( eeg->data!=NULL ){
	 for( c=0; c<MIN(eeg->nbchan,preview); c++ ){
		for( i=0; i<MIN(eeg->ntrials, preview); i++ ){
		  for( j=0; j<MIN(eeg->n, preview); j++ ){
			 fprintf(out, 
						"  [%i][%i][%i] = %f\n", c,i,j,eeg->data[c][i][j] );
		  }
		}
	 }
  } else {
	 fprintf(out, 
			 " data          = <NULL>\n");
  }  

  fprintf(out, " nmarkers      = \n");  
  if( eeg->nmarkers ){
	 for( i=0; i<MIN(eeg->ntrials, preview); i++ ){
		fprintf(out, 
				  "  [%i]  = %i\n", i, eeg->nmarkers[i] );
	 }
  } else {
	 fprintf(out, "  <NULL>\n");
  }
  fprintf(out, " markers       = \n");
  if( eeg->markers ){
	 for( i=0; i<MIN( eeg->ntrials, preview ); i++ ){
		for( j=0; j<MIN( eeg->nmarkers[i], preview ); j++ ){
		  fprintf(out, "  [%i][%i] = %i\n", i,j, eeg->markers[i][j] );
		}
	 }
  } else {
	 fprintf(out, "  <NULL>\n");
  }

  fprintf(out, " marker_labels = \n");
  if( eeg->marker_labels ){
	 for( i=0; i<MIN(eeg->ntrials, preview ); i++ ){
		for( j=0; j<MIN(eeg->nmarkers[i], preview ); j++ ){
		  fprintf(out, "  [%i][%i] = %s\n", i,j, eeg->marker_labels[i][j] );
		}
	 }
  } else {
	 fprintf(out,
			 "  <NULL>\n");
  }

}

/** deep-copy the complete EEG-struct and return
	 a freshly allocated one.
	 \param eeg to clone
	 \param flags an 'or' combination of EEG_CLONE_*
	 \return cloned struct
*/
EEG* eeg_clone( const EEG *eeg, int flags ){
  EEG *clone;
  int c, i, j;

  clone = eeg_init( eeg->nbchan, eeg->ntrials, eeg->n );
  if( eeg->filename ){
	 clone->filename = (char*)  malloc( strlen( eeg->filename )*sizeof(char) );
	 strcpy( clone->filename, eeg->filename );
  }
  if( eeg->comment ){
	 clone->comment = (char*)  malloc( strlen( eeg->comment )*sizeof(char) );
	 strcpy( clone->comment, eeg->comment );
  }
  clone->sampling_rate = eeg->sampling_rate;
  if( eeg->times ){
	 clone->times = (double*) malloc( eeg->n*sizeof(double) );
	 memcpy( clone->times, eeg->times, eeg->n*sizeof(double) );
  }

  if( !(flags & EEG_CLONE_NODATA) ){
	 for( c=0; c<eeg->nbchan; c++ ){
		for( i=0; i<eeg->ntrials; i++ ){
		  for( j=0; j<eeg->n; j++ ){
			 clone->data[c][i][j] = eeg->data[c][i][j];
		  }
		}
	 }
  }
  if( !(flags & EEG_CLONE_NOCHANINFO) ){
	 if( eeg->chaninfo ){
		clone->chaninfo = (ChannelInfo*) malloc( eeg->ntrials*sizeof(ChannelInfo) );
		for( i=0; i<eeg->ntrials; i++ ){
		  memcpy( &(clone->chaninfo[i]), &(eeg->chaninfo[i]), sizeof(ChannelInfo) );
		}
	 }
  }
  if( !(flags & EEG_CLONE_NOMARKERS) ){
	 if( eeg->nmarkers ){
		clone->nmarkers = (unsigned int*)malloc( eeg->ntrials*sizeof(unsigned int) );
		memcpy( clone->nmarkers, eeg->nmarkers, eeg->ntrials*sizeof(unsigned int) );
	 }

	 if( eeg->markers ){
		clone->markers = (unsigned int**)malloc( eeg->ntrials*sizeof(unsigned int*) );
		for( i=0; i<eeg->ntrials; i++ ){
		  clone->markers[i] =(unsigned int*)   malloc( eeg->nmarkers[i]*sizeof(unsigned int) );
		  memcpy( clone->markers[i], eeg->markers[i],  eeg->nmarkers[i]*sizeof(unsigned int) );
		}
	 }

	 if( eeg->marker_labels ){
		clone->marker_labels = (char ***)malloc( eeg->ntrials*sizeof(char**) );
		for( i=0; i<eeg->ntrials; i++ ){
		  clone->marker_labels[i] =(char**)malloc( eeg->nmarkers[i]*sizeof(char*) );
		  for( j=0; j<eeg->nmarkers[i]; j++ ){
			 clone->marker_labels[i][j] = (char*) malloc( strlen( eeg->marker_labels[i][j])*sizeof(char) );
			 strcpy( clone->marker_labels[i][j], eeg->marker_labels[i][j] );
		  }
		}
	 }
  }

  return clone;
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

int       eegdata_cmp_settings( const EEGdata *s1, const EEGdata *s2 ){
  if( s1->nbchan != s2->nbchan || s1->n != s2->n || s1->nmarkers != s2->nmarkers ){
	 return 1;
  } else {
	 return 0;
  }
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


void print_channelinfo( FILE* out, const ChannelInfo *c ){
  fprintf( out, "Channel No. %i/%i - '%s' at (%.2f, %.2f, %.2f)\n", 
			  c->num, c->num_chans, c->label, c->x, c->y, c->z );
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
EEGdata*         clone_eegdata( const EEGdata *source ){
  EEGdata *eeg;
  eeg = init_eegdata( source->nbchan, source->n, source->nmarkers );
  copy_similar_eegdata( eeg, source );
  return eeg;
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
