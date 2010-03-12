#include "eeg.h"
#include "helper.h"
#include "mathadd.h"


/** in-depth comparison of two eeg-structs;
	 \todo finish!
 */
bool eeg_cmp_depth( const EEG *eeg1, const EEG *eeg2 ){
  if( strcmp( eeg1->filename, eeg2->filename ) ){
	 dprintf("Filenames to not match: %s != %s\n", eeg1->filename,
				eeg2->filename );
	 return FALSE;
  }
  if( strcmp( eeg1->comment, eeg2->comment ) ){
	 dprintf("Comments to not match: %s != %s\n", eeg1->comment,
				eeg2->comment );
	 return FALSE;
  }
  if( eeg1->nbchan!=eeg2->nbchan ||
		eeg1->ntrials!=eeg2->ntrials ||
		eeg1->n!=eeg2->n ){
	 dprintf("Main params do not match: (%i,%i,%i) != (%i,%i,%i)\n",
				eeg1->nbchan, eeg2->nbchan,
				eeg1->ntrials, eeg2->ntrials,
				eeg1->n, eeg2->n );
	 return FALSE;
  }
  if( cmpdouble( eeg1->sampling_rate, eeg2->sampling_rate, 2 ) ){
	 dprintf("sampling rate is different: %f!=%f\n", eeg1->sampling_rate,
				eeg2->sampling_rate );
  }
  /* TODO */

  return TRUE;
}

EEG* eeg_init            ( int nbchan, int ntrials, int nsamples ){
  EEG *eeg;
  int c,i,j;

  dprintf("init\n");
  eeg = (EEG*)malloc( sizeof( EEG ) );
  eeg->nbchan = nbchan;
  eeg->ntrials= ntrials;
  eeg->n      = nsamples;
  eeg->sampling_rate = -1;
  
  dprintf("alloc data\n");
  eeg->data   = (double***) malloc( nbchan*sizeof( double** ));
  dprintf("done\n");
  for( c=0; c<nbchan; c++ ){
	 eeg->data[c] = (double**) malloc( ntrials*sizeof(double*) );
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
  int i, j;
  eeg->nmarkers= (unsigned int*) malloc( eeg->ntrials*sizeof(unsigned int) );
  eeg->markers = (unsigned int**)malloc( eeg->ntrials*sizeof(unsigned int*) );
  eeg->marker_labels=(char***)   malloc(eeg->ntrials*sizeof(char**) );
  for( i=0; i<eeg->ntrials; i++ ){
	 eeg->nmarkers[i] = nmarkers_per_trial;
	 eeg->markers[i] = (unsigned int*) malloc( nmarkers_per_trial*sizeof(unsigned int) );
	 eeg->marker_labels[i] = (char**)  malloc( nmarkers_per_trial*sizeof(char*) );
	 for( j=0; j<nmarkers_per_trial; j++ ){
		eeg->marker_labels[i][j] = (char*)malloc( MAX_LABEL_LENGTH*sizeof(char) );
		sprintf( eeg->marker_labels[i][j], "Marker %i", j );
	 }
  }
  return eeg;
}

/** append text to the comment-field in the EEG-struct.
	 \param eeg the input struct
	 \param comment the text to append
 */
void eeg_append_comment( EEG *eeg, const char *comment ){
  if( !eeg->comment ){
	 eeg->comment = (char*) malloc( (strlen(comment)+1)*sizeof(char) );
	 strcpy( eeg->comment, comment );
  } else {
	 eeg->comment = (char*) realloc( eeg->comment, (strlen(eeg->comment)+strlen(comment)+2)*sizeof(char) );
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
EEG* eeg_extract_channels( EEG* eeg, const int *channels, int nchannels, bool alloc ){
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
  bool moved;
  for( c=0; c<outeeg->nbchan; c++ ){
	 moved=FALSE;
	 if( !outeeg->data[c] ){	  /* move all remaining data */	
		for( i=c+1; i<outeeg->nbchan; i++ ){
		  outeeg->data[i-1] = outeeg->data[i];
		}
		moved=TRUE;
	 }

	 if( outeeg->chaninfo ){
		if( outeeg->chaninfo[c].num==-1 ){ /* move remaining chaninfos */
		  for( i=c+1; i<outeeg->nbchan; i++ ){
			 memcpy( &(outeeg->chaninfo[i-1]), &(outeeg->chaninfo[i]), 
						sizeof(ChannelInfo) );
		  }
		  moved=TRUE;
		}
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
				  \todo implement it!
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
EEG* eeg_extract_trials  ( EEG* eeg, const int *trials,   int ntrials,   bool alloc ){
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
  bool moved;
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

/** Pretty-print an EEG - struct to a string.
	 \todo implement this!
	 \param out the output string or ALLOC_IN_FCT
	 \param eeg the struct
	 \param preview number of items to preview in data
	 \return memory for output string
*/
char* eeg_sprint( char *out, const EEG *eeg, int preview ){
  
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
  dprintf("start init\n");
  clone = eeg_init( eeg->nbchan, eeg->ntrials, eeg->n );
  dprintf("...done\n");
  if( eeg->filename ){
	 clone->filename = (char*)  malloc( (strlen( eeg->filename )+1)*sizeof(char) );
	 strcpy( clone->filename, eeg->filename );
  }
  if( eeg->comment ){
	 clone->comment = (char*)  malloc( (strlen( eeg->comment )+1)*sizeof(char) );
	 strcpy( clone->comment, eeg->comment );
  }
  clone->sampling_rate = eeg->sampling_rate;
  if( eeg->times ){
	 clone->times = (double*) malloc( eeg->n*sizeof(double) );
	 memcpy( clone->times, eeg->times, eeg->n*sizeof(double) );
  }
  dprintf("stuff done\n");

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
			 clone->marker_labels[i][j] = (char*) malloc( (strlen( eeg->marker_labels[i][j])+1)*sizeof(char) );
			 strcpy( clone->marker_labels[i][j], eeg->marker_labels[i][j] );
		  }
		}
	 }
  }

  return clone;
}



void print_channelinfo( FILE* out, const ChannelInfo *c ){
  fprintf( out, "Channel No. %i/%i - '%s' at (%.2f, %.2f, %.2f)\n", 
			  c->num, c->num_chans, c->label, c->x, c->y, c->z );
}

