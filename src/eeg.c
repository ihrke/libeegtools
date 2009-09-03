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

EEG* eeg_extract_channels( EEG* eeg, int *channels, int nchannels, Boolean alloc );
EEG* eeg_extract_trials  ( EEG* eeg, int *trials,   int ntrials,   Boolean alloc );

void eeg_free( EEG *eeg ){

}

/** Pretty-print an EEG - struct.
	 \param out the output stream
	 \param eeg the struct
	 \param preview number of items to preview in data
*/
void eeg_print( FILE *out, EEG *eeg, int preview ){
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
	 fprintf( out, "   <NULL>");
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
EEG* eeg_clone( EEG *eeg ){

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
