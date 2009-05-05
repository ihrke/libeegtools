/**\file ml_padtw.c
 * \brief Matlab wrapper for PADTW
 *
 *   Compilation:
 *     MATLAB will need to see libGSL and libeegtools.
 *     Put these (or similar) lines to the bottom of
 *     ~/.matlab/<version>/mexopts.sh
 \code
 *          CFLAGS="$CFLAGS -I/home/ihrke/local/include"
 *	         LD="$LD -lm -lgsl -lgslcblas -lplot"
 *          LDFLAGS="$LDFLAGS -L/home/ihrke/local/lib"
 *    
 *    then: 
 *     mex -v ml_padtw.c -leegtools
 *   
 \endcode
 */
#include "mex.h"
#include <stdlib.h>
#include "warping.h"
#include "helper.h"


void progressbar_matlab( int flag, int num );
void settings_from_matlabstruct( SettingsHierarchicalDTW *s, const mxArray *in );


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  char *docstring = 
	 "This function is part of the EEGlab-Plugin for the removal of\n"
	 "temporal distortion. The function is called thus:\n"
	 "[padtw, markers] = ml_padtw( X, markers, settings );\n\n"
	 "Here, X is the input matrix, which is 3-dimensional, \n"
	 "(channel x sample x trial).\n"
	 "markers is a (trials x nmarkers) matrix, containing the time-markers\n"
	 "in sampling units, used for the regularization.\n"
	 "settings is a struct, containing the configuration for \n"
	 "the hierarchical averaging procedure:\n"
	 " -> fields:\n"
	 "    - regularize in {'none','gaussian_markers'}\n"
	 "    - sigma in [0, 1]\n"
	 "    - linkage in {single, complete, average}\n"
	 "    - pointdistance in {euclidean, euclidean_derivative, stft}\n"
	 "    - windowfct in { dirichlet, gaussian, hamming, hanning, kaiser }\n"
	 "    - theta has two entries (theta1,theta2) with theta in [0, 1]\n"
	 "    - corner_freqs in Hz, two entries (bottom, top cutoff)\n"
	 "    - samling_rate in Hz\n"
	 "    - winlength (integer)\n"
	 "    - N_freq (integer)\n"
	 "    - N_time (integer)\n"
	 "    - \n"
	 "\n"
	 "The return values are:\n"
	 " - padtw - the final average\n"
	 " - markers - the new markers after temporal averaging\n";
  char *msg;
  msg = (char*)mxCalloc( strlen(docstring)+500, sizeof(char) );

  /* check proper input and output */
  if(nrhs<2){
    sprintf(msg, "Need at least 2 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    sprintf(msg, "First and second Input must be double arrays.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if( nrhs>2 && !mxIsStruct(prhs[2]) ){
	 sprintf(msg, "Third Input must be a struct.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if( mxGetNumberOfDimensions(prhs[0])!=3 ){
    sprintf(msg, "Data must be three-dimensional, got \n%s\n", 
				mxGetNumberOfDimensions(prhs[0]), docstring);
    mexErrMsgTxt(msg);
  } else if( mxGetNumberOfDimensions(prhs[1])!=2 ){
    sprintf(msg, "Marker-array must be two-dimensional, got \n%s\n", 
				mxGetNumberOfDimensions(prhs[1]), docstring);
    mexErrMsgTxt(msg);
  }

  double *dptr, *curptr, *mptr;
  int index[3], idx, i;
  int ntrials, nsamples, nchannels, nmarkers;
  int trial, channel; 			  /* counter */
  EEGdata_trials *eeg;
  EEGdata *new;
  double **Delta;

  nchannels   = ((mxGetDimensions(prhs[0]))[0]);
  nsamples    = ((mxGetDimensions(prhs[0]))[1]);
  ntrials     = ((mxGetDimensions(prhs[0]))[2]);


  nmarkers    = ((mxGetDimensions(prhs[1]))[1]);
  mexPrintf( "ml_padtw: Input Data is (%i x %i x %i), with %i markers\n", nchannels, ntrials, nsamples, nmarkers );

  eeg = init_eegdata_trials( ntrials, nmarkers, nchannels, nsamples, NULL );

  /* -------------------------------------------
	  copy data and markers to struct for handling 
     ------------------------------------------- */
  dptr = mxGetPr( prhs[0] );
  mptr = mxGetPr( prhs[1] );
  for( trial=0; trial<ntrials; trial++ ){
	 for( channel=0; channel<nchannels; channel++ ){ /* copying data */
		for( i=0; i<nsamples; i++ ){
		  /* mwIndex mxCalcSingleSubscript(const mxArray *pm, mwSize nsubs, mwIndex *subs); */
		  index[0] = channel;
		  index[1] = i;
		  index[2] = trial;
		  idx = mxCalcSingleSubscript( prhs[0], 3, index );
		  /* fprintf(stderr, "trial=%i, chan=%i, idx=%i\n", trial, channel, idx ); */
		  eeg->data[trial]->d[channel][i]=dptr[idx];
		}
	 }

	 /* copying markers */
	 for( i=0; i<nmarkers; i++ ){
		index[0] = trial;
		index[1] = i;
		idx = mxCalcSingleSubscript( prhs[1], 2, index );
		eeg->data[trial]->markers[i] = (unsigned long)mptr[idx];
		eeg->markers[trial][i] = (unsigned long)mptr[idx];
	 }
  }
  

  SettingsHierarchicalDTW settings = init_dtw_hierarchical( eeg ); 
  settings.progress = progressbar_matlab; 
  if( nrhs>2 ){
	 settings_from_matlabstruct( &settings, prhs[2] ); /* parse struct */
  }

  /* -------------------------------------------
	  do the calculation 
     ------------------------------------------- */

  mexPrintf("ml_padtw: Calculating Delta...");
  Delta = eegtrials_distmatrix_channel( eeg, vectordist_euclidean, 0, ALLOC_IN_FCT); 
  mexPrintf("... Done\n");

  new = init_eegdata( nchannels, nsamples, nmarkers ); 

  print_settings_hierarchicaldtw( stderr, settings ); 
  if( !eegtrials_dtw_hierarchical( eeg, Delta, ntrials, new, settings ) ){
	 mexErrMsgTxt("Sorry, there was some error in 'eegtrials_dtw_hierarchical'\n");
  }


  /* -------------------------------------------
	  return to MATLAB
     ------------------------------------------- */

  plhs[0] = mxCreateDoubleMatrix(nchannels, nsamples, mxREAL); /* return data */
  plhs[1] = mxCreateDoubleMatrix(1, nmarkers, mxREAL); /* return markers */

  dptr = mxGetPr( plhs[0] );
  mptr = mxGetPr( plhs[1] );
  for( channel=0; channel<nchannels; channel++ ){ /* copying data */
	 for( i=0; i<nsamples; i++ ){
		/* mwIndex mxCalcSingleSubscript(const mxArray *pm, mwSize nsubs, mwIndex *subs); */
		index[0] = channel;
		index[1] = i;
		idx = mxCalcSingleSubscript( plhs[0], 2, index );
		/* fprintf(stderr, "trial=%i, chan=%i, idx=%i\n", trial, channel, idx ); */
		dptr[idx] = new->d[channel][i];
	 }
  } 
  /* copying markers */
  for( i=0; i<nmarkers; i++ ){
	 mptr[i] = new->markers[i];
  }

  
  
  /* -------------------------------------------
	  cleaning up
     ------------------------------------------- */
  free_eegdata_trials( eeg );
  free_eegdata( new );
  matrix_free( Delta, ntrials );

  return;
}
        

void settings_from_matlabstruct( SettingsHierarchicalDTW *s, const mxArray *in ){
  char buf[200];
  int nfields, nelements;
  double *ptr;
  double tmp;
  mxArray *a;
  const char *fname;
  int i;

  nfields = mxGetNumberOfFields( in );
  nelements = mxGetNumberOfElements( in ); 

  for( i=0; i< nfields; i++ ){
	 fname = mxGetFieldNameByNumber( in, i );
	 a   = mxGetFieldByNumber( in, 0, i );
	 ptr = mxGetPr( a );	 

	 if( !strcmp( fname, "sampling_rate" ) ){
		s->sampling_rate=mxGetScalar( a );
		mexPrintf( "ml_padtw: set sampling_rate to '%f'\n", s->sampling_rate );
	 } else if( !strcmp( fname, "sigma" ) ){
		s->sigma=mxGetScalar( a );
		mexPrintf( "ml_padtw: set sigma to '%f'\n", s->sigma );
	 } else if( !strcmp( fname, "theta" ) ){
		s->theta1=ptr[0];
		s->theta2=ptr[1];
		mexPrintf( "ml_padtw: set theta to (%f,%f)\n", s->theta1, s->theta2 );
	 } else if( !strcmp( fname, "corner_freqs" ) ){
		s->corner_freqs[0]=ptr[0];
		s->corner_freqs[1]=ptr[1];
		mexPrintf( "ml_padtw: set corner_freqs to (%f,%f)\n", s->corner_freqs[0], 
					s->corner_freqs[1] );
	 } else if( !strcmp( fname, "winlength" ) ){	
		s->winlength=(int)mxGetScalar( a );
		mexPrintf( "ml_padtw: set winlength to '%i'\n", s->winlength );
	 } else if( !strcmp( fname, "N_freq" ) ){	
		s->N_freq=(int)mxGetScalar( a );
		mexPrintf( "ml_padtw: set N_freq to '%i'\n", s->N_freq );
	 } else if( !strcmp( fname, "N_time" ) ){	
		s->N_time=(int)mxGetScalar( a );
		mexPrintf( "ml_padtw: set N_time to '%i'\n", s->N_time );
	 } else if( !strcmp( fname, "progress" ) ){	
		tmp=(int)mxGetScalar( a );
		if( tmp<1 ){
		  s->progress=NULL;
		  mexPrintf( "ml_padtw: disabling progress-bar\n");
		} else { 
		  s->progress=progressbar_matlab;
		  mexPrintf( "ml_padtw: enabling progress-bar\n");
		}
	 } else if( !strcmp( fname, "regularize" ) ){	
		mxGetString( a, buf, mxGetN(a)+1 );
		if( !strcasecmp( buf, "none" ) ){
		  s->regularize = NULL;
		  mexPrintf( "ml_padtw: set regularization to 'none'\n" );
		} else if( !strcasecmp( buf, "gaussian_markers" ) ){
		  s->regularize = eeg_regularization_gaussian_line;
		  mexPrintf( "ml_padtw: set regularization to 'gaussian_markers\n" );
		} else {
		  mexPrintf( "ml_padtw: sorry, don't know '%s'\n", buf );
		}
	 }  else if( !strcmp( fname, "linkage" ) ){	
		mxGetString( a, buf, mxGetN(a)+1 );
		if( !strcasecmp( buf, "single" ) ){
		  s->linkage = dgram_dist_singlelinkage;
		  mexPrintf( "ml_padtw: set linkage to 'single'\n" );
		} else if( !strcasecmp( buf, "complete" ) ){
		  s->linkage = dgram_dist_completelinkage;
		  mexPrintf( "ml_padtw: set linkage to 'complete'\n" );
		} else if( !strcasecmp( buf, "average" ) ){
		  s->linkage = dgram_dist_averagelinkage;
		  mexPrintf( "ml_padtw: set linkage to 'average'\n" );
		} else {
		  mexPrintf( "ml_padtw: sorry, don't know '%s'\n", buf );
		}
	 } else if( !strcmp( fname, "windowfct" ) ){	
		mxGetString( a, buf, mxGetN(a)+1 );
		if( !strcasecmp( buf, "dirichlet" ) ){
		  s->winfct = window_dirichlet;
		  mexPrintf( "ml_padtw: set winfct to 'dirichlet'\n" );
		} else if( !strcasecmp( buf, "gaussian" ) ){
		  s->winfct = window_gaussian;
		  mexPrintf( "ml_padtw: set winfct to 'gaussian'\n" );
		} else if( !strcasecmp( buf, "hamming" ) ){
		  s->winfct = window_hamming;
		  mexPrintf( "ml_padtw: set winfct to 'hamming'\n" );
		} else if( !strcasecmp( buf, "hanning" ) ){
		  s->winfct = window_hanning;
		  mexPrintf( "ml_padtw: set winfct to 'hanning'\n" );
		} else if( !strcasecmp( buf, "kaiser" ) ){
		  s->winfct = window_kaiser;
		  mexPrintf( "ml_padtw: set winfct to 'kaiser'\n" );
		} else {
		  mexPrintf( "ml_padtw: sorry, don't know '%s'\n", buf );
		}
	 } else if( !strcmp( fname, "pointdistance" ) ){	
		mxGetString( a, buf, mxGetN(a)+1 );
		if( !strcasecmp( buf, "euclidean" ) ){
		  s->pointdistance=eeg_distmatrix_euclidean_channel;
		  mexPrintf( "ml_padtw: set pointdistance to 'euclidean'\n" );
		} else if( !strcasecmp( buf, "euclidean_derivative" ) ){
		  s->pointdistance=eeg_distmatrix_euclidean_derivative_channel;
		  mexPrintf( "ml_padtw: set pointdistance to 'euclidean_derivative'\n" );
		} else if( !strcasecmp( buf, "stft" ) ){
		  s->pointdistance=eeg_distmatrix_stft_channel;
		  mexPrintf( "ml_padtw: set pointdistance to 'stft'\n" );
		} else {
		  mexPrintf( "ml_padtw: sorry, don't know '%s'\n", buf );
		}
	 }

  } /* for nfields */
}
/** calls MATLAB-function "waitbar()" to show a progress bar
	 from a mex-file.
	 \param flag one of PROGRESSBAR_INIT, PROGRESSBAR_CONTINUE_LONG, 
	             PROGRESSBAR_CONTINUE_SHORT, PROGRESSBAR_FINISH
	 \param num meaning depends on flag
 */
void progressbar_matlab( int flag, int num ){
  char buf[100];
  
  switch( flag ){
  case PROGRESSBAR_INIT:
	 mexEvalString("h=waitbar(0, \'Hierarchical Averaging\');");
	 progress_status.max_progress = num;
	 progress_status.cur_progress = 0;
	 break;
  case PROGRESSBAR_CONTINUE_LONG:
	 sprintf( buf, "waitbar(%f,h);", (double)num/(double)progress_status.max_progress );
	 mexEvalString( buf );
	 break;
  case PROGRESSBAR_CONTINUE_SHORT:
	 break;
  case PROGRESSBAR_FINISH:
	 mexEvalString("close(h);");
	 break;
  }
}
