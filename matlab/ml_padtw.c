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
 *  Doc: 
 *      http://www.mathworks.com/support/solutions/data/1-1BDU5.html
 *
 * nlhs (Type = int): This paramter represents the number of "left hand side" arguments. So in my example
 *               function call, nlhs = 2 (the outputs are z0 and z1).
 * plhs (Type = array of pointers to mxArrays): This parameter is the actual output arguments.  As we will see
 *               later, an mxArray is MATLAB's structure for holding data and each element in plhs holds an mxArray of data.
 * nrhs (Type = int): Similar to nlhs, this paramter holds the number of "right hand side" arguments.
 * prhs (Type = const array of pointers to mxArrays): This array hold all of the pointers to the mxArrays of input data
 * for instance, prhs[0] holds the mxArray containing x, prhs[1] holds
 * the mxArray containing y, etc). 
 */
#include "mex.h"
#include <stdlib.h>
#include "warping.h"

void progressbar_matlab( int flag, int num ){
  char buf[100];

  switch( flag ){
  case PROGRESSBAR_INIT:
	 mexEvalString("h=waitbar(0, \'Hierarchical Averaging\')");
	 progress_status.max_progress = num;
	 progress_status.cur_progress = 0;
	 break;
  case PROGRESSBAR_CONTINUE_LONG:
	 sprintf( buf, "waitbar(h,%f);", (double)num/(double)progress_status.max_progress );
	 mexEvalString( buf );
	 break;
  case PROGRESSBAR_CONTINUE_SHORT:
	 break;
  case PROGRESSBAR_FINISH:
	 mexEvalString("close(h)");
	 break;
  }

}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  char *docstring = 
	 "This function is part of the EEGlab-Plugin for the removal of\n"
	 "temporal distortion. The function is called thus:\n"
	 "[padtw] = ml_padtw( X, markers, settings );\n\n"
	 "Here, X is the input matrix, which is 3-dimensional, \n"
	 "(channel x trial x sample).\n"
	 "markers is a (nmarkers x trials) matrix, containing the time-markers\n"
	 "in sampling units, used for the regularization.\n"
	 "settings is a struct, containing the configuration for \n"
	 "the hierarchical averaging procedure:\n"
	 " -> fields:\n"
	 "    - regularization\n"
	 "    - sigma\n"
	 "    - linkage\n"
	 "    - distance\n";
  char msg[500];
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

  double *dptr, *curptr;
  int index[3], idx;
  int ntrials, nsamples, nchannels, nmarkers;
  int trial, channel; 			  /* counter */
  EEGdata_trials *eeg;

  nchannels   = ((mxGetDimensions(prhs[0]))[0]);
  ntrials     = ((mxGetDimensions(prhs[0]))[1]);
  nsamples    = ((mxGetDimensions(prhs[0]))[2]);

  nmarkers    = ((mxGetDimensions(prhs[1]))[0]);
  mexPrintf("Input Data is (%i x %i x %x), with %i markers\n", nchannels, ntrials, nsamples, nmarkers );

  eeg = init_eegdata_trials( ntrials, nmarkers, nchannels, nsamples, NULL );


  /* copy data and markers to struct for handing */
  dptr = mxGetPr( prhs[0] );
  for( trial=0; trial<ntrials; trial++ ){
	 /* copying data */
	 for( channel=0; channel<nchannels; channel++ ){
		/* mwIndex mxCalcSingleSubscript(const mxArray *pm, mwSize nsubs, mwIndex *subs); */
		index[0] = channel;
		index[1] = trial;
		index[2] = 0;
		idx = mxCalcSingleSubscript( prhs[0], 3, index );
		curptr = &(dptr[idx]);
		memcpy( eeg->data[trial]->d[channel], curptr, nsamples*sizeof(double) );
	 }
	 /* copying markers */	 
	 
  }
  
  mexPrintf("data: %f, %f, %f, %f, %f\n",
				eeg->data[0]->d[0][0], 
				eeg->data[1]->d[0][0], 
				eeg->data[0]->d[2][4], 
				eeg->data[3]->d[2][1], 
				eeg->data[2]->d[1][1] );
  
  
  /* EEGdata *new; */
  /* double **Delta; */
  /* int i, n; */
  
  /* /\* get data *\/				   */
  /* oprintf("Reading from '%s'\n", buf); */
  /* eeg=read_eegtrials_from_raw( buf ); */
  /* print_eegdata_trials(stderr, eeg); */
  /* N = eeg->ntrials; */

  /* Delta = eegtrials_distmatrix_channel( eeg, vectordist_euclidean, 0, ALLOC_IN_FCT); */
  /* matrix_print( Delta, N, N ); */

  /* fprintf(stderr, "eeg->n=%i\n", eeg->nsamples ); */
  /* n = eeg->nsamples; */
  /* SettingsHierarchicalDTW settings = init_dtw_hierarchical( eeg ); */
  /* new = init_eegdata( eeg->data[0]->nbchan, eeg->nsamples, eeg->nmarkers_per_trial ); */

  /* settings = init_dtw_hierarchical( eeg ); */
  /* settings.progress = progressbar_matlab; */
  /* print_settings_hierarchicaldtw( stderr, settings ); */
  /* eegtrials_dtw_hierarchical( eeg, Delta, N, new, settings ); */
  /* write_eegdata_ascii_file( "test.out", new ); */
  


  plhs[0] = mxCreateDoubleMatrix(1, 10, mxREAL); 

  return;
}
        
