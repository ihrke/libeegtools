/**
	\file distmatrix.c
	\brief Matlab wrapper for \ref matrix_distmatrix().

 It calculates the pointwise distance matrix between elements of the signal.

 The MATLAB-interface is documented in the corresponding M-file.

*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "distances.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	 char msg[2000];
	 const char *docstring ="";

	 /* check proper input and output */
	 if(nrhs<1){
		  sprintf(msg, ">> Need at least 1 input.\n%s\n", docstring);
		  mexErrMsgTxt(msg);
	 } else if(!mxIsDouble(prhs[0]) ){
		  sprintf(msg, ">> First Input must be double array.\n%s\n", docstring);
		  mexErrMsgTxt(msg);
	 }

	 /* get the metric */
	 const char *metric="euclidean";
	 if( nrhs>1 ){
		  metric=get_mex_string( prhs[1] );
	 }
	 VectorDistanceFunction f;

	 mexPrintf("Using metric '%s'\n", metric);
	 if( !strcmp( metric, "euclidean" ) ){
		  f = vectordist_euclidean;
		  mexPrintf("OK\n");
	 } else if ( !strcmp( metric, "euclidean_normalized" ) ){
		  f = vectordist_euclidean_normalized;
		  mexPrintf("OK\n");
	 } else if ( !strcmp( metric, "dtw" ) ){
		  f = vectordist_dtw;
		  mexPrintf("OK\n");
	 } else {
		  mexErrMsgTxt("not found! Exiting.\n");
	 }

	 /* getting the appropriate array representation */
	 Array *tmp = mex_mxarray_to_array( prhs[0] );
	 Array *s = array_fromptr2( DOUBLE, 2, tmp->data, tmp->size[0],
										  (tmp->ndim>1)?(tmp->size[1]):1 );
	 /* computation */
	 Array *d=matrix_distmatrix( f, s, NULL, NULL );

	 plhs[0]=mex_array_to_mxarray( d );

	 /* cleaning up */
	 array_free( s );
	 array_free( d );
	 array_free( tmp );

  return;
}

