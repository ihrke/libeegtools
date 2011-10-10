/**
	\file dtwadd.c
	\brief Matlab wrapper for dynamic time warping.

 This function adds two signals either directly computing
 the warping, or using a previously calculated warping-path.

 The MATLAB-interface is documented in the corresponding M-file.

*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "distances.h"
#include "warping.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
  char msg[2000];
  const char *docstring = "[ average ] = dtwadd( s1, s2 );\n"
	 "s1 and s2 are matrices or vectors giving the signals to be warped and\n"
	 "added. The signals must be n x N, where n is the sample and N the\n"
	 "dimensionality.\n\n"
	 "[ average ] = dtwadd( s1, s2, path );\n"
	 "s1 and s2 are matrices or vectors giving the signals to be\n"
	 "added. path is a previously calculated warping path, i.e.\n"
	 "a 2 x M integer array giving corresponding points in the\n"
	 "signals.\n";

  /* check proper input and output */
  if(nrhs<2){
	  sprintf(msg, ">> Need at least 2 inputs.\n%s\n", docstring);
	  mexErrMsgTxt(msg);
  } else if( !is_mex_matrix(prhs[0]) || !is_mex_matrix(prhs[1]) ){
	  sprintf(msg, ">> Inputs must be double-precision matrices.\n%s\n", docstring);
	  mexErrMsgTxt(msg);
  } else if( !mex_have_same_size(prhs[0], prhs[1])){
	  const mwSize *dims1=mxGetDimensions( prhs[0] );
	  const mwSize *dims2=mxGetDimensions( prhs[1] );
	  sprintf(msg, ">> Inputs must of the same size (%i x %i) vs. (%i x %i).\n%s\n",
			  dims1[0], dims1[1], dims2[0], dims2[1], docstring);
	  mexErrMsgTxt(msg);
  }

  Array *s1 = mex_mxarray_to_array( prhs[0] );
  Array *s2 = mex_mxarray_to_array( prhs[1] ); 
  Array *distmat, *path;
  int i;

  if( nrhs==2 ){ /* only the signals */
	 distmat=distmatrix_signaldist( vectordist_euclidean, s1, s2, NULL, NULL );
	 matrix_dtw_cumulate( distmat, FALSE, NULL );
	 path = matrix_dtw_backtrack( distmat );
	 array_free( distmat );
  } else { /* signals + path */
	 path = mex_int_mxarray_to_array( prhs[2] ); 
	 bool ispath;
	 warppath_CHECK( ispath, path);
	 if( !ispath ){
		 mexErrMsgTxt(">> 3rd input must be a warppath\n");
	 }
	 for( i=0; i<path->size[1]; i++ ){
		 if( array_INDEX2( path, uint, 0, i )>=s1->size[0] ||
			 array_INDEX2( path, uint, 1, i )>=s1->size[0] ||
			 array_INDEX2( path, uint, 1, i )<0 ||
			 array_INDEX2( path, uint, 1, i )<0 ){
			 sprintf( msg, ">> Invalid path: (%i,%i), max=%i\n", array_INDEX2( path, uint, 0, i ),
					  array_INDEX2( path, uint, 0, i ), s1->size[0]);
			 mexErrMsgTxt(msg );
		 }
	 }
  }

  Array *avg = dtw_add_signals( s1, s2, path, NULL );
  plhs[0]=mex_array_to_mxarray( avg );

  /* cleaning up */
  array_free( path );
  array_free( s1 );
  array_free( s2 );
  array_free( avg );

  return;
}
        
