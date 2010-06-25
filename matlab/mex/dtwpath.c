/**
	\file dtwpath.c
	\brief Matlab wrapper for dynamic time warping.

 The script can be used to calculate the warping-path for
 either a given distance matrix between individual points in
 two signals, or by giving the (potential multidimensional) signals
 directly.

 It returns the warping path (that is corresponding elements in 
 both series as a 2 x M integer-array).

 The MATLAB-interface is as follows:

 [ path ] = dtwpath( s1, s2 );
 s1 and s2 are matrices or vectors giving the signals to be warped.

 [ path ] = dtwpath( s1, s2, metric );
 s1 and s2 are matrices or vectors giving the signals to be warped.
 metric is a string giving the distance metric. \todo implement this!

 [ path ] = dtwpath( distmat );
 distmat is the n x n pointwise distance matrix between signal points.
*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "distances.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
  char msg[2000];
  const char *docstring = " [ path ] = dtwpath( s1, s2 );\n"
	 " s1 and s2 are matrices or vectors giving the signals to be warped.\n\n"
	 " [ path ] = dtwpath( s1, s2, metric );\n"
	 " s1 and s2 are matrices or vectors giving the signals to be warped.\n"
	 " metric is a string giving the distance metric. \todo implement this!\n\n"
	 " [ path ] = dtwpath( distmat );\n"
	 " distmat is the n x n pointwise distance matrix between signal points.\n\n";

  /* check proper input and output */
  if(nrhs<1){
    sprintf(msg, ">> Need at least 1 input.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!is_mex_matrix(prhs[0])){
    sprintf(msg, ">> First Input must be a double-precision matrix.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } 

  Array *distmat;
  if( nrhs==1 ){ /* distance matric */
	 distmat = mex_mxarray_to_array( prhs[0] );
  } else { /* signals + something */
	 Array *s1 = mex_mxarray_to_array( prhs[0] );
	 Array *s2 = mex_mxarray_to_array( prhs[1] ); 
	 distmat=distmatrix_signaldist( vectordist_euclidean, s1, s2, NULL, NULL );
	 array_free( s1 );
	 array_free( s2 );
  }

  matrix_dtw_cumulate( distmat, FALSE, NULL );
  Array *path = matrix_dtw_backtrack( distmat );

  plhs[0]=mex_int_array_to_mxarray( path );

  /* cleaning up */
  array_free( path );
  array_free( distmat );

  return;
}
        
