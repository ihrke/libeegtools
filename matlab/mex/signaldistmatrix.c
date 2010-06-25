/**
	\file signaldistancematrix.c
	\brief Matlab wrapper for \ref distmatrix_signaldist().

 It calculates the 
 pointwise distance matrix: \f$ d_{ij} = d( s1_i, s2_j ) \f$.


 The MATLAB-interface is as follows:

signaldistfunction Pointwise distance Matrix between two signals.
 [ distmat ] = signaldistfunction( signal1, signal2 );
 Calculate a euclidean distance matrix.

 [ distmat ] = signaldistfunction( signal1, signal2, metric, ... );
 Calculate a distance matrix with another metric. 
 Valid values for metric are:
 "euclidean", ...
 You can pass additional parameters for the distance function along.


*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "distances.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
  char msg[2000];
  const char *docstring = 
	 "signaldistfunction Pointwise distance Matrix between two signals.\n"
	 "[ distmat ] = signaldistfunction( signal1, signal2 );\n"
	 " Calculate a euclidean distance matrix.\n\n"
	 "[ distmat ] = signaldistfunction( signal1, signal2, metric, ... );\n"
	 " Calculate a distance matrix with another metric.\n"
	 " Valid values for metric are:\n"
	 " 'euclidean', ...\n"
	 " You can pass additional parameters for the distance function along.\n";
  
  /* check proper input and output */
  if(nrhs<2){
    sprintf(msg, ">> Need at least 2 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    sprintf(msg, ">> First and second Input must be double arrays.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } 

  /* getting the appropriate array representation */
  Array *s1 = mex_mxarray_to_array( prhs[0] );
  Array *s2 = mex_mxarray_to_array( prhs[1] );

  /* computation */
  Array *d=distmatrix_signaldist( vectordist_euclidean, s1, s2, NULL, NULL );
  plhs[0]=mex_array_to_mxarray( d ); 

  /* cleaning up */
  array_free( s1 );
  array_free( s2 );
  array_free( d );

  return;
}
        
