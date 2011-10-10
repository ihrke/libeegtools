/**
	\file crossrecplot.c
	\brief Matlab wrapper for Cross-Recurrence Plot computation.

	This script calculates the cross-recurrence plot of two (multi-dimensional)
	signals. See \ref recplot.h for details.

	The MATLAB-interface is documented in the corresponding M-file.
*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "recurrence_plot.h"
#include "optarg.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){ 
  char msg[2000];
  const char *docstring = "[ R ] = crossrecplot( s1, s2 );\n"
			" s1 and s2 are matrices or vectors giving the signals.\n\n"
			"[ R ] = crossrecplot( s1, s2, epsilon );\n"
			" s1 and s2 are matrices or vectors giving the signals, epsilon is\n"
			" the neighbourhood of the signals. The signals should be normalized before.\n\n"
			"[ R ] = crossrecplot( s1, s2, epsilon, neighbourhood_crit );\n"
			" s1 and s2 are matrices or vectors giving the signals, epsilon is\n"
			" the neighbourhood of the signals or a fixed amount of neighbours (FAN) that\n"
			" is used to calculate the recurrences. Whether an epsilon-ball or a FAN is\n"
			" used depends on neighbourhood_crit which is one of 'fan' or 'epsilon'.\n";

  /* check proper input and output */
  if(nrhs<2){
	 sprintf(msg, ">> Need at least 2 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if( !( is_mex_matrix(prhs[0]) || is_mex_vector(prhs[0])) ||
	     !( is_mex_matrix(prhs[1]) || is_mex_vector(prhs[1]))) {
		sprintf(msg, ">> Inputs must be double-precision matrices.\n%s\n", docstring);
		mexErrMsgTxt(msg);
  }

  Array *s1 = mex_mxarray_to_array( prhs[0] );
  Array *s2 = mex_mxarray_to_array( prhs[1] );

  int fan=0; /* boolean flag */
  double epsilon=(*((double*)array_max(s1)))/20.0;
  if( nrhs>=3 ) {
		if( nrhs>3 ){
			 char *method=get_mex_string( prhs[3] );
			 mexPrintf( "Method '%s' chosen...\n", method );
			 if( !strcmp( method, "fan") ){
				  fan=(int)get_mex_double( prhs[2] );
				  mexPrintf( " Using FAN with '%i' neighbours\n", fan );
			 }
		}
		epsilon=get_mex_double( prhs[2] );
  }
  OptArgList *opts=optarglist( "fan=int", fan );

  Array *R = recplot( s1, s2, NULL, epsilon, opts );

  plhs[0]=mex_array_to_mxarray( R );

  /* cleaning up */
  array_free( R );
  array_free( s1 );
  array_free( s2 );
  optarglist_free( opts );

  return;
}
        
