/**
	\file tdelayreconstruction.c
	\brief Constructing a time-delay reconstruction of scalar time-series.

Constructing a time-delay reconstruction of scalar time-series.

 The MATLAB-interface is documented in the corresponding M-file.

*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>

/* ----- LibEEGTools includes ----------*/
#include "array.h"
#include "nonlinear.h"
/* ----- /LibEEGTools includes ----------*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	 char msg[2000];
	 const char *docstring ="";

	 /* check proper input and output */
	 if(nrhs!=3){
		  sprintf(msg, ">> Need 3 inputs.\n%s\n", docstring);
		  mexErrMsgTxt(msg);
	 } else if(!mxIsDouble(prhs[0]) || !is_mex_vector(prhs[0]) ){
		  sprintf(msg, ">> First Input must be double vector.\n%s\n", docstring);
		  mexErrMsgTxt(msg);
	  } else if( !(mxIsScalar(prhs[1]) && mxIsScalar(prhs[2])) ){
		  sprintf(msg, ">> Second and Third Inputs must be scalars.\n%s\n", docstring);
		  mexErrMsgTxt(msg);
	  }

	 /* getting the appropriate array representation */
	 Array *s = mex_mxarray_to_array( prhs[0] );
	 int m=(int)get_mex_double( prhs[1]);
	 int tau=(int)get_mex_double( prhs[2]);
	 double *x=s->data;
	 int n=s->size[0];

	 /* computation */
	 TimeDelayReconstruction* td = tdelay_init ( m, tau, x, n );
	 Array *tda = tdelay_to_array( td );

	 /* conversion to MATLAB */
	 plhs[0]=mex_array_to_mxarray( tda );

	 /* cleaning up */
	 array_free( s );
	 tdelay_free( td );

  return;
}

