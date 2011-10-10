/**
	\file loscalc.c
	\brief Matlab wrapper for Line-of-synchronisation calculation.

	References:

	Marwan et al. Cross recurrence plot based synchronization
	of time series. Nonlinear Processes in Geophysics (2002) vol. 9 (3-4) pp. 325-331.

	Matthias Ihrke, Hecke Schrobsdorff and J. Michael Herrmann: Recurrence-Based
	Synchronization of Single Trials for EEG-Data Analysis. Lecture Notes on Computer
	Science 5788, Intelligent Data Engineering and Automated Learning - IDEAL 2009.
	118-125 doi:10.1007/978-3-642-04394-9

	The MATLAB-interface is documented in the corresponding M-file.
*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>
#include "array.h"
#include "recurrence_plot.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  char msg[2000];
  const char *docstring = "";

  /* check proper input and output */
  if(nrhs<1){
	 sprintf(msg, ">> Need at least 1 input.\n%s\n", docstring);
	 mexErrMsgTxt(msg);
  } else if( !is_mex_matrix(prhs[0]) ){
	 sprintf(msg, ">> Input must be double-precision matrix.\n%s\n", docstring);
	 mexErrMsgTxt(msg);
  }

  Array *R = mex_mxarray_to_array( prhs[0] );
  const char *method="ihrke2009";
  if( nrhs>=2 ){ /* method-string provided */
		method=get_mex_string( prhs[1] );
  }

  Array *los;
  mexPrintf("Using method '%s'\n", method);
  if( !strcmp( method, "ihrke2009") ){
		los = recplot_los_dtw_noise( R );
		mexPrintf("OK\n");
  } else if( !strcmp( method, "marwan" )){
		int dx=4,dy=4;
		if( nrhs>2 ){
			 dx=(int)get_mex_double( prhs[2]);
			 dy=dx;
		}
		if( nrhs>3 ){
			 dy=(int)get_mex_double( prhs[3]);
		}
		los = recplot_los_marwan( R, dx, dy );
		mexPrintf("dx=%i, dy=%i, OK\n", dx, dy);
  } else {
		mexErrMsgTxt("Method not found, exiting.\n");
  }

  plhs[0]=mex_int_array_to_mxarray( los );

  /* cleaning up */
  array_free( los );
  array_free( R );

  return;
}

