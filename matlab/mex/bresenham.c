/**
	\file bresenham.c
	\brief Matlab wrapper for \ref bresenham_linesegments().


 The MATLAB-interface is documented in the corresponding M-file.

*/

#include "mex.h"
#include "mex_utils.h"

#include <stdlib.h>

#define SCRIPTNAME "bresenham"

/* ----- LibEEGTools includes ----------*/
#include "array.h"
#include "imageproc.h"
/* ----- /LibEEGTools includes ----------*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	char msg[2000];
	char *docstring;
	docstring=get_mfile_as_string( SCRIPTNAME ".m" );

	/* check proper input and output */
	if(nrhs<1){
		 sprintf(msg, ">> Need exactly 1 input.\n%s\n", docstring);
		 mexErrMsgTxt(msg);
	} else if(!is_mex_matrix(prhs[0]) ){
		 sprintf(msg, ">> First Input must be matrix.\n%s\n", docstring);
		 mexErrMsgTxt(msg);
	 }
	const mwSize *mdims=mxGetDimensions( prhs[0] );
	if( mdims[0]!=2 ){
		 sprintf(msg, ">> Input must have first dimension 2.\n%s\n", docstring);
		 mexErrMsgTxt(msg);
	 }

	/* getting the appropriate array representation */
	Array *m = mex_mxarray_to_array( prhs[0] );
	array_typecast( m, INT );

	/* computation */
	Array *b=bresenham_linesegments( m );

	/* conversion to MATLAB */
	plhs[0]=mex_int_array_to_mxarray( b );

	/* cleaning up */
	array_free( m );
	array_free( b );

	if( docstring )
		free( docstring);
	return;
 }

