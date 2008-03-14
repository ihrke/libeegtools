/*
*  C Implementation: ml_testwarppath
*
* Description: 
*
*
* Author: Matthias Ihrke <mihrke@uni-goettingen.de>, (C) 2008
*
* Copyright: See COPYING file that comes with this distribution
*
*/
#include "mex.h"
#include <stdlib.h>
#include "denoising.h"  
#include "definitions.h"
#include "averaging.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	  /** [path] = ml_testwarppath((const double*)u, (double*)s)
	  	*/
	char msg[500];
	int J, K;
	double *u, *s;
	double *up, *sp;
	int i;
	WarpPath *path;

	/* check proper input and output */
	if(nrhs!=2)
		mexErrMsgTxt("Need 2 inputs");
	else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
		mexErrMsgTxt("First and second Input must be double array.");
	else if(mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1){
		sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.", mxGetM(prhs[0]), mxGetN(prhs[0]));
		mexErrMsgTxt(msg);
	} 
  
	J = mxGetN(prhs[0]);
	K = mxGetN(prhs[1]);

	dprintf("J=%i, K=%i\n", J, K);

	u = mxGetPr(prhs[0]);
	s = mxGetPr(prhs[1]);

	path = get_warppath2( u, J, s, K, 1.0, 1.0 );

	plhs[0] = mxCreateDoubleMatrix(1, J+K, mxREAL); /* hold u-path component */
	plhs[1] = mxCreateDoubleMatrix(1, J+K, mxREAL); /* hold s-path component */
 	up = mxGetPr(plhs[0]);
	sp = mxGetPr(plhs[1]);
	for(i=0; i<J+K; i++){
		up[i] = (double)path->upath[i];
		sp[i] = (double)path->spath[i];
	}
	free_warppath( path);
  
	return;	
}
