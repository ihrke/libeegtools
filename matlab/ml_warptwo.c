/*
*  C Implementation: ml_warptwo
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
	  /** [wavg] = ml_warptwo((const double*)u, (double*)s, int [zero sRu sRs])
		*/
	char msg[500];
	int n;
	double *u, *s, *tmp;
	double *avg;
	int i;
	int zero, sRu, sRs;

	/* check proper input and output */
	if(nrhs!=3)
		mexErrMsgTxt("Need 3 inputs");
	else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
		mexErrMsgTxt("First and second Input must be double array.");
	else if(mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1){
		sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.", mxGetM(prhs[0]), mxGetN(prhs[0]));
		mexErrMsgTxt(msg);
	} else if(mxGetN(prhs[2])!=3 || mxGetM(prhs[2])!=1){
		sprintf(msg, "3rd input must be 1x3, got %ix%i.", mxGetM(prhs[2]), mxGetN(prhs[2]));
		mexErrMsgTxt(msg);
	}
	
	n = mxGetN(prhs[0]);
	if(n!=mxGetN(prhs[1])){
		sprintf(msg, "1st and 2nd input must be equally long, got %i and %i.", n, mxGetN(prhs[1]));
		mexErrMsgTxt(msg);				
	}
	tmp = mxGetPr(prhs[2]);
	zero = (int)tmp[0];
	sRu  = (int)tmp[1];
	sRs  = (int)tmp[2];

	dprintf("n=%i, zero=%i, sRu=%i, sRs=%i\n", n, zero, sRu, sRs);

	u = mxGetPr(prhs[0]);
	s = mxGetPr(prhs[1]);

	plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL); 
	avg = mxGetPr(plhs[0]);
	
	avg = warpavg_two( u, s, n, zero, sRu, sRs, avg );
	  
	return;	
}
