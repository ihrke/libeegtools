/**\file ml_running_median.c
 * \brief Matlab wrapper for Running Median
 *
 *   Compilation:
 *     MATLAB will need to see libGSL, helper.o, denoising.o
 *     Put these (or similar) lines to the bottom of
 *     ~/.matlab/<version>/mexopts.sh
 \code
 *          CFLAGS="$CFLAGS -I/home/ihrke/local/include"
 *	    LD="$LD -lm -lgsl -lgslcblas -lplot"
 *          LDFLAGS="$LDFLAGS -L/home/ihrke/local/lib"
 *    
 *    then: 
 *     call 'make' in the directory and
 *     >> mex -v ml_timewarp.c denoising.o helper.o
 *     from MATLAB
 *   
 \endcode
 *  Doc: 
 *      http://www.mathworks.com/support/solutions/data/1-1BDU5.html
 *
 * nlhs (Type = int): This paramter represents the number of "left hand side" arguments. So in my example
 *               function call, nlhs = 2 (the outputs are z0 and z1).
 * plhs (Type = array of pointers to mxArrays): This parameter is the actual output arguments.  As we will see
 *               later, an mxArray is MATLAB's structure for holding data and each element in plhs holds an mxArray of data.
 * nrhs (Type = int): Similar to nlhs, this paramter holds the number of "right hand side" arguments.
 * prhs (Type = const array of pointers to mxArrays): This array hold all of the pointers to the mxArrays of input data
 * for instance, prhs[0] holds the mxArray containing x, prhs[1] holds
 * the mxArray containing y, etc). 
 */
#include "mex.h"
#include <stdlib.h>
#include "denoising.h"  
#include "averaging.h"  

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /* DOKU
   */
  char msg[500], *buf;
  int n, win;
  double *s, *fs;
  double(*metric)(double,double);

  char *docstring = "[filtered] = ml_weighted_running_median(x, win, [opt])\n"
    "x - vector to filter\n"
    "win - [int] window size for running median\n"
    "OPTS:\n"
    "  metric - values: 'euclidean'";
  /* check proper input and output */
  if(nrhs<2){
    sprintf(msg, "Need at least 2 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    sprintf(msg, "First and second Input must be double array.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1){
    sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.\n%s\n", 
	    mxGetM(prhs[0]), mxGetN(prhs[0]), docstring);
    mexErrMsgTxt(msg);
  } 
  
  n = mxGetN(prhs[0]);
  s = mxGetPr(prhs[0]);
  win = (int)mxGetPr(prhs[1])[0];
  metric=dist_euclidean;
  if(nrhs>2){
    buf = mxArrayToString(prhs[2]);
    if(!strcmp(buf, "euclidean"))
      metric = dist_euclidean;
    else{
      sprintf(msg, "Metric not known: %s\n", buf);
      mexErrMsgTxt(msg);
    }
    mxFree(buf);
  }

  /* printf("win=%i\n", n); */
  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  fs = mxGetPr(plhs[0]);
  fs = memcpy((void*)fs, (void*)s, n*sizeof(double));
  /*   printf("fs[0] = %f, fs[n-1]=%f\n", fs[0], fs[n-1]); */
  fs = weighted_running_median(fs, n, win, metric);

  return;
}
        
