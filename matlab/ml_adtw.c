/**\file ml_adtw.c
 * \brief Matlab wrapper for ADTW
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
  char msg[500];
  int n, sR1, sR2, zero;
  double *s1, *s2, *avg;

  char *docstring = "[avg] = ml_adtw(s1, s2, zero, sR)\n"
    "s1,s2 - signals to average\n"
    "zero  - stimulus onset in both signals\n"
    "sR    - 1x2 array of reaction times in both signals (sampling units)\n"
    "OPTS:\n";

  /* check proper input and output */
  if(nrhs<4){
    sprintf(msg, "Need at least 4 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    sprintf(msg, "First and second Input must be double array.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1 || mxGetM(prhs[2])!=1|| mxGetM(prhs[3])!=1){
    sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.\n%s\n", 
	    mxGetM(prhs[0]), mxGetN(prhs[0]), docstring);
    mexErrMsgTxt(msg);
  } else if(mxGetN(prhs[0])!=mxGetN(prhs[1])){
    sprintf(msg, "s1 and s2 must have same dimensions, got %i, %i.\n%s\n",
	    mxGetN(prhs[0]), mxGetN(prhs[1]), docstring);
    mexErrMsgTxt(msg);
  }

  n = mxGetN(prhs[0]);
  s1 = mxGetPr(prhs[0]);
  s2 = mxGetPr(prhs[1]);

  zero = (int)mxGetPr(prhs[2])[0];
  sR1  = (int)mxGetPr(prhs[3])[0];
  sR2  = (int)mxGetPr(prhs[3])[1];

  printf("zero=%i, sR1=%i, sR2=%i\n", zero, sR1, sR2); 
  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
  avg = mxGetPr(plhs[0]);

  ADTW_signal(s1, sR1, s2, sR2, zero, n, avg); 

  return;
}
        
