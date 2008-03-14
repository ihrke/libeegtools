/**\file ml_template.c
 * \brief Matlab wrapper for <INSERT>
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

  /* check proper input and output */
  if(nrhs!=2)
    mexErrMsgTxt("Need 2 inputs");
  else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]))
    mexErrMsgTxt("First and second Input must be double array.");
  else if(mxGetM(prhs[0])!=1 || mxGetM(prhs[1])!=1){
    sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.", mxGetM(prhs[0]), mxGetN(prhs[0]));
    mexErrMsgTxt(msg);
  } 
  
  /* Important functions.
     mxGetPr(prhs[0]);
     mxCalloc((int)MAX(J,K), sizeof(double));
     mxFree();
     mxCreateDoubleMatrix(1, J, mxREAL);
     mxSetPr(plhs[0], snew);
  */

  return;
}
        
