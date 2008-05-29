/**\file ml_padtw.c
 * \brief Matlab wrapper for PADTW
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
  int n, N, zero, i;
  int *sR;
  double **ui;
  double *d, *m;
  double *wavg;

  char *docstring = "[average] = ml_padtw(s, zero, sR, [opt])\n"
    "s    - Nxn matrix of EEG-data (N trials, n sampling points)\n"
    "zero - stimulus onset in sampling units\n"
    "sR   - N reaction times in sampling units\n"
    "OPTS:\n";

  /* check proper input and output */
  if(nrhs<3){
    sprintf(msg, "Need at least 3 inputs.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    sprintf(msg, "First and second Input must be double array.\n%s\n", docstring);
    mexErrMsgTxt(msg);
  } else if(mxGetM(prhs[1])!=1 || mxGetM(prhs[2])!=1){
    sprintf(msg, "Inputs must be 1-dim vectors, got %ix%i.\n%s\n", 
	    mxGetM(prhs[1]), mxGetN(prhs[1]), docstring);
    mexErrMsgTxt(msg);
  } 
  
  N = mxGetN(prhs[0]);
  n = mxGetM(prhs[0]);
  zero=(int)mxGetPr(prhs[1])[0];

  fprintf(stderr, "N=%i, n=%i, zero=%i\n", (int)N, (int)n, zero);
  fprintf(stderr, "markers(M=%i, N=%i)\n", mxGetM(prhs[2]), mxGetN(prhs[2]));

  sR = (int*)mxCalloc(N,sizeof(int));

  /* for convenience, create pointers to the correct segments */
  ui = (double**)mxCalloc(N, sizeof(double*));
  d = mxGetPr(prhs[0]);
  m = mxGetPr(prhs[2]);
  for(i=0; i<N; i++){
	  ui[i]=&(d[n*i]);
	  sR[i]=(int)m[i];
  }
  fprintf(stderr, "m[0]=%i, m[N]=%i\n", (int)m[0], (int)m[N-1]);
  
  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL); 
  wavg = mxGetPr(plhs[0]);

  wavg = PADTW((const double**)ui, N, n, zero, sR, wavg);

  return;
}
        
