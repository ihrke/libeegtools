/* Matlab wrapper for timewarping fct
 *   Compilation:
 *     MATLAB will need to see libGSL, helper.o, denoising.o
 *     Put these (or similar) lines to the bottom of ~/.matlab/<version>/mexopts.sh
 *          CFLAGS="$CFLAGS -I/home/ihrke/local/include"
 *	    LD="$LD -lm -lgsl -lgslcblas -lplot"
 *          LDFLAGS="$LDFLAGS -L/home/ihrke/local/lib"
 *    
 *    then: 
 *     call 'make' in the directory and
 *     >> mex -v ml_timewarp.c denoising.o helper.o
 *     from MATLAB
 *   
 *  Doc: 
 *      http://www.mathworks.com/support/solutions/data/1-1BDU5.html
 * 
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

/* #define DEBUG_LOCAL  */
/* #define PLOTTING_ON */

void matlab_plot(const double *times, const double *v, int n, const char *color){
#ifdef PLOTTING_ON
  mxArray *in[3];
  double *d1,*d2;
  char *d3;
  int dims[2]={1,255};
 
  in[0] = mxCreateDoubleMatrix(1,n,mxREAL);
  in[1] = mxCreateDoubleMatrix(1,n,mxREAL);
  in[2] = mxCreateCharArray(1, dims);
  d1 = mxGetPr(in[0]); d1 = memcpy(d1, times, n*sizeof(double));
  d2 = mxGetPr(in[1]); d2 = memcpy(d2, v, n*sizeof(double));
  d3 = mxGetPr(in[2]); d3 = memcpy(d3, color, strlen(color)*sizeof(char));
  mexCallMATLAB(0, NULL, 3, in, "plot");
  mxDestroyArray(in[0]);
  mxDestroyArray(in[1]);
  mxDestroyArray(in[2]);
  return;
#endif
} 

void  matlab_wait(){
#ifdef PLOTTING_ON
  mexCallMATLAB(0, NULL, 0, NULL, "waitforbuttonpress");
#endif
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /* [sig] = ml_warpavg((double*)times, (double**)s, (double*)rts, [
		(double)theta1, (double)theta2)
   *  -- times is a Nx1 array containing the real times in ms
   *     s is a Nxn matrix (N-number of trials, n-number of sampling points)
   *     rts is a Nx1 matrix with the reaction times in ms
   *  -- sig is the warpavg
   */
  char msg[500];
  /* parameters read from args */
  double theta1, theta2;
  double **ui, *wavg, *m, *d;
  int **markers, N,n, i;
  /* parse these later !! */
  theta1=1.0;
  theta2=1.0;

  /* check proper input and output */
  if(nrhs!=2){
    sprintf(msg, "[sig] = ml_warpavg((double**)s, (double*)markers)\n"
	    "*     s is a Nxn matrix (N-number of trials, n-number of sampling points)\n"
				 "*     markers is a Nx2 matrix (zero, sRu) for each trial"
	    "*  -- sig is the warpavg\n");
    mexErrMsgTxt(msg);
  } else if(!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1])){
    mexErrMsgTxt("All Inputs must be double arrays.");
  } 

  n = mxGetM(prhs[0]);
  N = mxGetN(prhs[0]);
  fprintf(stderr, "N=%i, n=%i\n", (int)N, (int)n);
  fprintf(stderr, "markers(M=%i, N=%i)\n", mxGetM(prhs[1]), mxGetN(prhs[1]));
  markers = (int**)mxCalloc(N,sizeof(int*));
  for(i=0; i<N; i++) markers[i] = (int*)mxCalloc(2,sizeof(int));
    
  /* for convenience, create pointers to the correct segments */
  ui = (double**)mxCalloc(N, sizeof(double*));
  d = mxGetPr(prhs[0]);
  m = mxGetPr(prhs[1]);
  for(i=0; i<N; i++){
	  ui[i]=&(d[n*i]);
	  markers[i][0]=(int)m[i];
	  markers[i][1]=(int)m[i+N];
  }
  fprintf(stderr, "m[0]=%i, m[N]=%i\n", (int)m[0], (int)m[N]);
  

  plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL); 
  wavg = mxGetPr(plhs[0]);
  wavg = warpaverage(ui, N, n, markers, wavg);
  
  mxFree(ui);	
  for(i=0; i<N; i++) mxFree(markers[i]);
  mxFree(markers);
  return;
}
        
