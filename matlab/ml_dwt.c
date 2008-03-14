/* Matlab wrapper for DWT
 *   Compilation:
 *     MATLAB will need to see libGSL, helper.o, denoising.o
 *     Put these (or similar) lines to the bottom of ~/.matlab/<version>/mexopts.sh
 *          CFLAGS="$CFLAGS -I/home/ihrke/local/include"
 *	    LD="$LD -lm -lgsl -lgslcblas -lplot"
 *          LDFLAGS="$LDFLAGS -L/home/ihrke/local/lib"
 *    
 *    then: 
 *     call 'make' in the directory and
 *     >> mex -v ml_denoise.c denoising.o helper.o
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

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  /* [coeff] = ml_dwt((double*)data, (int) level, [(int)sigext, (int)direction])
   *
   * where 
   *    level        - level of daubechies wavelet
   *    direction    - +1=forward (**), -1=backward
   *    sigext       - 0=zeros, 1=zerosr, 2=sym (**), 3=smooth
   *
   * (**) marks default values
  */
  char msg[500];
  int direction=1;
  int n, N, J,sigext=2;
  int daublevel=10, i;
  double *d, *dptr, *tmp;
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;

  double* (*sigext_fct)(double*, int, int); 

  /* check proper input and output */
  if(nrhs<2)
    mexErrMsgTxt("At least 2 input arguments required:\n [coeff] = ml_dwt((double*)data, (int) level, [(int)sigext])");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  else if(!mxIsDouble(prhs[0]))
    mexErrMsgTxt("First Input must be double array.");
  else if(mxGetM(prhs[0])!=1){
    sprintf(msg, "Inputs must be 1-dim vectors, got %i, %i.", mxGetM(prhs[0]), mxGetN(prhs[0]));
    mexErrMsgTxt(msg);
  } 
  
  daublevel = (int)(mxGetPr(prhs[1]))[0];
  if(nrhs>=3) sigext=(int)(mxGetPr(prhs[2]))[0];
  if(nrhs>=4) direction=(int)(mxGetPr(prhs[3]))[0];

  /* Error handling */
  if(abs(direction)!=1) {
    sprintf(msg, "Error: Wrong direction for dwt (+1 or -1)");
    mexErrMsgTxt(msg);    
  } else if(!(daublevel>=4 && daublevel<=20 && (daublevel%2)==0)){
    sprintf(msg, "Error: Wrong wavelet level '%i' (4<=level<=20, level even)", daublevel);
    mexErrMsgTxt(msg);
  }
  printf("Found: daublevel=%i, sigext=%i, direction=%i\n", daublevel, sigext, direction);
  
  switch(sigext){
  case 0: sigext_fct=sigext_zeros; break;
  case 1: sigext_fct=sigext_zerosr; break;
  case 2: sigext_fct=sigext_sym; break;
  case 3: sigext_fct=sigext_smooth; break;
  default:     mexErrMsgTxt("error sig");
  }
  
  plhs[0] = mxDuplicateArray(prhs[0]); 
  n = mxGetN(plhs[0]);
  d = mxGetPr(plhs[0]);
  J = (int)ceil(glog((double)n, 2));
  N = pow(2, J);

  if(n==N){
    tmp  = d;
    dptr = d;
  } else {
    tmp = mxMalloc(N*sizeof(double));
    for(i=0; i<n; i++)
      tmp[i]=d[i];

    dptr = sigext_fct(tmp, n, N); 
  }

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, daublevel);
  work = gsl_wavelet_workspace_alloc(N);

  gsl_wavelet_transform(w, tmp, 1, N, direction, work);
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  
  /*  set the output pointer to the output matrix */
  if(n!=N){
    for(i=0; i<n; i++) 
      d[i]=dptr[i]; 
  }
  
  
  return;
}
        
