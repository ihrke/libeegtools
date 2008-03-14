/**\file ml_denoise.c
 * \brief Matlab wrapper for denoising fct.
 *
 *   Compilation:
 *     MATLAB will need to see libGSL, helper.o, denoising.o
 \code
 *     Put these (or similar) lines to the bottom of ~/.matlab/<version>/mexopts.sh
 *          CFLAGS="$CFLAGS -I/home/ihrke/local/include"
 *	    LD="$LD -lm -lgsl -lgslcblas -lplot"
 *          LDFLAGS="$LDFLAGS -L/home/ihrke/local/lib"
 *    
 *    then: 
 *     call 'make' in the directory and
 *     >> mex -v ml_denoise.c denoising.o helper.o
 *     from MATLAB
 *   \endcode
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

/**
   \code
   [sig] = ml_denoise((double*)data 
 *          [, (int)L, (int)denoising, (int)thresholding, (int)sigext)]
 * if data is Nxn, the denoising is done for the N trials of size n each
 *
 * where 
 *    denoising    - 0=conventional, 1=ti, 2=sureshrink (**),
 *                   3=heuristic_sureshrink
 *    thresholding - 0=hard, 1=soft (**)
 *    sigext       - 0=zeros, 1=zerosr, 2=sym (**), 3=smooth
 *
 * (**) marks default values
 \endcode
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  char msg[500];
  int n, N, L=6, denoising=2, thresholding=1, sigext=2;
  int i;
  double *d;

  double (*den_fct)(const double*, int);
  double (*eta_fct)(double, double);
  double* (*sigext_fct)(double*, int, int); 

  /* check proper input and output */
  if(nrhs<1)
    mexErrMsgTxt("At least one input arguments required:\n"
		 "[sig] = ml_denoise((double*)data\n"
		 " [, (int)L, (int)denoising, (int)thresholding, (int)sigext)]");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  else if(!mxIsDouble(prhs[0]))
    mexErrMsgTxt("First Input must be double array.");
  
/*   sprintf(msg, "AAA got %i, %i.", mxGetM(prhs[0]), mxGetN(prhs[0])); */
/*   mexErrMsgTxt(msg); */

  if(nrhs>=2) L=(int)(mxGetPr(prhs[1]))[0];
  if(nrhs>=3) denoising=(int)(mxGetPr(prhs[2]))[0];
  if(nrhs>=4) thresholding=(int)(mxGetPr(prhs[3]))[0];
  if(nrhs>=5) sigext=(int)(mxGetPr(prhs[4]))[0];

  switch(denoising){
  case 0: den_fct=conventional_thresholding; break;
  case 1: den_fct=translation_invariant_thresholding; break;
  case 2: den_fct=sureshrink; break;
  case 3: den_fct=heuristic_sure; break;
  default:     mexErrMsgTxt("error den");
  }

  if(thresholding==0) eta_fct=eta_h;
  else eta_fct=eta_s;
  
  switch(sigext){
  case 0: sigext_fct=sigext_zeros; break;
  case 1: sigext_fct=sigext_zerosr; break;
  case 2: sigext_fct=sigext_sym; break;
  case 3: sigext_fct=sigext_smooth; break;
  default:     mexErrMsgTxt("error sig");
  }
  
  plhs[0] = mxDuplicateArray(prhs[0]); 
  n = mxGetM(plhs[0]);
  N = mxGetN(plhs[0]);
  d = mxGetPr(plhs[0]);

  printf("L=%i, den=%i, thr=%i, sig=%i, n=%i, N=%i\n", L, denoising, thresholding, sigext, n, N);

  for(i=0; i<N; i++)
    extend_and_denoise(&(d[i*n]), n, L, den_fct, eta_fct, sigext_fct);
  
  return;
}
        
