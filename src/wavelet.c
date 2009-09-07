/* **************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "wavelet.h"
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <stdarg.h>
#include <string.h> /* memcpy */
#include "helper.h"
#include "mathadd.h"


/** Wavelet estimation of single trial ERP's using Wang et al.'s
 * (2007) technique. 
 *
 * Formula \f[ \lambda = \sigma \sqrt{2\log_e{n\log_2{n}}} \f].
*/
double translation_invariant_thresholding(const double *data, int n){
  dprintf("Db: translation_invariant_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data, n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n*glog(n, 2)));
  return lambda;
}

/** Wavelet estimation of single trial ERP's using the technique
 * labelled "conventional" in Wang et al. (2007).
 *
 * Formula \f[ \lambda = \sigma \sqrt{2\log_e{n}} \f].
 */
double conventional_thresholding(const double *data, int n){
  dprintf("Db: conventional_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data,n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n));
  return lambda;
}
/** thresholding as proposed in the Matlab Wavelet-Toolbox.
    \code
       case 'heursure' 
         hthr = sqrt(2*log(n)); 
         eta = (norm(x).^2-n)/n; 
         crit = (log(n)/log(2))^(1.5)/sqrt(n); 
         if eta < crit 
             thr = hthr; 
         else 
             thr = min(thselect(x,'rigrsure'),hthr); 
         end 
	 \endcode
*/
double heuristic_sure(const double *data, int n){
  dprintf("Db: heuristic_sure\n");
  double hthr, etat, crit;
  hthr = sqrt(2*log(n));
  etat = (pow(vnorm(data, n, 2), 2)-n)/(double)n;
  crit = pow((log(n)/log(2)), 1.5)/sqrt(n);
  if(etat < crit)
    return hthr;
  else 
    return fmin(sureshrink(data, n), hthr);
}

/** Wavelet estimation of single trial ERP's using the procedure described
 * in Donoho et al., 1995. The implementation described in the outline
 * is used; */
double sureshrink(const double *data, int n){
  int i,k;
  double lambda, sigma, sure, suremin;
  double *tmp;
  dprintf("Db: sureshrink\n");

  tmp = (double*)malloc(n*sizeof(double));
  sigma  = mad(data, n)/0.6745;   

  for(i=0; i<n; i++)
    tmp[i] = fabs(data[i])/sigma;
  /*  tmp = memcpy(tmp, data, n*sizeof(double)); */
 
  /* compute the SURESHRINK threshold */
  suremin = DBL_MAX;
  qsort(tmp, n, sizeof(double), abscmp);
  
  lambda=0.0;
  for(k=0; k<n; k++){
    sure = n - 2*(k+1)+(n-k)*pow(fabs(tmp[k]), 2);
    for(i=0; i<k; i++)
      sure = sure + pow(fabs(tmp[i]), 2);
    if(sure<suremin){
      suremin = sure;
      lambda = fabs(tmp[k]);
    }
  }
  lambda = sigma * lambda;
  
  free(tmp);
  return lambda;
}

/** Generic Denoising: compute DWT of signal, call the thresholding
   function 'threshfct' for each resolution level and IDWT the signal */
int generic_denoising   ( double *data, int n, int L, 
								  ThresholdSelectionFunction threshfct,
								  ThresholdFunction etafct ) {
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  int j, J, k, offset;
  double lambda; /* population sd, threshold */
  dprintf("Db: generic_denoising\n");

  w = gsl_wavelet_alloc(gsl_wavelet_daubechies, 20);
  work = gsl_wavelet_workspace_alloc(n);

  gsl_wavelet_transform_forward(w, data, 1, n, work);

  /* -- thresholding here -- */
  J = (int)round(glog((double)n, 2));
  for(j=L; j<J; j++){ /* loop through levels */
    offset = (int)pow(2, j);
    lambda = (*threshfct)(&(data[offset]), offset);

    for(k=offset; k<2*offset; k++){ /* loop through coefficients */
      /* soft or hard thresholding */
      data[k]=(*etafct)(data[k], lambda);
    }
  }
  /* -- thresholding end -- */

  gsl_wavelet_transform_inverse(w, data, 1, n, work);
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  return 0;
}

/** Extend data to length 2^j using sigextfct and denoise it.
 * \see generic_denoising()
 */
int extend_and_denoise  ( double *data, int n, int L, 
								  ThresholdSelectionFunction threshfct,
								  ThresholdFunction etafct, 
								  SignalExtensionFunction sigextfct ){
  int J;
  double *tmp, *dptr;

  J = (int)round(glog((double)n, 2));
  J++;
  if(((int)pow(2, J))<n)
    dprintf(" Warning: Signal length %i is not a power of 2 -> extending signal to length %i\n",
	    n, (int)pow(2,J));
  tmp = (double*)malloc((int)pow(2,J) * sizeof(double));
  tmp = memcpy(tmp, data, n*sizeof(double));
  dptr=(*sigextfct)(tmp, n, (int)pow(2, J));
  generic_denoising(tmp, (int)pow(2,J), L, threshfct, etafct);
  
  data = memcpy(data, dptr, n*sizeof(double));

  free(tmp);
  return 0;
}

/** Implements soft-thresholding.
 * Formula: \f$ \eta_s(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_s(double d, double lambda){
  /* soft thresholding, formula (6) */
/*   dprintf("Db: eta_s\n"); */
  if(fabs(d)<=lambda){
    d=0.0;
  } else if(d>lambda){
    d=d-lambda;
  } else if(d<(-1)*lambda){
    d=d+lambda;
  }
  return d;
}
/** Implements hard-thresholding.
 * Formula: \f$ \eta_h(\lambda, w)=\f$
 * \ingroup thresholding
 */
double eta_h(double d, double lambda){
  /* hard thresholding */
/*   dprintf("Db: eta_h\n"); */
  if(fabs(d)<=lambda)
    d=0.0;
  return d;
}


/** Extend data to length 2^j using sigextfct and denoise it.
 * \see generic_denoising()
 * \see extend_and_denoise()
 * data is directly written into the eegdata-struct
 */
void eeg_wavelet_denoise( EEGdata *eeg, int L, 
								  ThresholdSelectionFunction threshfct,
								  ThresholdFunction etafct, 
								  SignalExtensionFunction sigextfct ){
  int c;
  
  for( c=0; c<eeg->nbchan; c++ ){
	 extend_and_denoise( eeg->d[c], eeg->n, L, threshfct, etafct, sigextfct );
  }
}
