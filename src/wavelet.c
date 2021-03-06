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
#include "eeg.h"


/** \brief Wavelet estimation of single trial ERP's using Wang et al.'s (2007) technique. 
	 
	 \ingroup grpwavelet
	 Formula \f[ \lambda = \sigma \sqrt{2\log_e{n\log_2{n}}} \f].
*/
double translation_invariant_thresholding(const double *data, int n){
  dprintf("Db: translation_invariant_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data, n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n*glog(n, 2)));
  return lambda;
}

/**\brief  Wavelet estimation of single trial ERP's using a conventional estimate.
 *
	 \ingroup grpwavelet
 * Formula \f[ \lambda = \sigma \sqrt{2\log_e{n}} \f].
 */
double conventional_thresholding(const double *data, int n){
  dprintf("Db: conventional_thresholding\n");
  double sigma, lambda; /* population sd, threshold */
  
  sigma  = mad(data,n)/0.6745; /* funny constant, eh? */
  lambda = sigma * sqrt(2*log(n));
  return lambda;
}

/** \brief thresholding as proposed in the Matlab Wavelet-Toolbox.

	 \ingroup grpwavelet
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

/**\brief procedure described in Donoho et al., 1995. 

 \ingroup grpwavelet
*/
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

/** \brief Generic Wavelet-Denoising.
																	
	 \ingroup grpwavelet
	 compute DWT of signal, call the thresholding
	 function 'threshfct' for each resolution level and IDWT the signal 
*/
int wavelet_denoise   ( double *data, int n, WaveletParameters P ){
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  int j, J, k, offset;
  double lambda; /* population sd, threshold */
  dprintf("Db: generic_denoising\n");
  
  w = gsl_wavelet_alloc( P.wavelet, P.vanishing_moments );
  work = gsl_wavelet_workspace_alloc( n );

  gsl_wavelet_transform_forward( w, data, 1, n, work );

  /* -- thresholding here -- */
  J = (int)round(glog((double)n, 2));
  for(j=P.first_thresholding_level; j<J; j++){ /* loop through levels */
    offset = (int)pow(2, j);
    lambda = (*(P.threshselfct))(&(data[offset]), offset);

    for(k=offset; k<2*offset; k++){ /* loop through coefficients */
      /* soft or hard thresholding */
      data[k]=(*(P.threshfct))(data[k], lambda);
    }
  }
  /* -- thresholding end -- */

  gsl_wavelet_transform_inverse(w, data, 1, n, work);
  gsl_wavelet_free(w);
  gsl_wavelet_workspace_free(work);
  return 0;
}

/** \brief set default parameters for a waveletParameters struct.

	 \ingroup grpwavelet
 */
WaveletParameters wavelet_init(){
  WaveletParameters P;
  P.first_thresholding_level = 5;
  P.wavelet = (gsl_wavelet_type*)gsl_wavelet_daubechies;
  P.vanishing_moments = 10;
  P.threshselfct = heuristic_sure;
  P.threshfct = soft_thresholding;
  P.sigextfct = sigext_smooth;
  return P;
}

/**\brief Extend data to length 2^j using sigextfct and denoise it.

	\ingroup grpwavelet
 * \see generic_denoising()
 */
int wavelet_extend_and_denoise  ( double *data, int n, WaveletParameters P ){
  int J;
  double *tmp, *dptr;

  J = (int)round(glog((double)n, 2));
  J++;
  if(((int)pow(2, J))<n)
    dprintf(" Warning: Signal length %i is not a power of 2 -> extending signal to length %i\n",
	    n, (int)pow(2,J));
  tmp = (double*)malloc((int)pow(2,J) * sizeof(double));
  tmp = memcpy(tmp, data, n*sizeof(double));
  dptr=(*(P.sigextfct))(tmp, n, (int)pow(2, J));
  wavelet_denoise(tmp, (int)pow(2,J), P );
  
  data = memcpy(data, dptr, n*sizeof(double));

  free(tmp);
  return 0;
}

/** \brief soft-thresholding.

	 \ingroup grpwavelet
 * Formula: \f$ \eta_s(\lambda, w)=\f$
 */
double soft_thresholding(double d, double lambda){
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
/** \brief hard-thresholding.

	 \ingroup grpwavelet
 * Formula: \f$ \eta_h(\lambda, w)=\f$
 */
double hard_thresholding(double d, double lambda){
  /* hard thresholding */
/*   dprintf("Db: eta_h\n"); */
  if(fabs(d)<=lambda)
    d=0.0;
  return d;
}


/**\brief Extend data to length 2^j using sigextfct and denoise it.

	\see grpwavelet
 * \see generic_denoising()
 * \see extend_and_denoise()
 * data is directly written into the eegdata-struct or copied;

 \ingroup grpeegproc

 \param eeg input
 \param P parameters
 \param alloc allocate memory, or not
 */
EEG* eeg_wavelet_denoise( EEG *eeg, WaveletParameters P, bool alloc ){
#ifdef FIXEEG
  int c, i;
  EEG *eeg_out;
  if( alloc ){
	 eeg_out = eeg_clone( eeg, EEG_CLONE_ALL );
  } else {
	 eeg_out = eeg;
  }

  for( c=0; c<eeg->nbchan; c++ ){
	 for( i=0; i<eeg->ntrials; i++ ){
		wavelet_extend_and_denoise( eeg->data[c][i], eeg->n, P );
	 }
  }
  return eeg_out;
#endif
}
