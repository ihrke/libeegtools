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

#include "denoising.h"

/* ---------------------------------------------------------------------------- 
   -- Denoising routines                                                     -- 
   ---------------------------------------------------------------------------- */
/** Simple Running Median filter.
 * Is applied to each channel separately.
 */
void eeg_filter_running_median(EEGdata *s, int win){
  int i;
  for(i=0; i<s->nbchan; i++){
    running_median(s->d[i], s->n, win);
  }
}

/** Simple Running Median filter.
 * \code d[i] = median{ d[i-win], ..., d[i+win] } \endcode
 * -- choice of time-window 'win' is crucial
 * use 2*l+1 where l is the size of to-be-preserved details in sampling units
 * Ref: Schelter et al, 2005, chap. 6 (6.2.1)
 */
double* running_median(double *d, int n, int win){
  double *tmp, *org;
  int i, awin, ewin;
  double med;

  org = (double*)malloc(n*sizeof(double));
  org = (double*)memcpy(org, d, n*sizeof(double));
  tmp = (double*)malloc((2*win+1)*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;
    
    tmp = (double*)memcpy(tmp, &(org[awin]), sizeof(double)*(ewin-awin));
   
    gsl_sort(tmp, 1, ewin-awin);
    med = gsl_stats_median_from_sorted_data(tmp, 1, ewin-awin);

    d[i] = med;
  }
  free(tmp);
  free(org);
  return d;
}

/** Simple Moving Average filter in the time-domain.
 * s[i] = mean{ s[i-win], ..., s[i+win] }
 */
double* moving_average(double *s, int n, int win){
  int i, awin, ewin;
  double m;
  double *org;

  org = (double*)malloc(n*sizeof(double));
  org = (double*)memcpy(org, s, n*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;
	 m = gsl_stats_mean(&(org[awin]), 1, ewin-awin);
	 s[i] = m;
  }
  free(org);
  return s;
}

/**  Weighted Running Median filter.
 * Is applied to each channel separately.
 */
void eeg_filter_weighted_running_median(EEGdata *s, int win){
  int i;
  for(i=0; i<s->nbchan; i++){
    weighted_running_median(s->d[i], s->n, win, pointdist_euclidean);
  }
}
/** Weighted Running Median filter.
 * Ref: Schelter et al, 2005, chap. 6 (6.2.1)
 */
double* weighted_running_median(double *d, int n, int win, 
				double(*dist)(double,double)){
  double *pptr, *wptr;
  int i, j;
  int awin, ewin;
  double med;
  
  pptr = (double*)malloc((2*win+1)*sizeof(double));
  wptr = (double*)malloc((2*win+1)*sizeof(double));

  for(i=0; i<n; i++){
	 awin = i-win;
	 ewin = i+win;
	 if(awin<0) awin=0;
	 if(ewin>n-1) ewin=n;

    memcpy(pptr, &(d[awin]), sizeof(double)*(ewin-awin));
    for(j=0; j<(ewin-awin); j++){
      wptr[j] = (*dist)(i, awin+j);
    }
/* 	 dprintf("i=%i, awin=%i, ewin=%i\n", i, awin, ewin); */
    med = weighted_median_from_unsorted(pptr, wptr, ewin-awin);
	 
    d[i] = med;
  }
  free(pptr);
  free(wptr);

  return d;
}

/** Weighted Median computation.
 * Formula:
 *  \f[ WM = y_{(k)}\qquad k = \max\left\{h:\sum_{i=h}^n w_{(i)} \ge \frac{1}{2}\sum_{i=1}^{n} w_i\right\} \f]
 * \param d - data
 * \param w - corresponding weights (>=0)
 */
double weighted_median_from_unsorted(const double *d, const double *w, int n){
  size_t *permut; 
  double ref=0.0, refh;
  int k;
  int i, h, idx;

  permut=(size_t*)malloc(n*sizeof(size_t));
  gsl_sort_index(permut, d, 1, n);  
/*   dprintf("i\tperm[i]\td[i]\td[permut[i]]\n"); */
/*   for(i=0; i<n; i++){ */
/* 	 dprintf("%i\t%i\t%.2f\t%.2f\n", i, permut[i], d[i], d[permut[i]]); */
/*   } */

  for(i=0; i<n; i++)
    ref+=w[i];
  ref /= 2.0;
  
  k=0;
  for(h=n-1; h>0; h--){
    refh = 0.0;
    for(i=h; i<n; i++){ /* works because w_i>=0 */
      refh += w[permut[i]];
    }
    if(refh>=ref){
      k=h;
      break;
    }
  }
  idx = permut[k];
  free(permut);
  return d[idx];
}


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

/** Filters EEG data using Fidlib from Jim Peters (http://uazu.net/fidlib/). 
	 Filtering is done in place.
	 \param eeg
	 \param sampling_rate in Hz
	 \param spec filter-specification as given in Fidlib
	 \code
	 The spec consists of a series of letters usually followed by the order
	 of the filter and then by any other parameters required, preceded by
	 slashes.  For example:
	 
	 LpBu4/20.4    Lowpass butterworth, 4th order, -3.01dB at 20.4Hz
	 BpBu2/3-4     Bandpass butterworth, 2nd order, from 3 to 4Hz
	 BpBu2/=3-4    Same filter, but adjusted exactly to the range given
	 BsRe/1000/10  Bandstop resonator, Q=1000, frequency 10Hz
	 
	 Note that filter frequencies should be specified in the same units as
	 the sampling rate, i.e. normally in Hz.  However, all filters work
	 internally with the ratio of freq/sampling_rate, so you can use a
	 sampling rate of 1.0 and frequency values of 0 to 0.5 (Nyquist) if you
	 prefer.

	 The auto-adjust option, selected by prefixing the frequency or frequency
	 range with '=', automatically adjusts the -3.01dB points to the given
	 frequencies with a kind of binary search.  
	 This might be useful with some lower-order filters where the
	 -3.01dB points don't come out quite right otherwise.
	 \endcode
	 Available Filters:
	 \code
	 BpRe/<value>/<freq>
    Bandpass resonator, Q=<value> (0 means Inf), frequency <freq>
	 BsRe/<value>/<freq>
    Bandstop resonator, Q=<value> (0 means Inf), frequency <freq>
	 ApRe/<value>/<freq>
    Allpass resonator, Q=<value> (0 means Inf), frequency <freq>
	 Pi/<freq>
    Proportional-integral filter, frequency <freq>
	 PiZ/<freq>
    Proportional-integral filter, matched z-transform, frequency <freq>
	 LpBe<order>/<freq>
    Lowpass Bessel filter, order <order>, -3.01dB frequency <freq>
	 HpBe<order>/<freq>
    Highpass Bessel filter, order <order>, -3.01dB frequency <freq>
	 BpBe<order>/<range>
    Bandpass Bessel filter, order <order>, -3.01dB frequencies <range>
	 BsBe<order>/<range>
    Bandstop Bessel filter, order <order>, -3.01dB frequencies <range>
	 LpBu<order>/<freq>
    Lowpass Butterworth filter, order <order>, -3.01dB frequency <freq>
	 HpBu<order>/<freq>
    Highpass Butterworth filter, order <order>, -3.01dB frequency <freq>
	 BpBu<order>/<range>
    Bandpass Butterworth filter, order <order>, -3.01dB frequencies <range>
	 BsBu<order>/<range>
    Bandstop Butterworth filter, order <order>, -3.01dB frequencies <range>
	 LpCh<order>/<value>/<freq>
    Lowpass Chebyshev filter, order <order>, passband ripple <value>dB,
	 -3.01dB frequency <freq>
	 HpCh<order>/<value>/<freq>
    Highpass Chebyshev filter, order <order>, passband ripple <value>dB,
	 -3.01dB frequency <freq>
	 BpCh<order>/<value>/<range>
    Bandpass Chebyshev filter, order <order>, passband ripple <value>dB,
	 -3.01dB frequencies <range>
	 BsCh<order>/<value>/<range>
    Bandstop Chebyshev filter, order <order>, passband ripple <value>dB,
	 -3.01dB frequencies <range>
	 LpBeZ<order>/<freq>
    Lowpass Bessel filter, matched z-transform, order <order>, -3.01dB
	 frequency <freq>
	 HpBeZ<order>/<freq>
    Highpass Bessel filter, matched z-transform, order <order>, -3.01dB
	 frequency <freq>
	 BpBeZ<order>/<range>
    Bandpass Bessel filter, matched z-transform, order <order>, -3.01dB
	 frequencies <range>
	 BsBeZ<order>/<range>
    Bandstop Bessel filter, matched z-transform, order <order>, -3.01dB
	 frequencies <range>
	 LpBuZ<order>/<freq>
    Lowpass Butterworth filter, matched z-transform, order <order>, -3.01dB
	 frequency <freq>
	 HpBuZ<order>/<freq>
    Highpass Butterworth filter, matched z-transform, order <order>, -3.01dB
	 frequency <freq>
	 BpBuZ<order>/<range>
    Bandpass Butterworth filter, matched z-transform, order <order>, -3.01dB
	 frequencies <range>
	 BsBuZ<order>/<range>
    Bandstop Butterworth filter, matched z-transform, order <order>, -3.01dB
	 frequencies <range>
	 LpChZ<order>/<value>/<freq>
    Lowpass Chebyshev filter, matched z-transform, order <order>, passband
	 ripple <value>dB, -3.01dB frequency <freq>
	 HpChZ<order>/<value>/<freq>
    Highpass Chebyshev filter, matched z-transform, order <order>, passband
	 ripple <value>dB, -3.01dB frequency <freq>
	 BpChZ<order>/<value>/<range>
    Bandpass Chebyshev filter, matched z-transform, order <order>, passband
	 ripple <value>dB, -3.01dB frequencies <range>
	 BsChZ<order>/<value>/<range>
    Bandstop Chebyshev filter, matched z-transform, order <order>, passband
	 ripple <value>dB, -3.01dB frequencies <range>
	 LpBuBe<order>/<value>/<freq>
    Lowpass Butterworth-Bessel <value>% cross, order <order>, -3.01dB
	 frequency <freq>
	 LpBq<optional-order>/<value>/<freq>
    Lowpass biquad filter, order <order>, Q=<value>, -3.01dB frequency <freq>
	 HpBq<optional-order>/<value>/<freq>
    Highpass biquad filter, order <order>, Q=<value>, -3.01dB frequency
	 <freq>
	 BpBq<optional-order>/<value>/<freq>
    Bandpass biquad filter, order <order>, Q=<value>, centre frequency <freq>
	 BsBq<optional-order>/<value>/<freq>
    Bandstop biquad filter, order <order>, Q=<value>, centre frequency <freq>
	 ApBq<optional-order>/<value>/<freq>
    Allpass biquad filter, order <order>, Q=<value>, centre frequency <freq>
	 PkBq<optional-order>/<value>/<value>/<freq>
    Peaking biquad filter, order <order>, Q=<value>, dBgain=<value>,
	 frequency <freq>
	 LsBq<optional-order>/<value>/<value>/<freq>
    Lowpass shelving biquad filter, S=<value>, dBgain=<value>, frequency
	 <freq>
	 HsBq<optional-order>/<value>/<value>/<freq>
    Highpass shelving biquad filter, S=<value>, dBgain=<value>, frequency
	 <freq>
	 LpBl/<freq>
    Lowpass Blackman window, -3.01dB frequency <freq>
	 LpHm/<freq>
    Lowpass Hamming window, -3.01dB frequency <freq>
	 LpHn/<freq>
    Lowpass Hann window, -3.01dB frequency <freq>
	 LpBa/<freq>
    Lowpass Bartlet (triangular) window, -3.01dB frequency <freq>

	 \endcode
*/
void eeg_filter_fidlib( EEGdata *eeg, double sampling_rate, const char *spec ){
  int len;
  int i,j;
  FidFilter *filt;
  char *filtspec, *orgspec;
  FidRun *run;
  FidFunc *funcp;
  void **fbuf;
  
  fbuf = (void**) malloc( eeg->nbchan*sizeof(void*));

  filtspec = (char*) malloc( (strlen(spec)+10)*sizeof(char));
  orgspec = filtspec;
  memset( filtspec, 0, (strlen(spec)+10)*sizeof(char));
  strcpy( filtspec, spec );
  fid_parse( sampling_rate, &filtspec, &filt);
  run= fid_run_new(filt, &funcp);

  len= fid_run_bufsize(run);
  for( i=0; i<eeg->nbchan; i++ ){
	 fbuf[i]= fid_run_newbuf(run);
  }
  dprintf("Buffer initialized\n");
  for( i=0; i<eeg->n; i++ ){
	 for( j=0; j<eeg->nbchan; j++ ){
		//		dprintf("Sample %i, channel %i\n", i, j );
		eeg->d[j][i] = funcp( fbuf[j], eeg->d[j][i] );
	 }
  }
  for( i=0; i<eeg->nbchan; i++ ){
	 fid_run_freebuf(fbuf[i]);
  }
  free( fbuf );
  fid_run_free(run);
  free( filt );
  free( orgspec ) ;
}


/* -------------------------------------------------------------------*/
/**\cond OBSOLETE */

/**\todo Not 100% certain that it works!

	reimplementation of Matlab's butter() function.
	 -# normalization of passband frequencies such that fs = 1
	 -# conversion to analog frequencies by
	    \f[
		  f_a = \frac{2}{T} \tan\left(\frac{fT}{2}\right)
		 \f]
		 where T is the sampling interval (here 1, because of normalization)
	 -# An analog Butterworth-Filter is designed in the 
 	    Laplace-Domain (S-Space)
	  - Filter coefficients are generated by (http://en.wikipedia.org/wiki/Butterworth_filter)
	  - order even:
	    \f$ i = 1 \ldots \frac {n}{2},
	        a_i = 2 \cos \frac{(2 i - 1) \pi}{2 n},
	        b_i = 1 \f$
	  - order odd 
	    \f$ i = 2 \ldots \frac{(n + 1)}{2},
	        a_1 = 1, b_1 = 0,
	        a_i = 2 \cos \frac{(i - 1) \pi}{n},
	        b_i = 1\f$
	  - This is a filter prototype
	 -# This prototype can be scaled/shifted to be a bandpass-filter
	    (http://en.wikipedia.org/wiki/Prototype_filter#Lowpass_to_bandpass)
	  - and \f[ \frac{i\omega}{\omega_c'} \to Q \left( \frac {i\omega}{\omega_0}+\frac {\omega_0}{i\omega} \right)\f]
	    where \f$ Q=\frac{\omega_0}{\Delta\omega}\f$ and \f$ \Delta\omega=\omega_2-\omega_1 \f$ and 
	    \f$ \omega_0=\sqrt{\omega_1\omega_2}\f$ 
	 -# It is converted to discrete Z-Domain (Z-space) using
	    the bilinear transform, such that
		 \f[
		 z = \frac{2+s}{2-s}
		 \f]

  \param order filter order
  \param sampling rate
  \param passband lower and upper cutoff-freq.
  \param bcoeff (2*order+1) double array for the (top) b-coefficients
  \param acoeff (2*order+1) double array for the (bottom) a-coefficients
 */
void butterworth_design_bandpass( int order, double sampling_rate, double passband[2],
											 double *bcoeff, double *acoeff ){
  double Bw, 						  /* bandwidth */
	 Cf;								  /* center freqeuncy */
  Complex *poles,					  /* poles in s-/z-plane */
	 *zeros;							  /* zeros in s-/z-plane */
  int i;
  Complex *coeff;					  /* dummy, since acoeff and bcoeff are real */

  dprintf("cutoff-frequencies (%f,%f) at %f Hz\n",
			 passband[0], passband[1], sampling_rate);

  /* step 0: normalize frequencies */
  passband[0] /= (sampling_rate);
  passband[1] /= (sampling_rate);

  dprintf( "passband=(%f,%f)\n",
			  passband[0], passband[1] );

  /* % step 1: get analog frequencies */
  passband[0] = 2.0*tan( PI*passband[0] );
  passband[1] = 2.0*tan( PI*passband[1] );

  Bw = passband[1] - passband[0];
  Cf = sqrt( passband[0]*passband[1] );

  dprintf( "Bw=%f, Cf=%f, passband=(%f,%f)\n",
			  Bw, Cf, passband[0], passband[1] );

  /* Get N-th order Butterworth analog lowpass prototype */
  /* It so happens that the poles of a Butterworth low-pass filter
	  with cut-off frequency wc are evenly-spaced around the
	  circumference of a half-circle of radius wc centered upon the
	  origin of the s-plane. The poles of a two-pole filter are at
	  +-45. Those of a four-pole filter are at +-22.5 and +-67.5
	  degrees.  */
  poles = (Complex*)malloc( 2*order*sizeof( Complex ) );
  zeros = (Complex*)malloc( 2*order*sizeof( Complex ) );

  if( ISODD(order) ){
  	 poles[order-1] = complex( -1.0, 0 );
  }
  for( i=0; i<order-1; i+=2 ){
  	 dprintf( "i=%i\n", i );
  	 poles[i]   = complex_exp ( complex( 0.0, (PI*(i+1))/(2*(double)order)+PI/2.0 ) );
  	 dprintf(" poles[%i] = %f + i %f\n", i, poles[i].re, poles[i].im );
  	 if( i+1<order ){
  		poles[i+1] = complex_conj( poles[i] );
  		dprintf(" poles[%i] = %f + i %f\n", i+1, poles[i+1].re, poles[i+1].im );
  	 }
  }

  /* http://en.wikipedia.org/wiki/Prototype_filter#Lowpass_to_bandpass */
  /* warp the filter to central freq and passband */
  Complex hba,m,temp;
  for (i=0; i < order; i++){
  	 hba = complex_mul_double(poles[i], 0.5*Bw);
  	 m = complex_div( complex( Cf, 0.0 ), hba );
  	 temp = complex_sqrt(complex_sub(complex(1.0, 0.0), complex_mul( m, m ) ) );
  	 poles[i] = complex_mul( hba, complex(temp.re+1.0,temp.im) );
  	 poles[order+i] = complex_mul( hba,  complex_sub( complex(1.0,0.0), temp ) );
  }

  /* Converting S-Plane poles to Z-Plane poles using bilinear transform */
  for( i=0; i<order*2; i++ ){
  	 dprintf("S-plane: poles[%i]=%f +i %f\n", i, poles[i].re, poles[i].im );
  	 poles[i] =( complex_bilinear_transform( poles[i] ) );
  	 dprintf("Z-plane: poles[%i]=%f +i %f\n", i, poles[i].re, poles[i].im );
  }

  /* compute the zeros which is just order x -1 and order x +1 */
  for( i=0; i<order; i++ ){
  	 zeros[i] = complex( 1.0, 0.0 );
  	 zeros[order+i] = complex( -1.0, 0.0 );
  }

  /* /\* expand the polynomials *\/ */
  /* /\* initialize *\/ */
  coeff = (Complex*) malloc( (2*order+1)*sizeof(Complex) ); 
  
  expand_polynomial_from_roots( poles, 2*order, coeff );
  /* /\* check computed coeffs of z^k are all real and put them into bcoeffs *\/ */
  vector_complex_to_real( coeff, bcoeff, 2*order+1 );
  
  /* expand_polynomial_from_roots( zeros, 2*order, coeff ); */
  vector_complex_to_real( coeff, acoeff, 2*order+1 );

  /* normalize */
  for (i = 0; i < 2*order+1; i++){
  	 bcoeff[i] = +(bcoeff[i] / acoeff[2*order+1]);
  	 acoeff[i] = -(acoeff[i] / acoeff[2*order+1]);
  }

  /* cleaning up */
  free( poles );
  free( coeff );
  free( zeros );
}

/** \todo
	 Directly apply the filter given in bcoeff/acoeff (polynomial coefficients of the transfer function)
	 to signal. For FIR-filters, set the recursive coefficients to 0.
	 That is, compute
	 \f[
	 y_n = \frac{1}{b_0}\left( a_0 x_n + \ldots + a_n x_0 - b_1 y_{n-1} -\ldots - b_n y_0 \right)
	 \f]
	 \param bcoeff recursive coefficients in array of length nb
	 \param acoeff input coefficients in array of length nb
	 \param signal to-be-filtered signal of length n
 */
void filter_direct( const double *bcoeff, int nb, 
						  const double *acoeff, int na, double *signal, int n ){
  int i;
  for( i=0; i<n; i++ ){
	 
  }
}

/** estimate kaiser's alpha-value and a good filter order M,
	 given stop-band attenuation in dB - A
	 http://ccrma.stanford.edu/~jos/sasp/case_you_don_t_have.html
	 \param passband
	 \param A  stop-band attenuation in dB 
	 \param samplingrate of the signal
	 \param alpha (output)
	 \param M (output)
*/
void estimate_kaiseralpha_bandpass( double passband[2], double A, double samplingrate,
												double *alpha, int *M ){
  double f1,f2;

  if( A>50 ){
	 *alpha = 0.1102*(A-8.7);
  } else if( A>21 ) {
	 *alpha = 0.5842*pow( (A-21), 0.4 ) + 0.07886*(A-21 );
  } else {
	 *alpha=0.0;
  }

  f1 = passband[0]/samplingrate;
  f2 = passband[1]/samplingrate;
  dprintf("f1=%f, f2=%f\n", f1, f2 );

  *M = (A-8)/( 2.285*2.0*PI*ABS( f2-f1 ) );
}

/** FIR-filter design using the window method.
	 http://ccrma.stanford.edu/~jos/sasp/Window_Method.html
	 
 */
double* fir_bandpass_design_window_method(  ){
  return NULL;
}

/**\endcond */
