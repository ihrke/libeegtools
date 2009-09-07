/* **************************************************************************
 *   Copyright (C) 2009 by Matthias Ihrke   *
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

#include "filter.h"
#include "distances.h"
#include "fidlib/fidlib.h"

/* ---------------------------------------------------------------------------- 
   -- Denoising routines                                                     -- 
   ---------------------------------------------------------------------------- */
/** Simple Running Median filter.
 * Is applied to each channel and trial separately.
 \param s input eeg
 \param win the window size in which to calculate the median
 \param alloc TRUE: a new struct is returned; FALSE: in-place modification
 \return filtered EEG
 */
EEG* eeg_filter_running_median(EEG *eeg, int win, Boolean alloc){
  EEG *s;
  int c,i;
  if( alloc ){
	 s = eeg_clone( eeg );
  } else {
	 s = eeg;
  }
  for(c=0; c<s->nbchan; c++){
	 for(i=0; i<s->ntrials; i++){
		running_median(s->data[c][i], s->n, win);
	 }
  }
  return s;
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
