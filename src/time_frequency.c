/***************************************************************************
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
#include "helper.h"
#include "mathadd.h"
#include "time_frequency.h"
#include "linalg.h"

/** \brief Easy interface for spectrogram calculation.

	 \param sig - data to be spectrogrammified (1D DOUBLE array)
	 \param spectgram - if NULL, own memory is alloated, else use this pointer.
	 \param optargs may contain:
	 - <tt>sample_frequency=double</tt> of the signal, default is \c 500 (should really be provided!)
	 - <tt>winfct=void*</tt> windowing function, default is \c window_hanning
	 - <tt>timepoints=int*</tt> compute the spectrogram at selected time-points; this array is N_time long
	 - <tt>winlength=int</tt> size of the window, default is <tt> MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 )</tt>
	 - <tt>winparam=double</tt> parameter for the window-function (e.g. sigma for gaussian); default=0.0
	 - <tt>N_freq=int</tt> number of frequency bins, default is \c winlength*4
	 - <tt>N_time=int</tt> number of time bins, default is \c n
	 - <tt>corner_freqs=double*</tt> (array with two double entries), default is \c (0.0,250.0)
 */
Spectrogram* spectrogram( const Array *sig, Spectrogram *spectgram, 
								  OptArgList *optargs ){
  double x;
  void *ptr;
  double sample_frequency;
  WindowFunction winfct; 
  int winlength;
  double winparam;
  int N_freq;
  int N_time;
  double corner_freqs[2];
  int *timepoints;

  /* check input */
  bool isvec;
  vector_CHECK( isvec, sig );
  if( !isvec ) return NULL;

  int n=sig->size[0];
  /* defaults */
  sample_frequency = 500.0;
  winfct = window_hanning;
  winlength =  MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 );
  winparam=0.0;
  N_freq = winlength*4;
  N_time = n;
  corner_freqs[0] = 0.0;
  corner_freqs[1] = 250.0;
  timepoints = NULL;
  bool N_time_set=FALSE;

  /* get params */
  optarg_PARSE_SCALAR( optargs, "sample_frequency", sample_frequency, double, x );
  optarg_PARSE_PTR( optargs, "winfct", winfct, WindowFunction, ptr );
  optarg_PARSE_SCALAR( optargs, "winlength", winlength, int, x );
  optarg_PARSE_SCALAR( optargs, "winparam", winparam, double, x );
  optarg_PARSE_SCALAR( optargs, "N_freq", N_freq, int, x );
  if( optarglist_has_key( optargs, "N_time" ) ){
	 x = optarglist_scalar_by_key( optargs, "N_time" );
	 if( !isnan( x ) ){
		N_time=(int)x;
		N_time_set=TRUE;
	 }
  }
  if( optarglist_has_key( optargs, "corner_freqs" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "corner_freqs" );
	 if( ptr ){
		corner_freqs[0] = ((double*)ptr)[0];
		corner_freqs[1] = ((double*)ptr)[1];
	 }
  }
  if( optarglist_has_key( optargs, "timepoints" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "timepoints" );
	 if( ptr && N_time_set ){
		timepoints = (int*)ptr;
	 }
  }

  dprintf("got sfreq=%f, winfct=%p, winlength=%i, winparam=%f, N_f=%i,"
			 "N_t=%i, cf=(%f,%f), timepoints=%p\n",
			 sample_frequency, winfct, winlength, winparam, N_freq, N_time, 
			 corner_freqs[0], corner_freqs[1], timepoints );

  Array *window=winfct( winlength, winparam );

  spectgram = spectrogram_stft( sig, sample_frequency, 
										  window, N_freq, N_time, corner_freqs, 
										  timepoints, spectgram );
  return spectgram;
}


/** \brief Calculate spectrogram of a signal (using STFT).

	 This function is inspired by (read: 'was shamelessly ripped of from') 
	 the TIME-FREQUENCY TOOLBOX by  Emmanuel Roy - Manuel DAVY 
	 (http://www-lagis.univ-lille1.fr/~davy/toolbox/Ctftbeng.html)
	 It's a bit more efficient and uses different storage formats.
	 It doesn't work on complex signals.
	 
	 THE ALGORITHM 

	 Given a signal to analyze in time and frequency, computes the 
	 Short Time Fourier Transform (STFT) :			                      
	 \f[
	 {STFT}(t,f) = \int x(s)h(t-s)e^{-2i\pi f s} ds
	 \f]

	 This function is complex valued. Its computation requires a window, which 
	 can be computed with one of the window_*() functions.

	 \param sig - data to be spectrogrammified (1D DOUBLE array)
	 \param sampling_rate - sampling frequency of signal
	 \param window - window for calculating the STFT (1D DOUBLE Array); get from
	                  window_*() functions
	 \param N_freq -  number of frequency bins = number of rows in the TFR matrix (the next
	                  power of 2 is chosen for N_freq)
	 \param N_time - number of cols in the TFR matrix        
	 \param corner_freqs - corner frequencies (lower, upper) for the returned
	                       spectrum at each time-sample in Hz; maximal would be
								  {0, srate/2} 
    \param timepoints compute the spectrogram at selected time-points; this array is N_time long
	 \param spectgram - if NULL, own memory is alloated, else use this pointer.
	 \return the spectrogram
 */

Spectrogram* spectrogram_stft(const Array* sig, double sampling_rate,
										const Array *Window, int N_freq, int N_time, 
										double corner_freqs[2], 
										int *timepoints, Spectrogram *spectgram){

  int            Nfft, timepoint, frequency, time;
  int            taumin, taumax, tau;
  int            half_Window_Length;
  double        *wind_sig_complex;		/* windowed signal */
  double         normh;
  double         inter;
  int            corner_freqs_idx[2]; /* nearest sample in discretized setting */
  int            N_freq_spect;	  /* number of frequencies in spectrgram */
  int            Window_Length;

  /* check for 1D and DOUBLE */
  bool isvec;
  vector_CHECK( isvec, sig );
  if( !isvec ) return NULL;
  vector_CHECK( isvec, Window );
  if( !isvec ) return NULL;
  
  int n=sig->size[0];
  Window_Length=Window->size[0];

  /* normalizing to [0, 100] */
  half_Window_Length = (Window_Length - 1) / 2;
  inter = 0.0;
  for (frequency = 0; frequency <Window_Length; frequency++){
	 inter += SQR (vec_IDX(Window,frequency));
  }
  normh = sqrt (inter);

  /* FFT-settings */
  Nfft = pow( 2, next_pow2( N_freq ));	  /* num of samples for FFT */
  if( corner_freqs[1]<=corner_freqs[0] || 
		corner_freqs[0]<0 || corner_freqs[1]>sampling_rate/2.0 ){
	 errprintf("corner_freqs funny: (%f,%f)\n", 
				  corner_freqs[0], corner_freqs[1] );
	 return NULL;
  }
  corner_freqs_idx[0] = (int)round( (corner_freqs[0]/(sampling_rate/2.0)) * (N_freq-1) );
  corner_freqs_idx[1] = (int)round( (corner_freqs[1]/(sampling_rate/2.0)) * (N_freq-1) );
  dprintf(" Corner-freqs=(%.2f,%.2f) in idx=(%i,%i)\n", corner_freqs[0], corner_freqs[1],
			 corner_freqs_idx[0], corner_freqs_idx[1] );
  dprintf(" Returning actual frequencies: (%.3f, %.3f)\n", 
			 (corner_freqs_idx[0]/(double)N_freq)*sampling_rate/2.0, 
			 (corner_freqs_idx[1]/(double)N_freq)*sampling_rate/2.0 ); 
  N_freq_spect = corner_freqs_idx[1]-corner_freqs_idx[0];
  

  /*                memory allocation for the windowed signal           */
  wind_sig_complex = (double *) calloc (2*Nfft, sizeof (double));

  if( spectgram==NULL ){
	 spectgram = spectrogram_init( N_freq_spect, N_time );
  }
  spectgram->low_corner_freq = corner_freqs[0];
  spectgram->up_corner_freq  = corner_freqs[1];

  /* computation of the fft for the current windowed signal */
  for (timepoint = 0; timepoint < N_time; timepoint++) {
	 /* initialization of the intermediary vector */
	 for (frequency = 0; frequency < 2*Nfft; frequency++){
		wind_sig_complex[frequency] = 0.0;
	 }

	 /* current time to compute the stft */
	 if( timepoints ){
		time = timepoints[timepoint];
	 } else {
		time = timepoint*(int)(n/N_time);
	 }
	 /* dprintf( "Spectrum at sample '%i'\n", time ); */

	 /* the signal is multipied by the window between the instants
		 time-taumin and time+taumax 
	    when the window is wider than the number of desired frequencies (tfr.N_freq),
		 the range is limited to tfr.N_freq */
	 taumin = MIN (N_freq / 2, half_Window_Length);
	 taumin = MIN (taumin, time);

	 taumax = MIN ((N_freq / 2 - 1), half_Window_Length);
	 taumax = MIN (taumax, (n - time - 1));

	 /* dprintf(" tau = (%i, %i)\n", taumin, taumax); */
	 /* The signal is windowed around the current time */
	 frequency=0;
	 for (tau = -taumin; tau <= taumax; tau++){
		wind_sig_complex[2*frequency  ] = vec_IDX(sig,time + tau) *
		  vec_IDX( Window, half_Window_Length + tau ) / normh;
		wind_sig_complex[2*frequency+1] = 0.0;
		frequency++;
	 }
	 /* fft of the windowed signal */
	 fft( wind_sig_complex, Nfft, 1 );  
	 /* gsl_fft_complex_radix2_forward ( wind_sig_complex, 1, Nfft );  */

	 /* the first half of the fft is put in the stft matrix  */
	 for (frequency = corner_freqs_idx[0]; frequency < corner_freqs_idx[1]; frequency++){
		spectgram->sgram[timepoint][frequency].re = wind_sig_complex[2*frequency  ];
		spectgram->sgram[timepoint][frequency].im = wind_sig_complex[2*frequency+1];
	 }
  } /* for t */

 /*--------------------------------------------------------------------*/
  free (wind_sig_complex);

  return spectgram;
}


/** \brief Calculate the powerspectrum for a Spectrogram. 

	 Returns a s->N_time x s->N_freq matrix the absolute
	 value squared of the complex numbers in spectrogram.

	 \param s - the spectrogram
	 \return the powerspectrum for all time-points
 */
Array* spectrogram_powerspectrum( const Spectrogram *s ){
  int i,j;
  Array *p = array_new2( DOUBLE, 2, s->N_time, s->N_freq );
  
  for( i=0; i<s->N_time; i++ ){
	 for( j=0; j<s->N_freq; j++ ){
		mat_IDX( p, i,j) = SQR( complex_abs( s->sgram[i][j] ) );
	 }
  }
  return p;
}

/**\brief Initialize Spectrogram.

	\param N_freq resolution (number of points) in frequency
	\param N_time resolution in time
	\return a freshly allocated Spectrogram struct
 */
Spectrogram* spectrogram_init( int N_freq, int N_time ){
  Spectrogram *s;
  int i;
  
  dprintf("N_freq=%i, N_time=%i\n", N_freq, N_time );
  s = (Spectrogram*) malloc( sizeof(Spectrogram) );
  s->samplingrate = 500.0;
  s->N_freq = N_freq;
  s->N_time = N_time;
  s->low_corner_freq=0;
  s->up_corner_freq=0;
  s->sgram = (Complex**) malloc( N_time*sizeof(Complex*) );
  for( i=0; i<N_time; i++ ){
	 s->sgram[i] = (Complex*) malloc( N_freq*sizeof(Complex) );
  }
  return s;
}

/**\brief Free memory associated with Spectrogram.
 */
void spectrogram_free( Spectrogram *s ){
  int i;

  for( i=0; i<s->N_time;  i++ ){
	 free( s->sgram[i] );
  }
  free( s->sgram );
  free( s );
}

/*----------------- WINDOW-FUNCTIONS ------------------ */

/** \brief Dirichlet (Rectangular) window.

	 \f[
	 w(x) = 1;
	 \f]
	 \param n - ODD number for window (else it is incremented by one)
	 \param noparam ignored (just for comparability to \ref WindowFunction)
	 \return 1D DOUBLE array of size n
 */
Array* window_dirichlet(  int n, double noparam ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  Array *window = array_new2( DOUBLE, 1, n );
  for( i=0; i<n; i++ ){
	 vec_IDX( window,i ) = 1.0;
  }

  return window;
}

/** \brief Gaussian window.

	 \f[
	 w(n) = e^{-\frac{1}{2}\left( \frac{n-(N-1)/2}{\sigma(N-1)/2}\right)^2 }
	 \f]
	 with \f$ \sigma\le 0.5\f$

	 \param n - ODD number for window (else it is incremented by one)	
	 \param sigma - gaussian std
	 \return 1D DOUBLE array of size n
 */
Array* window_gaussian( int n, double sigma ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }

  Array *window = array_new2( DOUBLE, 1, n );
  for( i=0; i<n; i++ ){
	 vec_IDX(window,i) = exp(-0.5*SQR( (i-(n-1)/2.0)/(double)(sigma*(n-1)/2.0)) );
  }

  return window;
}

/** \brief Hamming window.

	 \f[
	 w(n) = 0.53836 - 0.46164 \cos\left( \frac{2\pi n}{N-1}\right)
	 \f]
	 \param n - ODD number for window (else it is incremented by one)	
	 \param noparam ignored (just for comparability to \ref WindowFunction)
	 \return 1D DOUBLE array of size n
*/
Array* window_hamming( int n, double noparam ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }

  Array *window = array_new2( DOUBLE, 1, n );
  for( i=0; i<n; i++ ){
	 vec_IDX(window,i) =  0.53836 - 0.46164*cos( (2* PI *i)/(double)(n-1) );
  }

  return window;
}
/** \brief Hanning (Hann) window.

	 \f[
	 w(n) = 0.5\left(1- \cos\left( \frac{2\pi n}{N-1}\right)\right)
	 \f]

	 \param n - ODD number for window (else it is incremented by one)	
	 \param noparam ignored (just for comparability to \ref WindowFunction)
	 \return 1D DOUBLE array of size n
*/
Array* window_hanning( int n, double noparam ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  Array *window = array_new2( DOUBLE, 1, n );
  for( i=0; i<n; i++ ){
	 vec_IDX( window,i) =  0.5*(1 - cos( (2* PI *i)/(double)(n-1) ));
  }

  return window;
}


/** \brief Kaiser Window.

	 http://en.wikipedia.org/wiki/Kaiser_window

	 \param n - ODD number for window (else it is incremented by one)	
	 \param alpha - parameter for window's steepness
	 \return 1D DOUBLE array of size n
*/
Array* window_kaiser( int n, double alpha ){
  int i;
  double bessel_alpha;

  Array *window = array_new2( DOUBLE, 1, n );

  bessel_alpha=gsl_sf_bessel_I0( alpha );
  for( i=0; i<n-1; i++ ){
	 vec_IDX(window,i) = gsl_sf_bessel_I0( alpha*sqrt(1 - SQR( (2*i)/(double)(n-1) - 1  ) ) ) / bessel_alpha;
  }
  vec_IDX(window,n-1)=0;

  return window;
}

