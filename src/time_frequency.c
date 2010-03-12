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

#include "time_frequency.h"

/** Easy interface for spectrogram calculation.
	 Calculates also the power spectrum.

	 \param s - data to be spectrogrammified (real function)
	 \param n - length(s)
	 \param spectgram - if NULL, own memory is alloated, else use this pointer.
	 \param optargs may contain:
	 - <tt>sample_frequency=double</tt> of the signal, default is \c 500 (should really be provided!)
	 - <tt>winfct=void*</tt> windowing function, default is \c window_hanning
	 - <tt>timepoints=int*</tt> compute the spectrogram at selected time-points; this array is N_time long
	 - <tt>winlength=int</tt> size of the window, default is <tt> MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 )</tt>
	 - <tt>N_freq=int</tt> number of frequency bins, default is \c winlength*4
	 - <tt>N_time=int</tt> number of time bins, default is \c n
	 - <tt>corner_freqs=double*</tt> (array with two double entries), default is \c (0.0,250.0)
 */
Spectrogram* spectrogram( const double *s, int n, Spectrogram *spectgram, 
								  OptArgList *optargs ){
  double x;
  void *ptr;
  double sample_frequency;
  WindowFunction winfct; 
  int winlength;
  int N_freq;
  int N_time;
  double corner_freqs[2];
  int *timepoints;

  /* defaults */
  sample_frequency = 500.0;
  winfct = window_hanning;
  winlength =  MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 );
  N_freq = winlength*4;
  N_time = n;
  corner_freqs[0] = 0.0;
  corner_freqs[1] = 250.0;
  timepoints = NULL;
  bool N_time_set=FALSE;

  /* get params */
  if( optarglist_has_key( optargs, "sample_frequency" ) ){
	 x = optarglist_scalar_by_key( optargs, "sample_frequency" );
	 if( !isnan( x ) ) sample_frequency=x;
  }
  if( optarglist_has_key( optargs, "winfct" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "winfct" );
	 if( ptr ) winfct = (WindowFunction)ptr;
  }
  if( optarglist_has_key( optargs, "winlength" ) ){
	 x = optarglist_scalar_by_key( optargs, "winlength" );
	 if( !isnan( x ) ) winlength=(int)x;
  }
  if( optarglist_has_key( optargs, "N_freq" ) ){
	 x = optarglist_scalar_by_key( optargs, "N_freq" );
	 if( !isnan( x ) ) N_freq=(int)x;
  }
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


  dprintf("got sfreq=%f, winfct=%p, winlength=%i, N_f=%i, N_t=%i, cf=(%f,%f), timepoints=%p\n",
			 sample_frequency, winfct, winlength, N_freq, N_time, corner_freqs[0], corner_freqs[1],
			 timepoints );

  double *window=winfct( ALLOC_IN_FCT, winlength );

  spectgram = spectrogram_stft( s, n, sample_frequency, 
										  window, winlength, N_freq,
										  N_time, corner_freqs, timepoints,
										  spectgram );

  spectrogram_compute_powerspectrum( spectgram );
  return spectgram;
}


/** This function is inspired by (read: 'was shamelessly ripped of from') 
	 the TIME-FREQUENCY TOOLBOX by  Emmanuel Roy - Manuel DAVY 
	 (http://www-lagis.univ-lille1.fr/~davy/toolbox/Ctftbeng.html)
	 It's a bit more efficient and uses different storage formats.
	 Also it doesn't work with arbitrary time-points for the TFR and not on
	 complex signals.
	 
	 THE ALGORITHM 
	 Given a signal to analyze in time and frequency, computes the 
	 Short Time Fourier Transform (STFT) :			                      
	 \f[
	 {STFT}(t,f) = \int x(s)h(t-s)e^{-2i\pi f s} ds
	 \f]

	 This function is complex valued. Its computation requires a window, which 
	 can be computed with one of the window_*() functions.
	 \param s - data to be spectrogrammified (real function)
	 \param n - number of samples in s
	 \param sampling_rate - sampling frequency of signal
	 \param N_freq -  number of frequency bins = number of rows in the TFR matrix
	 \param N_time - number of cols in the TFR matrix        
	 \param corner_freqs - corner frequencies (lower, upper) for the returned
	                       spectrum at each time-sample in Hz; maximal would be
								  {0, srate/2} 
    \param timepoints compute the spectrogram at selected time-points; this array is N_time long
	 \param spectgram - if NULL, own memory is alloated, else use this pointer.
 */

Spectrogram* spectrogram_stft(const double* s, int n, double sampling_rate,
										const double *Window, int Window_Length,
										int N_freq, int N_time, double corner_freqs[2], 
										int *timepoints, Spectrogram *spectgram){

  int            Nfft, timepoint, frequency, time;
  int            taumin, taumax, tau;
  int            half_Window_Length;
  double        *wind_sig_complex;		/* windowed signal */
  double         normh;
  double         inter;
  int            corner_freqs_idx[2]; /* nearest sample in discretized setting */
  int            N_freq_spect;	  /* number of frequencies in spectrgram */

  dprintf("n=%i, sf=%f, wl=%i, Nf=%i, Nt=%i\n", n, sampling_rate, Window_Length, N_freq, N_time);
 /*--------------------------------------------------------------------*/
 /*                   checks that the window length is odd             */
 /*--------------------------------------------------------------------*/
  if ( !ISODD(Window_Length) ){
	 errprintf ("The window Length must be an ODD number, got %i\n", Window_Length);
	 return NULL;
  }
  
  /* normalizing to [0, 100] */
  half_Window_Length = (Window_Length - 1) / 2;
  inter = 0.0;
  for (frequency = 0; frequency <Window_Length; frequency++){
	 inter += SQR (Window[frequency]);
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
  
  /*--------------------------------------------------------------------*/
  /*                memory allocation for the windowed signal           */
  /*--------------------------------------------------------------------*/
  wind_sig_complex = (double *) calloc (2*Nfft, sizeof (double));

  if( spectgram==NULL ){
	 spectgram = spectrogram_init( N_freq_spect, N_time );
	 dprintf(" spectgram=%p, real[0]=%f\n", spectgram, spectgram->sgram[0][0].re);
  }
  spectgram->low_corner_freq = corner_freqs[0];
  spectgram->up_corner_freq  = corner_freqs[1];

  /*--------------------------------------------------------------------*/
  /*      computation of the fft for the current windowed signal        */
  /*--------------------------------------------------------------------*/
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
		wind_sig_complex[2*frequency  ] = s[time + tau] * Window[half_Window_Length + tau] / normh;
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
 /*                free the memory used in this program                */
 /*--------------------------------------------------------------------*/
  free (wind_sig_complex);
  dprintf("done\n");

  return spectgram;
}


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
  s->has_power_spectrum=0;
  s->sgram = (Complex**) malloc( N_time*sizeof(Complex*) );
  s->powerspect = (double**) malloc( N_time*sizeof(double*) );
  for( i=0; i<N_time; i++ ){
	 s->sgram[i] = (Complex*) malloc( N_freq*sizeof(Complex) );
	 s->powerspect[i] = vector_init( NULL, N_freq, -1.0 ); /*(double*) malloc( N_freq*sizeof(double) );*/
  }
  return s;
}

void spectrogram_free( Spectrogram *s ){
  int i;

  for( i=0; i<s->N_time;  i++ ){
	 free( s->sgram[i] );
	 free( s->powerspect[i] );
  }
  free( s->sgram );
  free( s->powerspect );
  free( s );
}


/** fill window with data for a Dirichlet (Rectangular) window.
	 \f[
	 w(x) = 1;
	 \f]
	 \param window NULL (own memory allocation) or own pointer 
	 \param n - ODD number for window
 */
double* window_dirichlet( double *window, int n ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  if( window==NULL ){
	 window = (double*) malloc( n*sizeof(double) );
  }
  for( i=0; i<n; i++ ){
	 window[i] = 1.0;
  }

  return window;
}

/** fill window with data for a Gaussian window.
	 \f[
	 w(n) = e^{-\frac{1}{2}\left( \frac{n-(N-1)/2}{\sigma(N-1)/2}\right)^2 }
	 \f]
	 with \f$ \sigma\le 0.5\f$
	 \param window NULL (own memory allocation) or own pointer 
	 \param n - ODD number for window
	 \param sigma - gaussian std
 */
double* window_gaussian( double *window, int n, double sigma ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  if( window==NULL ){
	 window = (double*) malloc( n*sizeof(double) );
  }
  for( i=0; i<n; i++ ){
	 window[i] = exp(-0.5*SQR( (i-(n-1)/2.0)/(double)(sigma*(n-1)/2.0)) );
  }

  return window;
}

/** fill window with data for a Hamming window.
	 \f[
	 w(n) = 0.53836 - 0.46164 \cos\left( \frac{2\pi n}{N-1}\right)
	 \f]
	 \param window NULL (own memory allocation) or own pointer 
	 \param n - ODD number for window
*/
double* window_hamming( double *window, int n ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  if( window==NULL ){
	 window = (double*) malloc( n*sizeof(double) );
  }
  for( i=0; i<n; i++ ){
	 window[i] =  0.53836 - 0.46164*cos( (2* PI *i)/(double)(n-1) );
  }

  return window;
}
/** fill window with data for a Hanning (Hann) window.
	 \f[
	 w(n) = 0.5\left(1- \cos\left( \frac{2\pi n}{N-1}\right)\right)
	 \f]
	 \param window NULL (own memory allocation) or own pointer 
	 \param n - ODD number for window
*/
double* window_hanning( double *window, int n ){
  int i;

  if( !ISODD( n ) ){
	 errprintf("window size must be ODD! Adding 1, n=%i!\n", n+1);
	 n++;
  }
  if( window==NULL ){
	 window = (double*) malloc( n*sizeof(double) );
  }
  for( i=0; i<n; i++ ){
	 window[i] =  0.5*(1 - cos( (2* PI *i)/(double)(n-1) ));
  }

  return window;
}


/** fill window with data for a Kaiser Window
	 http://en.wikipedia.org/wiki/Kaiser_window

	 \param window NULL (own memory allocation) or own pointer 
	 \param n - ODD number for window
	 \param alpha - parameter for window's steepness
*/
double* window_kaiser( double *window, int n, double alpha ){
  int i;
  double bessel_alpha;

  if( window==NULL ){
	 window = (double*) malloc( n*sizeof(double) );
  }

  bessel_alpha=gsl_sf_bessel_I0( alpha );
  dprintf("alpha=%f, b_alpha=%f\n", alpha, bessel_alpha );
  for( i=0; i<n-1; i++ ){
	 window[i] = gsl_sf_bessel_I0( alpha*sqrt(1 - SQR( (2*i)/(double)(n-1) - 1  ) ) ) / bessel_alpha;
	 dprintf("w[%i]=%f\n", i, window[i]);
  }
  window[n-1]=0;

  return window;
}

/** Fills field 'powerspect' from a Spectrogram struct with |x|^2 of 
	 the complex numbers in spectrogram.
 */
void         spectrogram_compute_powerspectrum( Spectrogram *s ){
  int i,j;
  for( i=0; i<s->N_time; i++ ){
	 for( j=0; j<s->N_freq; j++ ){
		s->powerspect[i][j] = SQR( complex_abs( s->sgram[i][j] ) );
	 }
  }
  s->has_power_spectrum=1;
}
