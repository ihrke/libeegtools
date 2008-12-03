#include "time_frequency.h"

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
	 \param sample_freq - sampling frequency of signal
	 \param N_freq -  number of frequency bins = number of rows in the TFR matrix
	 \param N_time - number of cols in the TFR matrix        
	 \param spectgram - if NULL, own memory is alloated, else use this pointer.
 */

Spectrogram* spectrogram_stft(const double* s, int n, double sample_freq,
										const double *Window, int Window_Length,
										int N_freq, int N_time,
										Spectrogram *spectgram){

  int            Nfft, timepoint, frequency, time;
  int            taumin, taumax, tau;
  int            half_Window_Length;
  double        *wind_sig_complex;		/* windowed signal */
  double         normh;
  double         inter;

  dprintf("n=%i, sf=%f, wl=%i, Nf=%i, Nt=%i\n", n, sample_freq, Window_Length, N_freq, N_time);
 /*--------------------------------------------------------------------*/
 /*                   checks that the window length is odd             */
 /*--------------------------------------------------------------------*/
  if ( !ISODD(Window_Length) ){
	 errprintf ("The window Length must be an ODD number, got %i\n", Window_Length);
	 return NULL;
  }
  
  half_Window_Length = (Window_Length - 1) / 2;
  inter = 0.0;
  for (frequency = 0; frequency <Window_Length; frequency++){
	 inter = inter + SQR (Window[frequency]);
  }
  normh = sqrt (inter);

  Nfft = next_pow2( N_freq );
  
  /*--------------------------------------------------------------------*/
  /*                memory allocation for the windowed signal           */
  /*--------------------------------------------------------------------*/
  wind_sig_complex = (double *) calloc (2*Nfft, sizeof (double));

  if( spectgram==NULL ){
	 spectgram = init_spectrogram( N_freq, N_time );
	 dprintf(" spectgram=%p, real[0]=%f\n", spectgram, spectgram->real[0][0]);
  }


  /*--------------------------------------------------------------------*/
  /*      computation of the fft for the current windowed signal        */
  /*--------------------------------------------------------------------*/
  for (timepoint = 0; timepoint < N_time; timepoint++) {
	 /* initialization of the intermediary vectors */
	 for (frequency = 0; frequency < 2*Nfft; frequency++){
		wind_sig_complex[frequency] = 0.0;
	 }

	 /* current time to compute the stft */
	 /* time = ((int) tfr.time_instants[column]) - 1;*/
	 time = timepoint*(int)(n/N_time);
	 dprintf( "Spectrum at sample '%i'\n", time );

	 /* the signal is multipied by the window between the instants
		 time-taumin and time+taumax 
	    when the window is wider than the number of desired frequencies (tfr.N_freq),
		 the range is limited to tfr.N_freq */
	 taumin = MIN (N_freq / 2, half_Window_Length);
	 taumin = MIN (taumin, time);

	 taumax = MIN ((N_freq / 2 - 1), half_Window_Length);
	 taumax = MIN (taumax, (n - time - 1));

	 /* The signal is windowed around the current time */
	 for (tau = -taumin; tau <= taumax; tau++){
		frequency = iremainder( (N_freq+tau), N_freq ) ;
		wind_sig_complex[2*frequency  ] = s[time + tau] * Window[half_Window_Length + tau] / normh;
		wind_sig_complex[2*frequency+1] = 0.0;
	 }
	 /* fft of the windowed signal */
	 fft ( wind_sig_complex, Nfft, 1 );

	 /* the first half of the fft is put in the stft matrix  */
	 for (frequency = 0; frequency < N_freq; frequency++){
		spectgram->real[timepoint][frequency] = wind_sig_complex[2*frequency  ];
		spectgram->imag[timepoint][frequency] = wind_sig_complex[2*frequency+1];
	 }
  } /* for t */

 /*--------------------------------------------------------------------*/
 /*                free the memory used in this program                */
 /*--------------------------------------------------------------------*/
  free (wind_sig_complex);

  return spectgram;
}


Spectrogram* init_spectrogram( int N_freq, int N_time ){
  Spectrogram *s;
  int i;
  
  dprintf("N_freq=%i, N_time=%i\n", N_freq, N_time );
  s = (Spectrogram*) malloc( sizeof(Spectrogram) );
  s->N_freq = N_freq;
  s->N_time = N_time;
  s->has_power_spectrum=0;
  s->real = (double**) malloc( N_time*sizeof(double*) );
  s->imag = (double**) malloc( N_time*sizeof(double*) );
  s->powerspect = (double**) malloc( N_time*sizeof(double*) );
  for( i=0; i<N_time; i++ ){
	 s->real[i] = (double*) malloc( N_freq*sizeof(double) );
	 s->imag[i] = (double*) malloc( N_freq*sizeof(double) );
	 s->powerspect[i] = vector_init( NULL, N_freq, -1.0 ); /*(double*) malloc( N_freq*sizeof(double) );*/
  }
  return s;
}

void free_spectrogram( Spectrogram *s ){
  int i;

  for( i=0; i<s->N_time;  i++ ){
	 free( s->real[i] );
	 free( s->imag[i] );
	 free( s->powerspect[i] );
  }
  free( s->real );
  free( s->imag );
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
/** Fills field 'powerspect' from a Spectrogram struct with |x|^2 of 
	 the complex numbers in spectrogram.
 */
void         spectrogram_compute_powerspectrum( Spectrogram *s ){
  int i,j;
  for( i=0; i<s->N_time; i++ ){
	 for( j=0; j<s->N_freq; j++ ){
		s->powerspect[i][j] = SQR( sqrt( SQR(s->real[i][j])+SQR(s->imag[i][j]) ) );
	 }
  }
}
