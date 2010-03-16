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
/**\file time_frequency.h
	\brief \ref status_stable Time-Frequency represenation of time-series.
*/

#ifndef TIME_FREQUENCY_H
# define TIME_FREQUENCY_H

#include "definitions.h"
#include "mathadd.h"
#include "optarg.h"
#include <math.h>
/* #include <gsl/gsl_fft_complex.h> */

#ifdef __cplusplus
extern "C" {
#endif

  /** \addtogroup spectgram 
		\{ */
  typedef struct _Spectrogram  {
    int            N_freq;		 /**< number of freq bins in the TFR matrix */
    int            N_time;	/**< number of time_bins in the TFR matrix */
	 double         samplingrate; /**< sampling rate with which the signal was sampled */
	 double         low_corner_freq; /**< lower corner frequency */
	 double         up_corner_freq;	/**< upper corner frequency */
	 Complex        **sgram; /**< the spectrogram data */
	 double        **powerspect;
	 int           has_power_spectrum; /**< true if power spectrum is filled */
  } Spectrogram;
  
  /** \addtogroup spectspect
		\{ */
  Spectrogram* spectrogram_stft(const double* s, int n, double sampling_rate,
										  const double *Window, int Window_Length,
										  int N_freq, int N_time, double corner_freqs[2],
										  int *timepoints, Spectrogram *spectgram);

  Spectrogram* spectrogram( const double *s, int n, Spectrogram *spectgram, 
									 OptArgList *optargs );
  /** \} */

  /** \ingroup specthelp
		\{ */
  Spectrogram* spectrogram_init( int N_freq, int N_time );
  void         spectrogram_free( Spectrogram *s );
  void         spectrogram_compute_powerspectrum( Spectrogram *s );
  /** \} */

  /** \addtogroup windows
	* \{ */
  double* window_dirichlet( double *window, int n );
  double* window_gaussian ( double *window, int n, double sigma );
  double* window_hamming  ( double *window, int n );
  double* window_hanning  ( double *window, int n );
  double* window_kaiser   ( double *window, int n, double alpha );
  /** \} */

  /** \} */
#ifdef __cplusplus
}
#endif



#endif /* TIME_FREQUENCY_H */
