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

#ifndef TIME_FREQUENCY_H
# define TIME_FREQUENCY_H

#include "definitions.h"
#include "mathadd.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  /** \defgroup spectgram Spectrogram (Time-Frequency Representations)
		\{ */
  typedef struct _Spectrogram  {
    int            N_freq;		 /** number of freq bins in the TFR matrix */
    int            N_time;	/** number of time_bins in the TFR matrix */
	 double         samplingrate; /** sampling rate with which the signal was sampled */
    double        **real;	/** real part of the TFR, N_time x N_freq */
    double        **imag;	/** imaginary part of the TFR, N_time x N_freq*/
	 double        **powerspect;
	 int           has_power_spectrum; /** true if power spectrum is filled */
  } Spectrogram;
  
  /** \defgroup spectspect Spectrogram Algorithms (TFR)
		\{ */
  Spectrogram* spectrogram_stft(const double* s, int n, double sample_freq,
										  const double *Window, int Window_Length,
										  int N_freq, int N_time,
										  Spectrogram *spectgram);
  /** \} */

  /** \defgroup specthelp Spectrogram Convenience
		\{ */
  Spectrogram* init_spectrogram( int N_freq, int N_time );
  void         free_spectrogram( Spectrogram *s );
  void         spectrogram_compute_powerspectrum( Spectrogram *s );
  /** \} */

  /** \defgroup windows Windowing Functions
	* \{ */
  double* window_dirichlet( double *window, int n );
  double* window_gaussian ( double *window, int n, double sigma );
  double* window_hamming  ( double *window, int n );
  /** \} */

  /** \} */
#ifdef __cplusplus
}
#endif



#endif /* TIME_FREQUENCY_H */