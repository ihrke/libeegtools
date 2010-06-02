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
	\brief \ref status_unstable Time-Frequency representation of time-series.

	A time-frequency represenation (TFR) of a signal displays the 
	varying frequency content of a signal over time.

	\image html spectrogram.jpg

	\todo add some docu here

	\note the test-functions for this functionality are not excessive.
*/

#ifndef TIME_FREQUENCY_H
# define TIME_FREQUENCY_H

#include "definitions.h"
#include "complex.h"
#include "array.h"
#include "optarg.h"
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

  /** \brief A Window Function.
		\param n is the number of points in the window
		\param param is a parameter depending on the choice of the
		         window (e.g. sigma for gaussian)
	*/
  typedef Array* (*WindowFunction)(int,double);

  /** \brief Complex spectrogram (Time-Frequency-Representation) of a signal.
	*/
  typedef struct {
    int            N_freq;		 /**< number of freq bins in the TFR matrix */
    int            N_time;	/**< number of time_bins in the TFR matrix */
	 double         samplingrate; /**< sampling rate with which the signal was sampled */
	 double         low_corner_freq; /**< lower corner frequency */
	 double         up_corner_freq;	/**< upper corner frequency */
	 Complex        **sgram; /**< the spectrogram data (N_time x N_freq) */
  } Spectrogram;
  
  Spectrogram* spectrogram_stft(const Array* sig, double sampling_rate,
										  const Array *Window, int N_freq, int N_time, 
										  double corner_freqs[2], 
										  int *timepoints, Spectrogram *spectgram);
  Spectrogram* spectrogram( const Array *sig, Spectrogram *spectgram, 
									 OptArgList *optargs );
  Array* spectrogram_powerspectrum( const Spectrogram *spect );

  Spectrogram* spectrogram_init( int N_freq, int N_time );
  void         spectrogram_free( Spectrogram *s );


  /*------------------- window functions --------------- */
  Array* window_dirichlet( int n, double noparam );
  Array* window_gaussian ( int n, double sigma );
  Array* window_hamming  ( int n, double noparam );
  Array* window_hanning  ( int n, double noparam );
  Array* window_kaiser   ( int n, double alpha );

#ifdef __cplusplus
}
#endif



#endif /* TIME_FREQUENCY_H */
