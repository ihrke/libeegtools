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

/**\file wavelet.h
 * \brief \ref status_stable Contains functions used for wavelet-based signal denoising.

\section secwavelet Wavelet-based signal denoising.

\subsection secwavtheory Wavelet Theory

Similar to the Fourier transform, which decomposes a function into
bases of sines and cosines, the wavelet transform can be used to
obtain a different representation of a function. Unlike the Fourier
transform however, the wavelet transform is not limited to specific
functions, but rather decomposes a function into bases of scaled and
shifted functions of a mother function, commonly referred to as the
mother wavelet. Due to the fact that the used mother wavelets are well
localized in frequency and time (unlike the sines and cosines used in
the Fourier transform), a much more compact and sparse representation
of the function can be achieved using wavelets.

Technically, the wavelets are scaled and shifted versions of a
function \f$\psi\f$
\f[
\psi_{a,b}(t) = \frac{1}{\sqrt{a}}\psi\left(\frac{t-b}{a}\right)
\f]
and the continuous wavelet transform (CWT) is characterised by
\f[
  \mathcal{W} s(a,b) = \int_{-\infty}^{+\infty}s(t)\psi_{a,b}(t)dt
\f]
where s is the signal to be decomposed. Restrictions put on the
mother wavelet involve its square and absolute integrability as well
as the requirement that the function integrates to 0, thus ensuring
a tight localisation in time and frequency space.

In order to overcome the redundancy present in this reconstruction of
the function, as well as to make the transform applicable to practical
problems, the discrete wavelet transform (DWT) has been developed:
\f[
  \mathcal{W} s(j,k) = \frac{1}{\sqrt{2^j}}\int_{-\infty}^{+\infty}s(t)\psi\left(\frac{t-k}{2^j}\right)dt
\f]
where j and k are integers.

It has been shown [Daubechies1992] that a function can be
reproduced accurately, if and only if the energy of the wavelet
coefficients lies between two positive bounds. It can be shown, that
the discrete wavelet transform for a discrete signal is equivalent to
passing the signal through a filter bank, filtering and
downsampling the signal in successive steps to obtain wavelet
coefficients at different resolution levels [Jansen2001].
This procedure can be implemented in an efficient way (O(n)) using
a pyramidal algorithm, the fast discrete wavelet transform (FDWT) as
described in [Press1992].

\subsection secwavfilter Wavelet-based Denoising

test test test

*/

#ifndef WAVELET_H
# define WAVELET_H

#include "definitions.h"
#include <gsl/gsl_wavelet.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "eeg.h"

#ifdef __cplusplus
extern "C" {
#endif
  /*-----------------------------------------------------------
	 - Function Pointers
	 ---------------------------------------------------------*/
  /**\brief function for threshold selection.

	  \ingroup grpwavelet
  */
  typedef double  (*ThresholdSelectionFunction)   (const double*,int);
  /**\brief function for thresholding.

	  \ingroup grpwavelet
  */
  typedef double  (*ThresholdFunction)            (double,double);
  /**\brief function for extending the signal

	  \ingroup grpwavelet
  */
  typedef double* (*SignalExtensionFunction)      (double*,int,int);

  /**\brief parameters for wavelet-filtering.

	  \ingroup grpwavelet
  */
  typedef struct {
	 int first_thresholding_level; /**< level to start thresholding on */
	 gsl_wavelet_type *wavelet;         /**< GSL-specification of wavelet */
	 int vanishing_moments;        /**< number of vanishing moments of the wavelet */
	 ThresholdSelectionFunction threshselfct; /**< function to select the threshold */
	 ThresholdFunction threshfct;             /**< hard/soft thresholding */
	 SignalExtensionFunction sigextfct;       /**< how to extend the signal to 2^K */
  } WaveletParameters;

  /* ---------------------------------------------------------------------------- 
	  -- Wavelet-based Denoising routines                                       -- 
	  ---------------------------------------------------------------------------- */
  double soft_thresholding( double d, double lambda );
  double hard_thresholding( double d, double lambda );

  double translation_invariant_thresholding( const double *data, int n );
  double conventional_thresholding         ( const double *data, int n );
  double sureshrink                        ( const double *data, int n );
  double heuristic_sure                    ( const double *data, int n );

  WaveletParameters wavelet_init();
  int  wavelet_denoise           ( double *data, int n, WaveletParameters P );
  int  wavelet_extend_and_denoise( double *data, int n, WaveletParameters P );
  EEG* eeg_wavelet_denoise       ( EEG *eeg, WaveletParameters P, bool alloc );

#ifdef __cplusplus
}
#endif
#endif /* - WAVELET_H - */
