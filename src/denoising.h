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

/**\file denoising.h
 * \brief Contains functions used for signal denoising.
 *
 * Especially wavelet-based denoising. 
 */

#ifndef DENOISING_H
# define DENOISING_H

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
#include "distances.h"


#ifdef __cplusplus
extern "C" {
#endif
  /* ------------------------------ 
	  -- Robust filtering methods  --
	  ------------------------------ */

  /** \ingroup robust_filtering
		\{
	*/
  void    eeg_filter_running_median(EEGdata *s, int win);
  double* running_median           (double *d, int n, int win);

  void    eeg_filter_weighted_running_median(EEGdata *s, int win);
  double* weighted_running_median           (double *d, int n, int win, 
															PointDistanceFunction dist);

  double  weighted_median_from_unsorted(const double *d, const double *w, int n);


  /** \} */

  /* ------------------------------ 
	  -- Other filtering methods  --
	  ------------------------------ */

  /** \ingroup other_filtering
		\{
	*/
  double* moving_average(double *s, int n, int win);
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Wavelet-based Denoising routines                                                     -- 
	  ---------------------------------------------------------------------------- */
  /** \ingroup thresholding
		\{
  */
  double eta_s(double d, double lambda);
  double eta_h(double d, double lambda);
  /** \} */

  /** \ingroup wavelet 
		\{
  */
  int generic_denoising   ( double *data, int n, int L, 
									 ThresholdSelectionFunction threshfct,
									 ThresholdFunction etafct );
  int extend_and_denoise  ( double *data, int n, int L, 
									 ThresholdSelectionFunction threshfct,
									 ThresholdFunction etafct, 
									 SignalExtensionFunction sigextfct );
  void eeg_wavelet_denoise( EEGdata *eeg, int L, 
									 ThresholdSelectionFunction threshfct,
									 ThresholdFunction etafct, 
									 SignalExtensionFunction sigextfct );
  /**\} */

  /** \ingroup select_thresh
		\{
	*/
  double translation_invariant_thresholding( const double *data, int n );
  double conventional_thresholding         ( const double *data, int n );
  double sureshrink                        ( const double *data, int n );
  double heuristic_sure                    ( const double *data, int n );
  /** \} */

 /* ------------------------------ 
	  -- Bandpass filter design  --
	  ------------------------------ */
  /** \ingroup filter_design
	*\{
	*/
  void butterworth_design_bandpass( int order, double sampling_rate, double passband[2],
												double *bcoeff, double *acoeff );
  /** \} */


#ifdef __cplusplus
}
#endif
#endif /* - DENOISING_H - */
