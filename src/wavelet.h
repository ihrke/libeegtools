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
 *
 */

#ifndef WAVELET_H
# define WAVELET_H

#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* ---------------------------------------------------------------------------- 
	  -- Wavelet-based Denoising routines                                       -- 
	  ---------------------------------------------------------------------------- */
  /** \ingroup thresholding
		\{
  */
  double soft_thresholding( double d, double lambda );
  double hard_thresholding( double d, double lambda );
  /** \} */

  /** \ingroup select_thresh
		\{
	*/
  double translation_invariant_thresholding( const double *data, int n );
  double conventional_thresholding         ( const double *data, int n );
  double sureshrink                        ( const double *data, int n );
  double heuristic_sure                    ( const double *data, int n );
  /** \} */

  /** \ingroup wavelet 
		\{
  */
  WaveletParameters wavelet_init();
  int  wavelet_denoise           ( double *data, int n, WaveletParameters P );
  int  wavelet_extend_and_denoise( double *data, int n, WaveletParameters P );
  EEG* eeg_wavelet_denoise       ( EEG *eeg, WaveletParameters P, bool alloc );
  /**\} */


#ifdef __cplusplus
}
#endif
#endif /* - WAVELET_H - */
