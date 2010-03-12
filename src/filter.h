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

/**\file filter.h
 * \brief Contains functions used for frequency-domain signal filtering.
 *
 */

#ifndef FILTER_H
# define FILTER_H

#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* ------------------------------ 
	  -- Robust filtering methods  --
	  ------------------------------ */
  /** \ingroup robust_filtering\{	*/ 
  double* running_median         ( double *d, int n, int win);
  double* weighted_running_median( double *d, int n, int win, 
											  PointDistanceFunction dist);

  EEG*    eeg_filter_running_median( EEG *eeg, int win, bool alloc );
  EEG*    eeg_filter_weighted_running_median( EEG *eeg, int win, bool alloc );
  /** \} */

  /* ------------------------------ 
	  -- Other filtering methods  --
	  ------------------------------ */

  /** \ingroup other_filtering\{	*/
  double* moving_average(double *s, int n, int win);
  /** \} */

 /* ------------------------------ 
	  -- Bandpass filter  --
	  ------------------------------ */ 
 /** \ingroup frequence_filter \{ */
  EEG* eeg_filter_fidlib( EEG *eeg, const char *spec, bool alloc );
  /** \} */
  
#ifdef __cplusplus
}
#endif


#endif
