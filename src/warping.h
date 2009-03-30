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

/**\file warping.h
 * Warping functions.
 
 Generally, to warp EEG-signals with a hierarchical method, you need:
  -# distance matrix Delta between trials (use clusterdist_*() functions)
  -# a pointwise distance matrix d between 2 trials (use eeg_distmatrix_*())
  -# optionally a regularization matrix G (multiply with d)
  -# the warp-path P obtained from d (or d*G)
  -# an averaging scheme to put together s1 and s2 using P.
 */
#ifndef WARPING_H
#define WARPING_H
#include "mathadd.h"
#include "definitions.h"
#include "clustering.h"
#include "time_frequency.h"
#include "regularization.h"
 
#ifdef __cplusplus
extern "C" {
#endif


  /** \addtogroup distance
	*\ingroup warping
	*\{
	*/
  double** signaldist_euclidean( const double *s1, int n1, const double *s2, int n2, 
											double **d, void *userdata );
  double** signaldist_euclidean_derivative( const double *s1, int n1, const double *s2, int n2, 
														  double **d, void *userdata );
  double** signaldist_stft( const double *s1, int n1, const double *s2, int n2, 
									 double **d, void *userdata );

  /* wrapper for eeg */
  double** eeg_distmatrix_euclidean_channel( const EEGdata *s1, const EEGdata *s2, 
															int channel, double **d, void* params );
  double** eeg_distmatrix_euclidean_derivative_channel( const EEGdata *s1, const EEGdata *s2, 
																		  int channel, double **d, void *params ); 
  double** eeg_distmatrix_stft_channel( const EEGdata *s1,const  EEGdata *s2, 
													 int channel, double **d, void *params );
  /** \} */


  /* ---------------------------------------------------------------------------- 
	  -- Timewarping                                                            -- 
	  ---------------------------------------------------------------------------- */
  /** \addtogroup timewarp
	*\ingroup warping
	*\{
	*/


  WarpPath* DTW_path_from_square_distmatrix(const double **d, int n, WarpPath *P);

  double**  DTW_build_cumdistmatrix(const double *u, int J, const double *s, int K, 
												double theta1, double theta2, double **d);
  double**  DTW_build_restricted_cumdistmatrix(const double *u, int J, 
															  const double *s, int K, 
															  double R, double **d);
  
  WarpPath* DTW_get_warppath2(const double *u, int J, const double *s, int K,	
										double theta1, double theta2, double *Djk);
  double    DTW_get_warpdistance(const double *u, int J, const double *s, int K,
											double theta1, double theta2);
  WarpPath* DTW_path_from_cumdistmatrix(const double **d, int J, int K, WarpPath *P);

  double    DTW_distance_between_paths(const WarpPath *P1, const WarpPath *P2);
  /** \} */

  /** \addtogroup warpaveraging
	*\ingroup warping
	*\{
	*/  
  EEGdata*  eeg_ADTW_from_path(const EEGdata *s1, const EEGdata *s2, 
										 EEGdata *target, int channel, const WarpPath *P);  
  double*   ADTW (const double *s1, int n1, const double *s2, int n2, double *avg);
  double*   ADTW_from_path(const double *u, int J, const double *s, int K, 
									const WarpPath *P, double *avg);
  double*   DTW_add_signals_by_path(const double *s1, int n1,
												const double *s2, int n2, const WarpPath *P, double *avg,
												const double weights[2]);
  EEGdata*  eeg_DTW_add_signals_by_path(const EEGdata *s1, const EEGdata *s2, EEGdata *target, 
													 int channel, const WarpPath *P, const double weights[2]);

  EEGdata* eegtrials_PADTW( EEGdata_trials *eeg_in, const double **distmatrix, 
									 int N,  EEGdata *out, SettingsPADTW settings );

  SettingsPADTW init_PADTW();
  /** \} */

  /** \addtogroup otherwarp
	*\ingroup warping
	*\{
	*/  
  EEGdata* eegtrials_gibbons( const EEGdata_trials *eeg_in, EEGdata *target, 
										int stimulus_marker, int response_marker, double k );
  /** \} */

#ifdef __cplusplus
}
#endif

#endif
