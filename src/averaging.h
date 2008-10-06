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

/**\file averaging.h
 * Averaging functions.
 * \defgroup averaging Averaging functions
 *\{
 *    \defgroup timewarp Timewarping-functions
 *    \defgroup otheravg Other Averaging functions
 *\}
 */
#ifndef AVERAGING_H
#define AVERAGING_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif
/* ---------------------------------------------------------------------------- 
   -- Timewarping                                                            -- 
   ---------------------------------------------------------------------------- */
/** \addtogroup timewarp
 *\ingroup averaging
 *\{
 */

  double**  DTW_build_cumdistmatrix(const double *u, int J, const double *s, int K, 
												double theta1, double theta2, double **d);
  double**  DTW_build_restricted_cumdistmatrix(const double *u, int J, 
															  const double *s, int K, 
															  double R, double **d);
  double    DTW_get_warppath(const double *u, int J, const double *s, int K,
									  double theta1, double theta2, int* path);
  WarpPath* DTW_get_warppath2(const double *u, int J, const double *s, int K,	
										double theta1, double theta2, double *Djk);
  double    DTW_get_warpdistance(const double *u, int J, const double *s, int K,
											double theta1, double theta2);
  WarpPath* DTW_path_from_cumdistmatrix(const double **d, int J, int K, WarpPath *P);
  void      DTW_cumulate_distmatrix(double **d, int J, int K);
  double**  DTW_build_distmatrix(const double *u, int J, const double *s, int K, 
											double theta1, double theta2, double **d);
  void      DTW_markers_to_distmatrix(double **d, int J, int K, 
												  const unsigned long **markers, 
												  unsigned nmarkers);
  WarpPath* DTW_warppath_with_markers(const double *u, int J, const double *s, int K,
												  double theta1, double theta2, 
												  const unsigned long **markers, unsigned nmarkers);
  double    DTW_get_warpdistance_markers(const double *u, int J, 
													  const double *s, int K,
													  double theta1, double theta2, 
													  const unsigned long **markers, 
													  unsigned nmarkers);
  double    DTW_distance_between_paths(const WarpPath *P1, const WarpPath *P2);
  WarpPath* eeg_DTW_get_paths_by_markers( const EEGdata *s1, const EEGdata *s2, 
														int channel, double theta, WarpPath *P );

  double*   ADTW (const double *s1, int n1, const double *s2, int n2, double *avg);
  double*   ADTW_signal(const double *s1, int sR1, 
								const double *s2, int sR2, 
								int zero, int n, double *avg);
  double*   ADTW_from_path(const double *u, int J, const double *s, int K, 
									const WarpPath *P, double *avg);
  void      eeg_ADTW_markers_channel(const EEGdata *s1, const EEGdata *s2, 
												 EEGdata *target, int channel, double theta);
  void      eeg_ADTW_markers(const EEGdata *s1, const EEGdata *s2, EEGdata *target, double theta);
  EEGdata*  eeg_ADTW_from_path(const EEGdata *s1, const EEGdata *s2, 
										 EEGdata *target, int channel, const WarpPath *P);
  double*   PADTW(const double **s, int N, int n, int zero, int *sR, double *wa);
/** \} */


/** \addtogroup otheravg
 *\ingroup averaging
 *\{
 */
  double* simple_average_2v(const double *s1, const double *s2, int n, double *avg);
  double* simple_average_nv(const double **s, int N, int n, double *avg);
  double* alternate_average_nv(const double **s, int N, int n, double *avg);

  EEGdata* eegtrials_simple_average( EEGdata_trials *eeg, EEGdata *avg);
  double*  eegtrials_simple_average_channel( EEGdata_trials *eeg, int channel, double *avg );

/** \} */
#ifdef __cplusplus
}
#endif

#endif
