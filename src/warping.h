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
	\brief \ref status_stable Warping functions.

 */
#ifndef WARPING_H
#define WARPING_H
#include "mathadd.h"
#include "definitions.h"
#include "clustering.h"
#include "time_frequency.h"
#include "regularization.h"
#include "distances.h"
 
#ifdef __cplusplus
extern "C" {
#endif

  /* ---------------------------------------------------------------------------- 
	  -- Timewarping                                                            -- 
	  ---------------------------------------------------------------------------- */
  /** \addtogroup dtw
		Dynamic Time-Warping is a method to account for temporal distortion when 
		comparing two signals.
		It is done by finding a function that minimizes the sum 
		of entries through a distance matrix, such that
		\f[ 
		\mbox{argmin}_\phi \int d(s_1(t_1),s_2(\phi(t_2)))
		\f]
		for a given pointwise distance between two signals. Such distances can be 
		computed using the functions in distances.h.

		The minimization is done by cumulating the matrix
		\f[
		D_{jk} = \mathbf{d}_{jk}+min{\{D_{j,k-1}, D_{j-1,k}, D_{j-1, k-1}\}}
		\f]
		and backtracking via the minimum of the three neighboring entries (down,
		down-right, right) from \f$D_{J,K}\f$ to \f$D_{1,1}\f$.
		Here, the functions dtw_cumulate_matrix() and dtw_backtrack() do that.

		Finally, the signals need to be mapped to one another to get time-amplitude 
		averaging with 
		\f[
		s'(t)=\frac{\omega_1 s_1(p_1^t) + \omega_2 s_2(p_2^t)}{\omega_1 +
		\omega_2}
		\f]
		using weights \f$\omega_1,\omega_2\f$. Here, you can use 
		warp_add_signals_by_path().

		\todo add a high-level interface where a dtw-struct is used
		\{
	*/

  void       dtw_regularize_matrix( double **d, const double **R, int M, int N );
  void       dtw_cumulate_matrix  ( double **d, int M, int N, OptArgList *opts );
  /* void dtw_cumulate_matrix_chiba_band( double **d, int M, int N, double restriction );*/
  WarpPath*  dtw_backtrack        ( const double **d, int M, int N, WarpPath *P );
  

  /** \} */
  
  /** \addtogroup warpaveraging
	*\{
	*/  
  double*   warp_add_signals_by_path(const double *s1, int n1,
												 const double *s2, int n2, 
												 const WarpPath *P, double *avg,
												 const double weights[2]);
  int*      warp_adjust_time_markers(const int *m1, const int *m2, 
												 int nmarkers, int *outmarkers,
												 const double weights[2] );
  /** \} */

  /** \addtogroup hierarchical		
		Generally, to warp EEG-signals with a hierarchical method, you need:
		-# distance matrix Delta between trials (use vectordist_distmatrix())
		-# a pointwise distance matrix d between 2 trials (use eeg_distmatrix_*())
		-# optionally a regularization matrix G (multiply with d)
		-# the warp-path P obtained from d (or d*G)
		-# an averaging scheme to put together s1 and s2 using P.
		
		The distance matrix for between-trials must be provided separately.
		\{
  */

  
  EEG* eeg_dtw_hierarchical( EEG *eeg_in, const double **distmatrix, 
									  EEG *out, OptArgList *optargs );

  /** \} */
  

  /** \addtogroup structwarppath
		These are convenience functions for handling the WarpPath struct.
		\{
  */
  WarpPath* init_warppath ( WarpPath *path, int n1, int n2 );
  void      free_warppath ( WarpPath *p );
  void      reset_warppath( WarpPath *P, int n1, int n2 );  
  void      print_warppath( FILE *out, WarpPath *P );
  /** \} */

  /** \addtogroup otherwarp
	*\{
	*/  
  EEG* eeg_gibbons( EEG *eeg, int stimulus_marker, int response_marker, double k );
  /** \} */


#ifdef __cplusplus
}
#endif

#endif
