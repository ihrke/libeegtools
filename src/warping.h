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

	Example for warping 2D data:

	\image html 2dwarping.jpg

	\section dtw Dynamic Time Warping

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

*/


#ifndef WARPING_H
#define WARPING_H
#include "mathadd.h"
#include "definitions.h"
#include "time_frequency.h"
#include "regularization.h"
#include "eeg.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*-----------------------------------------------------------
	 -  slope constraint for DTW
	 ---------------------------------------------------------*/
  /**\brief Severity of the slope constraint in the Dynamic Time-Warping Alg.
	*/
  typedef enum {
	 SLOPE_CONSTRAINT_NONE=0,
	 SLOPE_CONSTRAINT_LAX,
	 SLOPE_CONSTRAINT_MEDIUM,
	 SLOPE_CONSTRAINT_SEVERE
  } SlopeConstraint;

  /*-----------------------------------------------------------
	 - WARPING -
	 ---------------------------------------------------------*/

  typedef struct{
	 int *t1; /**< time-scale of first signal */
	 int *t2; /**< time-scale of second signal */
	 int n1; /**< length of signal 1 */
	 int n2; /**< length of signal 2 */
	 int n; /**< length(t1)==length(t2) */
  } WarpPath;
	

  /** \brief is w a warppath?
		
		Usage:
		\code  
		bool ispath;
		warppath_CHECK( ispath, X );
		if( !ispath ) return NULL;
		\endcode
		\param flag (output) (bool) set by macro
		\param w (input) Array* to check
  */
#define warppath_CHECK( flag, w )													\
  if(!( (w)->ndim==2 && (w)->dtype==UINT )){									\
	 char *dts="";																		\
	 array_DTYPESTRING( dts, w->dtype );										\
	 errprintf("not a warppath, ndim=%i, dtype=%s\n", w->ndim, dts );	\
	 flag=FALSE;																		\
  } else { flag=TRUE; }																
  
  /* ---------------------------------------------------------------------------- 
	  -- Timewarping                                                            -- 
	  ---------------------------------------------------------------------------- */

  Array*  matrix_dtw_cumulate ( Array *mat, bool alloc, OptArgList *optargs );
  Array*  matrix_dtw_backtrack ( const Array *d );

  Array*  dtw_add_signals( const Array *s1, const Array *s2, const Array *path, OptArgList *opts );

  
  EEG* eeg_dtw_hierarchical( EEG *eeg_in, const double **distmatrix, 
									  EEG *out, OptArgList *optargs );

  EEG* eeg_gibbons( EEG *eeg, int stimulus_marker, int response_marker, double k );


  /******************************************************************/
  /** GOING TO BE OBSOLETE */
  /******************************************************************/

  void       dtw_cumulate_matrix  ( double **d, int M, int N, OptArgList *opts );
  /* void dtw_cumulate_matrix_chiba_band( double **d, int M, int N, double restriction );*/
  WarpPath*  dtw_backtrack        ( const double **d, int M, int N, WarpPath *P );

   WarpPath* init_warppath ( WarpPath *path, int n1, int n2 );
   void      free_warppath ( WarpPath *p );
   void      reset_warppath( WarpPath *P, int n1, int n2 );
   void      print_warppath( FILE *out, WarpPath *P );
 
#ifdef __cplusplus
}
#endif

#endif
