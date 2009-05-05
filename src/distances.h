/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke                                  *
 *   mihrke@uni-goettingen.de                                              *
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

/**\file distances.h
 * Distances between: points, signals, trials.
 */
#ifndef DISTANCES_H
# define DISTANCES_H
#include "clustering.h"
#include "definitions.h"
#include "recurrence_plot.h"
#include "nonlinear.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\addtogroup signaldist
	  These functions (signaldist_*()) calculate pointwise differences
	  between individual points in the signals, returning the matrix
	  \f[
	  d(t_1, t_2) = ||s_1(t_1) - s_2(t2)||
	  \f]
	  for two signals \f$ s_1(t), s_2(t)\f$ and some metric \f$||\cdot ||\f$.
	*\{
	*/  
  double** signaldist_euclidean( const double *s1, int n1, const double *s2, int n2, 
											double **d, void *userdata );
  double** signaldist_euclidean_derivative( const double *s1, int n1, const double *s2, int n2, 
														  double **d, void *userdata );
  double** signaldist_stft( const double *s1, int n1, const double *s2, int n2, 
									 double **d, void *userdata );
  /*\}*/

  
  /** \addtogroup eegdist
		These functions are wrappers for the signaldist_*() functions from group 
		\ref signaldist which take EEGData structs instead of raw signals.
	*\{
	*/
  double** eeg_distmatrix_euclidean_channel( const EEGdata *s1, const EEGdata *s2, 
															int channel, double **d, void* params );
  double** eeg_distmatrix_euclidean_derivative_channel( const EEGdata *s1, const EEGdata *s2, 
																		  int channel, double **d, void *params ); 
  double** eeg_distmatrix_stft_channel( const EEGdata *s1,const  EEGdata *s2, 
													 int channel, double **d, void *params );
  double** eeg_distmatrix_recplot_losdtwnoise_channel( const EEGdata *s1,const  EEGdata *s2, 
																		 int channel, double **d, void *params );  
  /** \} */

  /** \addtogroup vectordist 
		These distances (functions starting with vectordist_*()) are between 
		two vectors of the same size, i.e. they calculate
		\f[
		||\vec{x}_i - \vec{x}_j||
		\f]
		for two vectors and some metric. vectordist_distmatrix() applies one of 
		these functions to all pairs of columns in the matrix \f$\mathbf{X}\f$, 
		yielding a distance matrix
		\f[
		D_{ij} =	||\vec{x}_i - \vec{x}_j||
		\f]
		\{
  */  
  double** vectordist_distmatrix          ( VectorDistanceFunction f, const double **X, 
														  int n, int p, double **D, 
														  ProgressBarFunction progress, void* userdata );
  double   vectordist_euclidean           ( double *x1, double *x2, int p, void *userdata );
  double   vectordist_euclidean_normalized( double *x1, double *x2, int p, void *userdata );
  double   vectordist_dtw                 ( double *x1, double *x2, int p, void *userdata );
  double   vectordist_regularized_dtw     ( double *x1, double *x2, int p, void *userdata );

  double** eegtrials_distmatrix_channel( EEGdata_trials *eeg, VectorDistanceFunction f, 
													  int channel, double **d, void *userdata );
  /** \} */

  /** \addtogroup warppathdist 
		These distances (functions starting with pathdist_*()) are between 
		two WarpPath structs.
		\{
  */  
  double   pathdist_euclidean_dt(WarpPath *p1, WarpPath *p2); 
  /** \} */

  /** \addtogroup pointdist 
		These distances (functions starting with pointdist_*()) are between 
		two numbers.
		\{
  */  
  double   pointdist_euclidean(double x, double y); 
  /** \} */

  /** \addtogroup otherdist
		Other distance functions.
		\{
  */  
  double   dist_point_line(double *p, double *x, double *y); 
  /** \} */
  
#ifdef __cplusplus
}
#endif


#endif
