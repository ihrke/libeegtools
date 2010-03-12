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

/**\file nonlinear.h
 \brief Functions using nonlinear systems-theory.

 Phase-Space reconstruction is done as in
 \f[
 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
 \f]
 where \f$ \tau = \nu \Delta t\f$ is a multiple of the sampling step.
	
 */
#ifndef NONLINEAR_H
# define NONLINEAR_H
#include "mathadd.h"
#include "distances.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\ingroup nonlinear
	*\{
	*/ 

  /* parameter estimation */
  int         phspace_estimate_timelag_mutual( PhaseSpace *p, long partitions, long corrlength, double *mutual );
  int         phspace_estimate_timelag_autocorr( PhaseSpace *p );
  int         phspace_estimate_dimension( double Rtol, int m_start, int m_end );
  double      phspace_fnn_ratio( PhaseSpace *p, double Rtol, double Atol );
  double      phspace_attractor_size( PhaseSpace *p );

  /* struct handling */
  PhaseSpace* phspace_init ( int m, int tau, double *x, int n );
  void        phspace_free ( PhaseSpace *p );
  void        phspace_print( FILE *out, PhaseSpace *p);

  /* indexing */
  double      phspace_index_ij( PhaseSpace *p, int i, int j );
  void        phspace_index_i ( PhaseSpace *p, int i, double *x);
  void        phspace_index_j ( PhaseSpace *p, int j, double *x);

  /* nonlinear prediction */
  double phspace_predict_simple( PhaseSpace *p, double *sample, int npredict, double epsilon );
  double phspace_simple_nonlinear_prediction_error( PhaseSpace *reference, double *y, int yn, 
																	 int npredict, double epsilon );
  /*\}*/

  double** eeg_nonlinear_prediction_error( const EEG *eeg, int embedding_dim, int time_lag,
														 int npredict, double epsilon, 
														 double** output, OptArgList *optargs );
#ifdef __cplusplus
}
#endif


#endif
