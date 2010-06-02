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
 \brief \ref status_unstable Functions using nonlinear systems-theory.

 In general, the phase space is represented through a 
 two-dimensional DOUBLE-array (a matrix). This can e.g. be all
 electrodes over time or any tranformation of the raw data.

 \section tdelay Time Delay Reconstruction

 Time-Delay reconstruction is done as in
 \f[
 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
 \f]
 where \f$ \tau = \nu \Delta t\f$ is a multiple of the sampling step.

 There is a special structure for time-delay reconstruction to avoid the
 space-overhead. This structure saves only the 1D original time-series
 and calculates the other dimensions on the fly.

 The disadvantage of this approach is that there need to be special functions
 to deal with these structs. Some of the functions are implemented both for
 general phase-spaces (full matrices) and the TimeDelayReconstruction struct. 
 Others are not. You can always convert from a TimeDelayReconstruction-struct to
 a Matrix using tdelay_to_array()
 \code
 double *data=getdata(100);
 TimeDelayReconstruction *p=tdelay_init( 4, 10, data, 100 );
 Array *space=tdelay_to_array( p );
 \endcode

 Functions that are designed to be used with the TimeDelayReconstruction begin 
 with the prefix tdelay_*.

 Handling of the struct is done via 
 \code
  TimeDelayReconstruction* tdelay_init ( int m, int tau, double *x, int n );
  void        tdelay_free ( TimeDelayReconstruction *p );
  void        tdelay_print( FILE *out, TimeDelayReconstruction *p);
 \endcode
 and it can be indexed with
 \code
  double      tdelay_index_ij( TimeDelayReconstruction *p, int i, int j );
  void        tdelay_index_i ( TimeDelayReconstruction *p, int i, double *x);
  void        tdelay_index_j ( TimeDelayReconstruction *p, int j, double *x);
 \endcode

 */
#ifndef NONLINEAR_H
# define NONLINEAR_H
#include "mathadd.h"
#include "distances.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*-----------------------------------------------------------
	 - Time-Delay Embedding
	 ---------------------------------------------------------*/

  typedef struct {
	 int m;     /**< embedding dimension */
	 int tau;   /**< time-lag for embedding (>0, multiple of sampling_step) */
	 double *x; /**< pointer to data */
	 int xn;    /**< len(x) */
  } TimeDelayReconstruction;

  Array* tdelay_to_array( TimeDelayReconstruction *p );

  /* parameter estimation */
  int         tdelay_estimate_timelag_mutual( TimeDelayReconstruction *p, 
															 long partitions, long corrlength, double *mutual );
  int         tdelay_estimate_timelag_autocorr( TimeDelayReconstruction *p );
  int         tdelay_estimate_dimension( double Rtol, int m_start, int m_end );
  double      tdelay_fnn_ratio( TimeDelayReconstruction *p, double Rtol, double Atol );
  double      tdelay_attractor_size( TimeDelayReconstruction *p );

  /* struct handling */
  TimeDelayReconstruction* tdelay_init ( int m, int tau, double *x, int n );
  void        tdelay_free ( TimeDelayReconstruction *p );
  void        tdelay_print( FILE *out, TimeDelayReconstruction *p);

  /* indexing */
  double      tdelay_index_ij( TimeDelayReconstruction *p, int i, int j );
  void        tdelay_index_i ( TimeDelayReconstruction *p, int i, double *x);
  void        tdelay_index_j ( TimeDelayReconstruction *p, int j, double *x);

  /* nonlinear prediction */
  double tdelay_predict_simple( TimeDelayReconstruction *p, double *sample, int npredict, double epsilon );
  double tdelay_simple_nonlinear_prediction_error( TimeDelayReconstruction *reference, double *y, int yn, 
																	 int npredict, double epsilon );

  double** eeg_nonlinear_prediction_error( const EEG *eeg, int embedding_dim, int time_lag,
														 int npredict, double epsilon, 
														 double** output, OptArgList *optargs );
#ifdef __cplusplus
}
#endif


#endif
