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
 * \brief \ref status_unstable Distances between: points, signals, trials.
 
 \section vectdist Vector Distance Functions
 The distances (functions starting with vectordist_*()) are between 
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
*/

#ifndef DISTANCES_H
#define DISTANCES_H

#include "array.h"
#include "linalg.h"
#include "clustering.h"
#include "definitions.h"
#include "recurrence_plot.h"
#include "nonlinear.h"
#include "warping.h"
#include "optarg.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  typedef double**(*PointwiseDistanceFunction)    (double*,int,double*,int,double**,OptArgList*);
  typedef double  (*VectorDistanceFunction)       (const double*,const double*,int,OptArgList*);
 
  Array* distmatrix_signaldist( VectorDistanceFunction f, const Array *s1,
										  const Array *s2, Array *out, OptArgList *optargs );
  Array* matrix_distmatrix( VectorDistanceFunction f, 
									 const Array *X, Array *D, 
									 OptArgList *optargs );

  double** vectordist_distmatrix          ( VectorDistanceFunction f, const double **X, 
														  int n, int p, double **D, 
														  ProgressBarFunction progress, 
														  OptArgList* optargs );

  double   vectordist_euclidean           ( const double *x1, const double *x2, 
														  int p, OptArgList *optargs );
  double   vectordist_euclidean_normalized( const double *x1, const double *x2, 
														  int p, OptArgList *optargs );
  double   vectordist_dtw                 ( const double *x1, const double *x2, 
														  int p, OptArgList *optargs );
  double   vectordist_regularized_dtw     ( const double *x1, const double *x2, 
														  int p, OptArgList *optargs );

  double** eeg_distmatrix( EEG *eeg, VectorDistanceFunction f, 
									double **d, OptArgList *optargs );
 
  /*	These distances (functions starting with pathdist_*()) are between 
		two WarpPath structs.
  */  
  double   pathdist_euclidean_dt(WarpPath *p1, WarpPath *p2); 

  /* Other distance functions.
  */  
  double   dist_point_line(double *p, double *x, double *y); 
  

  /***********************************************************
    GOING TO BE OBSOLETE 
  ***********************************************************/
 
  double** signaldist_euclidean( double *s1, int n1, double *s2, int n2, 
											double **d, 
											OptArgList *optargs );
  double** signaldist_euclidean_derivative( double *s1, int n1, double *s2, int n2, 
														  double **d, 
														  OptArgList *optargs );
  double** signaldist_stft( double *s1, int n1, double *s2, int n2, 
									 double **d, 
									 OptArgList *optargs );

 
#ifdef __cplusplus
}
#endif


#endif
