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

/** \file mathadd.h
	 \brief \ref status_unstable Additional math-functions.
*/

#ifndef MATH_ADD_H
# define MATH_ADD_H

#ifndef _ISOC99_SOURCE
#define _ISOC99_SOURCE
#endif

#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "helper.h"
#include "definitions.h"

/** \ingroup helpermath*/
#define MAX(a,b) ((a) > (b) ? (a):(b))
/** \ingroup helpermath*/
#define MIN(a,b) ((a) < (b) ? (a):(b))
/** \ingroup helpermath*/
#define SQR(a) ((a)*(a))
/** \ingroup helpermath*/
#define ABS(a) ( ((a)<0) ? (-1*(a)) : (a) )
/** \ingroup helpermath*/
#define ISODD(x)        ((((x)%2)==0)? (0) : (1))

# define PI           3.14159265358979323846  /* pi */

#ifdef __cplusplus
extern "C" {
#endif
  /* -------- Math ----------- */  

  /* ---------------------------------------------------------------------------- 
	  -- Complex Arithmetic
	  ---------------------------------------------------------------------------- */ 
  /** \ingroup complex 
		\{
  */
  typedef struct {
	 double re;
	 double im;
  } Complex;
  
  Complex complex     ( double re, double im );
  Complex complex_add ( Complex a, Complex b );
  Complex complex_add_dbl ( Complex a, double b );  
  Complex complex_sub ( Complex a, Complex b );
  Complex complex_mul ( Complex a, Complex b );
  Complex complex_mul_double( Complex a, double b );
  double  complex_abs ( Complex a );
  Complex complex_exp ( Complex a );
  Complex complex_conj( Complex a );
  Complex complex_neg ( Complex a );
  Complex complex_div ( Complex a, Complex b);
  Complex complex_sqrt( Complex x );
  Complex complex_bilinear_transform(Complex pz);
  /** \} */

  /**\ingroup helpermath
	*\{*/
  double  glog(double v, int b);  
  double  mad(const double *data, int n); 
  int     abscmp(const void *p1, const void *p2);
  double  vnorm(const double *v, int n, int p);
  int     sgn(int x);
  double  maxel(double *v, int n);
  int     maxeli(int *v, int n);
  int     closest_index(const double *v, int n, double c);
  double* sampled_line(double *ntimes, int n, double start, double end);
  double* lininterp(const double *x1, const double *y1, int n1, 
						  const double *x2,       double *y2, int n2);
  int*    linspace(int first, int last);
  double* linspace_dbl(double first, double last, double step, double *v, int *n);
  double  weighted_median_from_unsorted(const double *d, const double *w, int n);
  int cmpdouble(double d1, double d2, int precision);

  int     bresenham_howmany_points( int xstart,int ystart,int xend,int yend );
  int*    bresenham(int xstart,int ystart,int xend,int yend, int *points);

  void    swap2i(int *v1, int *v2);
  void    swap2d(double *v1, double *v2);

  double* flip_array( double *v, int n );

  double  gaussian( double x, double sigma, double mu );

  double** disttransform_deadreckoning(int **I, int X, int Y, double **d);

  int     next_pow2( int n );
  int     iremainder( double x, double y);
  Complex* expand_polynomial_from_roots( const Complex *roots, int n, Complex *coeffs );

  /* ---------------------------------------------------------------------------- 
	-- Fourier methods
	---------------------------------------------------------------------------- */
  void fft(double *data,  unsigned long nn, int isign);
  
  /**\}*/
  
/* ---------------------------------------------------------------------------- 
   -- Signal extension routines                                              -- 
   ---------------------------------------------------------------------------- */
  /** \ingroup sigext
	\{
  */
 
  double* sigext_zeros(double *data, int ns, int n);
  double* sigext_zerosr(double *data, int ns, int n);
  double* sigext_sym(double *data, int ns, int n);
  double* sigext_smooth(double *data, int ns, int n);
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Interpolation
	  ---------------------------------------------------------------------------- */ 
  /** \ingroup interp
	\{
  */
  double  drawsample_nearest_neighbour( const double *v, int n, double x );
  double* resample_linear( const double *s, int n, int newn, double *news );
  double* resample_nearest_neighbour( const double *s, int n, int newn, double *news );
  double* resample_gsl( const double *s, int n, int newn, double *news, const gsl_interp_type *method );
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Merit Measures                                                         -- 
	  ---------------------------------------------------------------------------- */
  double rmse(const double *r, const double *d, int n);
  double snr (const double *r, const double *d, int n);

  /* ---------------------------------------------------------------------------- 
	  -- dblp ops                                                             -- 
	  ---------------------------------------------------------------------------- */  
  /** \ingroup dblpops
		\{
  */
  double  dblp_min( double *v, int n, int *idx );
  double  dblp_max( double *v, int n, int *idx );
  double* dblp_init( double *v, int n,  double val );
  void    dblp_minus_scalar( double *v, int n, double val );
  double* dblp_complex_to_real( const Complex *vc, double *vr, int n );
  double  dblp_euclidean_distance( const double *v1, const double *v2, int n );
  void    dblp_shuffle_int( int *permut, int n );
  void    dblp_print( double *v, int n );
  void    dblp_print_int( int *v, int n );
  double  dblp_mean( double *v, int n );
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Dblpp ops                                                             -- 
	  ---------------------------------------------------------------------------- */ 
  /** \ingroup dblppops
		\{
  */
  double** dblpp_delrow(double **m, int N, int n, int row);
  double** dblpp_delcol(double **m, int N, int n, int col);
  double   dblpp_min(const double **m, int N, int n, int *i1, int *i2);
  double   dblpp_max(const double **m, int N, int n, int *i1, int *i2);
  void     dblpp_print(const double **m, int N, int n);
  double** dblpp_init(int N, int M);
  int**    dblpp_init_int(int N, int M);
  void     dblpp_divide_scalar(double **m, int N, int n, double s);
  void     dblpp_add_scalar(double **m, int N, int n, double s);
  void     dblpp_mul_scalar(double **m, int N, int n, double s);
  void     dblpp_normalize_by_max( double **m, int M, int N );
  void     dblpp_add_dblpp(double **m1, const double **m2, int N, int n);
  void     dblpp_sub_dblpp(double **dest, const double **src, int N, int n);
  void     dblpp_dottimes_dblpp( double **m1, const double **m2, int N, int M );
  void     dblpp_copy( const double **src, double **dest, int N, int M );
  void     scalar_minus_dblpp( double scalar, double **m, int N, int M );
  double** dblpp_rand( double **m, int N, int M, double lower, double upper );
  void     dblpp_free(double **m, int N);
  /**\}*/

#ifdef __cplusplus
}
#endif

#endif
