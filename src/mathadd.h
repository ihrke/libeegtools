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

#ifndef MATH_ADD_H
# define MATH_ADD_H
#define _ISOC99_SOURCE
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>
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
  /**\addtogroup helpermath
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
  void    bresenham(int xstart,int ystart,int xend,int yend, int *points);

  void    swap2i(int *v1, int *v2);
  void    swap2d(double *v1, double *v2);

  double* flip_array( double *v, int n );

  int     next_pow2( int n );
  int     iremainder( double x, double y);
/* ---------------------------------------------------------------------------- 
	-- Fourier methods
	---------------------------------------------------------------------------- */
  void fft(double *data,  unsigned long nn, int isign);
  

/* ---------------------------------------------------------------------------- 
   -- Signal extension routines                                              -- 
   ---------------------------------------------------------------------------- */
  /** \defgroup sigext Signal Extension Functions 
	\{
  */
  /**
	* Signal extension schemes to extend signal of length n to length 2^j
	*    with 2^j being the closest power of 2 to n.
	* The extension functions return a pointer to the former data[0],
	* because this is where the unextended signal began;\n
	* Example:
	* \code
	* sigext([1 2 3 - - - -]) -> [0 0 1 2 3 0 0] 
	*                                 ^ ptr
	* \endcode
	*  Assumptions (not for full generality!):
	* -# ns <= n
	* -# n <= 2*ns
	*
	*/
  double* sigext_zeros(double *data, int ns, int n);
  double* sigext_zerosr(double *data, int ns, int n);
  double* sigext_sym(double *data, int ns, int n);
  double* sigext_smooth(double *data, int ns, int n);
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Interpolation
	  ---------------------------------------------------------------------------- */
  double  drawsample_linear( const double *v, int n, double x );
  double* resample_linear( const double *s, int n, int newn, double *news );
  double  drawsample_nearest_neighbour( const double *v, int n, double x );
  double* resample_nearest_neighbour( const double *s, int n, int newn, double *news );

  /* ---------------------------------------------------------------------------- 
	  -- Merit Measures                                                         -- 
	  ---------------------------------------------------------------------------- */
  double rmse(const double *r, const double *d, int n);
  double snr (const double *r, const double *d, int n);

  /* ---------------------------------------------------------------------------- 
	  -- vector ops                                                             -- 
	  ---------------------------------------------------------------------------- */  
  /** \defgroup vectorops Vector Operations
		\{
  */
  double* vector_init( double *v, int n,  double val );
  void    vector_minus_scalar( double *v, int n, double val );
  /** \} */

  /* ---------------------------------------------------------------------------- 
	  -- Matrix ops                                                             -- 
	  ---------------------------------------------------------------------------- */
  double** matrix_delrow(double **m, int N, int n, int row);
  double** matrix_delcol(double **m, int N, int n, int col);
  double   matrix_min(const double **m, int N, int n, int *i1, int *i2);
  void     matrix_print(double **m, int N, int n);
  double** matrix_init(int N, int M);
  void     matrix_divide_scalar(double **m, int N, int n, double s);
  void     matrix_add_matrix(double **m1, const double **m2, int N, int n);
  void     matrix_free(double **m, int N);
  /**\}*/

#ifdef __cplusplus
}
#endif

#endif
