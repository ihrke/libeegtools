/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke                                  *
 *   ihrke@nld.ds.mpg.de
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
 *   aint with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**\file linalg.h
 \brief \ref status_inprogress Linear Algebra.

 \section matrices Vectors and Matrices
 1D and 2D Array structs of DType DOUBLE are used to represent vectors and matrices, respectively. 
 I.e. that each function starting with matrix_*( Array *m ) accepts only 2D Arrays of DType DOUBLE 
 functions starting with vector_*( Array *m ) accept only 1D Arrays of DType DOUBLE.
 This is checked with the Macros matrix_CHECK() and vector_CHECK().
	
 */
#ifndef LINALG_H
# define LINALG_H
#include "definitions.h"
#include "array.h"
#include <gsl/gsl_matrix.h>

/** \brief is m a matrix?

	 Usage:
	 \code  
	 bool ismatrix;
	 matrix_CHECK( ismatrix, X );
	 if( !ismatrix ) return NULL;
	 \endcode
	 \param flag (output) (bool) set by macro
	 \param m (input) Array* to check
*/
#define matrix_CHECK( flag, m )													\
  if(!( (m)->ndim==2 && (m)->dtype==DOUBLE )){								\
	 char *dts="";																		\
	 array_DTYPESTRING( dts, m->dtype );										\
	 errprintf("not a matrix, ndim=%i, dtype=%s\n", m->ndim, dts );	\
	 flag=FALSE;																		\
  } else { flag=TRUE; }	
															
/** \brief is m a square matrix?

	 Usage:
	 \code  
	 bool ismatrix;
	 matrix_CHECKSQR( ismatrix, X );
	 if( !ismatrix ) return NULL;
	 \endcode
	 \param flag (output) (bool) set by macro
	 \param m (input) Array* to check
*/
#define matrix_CHECKSQR( flag, m )													\
  if(!( (m)->ndim==2 && (m)->dtype==DOUBLE )){									\
	 char *dts="";																			\
	 array_DTYPESTRING( dts, (m)->dtype );											\
	 errprintf("not a matrix, ndim=%i, dtype=%s\n", (m)->ndim, dts );		\
	 flag=FALSE;																			\
  } else if( (m)->size[0]!=(m)->size[1] ) {										\
	 errprintf("Matrix is not square, got (%i,%i)\n",(m)->size[0],(m)->size[1]); \
	 flag=FALSE;																			\
  } else { flag=TRUE; }																

/** \brief index matrix at i,j.
	 \param m matrix
	 \param i,j indices
 */
#define mat_IDX( m, i, j )								\
  array_INDEX2( m, double, i, j )

/** \brief is v a vector?

	 Usage:
	 \code
	 bool isvector;
	 vector_CHECK( isvector, X );
	 if( !isvector ) return NULL;
	 \endcode
	 \param flag (output) (bool) set by macro
	 \param v (input) Array* to check
*/
#define vector_CHECK( flag, v )													\
  if(!( (v)->ndim==1 && (v)->dtype==DOUBLE ) ){								\
	 char *dts="";																			\
	 array_DTYPESTRING( dts, v->dtype );										\
	 errprintf("not a vector, ndim=%i, dtype=%s\n", v->ndim, dts );	\
	 flag=FALSE;																		\
  } else { flag=TRUE; }																
	 

#ifdef __cplusplus
extern "C" {
#endif

  /* -------------- FUNCTIONS ---------------- */
  gsl_matrix* matrix_to_gsl( Array *in, bool alloc );

  Array* matrix_get_col( Array *m, int col );
  Array* matrix_get_row( Array *m, int row, bool alloc );
  
  Array* matrix_mult( const Array *m1, const Array *m2 );
  Array* matrix_transpose( Array *m, bool alloc );

  Array* matrix_mean( Array *a, int dim );
  Array* matrix_pca( Array *X, Array **var, bool alloc );


#ifdef __cplusplus
}
#endif

#endif /* LINALG_H */
