/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke   *
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
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "linalg.h"
#include "helper.h"
#include <string.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

/** \brief convert 2D- Array struct to gsl_matrix.

	 This function works also for non-double Array's,
	 but the input array MUST be of DType DOUBLE if alloc
	 is FALSE. If alloc is TRUE, the type is casted
	 to double during the copy-process.
	 
	 \note the overhead for alloc=FALSE is very small,
	    for alloc=TRUE MxN copy processes are necessary.
	 \warning use alloc=FALSE only if you have a double-Array
	 \note  if you free with gsl_matrix_free, the memory is 
	       *NOT* free'd but still belongs to the Array!
			 
	 Example:
	 \code
	 gsl_matrix *g;
	 Array *a = array_new_dummy( DOUBLE, 2, 10, 11 );

	 g = matrix_to_gsl( a, FALSES );
	 gsl_matrix_free( g ); //does not free a->data (!)

	 array_free( a ); // does free a->data
	 \endcode

	 \param in the input 2D-array
	 \param alloc if true, copy the memory pointed to by in;
	            if false, simply put the pointer into the 
					gsl_matrix struct (bad stuff happens if
					it is not DType DOUBLE);
	 \return the gsl_matrix
 */
gsl_matrix* matrix_to_gsl( Array *in, bool alloc ){
  gsl_matrix *out;
  if( in->ndim != 2 ){
	 errprintf("in is not a matrix, ndim=%i\n", in->ndim );
  } 
  if( alloc ){
	 out = gsl_matrix_alloc( (size_t)in->size[0], (size_t)in->size[1] );
	 if( in->dtype==DOUBLE ){ /* faster */
		memcpy( out->data, in->data, in->nbytes );
	 } else { /* slow */
		int i,j;
		double el;
		for( i=0; i<in->size[0]; i++ ){
		  for( j=0; j<in->size[1]; j++ ){
			 array_dtype_to_double( &el, array_index2( in, i, j ), in->dtype );
			 gsl_matrix_set( out, i, j, el );
		  }
		}
	 }
  } else { /* no mallocing */
	 MALLOC( out, 1, gsl_matrix );
	 out->size1=(size_t)in->size[0];
	 out->size2=(size_t)in->size[1];
	 out->tda=(size_t)in->size[1];
	 out->data=(double*)in->data;
	 out->owner=0; /* GSL is not the owner of the data */
	 out->block=NULL;
  }
  
  return out;
}

/** \brief calculate the mean of matrix a along dimension dim.
	 \param a the matrix
	 \param dim along which dimension?
	 \return vector of size a->size[1-dim]
 */
Array* matrix_mean( Array *a, int dim ){
  Array *v;
  int i,j;
  bool ismatrix;
  matrix_CHECK( ismatrix, a );
  if( !ismatrix ) return NULL;
  if( dim!=0 && dim!=1 ){
	 errprintf("wrong dimension dim=%i, must be in {0,1}\n", dim );
	 return NULL;
  }
  double *p;
  v = array_new2( DOUBLE, 1, a->size[1-dim] );
  for( i=0; i<a->size[1-dim]; i++ ){
	 p=&array_INDEX1( v, double, i );
	 *p= 0.0;
	 for( j=0; j<a->size[dim]; j++ ){
		array_INDEX1( v, double, i ) += 
		  array_INDEX2( a, double, (dim==1)?i:j,(dim==1)?j:i );
	 }
	 array_INDEX1( v, double, i ) /= (double)a->size[dim];
  }
  return v;
}

/** \brief get a row of a matrix as vector.

	 \note if alloc is TRUE, the returned vector is 
	 independant. If alloc is FALSE, the memory 
	 is shared between m and the row vector but
	 is owned by the matrix (an array_free(vec) call 
	 does not free the memory, but array_free(m) does).

	 \param m the matrix
	 \param row index to the row
	 \param alloc allocate memory for the row-data, or 
	           only pass the pointer to the new Array
 */
Array* matrix_get_row( Array *m, int row, bool alloc ){
  Array *rowvec;
  bool ismatrix;
  matrix_CHECK( ismatrix, m );
  if( !ismatrix ) return NULL;  
  if( row<0 || row>=m->size[0] ){
	 errprintf("Index out of bounds: %i/%i\n", row, m->size[0] );
	 return NULL;
  }

  if( !alloc ){
	 rowvec = array_fromptr2( DOUBLE, 1, array_INDEXMEM2( m, row, 0 ), 
									  m->size[1] );
  } else {
	 rowvec = array_new2( DOUBLE, 1, m->size[1] );
	 memcpy( rowvec->data, array_INDEXMEM2( m, row, 0 ), 
				m->size[1]*m->dtype_size );
  }
  return rowvec;
}

/** \brief get a column of a matrix as vector.

	 \note memory copy is necessary, because the 
	       memory of a column is not aligned.
			 It is therefore slower than calling
			 matrix_get_row() with alloc=FALSE.
			 If you need a lot of columns, you might
			 want consider calling matrix_transpose() first
			 and matrix_get_row() afterwards.

	 \param m the matrix
	 \param col index to the column
	 \return a column vector
 */
Array* matrix_get_col( Array *m, int col ){
  Array *colvec;
  int i;
  bool ismatrix;
  matrix_CHECK( ismatrix, m );
  if( !ismatrix ) return NULL;
  if( col<0 || col>=m->size[1] ){
	 errprintf("Index out of bounds: %i/%i\n", col, m->size[1] );
	 return NULL;
  }
  colvec = array_new2( DOUBLE, 1, m->size[0] );
  for( i=0; i<m->size[0]; i++ ){
	 array_INDEX1( colvec, double, i )=array_INDEX2(m, double, i, col );
  }
  return colvec;
}

/** \brief Matrix Multiplication.

	 m1 and m2 must be matrices such that m1->size[1]==m2->size[0].
	 \param m1,m2 the matrices
	 \return new matrix of size m1->size[1], m2->size[0]
 */
Array* matrix_mult( const Array *m1, const Array *m2 ){
  Array *out=NULL;
  int i,j,k;
  bool ismatrix;
  matrix_CHECK( ismatrix, m1 ); if( !ismatrix ) return NULL;
  matrix_CHECK( ismatrix, m2 ); if( !ismatrix ) return NULL;
  if( m1->size[1] != m2->size[0] ){
	 errprintf("Mismatch in size, cannot calculate matrix product: %i!=%i\n",
				  m1->size[1], m2->size[0] );
	 return NULL;
  }

  out = array_new2( DOUBLE, 2, m1->size[0], m2->size[1] );

  double *field;
  for( i=0; i<out->size[0]; i++ ){
	 for( j=0; j<out->size[1]; j++ ){
		field=(double*)array_INDEXMEM2( out, i, j );
		*field=0.0;
		for( k=0; k<m1->size[1]; k++ ){
		  *field += array_INDEX2( m1, double, i, k )*
			 array_INDEX2( m2, double, k, j );
		}
	 }
  }
  
  return out;
}
/** \brief Matrix transpose.
	 
	 \todo implement a faster way

	 \param m matrix
	 \param alloc if TRUE, return a new matrix that is the transpose of m;
	       else do in-place transpose
 */
Array* matrix_transpose( Array *m, bool alloc ){
  int i,j;
  bool ismatrix;
  matrix_CHECK( ismatrix, m );
  if( !ismatrix ) return NULL;

  Array *out;
  out = array_new2( DOUBLE, 2, m->size[1], m->size[0] );
  for( i=0; i<m->size[0]; i++ ){
	 for( j=0; j<m->size[1]; j++ ){
		array_INDEX2( out, double, j,i) = array_INDEX2( m, double, i,j);
	 }
  }
  if( !alloc ){
	 memcpy( m->data, out->data, out->nbytes );
	 array_free( out );
	 out = m;
	 double tmp=out->size[0];
	 out->size[0]=out->size[1];
	 out->size[1]=tmp;
  }
  
  return out;
}

/** \brief Principal Components analysis.

	 This implementation uses SVD to calculate the PCA.
	 Example:
	 \code
	 Array *X = get_data_from_somewhere();
	 Array *var, *X_pca;
	 matrix_pca( X, &var, FALSE ); // overwrite X
	 X_pca = matrix_pca( X, &var, TRUE ); // do not touch X
	 \endcode
	 \param X a 2D array observations x variables containing the data
	 \param var output: vector of eigenvalues of XX^T in decreasing order; if you
 	         pass NULL, it is ignored; 
	 \param alloc if true, new memory is allocated and returned. Else
	            X is overwritten.
	 \return pointer to the PCA'd matrix 
*/
Array* matrix_pca( Array *X, Array **var, bool alloc ){
  Array *out=NULL;
  int i,j;
  bool ismatrix;
  matrix_CHECK( ismatrix, X );
  if( !ismatrix ) return NULL;

  int N,K; /* N observations, K variables */
  N = X->size[0];
  K = X->size[1];

  Array *tmp = array_copy( X, TRUE );

  /* subtract mean from observations */
  Array *mean=matrix_mean( X, 0 ); 
  for( i=0; i<N; i++ ){ 
	 for( j=0; j<K; j++ ){
		array_INDEX2( tmp, double, i, j ) -=
		  array_INDEX1( mean, double, j );
	 }
  }

  array_scale( tmp, 1.0/sqrt((double) N-1 ) );

  gsl_matrix *A=matrix_to_gsl( tmp, TRUE ); /* copy */
  gsl_matrix *V=gsl_matrix_alloc( K, K );
  gsl_vector *S=gsl_vector_alloc( K );
  gsl_vector *workspace=gsl_vector_alloc( K );

  /* A->U, V->V, S->S */
  gsl_linalg_SV_decomp( A, V, S, workspace);
  gsl_matrix_transpose( V );

  if( var ){
	 (*var)=array_fromptr2( DOUBLE, 1, S->data, S->size );
	 S->owner=0; /* transfer ownership to array */
	 (*var)->free_data=1;
	 for( i=0; i<array_NUMEL( *var ); i++ ){ 
		array_INDEX1( *var, double, i ) = SQR( array_INDEX1( *var, double, i ) );
	 }
  }

  Array *Vp=array_fromptr2( DOUBLE, 2, V->data, V->size1, V->size2 );
  matrix_transpose( tmp, FALSE );
  out=matrix_mult( Vp, tmp ); /* PCA'd data */
  matrix_transpose( out, FALSE );

  if( out->size[0]!=X->size[0] || out->size[1]!=X->size[1] ){
	 errprintf("Input/Output dimension mismatch: (%i,%i) vs. (%i,%i)\n",
				  X->size[0], X->size[1], out->size[0], out->size[1] );
  }

  if( !alloc ){ /* write back out->X */
	 memcpy( X->data, out->data, out->nbytes );
	 array_free( out );
	 out = X;
  }

  /* clean up */
  gsl_matrix_free( A );
  gsl_matrix_free( V );
  gsl_vector_free( S );
  gsl_vector_free( workspace );
  array_free( mean );
  array_free( Vp );
  array_free( tmp );

  return out;
}
