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

/**\file array.h
 \brief \ref status_inprogress Numerical Arrays.

 This is an efficient implementation of a numerical array of
 - arbitrary dimension and
 - arbitrary type.
 It wraps a thin layer around a one-dimensional (void *) for 
 efficiency.

 The library provides some convenient functionality:
 - \ref indexing
 - \ref slicedesc
 - printing

 \section arb Arbitrary Dimensions and Types
 \note in applications you almost always know the type/dimensionality
 of your array. This section is for you only if you want to write
 functions that operate on arrays of arbitrary dimension and/or type.

 In principle, an Array struct can hold n-dimensional arrays of any type 
 implemented in DType. All functions starting with the prefix array_*() 
 operate on arbitrary dimensions and types.

 \subsection arbtyp Arbitrary types
 Working with arbitray types requires some sophistication if it is conducted
 at runtime. 

 To \b read  elements in the array where you don't know the type, you can 
 convert the element to double using array_dtype_to_double()
 \code
 // void   array_dtype_to_double( double *out, void *mem, DType dt );
 void *mem = array_index( a, 3, 4 ); // access element (3,4)
 double el;
 array_dtype_to_double( &el, mem, a->dtype );
 // work with el
 \endcode

 To \b write elements in the array of unknown type, you can use the macro 
 array_MEMSET()
 \code
 // #define array_MEMSET( mem, dtype, val )
 void *mem = array_index( a, 3, 4 ); // access element (3,4)
 array_MEMSET( mem, a->dtype, 3 ); // set element to 3 (value is casted)
 \endcode 

 \subsection arbdim Arbitrary dimension
 To write functions that operate on arbitrary dimensions is either very simple
 (if the operation you want to apply is the same for all elements) or
 very complicated.

 If you want to apply a per-element operation on all elements in the array,
 you can simply loop over them:
 \code
 Array *a=get_array_from_somewhere();
 long i, n = array_NUMEL( a ); // number of elements
 double *mem;
 for( i=0; i<n; i++ ){
    mem=&array_INDEX1( a, double, i ); // if you know the type
	 // ...
	 // operate on *mem
	 // ...
 }
 \endcode

 \subsection vec Vectors and Matrices
 However, for numerical work, it is in most cases sufficient to
 operate on floating point numbers. Therefore, 1D and 2D Array structs
 of DType DOUBLE are used to represent vectors and matrices,
 respectively.  I.e. that each function starting with matrix_*( Array
 *m ) accepts only 2D Arrays of DType DOUBLE functions starting with
 vector_*( Array *m ) accept only 1D Arrays of DType DOUBLE.  This is
 checked with the Macros matrix_CHECK() and vector_CHECK(). See also \ref linalg.h .
	
 \section indexing Convenient Indexing
 Indexing the n-dimensional array can be done via two routes:
 - a very fast but dangerous (not type-safe, not boundary checks) Macro-based approach
    - this is only available for 1, 2 and 3 dimensions
 - a slower but safer approach based on functions
 <ul>
 <li>this is implemented for arbitray dimensions
 <li>it's not all that slow, it has O(ndim) complexity plus some overhead for the
	   function call
 </ul>

 \section slicedesc Slicing Array Objects
 It is often necessary to get parts of arrays for computation.
 Inspired by slicing in MATLAB, a string representing the slice is
 composed as follows:
 - dimensions are separated by commas
 - each dimension can be indexed by:
 <ul>
 <li> 'n' -> a single integer; dimension is deleted 
 <li> ':' -> all elements in this dimension
 <li> 'n1-n2' ->  (n1 and n2 are integers) elements between n1 and n2 (including n1 and n2) are taken
 <li> '[n1 n2 ...]' -> only elements [n1 n2 ...] are taken from the array in this dimension
 </ul>

\note indexing assumes c-like indices (i.e. from 0 to n-1)
\note you need to provide all dimensions 

Examples:

\code
 // 3D-array 
 Array *a = array_new( DOUBLE, 3, 10, 11, 12 );

 // Here b is a new 3D-array containing all elements in the first dimension,
 // the first 4 in the second and elements 1,2,3 in the third. 
 Array *b = array_slice( a, ":,0-3,[1 2 3]");

 // c is a 2D array 
 Array *c = array_slice( a, ":,:,1");

 // d is a 1D array 
 Array *d = array_slice( a, ":,3,1");
\endcode


 \section owner Data Ownership
 Depending on how you initialize and Array, the struct takes over
 responsibility to free the array->data pointer or not. 
 This is implemented using the array->free_data flag. Some functions take a 
 'alloc'-flag which indicates whether memory is shared or not. Usually, the 
 'parent'-array has the resonsibility then.
 \note in practice, you hardly have to care about data-ownership. If you simply
       call array_free() on all arrays in your code you should be fine.

 */

#ifndef ARRAY_H
# define ARRAY_H

#include "definitions.h"
#include "mathadd.h"

#ifdef __cplusplus
extern "C" {
#endif

#define DIM_ROWS 0
#define DIM_COLS 1

  /** \def array_SIZEOF_DTYPE( size, dtype )
		\brief get size of datatype.
		
		fill variable size with the size of the datatype dtype in
		bytes. E.g.
		\code
		array_SIZEOF_DTYPE( array->dtype_size, dtype );
		\endcode 
		\param size (output) uint to hold the size
		\param dtype (input) the datatype
	*/
#define array_SIZEOF_DTYPE( size, dtype )										\
  switch( dtype ) {																	\
  case CHAR:																			\
	size=sizeof(char);break;														\
  case UINT:																			\
	size=sizeof(uint);break;														\
  case INT:																				\
	size=sizeof(int);break;															\
  case LONG:																			\
	size=sizeof(long);break;														\
  case ULONG:																			\
	size=sizeof(unsigned long); break;											\
  case FLOAT:																			\
	size=sizeof(float); break;														\
  case DOUBLE:																			\
	size=sizeof(double); break;													\
  default:																				\
	 warnprintf("Do not know this datatype: %i\n", dtype );				\
  }
   /** \def array_DTYPESTRING(dtstring,dtype)
		 \brief convert datatype to a string representation.
		 \param dtstring (output), char* to hold the string
		 \param dtype (input) the datatype
	*/
#define array_DTYPESTRING(dtstring,dtype)										\
  switch( dtype ){																	\
  case CHAR:																			\
	 dtstring="char";break;															\
  case UINT:																			\
	 dtstring="uint";break;															\
  case INT:																				\
	 dtstring="int";break;															\
  case LONG:																			\
	 dtstring="long";break;															\
  case ULONG:																			\
	 dtstring="unsigned long"; break;											\
  case FLOAT:																			\
	 dtstring="float"; break;														\
  case DOUBLE:																			\
	 dtstring="double"; break;														\
  default:																				\
	 warnprintf("Do not know this datatype: %i\n", dtype );				\
  }
 
  /** \brief get number of elements in array.
		\param a Array*
	*/
#define array_NUMEL( a )								\
  ( ((a)->nbytes)/((a)->dtype_size) )

  /** \brief Print an array-datatype.
		
		\param out FILE* output stream
		\param dtype one of DType
		\param mem void* memory containing the to-be-printed value
	*/
#define array_DTYPEPRINT( out, dtype, mem )									\
  switch( dtype ){																	\
  case CHAR:																			\
	 fprintf( out, "%c", *(char*)(mem) ); break;								\
  case UINT:																			\
	 fprintf( out, "%i", *(uint*)(mem) ); break;								\
  case INT:																				\
	 fprintf( out, "%i", *(int*)(mem) ); break;								\
  case LONG:																			\
	 fprintf( out, "%li", *(long*)(mem) ); break;							\
  case ULONG:																			\
	 fprintf( out, "%li", *(ulong*)(mem) ); break;							\
  case FLOAT:																			\
	 fprintf( out, "%.2f", *(float*)(mem) ); break;							\
  case DOUBLE:																			\
	 fprintf( out, "%.2f", *(double*)(mem) ); break;							\
  default:																				\
	 warnprintf("Do not know this datatype: %i\n", dtype );				\
  }
  
  /** \brief set the memory mem to value of type dtype.

		\param (out) void * mem
		\param dtype one of DType
		\param val int,float,double...
	*/
#define array_MEMSET( mem, dtype, val )											\
  switch( dtype ){																		\
  case CHAR:																				\
	 *((char*)(mem)) = (char)val;break;												\
  case UINT:																				\
	 *((uint*)(mem)) = (uint)val;break;												\
  case INT:																					\
	 *((int*)(mem)) = (int)val;break;													\
  case LONG:																				\
	 *((long*)(mem)) = (long)val;break;												\
  case ULONG:																				\
	 *((ulong*)(mem)) = (ulong)val;break;												\
  case FLOAT:																				\
	 *((float*)(mem)) = (float)val;break;												\
  case DOUBLE:																				\
	 *((double*)(mem)) = (double)val;break;											\
  default:																					\
	 warnprintf("Do not know this datatype: %i\n", dtype );					\
  }

/** \def array_INDEX1( array, dtype, i1 )
	 \brief fast 1D array indexing typecasting to C-type.

	 index the array and assume only one dimension (not checked).
	 These indexing macros are not safe but as fast as you can
	 get. For more convenience, you can use the dimension independant
	 array_index() and array_index2() functions.
	 
	 Example:
	 \code
	 double a;
	 a = array_INDEX1( a, double, 10 ); // accesses a[10] and converts to double
	 \endcode

	 \param array an Array object
	 \param dtype a valid c-data type (int, float, double, ... )
	 \param i1,...,iN  the index for each of the dimensions
*/
#define array_INDEX1( array, dtype, i1 )							\
  (*((dtype*)( ((array)->data) +										\
					((i1)*((array)->dtype_size)) )))
/** \brief fast 2D array indexing typecasting to C-type.
*/
#define array_INDEX2( array, dtype, i1, i2 )										\
  (*((dtype*)( ((array)->data)+															\
					((i1)*((array)->size[1])*((array)->dtype_size))+					\
					((i2)*((array)->dtype_size)) )))
/** \brief fast 3D array indexing typecasting to C-type.
*/
#define array_INDEX3( array, dtype, i1, i2, i3 )								\
  (*((dtype*)( ((array)->data)+															\
					((i1)*((array)->size[1])*((array)->size[2])*((array)->dtype_size))+ \
					((i2)*((array)->size[2])*((array)->dtype_size))+					\
					((i3)*((array)->dtype_size)) )))

/** \def array_INDEXMEM1( array, i1 )
	 \brief fast 1D array indexing returning a void pointer.

	 index the array and assume only one dimension (not checked).
	 These indexing macros are not safe but as fast as you can
	 get. For more convenience, you can use the dimension independant
	 array_index() and array_index2() functions.
	 
	 Example:
	 \code
	 void *a;
	 a = array_INDEXMEM1( a, 10 ); // accesses a[10] of any type
	 \endcode

	 \param array an Array object
	 \param i1,...,iN  the index for each of the dimensions
*/
#define array_INDEXMEM1( array, i1 )									\
  ( ((array)->data) +														\
	 ((i1)*((array)->dtype_size)) )
/** \brief fast 2D array indexing returning a void pointer.
*/
#define array_INDEXMEM2( array, i1, i2 )											\
  ( ((array)->data)+																		\
	 ((i1)*((array)->size[1])*((array)->dtype_size))+								\
	 ((i2)*((array)->dtype_size)) )
/** \brief fast 3D array indexing returning a void pointer.
*/
#define array_INDEXMEM3( array, i1, i2, i3 )										\
  ( ((array)->data)+																		\
	 ((i1)*((array)->size[1])*((array)->size[2])*((array)->dtype_size))+			\
	 ((i2)*((array)->size[2])*((array)->dtype_size))+								\
	 ((i3)*((array)->dtype_size)) )

  /** \brief Data-types for Array struct. */
  typedef enum {
	 CHAR,
	 UINT,
	 INT,
	 LONG,
	 ULONG,
	 FLOAT,
	 DOUBLE
  } DType;

  /** \brief Numerical Array struct. */
  typedef struct {
	 void *data;       /**< data in C-storage format */
	 DType dtype;      /**< datatype of data */
	 uint  dtype_size; /**< sizeof(dtype) in bytes */
	 uint  ndim;       /**< number of dimensions */
	 ulong nbytes;     /**< number of bytes in array in total */
	 uint  *size;      /**< number of elements per dimension */
	 bool  free_data;  /**< should data be free'd ? */
  } Array;

  /* -------------- FUNCTIONS ---------------- */
  Array *array_new ( DType dtype, uint ndim, const uint *dims );
  Array *array_new2( DType dtype, uint ndim, ... );
  Array *array_fromptr ( DType dtype, uint ndim, void *data, const uint *size );
  Array *array_fromptr2( DType dtype, uint ndim, void *data, ... );
  Array *array_new_dummy( DType dtype, uint ndim, ... );
  Array *array_randunif( unsigned long seed, uint ndim, ... );
  bool   array_comparable( const Array *a, const Array *b );

  Array* array_copy( const Array *in, bool allocdata );
  void   array_free( Array *a );

  void   array_shuffle( Array *a, unsigned long seed );

  void   array_typecast( Array *a, DType target_type );

  void*  array_index ( const Array *a, uint *idx );
  void*  array_index2( const Array *a, ... );
  Array* array_slice ( const Array *a, const char *slicedesc );

  Array* array_concatenate( const Array *a, const Array *b, int dim );

  int    array_scale (Array * a, double x);

  void   array_reverse( Array *a );

  int    array_dimred( Array *a );

  void*  array_max( const Array *a );
  void*  array_min( const Array *a );


  Array* array_convert_rowcolmajor( Array *a, bool alloc);
  void   array_calc_colindex( ulong idx, const uint *size, uint nsize, uint *index );
  void   array_calc_rowindex( ulong idx, const uint *size, uint nsize, uint *index );

  void   array_print ( Array *a, uint nel_per_dim, FILE *out );
  void   array_dtype_to_double( double *out, void *mem, DType dt );

#ifdef __cplusplus
}
#endif

#endif /* ARRAY_H */
