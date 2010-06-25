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

#include "array.h" 
#include "helper.h" 
#include <stdarg.h> 
#include <string.h>
#include <gsl/gsl_rng.h>
#include <time.h>

/** \brief Get maximum element from array.
	 
	 You can easily find the index of this element by doing
	 \code
	 Array *a = get_array_from_somewhere();
	 long idx = (array_max(a)-(a->data))/a->dtype_size;
	 // if you need an index-tuple, you can do
	 uint idxtuple[num_dimension];
	 array_calc_rowindex( idx, a->size, a->ndim, idxtuple );
	 \endcode

	 \note this function is a bit slow, as it converts everything
	      to double before comparing. If you need it fast, you
			should consider writing a specialised function.

	 \param a the array
	 \return pointer to the maximum element
 */
void*  array_max( const Array *a ){
  ulong i;
  double maxel=DBL_MIN;
  double tmp;
  void *rmaxel=NULL, *mem=NULL;
  for( i=0; i<array_NUMEL(a); i++ ){
	 mem = array_INDEXMEM1(a,i);
	 array_dtype_to_double( &tmp, mem, a->dtype );
	 if( tmp>maxel ){
		maxel=tmp;
		rmaxel=mem;
	 }
  }
  return rmaxel;
}

/** \brief Get minimum element from array.
	 
	 You can easily find the index of this element by doing
	 \code
	 Array *a = get_array_from_somewhere();
	 long idx = (array_min(a)-(a->data))/a->dtype_size;
	 // if you need an index-tuple, you can do
	 uint idxtuple[num_dimension];
	 array_calc_rowindex( idx, a->size, a->ndim, idxtuple );
	 \endcode

	 \note this function is a bit slow, as it converts everything
	      to double before comparing. If you need it fast, you
			should consider writing a specialised function.

	 \param a the array
	 \return pointer to the minimum element
 */
void*  array_min( const Array *a ){
  ulong i;
  double minel=DBL_MAX;
  double tmp;
  void *rminel=NULL, *mem=NULL;
  for( i=0; i<array_NUMEL(a); i++ ){
	 mem = array_INDEXMEM1(a,i);
	 array_dtype_to_double( &tmp, mem, a->dtype );
	 if( tmp<minel ){
		minel=tmp;
		rminel=mem;
	 }
  }
  return rminel;
}

/** \brief Array-typecast.

	 Typecast the whole array.

	 Warning: 
	 - in case that sizeof(target_type)<sizeof(current_type), you lose precision
	 - in case that sizeof(target_type)!=sizeof(current_type), it is slow because
	   the whole data is copied
	 
	 \param a the array
	 \param target_type the target type
 */
void   array_typecast( Array *a, DType target_type ){
  int srcsize=0,destsize=0;
  ulong nel;
  ulong i;
  array_SIZEOF_DTYPE( srcsize, a->dtype );
  array_SIZEOF_DTYPE( destsize, target_type );
  nel = array_NUMEL( a );
  
  double tmpel;
  void *mem;
  
  if( srcsize==destsize ){
	 a->dtype=target_type;
	 return;
  } 

  /* cast everything to double and back */
  void *newdata=malloc( destsize*nel );

  for( i=0; i<nel; i++ ){
	 mem=a->data+(i*srcsize);
	 array_dtype_to_double( &tmpel, mem, a->dtype );
	 array_MEMSET( newdata+(i*destsize), target_type, tmpel );
  }

  free(a->data);
  a->data=newdata;
  a->dtype=target_type;
  a->dtype_size=destsize;
  a->nbytes=destsize*nel;
}


/** \brief concatenate two arrays.

	 The two arrays must have the same number of elements in all dimensions 
	 except in the concatenated dimension.

	 \todo currently works only for 1D and 2D. implement nD if necessary.

	 Example:
	 \verbatim
	 [ 1 2     [ 7 8 
	   3 4       9 10
      5 6 ]    11 12] 

		=> [ 1 2
		     3 4
           5 6
           7 8
           9 10
			  11 12 ]
	 \endverbatim

	 \param a,b the two arrays to concatenate (must be 1D/2D arrays of arbitrary
	      type); if one of them is NULL, return a copy of the other one
	 \param dim the dimension along which they are to be concatenated (0-rows, 1-columns)
	 \return the concatenation of a and b
 */
Array* array_concatenate( const Array *a, const Array *b, int dim ){
  Array *out;
  if( a==NULL && b==NULL ){
	 warnprintf("concatenating two NULL-arrays\n");
	 return NULL;
  }
  if( a==NULL ){
	 out=array_copy( b, TRUE ); return out;
  }
  if( b==NULL ){
	 out=array_copy( a, TRUE ); return out;
  }
  if( dim!=0 && dim!=1 ){
	 errprintf("can only concatenate rows (0) or columns (1), got %i\n",dim);
	 return NULL;
  }
  if( a->dtype!=b->dtype ){
	 errprintf("Arrays must be of same data-type\n");
	 return NULL;
  }
  if( a->ndim!=b->ndim ){
	 errprintf("Arrays must be of same dimensionality\n");
	 return NULL;
  }
  if( a->ndim>2 || b->ndim>2 ){
	 errprintf("Arrays must be 1D or 2D\n");
	 return NULL;
  }
  if( a->size[1-dim]!=b->size[1-dim] ){
	 errprintf("Arrys must be of same dimension in dim %i, have %i vs. %i\n",
				  1-dim, a->size[1-dim], b->size[1-dim] );
	 return NULL;
  }

  Array *sa,*sb;					  /* wrap 1D arrays */
  sa = array_fromptr2( a->dtype, 2, a->data, a->size[0], (a->ndim>1)?(a->size[1]):1 );
  sb = array_fromptr2( b->dtype, 2, b->data, b->size[0], (b->ndim>1)?(b->size[1]):1 );

  int outr, outc;
  if( dim==0 ){
	 outr=sa->size[0]+sb->size[0];
	 outc=sa->size[1];
  } else {
	 outr=sa->size[0];
	 outc=sa->size[1]+sb->size[1];
  }
  out = array_new2( a->dtype, 2, outr, outc );

  if( dim==0 ){
	 memcpy( out->data, sa->data, sa->nbytes );
	 memcpy( out->data+sa->nbytes, sb->data, sb->nbytes );
  } else {
	 int i;
	 for( i=0; i<a->size[0]; i++ ){
		memcpy( array_INDEXMEM2( out, i, 0), array_INDEXMEM2( a, i, 0 ), a->dtype_size*a->size[1] );
		memcpy( array_INDEXMEM2( out, i, a->size[1]), 
				  array_INDEXMEM2( b, i, 0 ), b->dtype_size*b->size[1] );
	 }
  }

  array_free( sa );
  array_free( sb );
  
  array_dimred( out );
  return out;
}


/** \brief calculate the row-major index-tuple given the offset from element 0.

	 It is assumed, that the data is in row-major format.

	 \param offset 
	 \param size the dimensions of the array
	 \param nsize number of dimensions
	 \param index (output) the index tuple
 */
void array_calc_rowindex( ulong offset, const uint *size, uint nsize, uint *index ){
  ulong prod=1;
  ulong accounted=0;
  int i;
  if( nsize<=0 ){
	 errprintf("invalid dimensionality=%i\n", nsize);
	 return;
  }
  for( i=0; i<nsize; i++ ) prod*=size[i];

  for( i=0; i<nsize; i++ ){
	 index[i] = (offset-accounted)/(prod/size[i]);

	 prod /= size[i];
	 accounted += index[i]*prod;
  }
}

/** \brief calculate the column-major index-tuple given the offset from element 0.

	 It is assumed, that the data is in row-major format.
	 
	 \param offset 
	 \param size the dimensions of the array
	 \param nsize number of dimensions
	 \param index (output) the index tuple
 */
void array_calc_colindex( ulong offset, const uint *size, uint nsize, uint *index ){
  ulong prod=1;
  ulong accounted=0;
  int i;
  if( nsize<=0 ){
  	 errprintf("invalid dimensionality=%i\n", nsize);
  	 return;
  }
  for( i=0; i<nsize; i++ ) prod*=size[i];

  for( i=0; i<nsize; i++ ){
  	 index[nsize-i-1] = (offset-accounted)/(prod/size[nsize-i-1]);
	 
  	 prod /= size[nsize-i-1];
  	 accounted += index[nsize-i-1]*prod;
  }
}

/** \brief convert a row-major array to col-major.
	 
	 This is a d-dimensional matrix-transpose.

	 \param a the row-major array
	 \param alloc if TRUE, allocate new memory, else overwrite a
	 \return a new col-major array
*/
Array* array_convert_rowcolmajor( Array *a, bool alloc){
  ulong i;
  uint *idx=(uint*)malloc( a->ndim*sizeof(uint) );
  Array *b = array_copy( a, TRUE );

  for( i=0; i<array_NUMEL( a ); i++ ){
	 array_calc_colindex( i, a->size, a->ndim, idx );
	 memcpy( array_INDEXMEM1( b, i ), array_index( a, idx ), a->dtype_size );
  }
 
  if( !alloc ){
	 memcpy( a->data, b->data, b->nbytes );
	 array_free( b );
	 b=a;
  }

  return b;
}

/** \brief Reverse order of elements in a (in-place).
	 
	 Dimensionality is not considered. I.e., the function
	 loops over all elements such as they are stored in memory.
	 
	 \param a the array
*/
void   array_reverse( Array *a ){
  int i;
  void *loc1,*loc2;
  int n = array_NUMEL(a);
  void *tmp;
  tmp = malloc( a->dtype_size );
  for( i=0; i<n/2; i++ ){
	 loc1=array_INDEXMEM1( a, i );
	 loc2=array_INDEXMEM1( a, n-i-1 );
	 memcpy( tmp, loc1, a->dtype_size );
	 memcpy( loc1, loc2, a->dtype_size );
	 memcpy( loc2, tmp, a->dtype_size );
  }
  free( tmp );
}

/** \brief Compare two array's dimensions and datatype.

	 \param a,b the two arrays to be tested
	 \return TRUE if both arrays have same dimensionality, FALSE otherwise.
*/
bool array_comparable( const Array *a, const Array *b ){
  int i;
  if( a->dtype!=b->dtype ){															
	 char *dts1="", *dts2="";															
	 array_DTYPESTRING( dts1, a->dtype );											
	 array_DTYPESTRING( dts2, b->dtype );											
	 errprintf("arrays do not have the same datatype: '%s' vs. '%s'\n", dts1, dts2 ); 
	 return FALSE;
  } 
  if( a->ndim != b->ndim ){													
	 errprintf("arrays do not have the same dimensionality, %i vs. %i\n", a->ndim, b->ndim ); 
	 return FALSE;
  } 
  for( i=0; i<a->ndim; i++ ){
	 if( a->size[i] != b->size[i] ){
		errprintf("Arrays differ in dimension %i: %i vs. %i elements\n",
					 i,  a->size[i], b->size[i] );
		return FALSE;
	 }
  }
  return TRUE;
}
  
/** \brief makes a copy of an array.

	 The flag allocdata determines, whether the memory is copied. 
	 \param in input
	 \param allocdata allocate own data, or share 
	       data between arrays. Shared memory belongs to
			 to the passed array.
	 \return the copy 
 */
Array* array_copy( const Array *in, bool allocdata ){
  Array *out;
  if( allocdata ){
	 out = array_new( in->dtype, in->ndim, in->size );
	 memcpy( out->data, in->data, in->nbytes );
  } else {
	 out = array_fromptr( in->dtype, in->ndim, in->data, in->size );
  }
  return out;
}

/** \brief Delete unnecessary dimensions (of size 1).
	 
	 Dimensions of size 1 are removed from the size-array of the
	 struct.

	 \param a in-place reduction
	 \return 0 on success, other on failure
 */
int array_dimred( Array *a ){
  int i,j;
  int nd=0;
  uint *nsize;
  for( i=0; i<a->ndim; i++ ){
	 if( a->size[i]!=1 )
		nd++;
  }
  if( nd==a->ndim ) return 0;

  MALLOC( nsize, nd, uint );
  j=0;
  for( i=0; i<a->ndim; i++ ){
	 if( a->size[i]!=1 )
		nsize[j++] = a->size[i];
  }
  free( a->size );
  a->size=nsize;
  a->ndim=nd;

  return 0;
}

/** \brief multiply all entries in array with x.
	 \param a the input/output array
	 \param x the scalar to multiply a with
 */
int array_scale (Array * a, double x){
  long i;
  double r;
  for( i=0; i<a->nbytes/a->dtype_size; i++ ){
	 array_dtype_to_double( &r, a->data+i*a->dtype_size, a->dtype );
	 array_MEMSET( a->data+i*a->dtype_size, a->dtype, x*r );
  }
  return 0;
}

/**\cond PRIVATE */
/* helper for array-slicing */
uint fill_buffer( void *data, void *buf, int nd, uint **ind, uint *size, uint dtsize, ulong *nbelow ){
  int i;
  int writ=0;
  dprintf("nd=%i,d=%p,b=%p\n",nd, data, buf);
  if( nd==1 ){ /* base case */
	 for( i=0; i<*size; i++ ){
		memcpy( buf+i*dtsize, data+(*ind)[i]*dtsize, dtsize );
		writ += dtsize;
	 }
	 return writ;
  } else { /* recurse into dimensions */
	 for( i=0; i<*size; i++ ){
		writ += fill_buffer( data+((*ind)[i])*(*nbelow)*dtsize, buf+writ, nd-1, ind+1, 
									size+1, dtsize, nbelow+1 );
	 }
	 return writ;
  }
}
/**\endcond */


/**\cond PRIVATE */
/* helper for slicing */
void parse_slicedesc( const Array *a, const char *slicedesc, uint *size, uint **ind ){ 
  int i,j;
  char *desc, *cptr1, *cptr2, *orgdesc;

  desc = strdup( slicedesc );

  orgdesc=desc;
  for( i=0; i<a->ndim; i++ ){
	 cptr1=strchr( desc, ',' );
	 if( cptr1 )
		*cptr1='\0';
	 if( strchr( desc,':' ) ){ /* keep dimension */
		size[i] = a->size[i];
		MALLOC( ind[i], size[i], uint );
		for( j=0; j<size[i]; j++ )
		  ind[i][j]=j;
		dprintf("Keep original dimenson\n");
	 } else if( (cptr2=strchr( desc,'[' ))!=NULL ){ /* array indexing */
		MALLOC( ind[i], a->size[i], uint ); /* 'worst' case */
		j=0;
		size[i]=0;
		cptr2++;
		while( *cptr2 != ']' ){
		  ind[i][j]=atoi( cptr2 );
		  dprintf("Found index '%i'\n", ind[i][j] );
		  size[i]++;
		  j++;
		  cptr2=strchr( cptr2+1, ' ' );
		  if(!cptr2) break;
		}
	 } else if( (cptr2=strchr( desc, '-' ))!=NULL ){ /* range indexing */
		int a1,a2;
		a1=atoi( desc );
		a2=atoi( cptr2+1 );
		dprintf("range %i-%i\n", a1, a2 );
		size[i]=(a2-a1)+1;
		MALLOC( ind[i], size[i], uint );
		for( j=0; j<(a2-a1)+1; j++ ){
		  dprintf("j=%i\n", j);
		  ind[i][j] = a1+j;
		}
	 } else { /* drop dimension using index */
		size[i] = 1;
		MALLOC( ind[i], 1, uint );
		ind[i][0]=atoi( desc );
		dprintf("Found index '%i'\n", ind[i][0] );
	 }
	 desc=cptr1+1;
  }
  free( orgdesc );
}
/**\endcond */

/** \brief Extract sub-arrays from array.

	 This is a function for creating array slices from
	 an existing array. I.e. you can extract sub-arrays
	 similar to interpreted languages like python or MATLAB. 
	 
	 See \ref slicedesc for the format of the slice string.

	 \param a the array
	 \param slicedesc the description of the slice as described in \ref slicedesc
	 \return a freshly allocated array containing the slice
 */
Array* array_slice( const Array *a, const char *slicedesc ){
  Array *b;
  uint *size, **ind;
  int i;

  if( strcount( slicedesc, ',' )!=a->ndim-1 ){
	 errprintf("Slice Description does not contain enough dimensions (need %i)\n",a->ndim );
	 return NULL;
  }
  MALLOC( size, a->ndim, uint );
  MALLOC( ind,  a->ndim, uint*);

  /* parse description */
  parse_slicedesc( a, slicedesc, size, ind );

  /* prepare array slice */
  long bufn=0;
  ulong *nbelow; /* number of elements below this dimension */
  MALLOC( nbelow, a->ndim, ulong );
  nbelow[0]=a->nbytes/a->dtype_size; // total number of elements
  for( i=1; i<a->ndim; i++ ){
	 nbelow[i] = nbelow[i-1]/a->size[i-1];
  }
  dprintf("a->ndim=%i\n", a->ndim );
  b = array_new( a->dtype, a->ndim, size );  
  bufn=fill_buffer( a->data, b->data, a->ndim, ind, size, a->dtype_size, nbelow+1 );
  if( bufn!=b->nbytes ){
	 errprintf("Somethings wrong, wrote %li bytes but needed to write %li\n",
				  bufn, b->nbytes );
  }

  /* remove 1-element dimensions */
  array_dimred( b );

  /* cleaning up */
  free( nbelow );
  free( size );
  for( i=0; i<a->ndim; i++ ){
	 free( ind[i] );
  }
  free( ind );

  return b;
}

/**\cond PRIVATE */
/* helper for data printing */
void dump_data( FILE *out, DType dt, uint dtsize, uint nd, uint nel, void *data, 
					 uint *size, ulong *nbelow ){
  int i, len=nel;
  if( nel<=0 || nel>=*size )
	 len=*size;
  if( nd==1 ){ /* base case */
	 fprintf( out, "[ ");
	 for( i=0; i<len; i++ ){
		array_DTYPEPRINT( out, dt, data+i*dtsize );
		fprintf( out, ", " );
	 }
	 if( len<*size )
		fprintf( out, "... " );
	 fprintf( out, "]\n" );
  } else { /* recurse */
	 fprintf( out, "[ " );
	 for( i=0; i<len; i++ ){
		dump_data( out, dt, dtsize, nd-1, nel, data+i*(*nbelow)*dtsize, size+1, nbelow+1 );
	 }
	 if( len<*size )
		fprintf( out, "...\n" );
	 fprintf( out, "]\n" );
  }
}
/**\endcond */

/** \brief print an array.

	 Print to an output stream, nel_per_dim elements per 
	 dimension. Passing -1 for nel_per_dim prints everything.
	 \param a the array
	 \param nel_per_dim number of elements to print per dimension
	 \param out the output file/stream
 */
void array_print( Array *a, uint nel_per_dim, FILE *out ){ 
  int i;
  ulong *nbelow; /* number of elements below this dimension */
  MALLOC( nbelow, a->ndim, ulong );
  nbelow[0]=a->nbytes/a->dtype_size; // total number of elements
  for( i=1; i<a->ndim; i++ ){
	 nbelow[i] = nbelow[i-1]/a->size[i-1];
  }
  char *dt=NULL;
  array_DTYPESTRING(dt, a->dtype);
  fprintf( out, "array(%s, ", dt ); 
  for( i=0; i<a->ndim; i++ ){
	 fprintf( out, "%i, ", a->size[i] );
  }
  fprintf( out, "\b\b):\n" );
  dump_data( out, a->dtype, a->dtype_size, a->ndim, nel_per_dim, a->data, a->size, nbelow+1 );
  fprintf( out, "\n");

  free( nbelow );
}

/** \brief Index the array. 

	 You need to provide as many arguments 
	 as there are dimensions in the array.
	 A pointer to the corresponding memory in the array is returned.
	 To use it you need to cast it:
	 \code
	 uint idx[3]={1,2,3};
	 float a = *(float*)array_index( arr, idx );
	 \endcode
	 \param a the array
	 \param idx array of size a->ndim giving the indices
	 \return pointer to the memory at the specific location
 */
void*  array_index( const Array *a, uint *idx ){
  int i;
  uint index=0;
  uint cumdim;

  if( a->nbytes % a->dtype_size!=0 ){
	 errprintf("Corrupt array, nbytes not divisible by dtype_size\n");
  }
  cumdim = a->nbytes/a->dtype_size; // total number of elements
  for( i=0; i<a->ndim; i++ ){
 	 /* if( idx[i]>=a->size[i] ){  /\* check array bounds *\/ */
	 /* 	errprintf("Array out of bounds at dim=%i: i=%i, size=%i\n", */
	 /* 				 i, idx[i], a->size[i] ); */
	 /* 	return NULL; */
	 /* } */
	 cumdim /= a->size[i];
	 index += idx[i]*cumdim*a->dtype_size;
  }
  return a->data+index;
}

/** \brief Index the array. 

	 You need to provide as many arguments 
	 as there are dimensions in the array.
	 A pointer to the corresponding memory in the array is returned.
	 To use it you need to cast it:
	 \code
	 float a = *(float*)array_index2( arr, 1, 2 );
	 \endcode
	 \param a the array
	 \param ... a->ndim integers
	 \return pointer to the memory at the specific location
 */
void*  array_index2( const Array *a, ... ){
  va_list ap;
  uint *idx;
  int i;
  void *r;
  MALLOC( idx, a->ndim, uint );

  va_start (ap, a ); 
  for( i=0; i<a->ndim; i++ ){
	 idx[i] = (uint)va_arg( ap, uint );
  }
  va_end (ap);                  /* Clean up. */
  r = array_index( a, idx );
  free( idx );

  return r;
}

/** \brief Initialize new array struct. 

	 Memory for the data is allocated
	 and set to zero.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param dims number of elements in each of the dimensions
	 \return array
*/
Array *array_new ( DType dtype, uint ndim, const uint *dims ){
  Array *a; 
  int i;
  MALLOC( a, 1, Array );
  a->dtype=dtype;
  array_SIZEOF_DTYPE( a->dtype_size, dtype );
  a->ndim=ndim;  
  a->free_data=TRUE;
  MALLOC( a->size, ndim, uint );

  /* get size of dimensions */
  a->nbytes=1;
  for( i=0; i<ndim; i++ ){
	 a->size[i] = dims[i];
	 a->nbytes *= a->size[i];
  }
  a->nbytes *= a->dtype_size;

  /* allocate memory for the data */
  a->data = (void*)malloc( a->nbytes );
  if( !a->data ){
	 errprintf("Failed to allocate memory!\n");
  } else {
	 memset( a->data, 0, a->nbytes );
  }

  return a;
}

/** \brief Initialize new array struct. 

	 Memory for the data is allocated
	 and set to zero.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_new2( DType dtype, uint ndim, ... ){ 
  Array *a;
  va_list ap;
  uint *size;
  int i;
  MALLOC( size, ndim, uint );

  /* get size of dimensions */
  va_start (ap, ndim ); 
  for( i=0; i<ndim; i++ ){
	 size[i] = (uint)va_arg( ap, uint );
  }
  va_end (ap);                  /* Clean up. */

  a = array_new( dtype, ndim, size );

  free( size );
  return a;
}

/** \brief create array that is filled with random values 
	     from [0, ..., 1].
	 
	 Mainly used for debugging stuff.
	 \param seed to use for the random numbers; if seed=0, the current
	        time is used
	 \param ndim number of dimensions
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_randunif( unsigned long seed, uint ndim, ... ){
  va_list ap;
  Array *a;
  long i;
  int n; 
  uint *size;
  MALLOC( size, ndim, uint );
  
  /* get size of dimensions */
  va_start (ap, ndim ); 
  for( i=0; i<ndim; i++ ){
	 size[i] = (uint)va_arg( ap, uint );
  }
  va_end (ap);                  /* Clean up. */

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  if( seed==0 )
	 seed=(unsigned long)time(NULL);
  gsl_rng_set (r, seed );
  a = array_new( DOUBLE, ndim, size );
  n = a->nbytes/a->dtype_size;
  for( i=0; i<n; i++ ){
	 array_INDEX1( a, double, i )= gsl_rng_uniform (r);
  }
     
  gsl_rng_free (r);
  free( size );
  return a;
}

/** \brief shuffle the entries of an array.
	 \param a the array (arbitrary type and dimensionality)
	 \param seed if 0, use time(NULL), else the seed
 */
void   array_shuffle( Array *a, unsigned long seed ){
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  if( seed==0 )
	 seed=(unsigned long)time(NULL);
  gsl_rng_set (r, seed );

  gsl_ran_shuffle( r, a->data, array_NUMEL(a), a->dtype_size );

  gsl_rng_free (r);
}


/** \brief create array that is filled with values from 
	        1 to n over all dimensions.
	 
	 Mainly used for debugging stuff.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_new_dummy( DType dtype, uint ndim, ... ){
  va_list ap;
  Array *a;
  long i;
  int n; 
  uint *size;
  MALLOC( size, ndim, uint );
  
  /* get size of dimensions */
  va_start (ap, ndim ); 
  for( i=0; i<ndim; i++ ){
	 size[i] = (uint)va_arg( ap, uint );
  }
  va_end (ap);                  /* Clean up. */

  a = array_new( dtype, ndim, size );
  n = a->nbytes/a->dtype_size;
  for( i=0; i<n; i++ ){
	 array_MEMSET( a->data+i*a->dtype_size, a->dtype, i );
  }

  free( size );
  return a;
}

/** \brief Initialize new array struct. 

	 Memory for the data is not allocated 
	 but set to data. 
	 \warning the  array does NOT take over the responsibility 
	 for the memory. I.e. if you call array_free() it is NOT free'd.
	 If you want this to happen, you need to set array->free_data 
	 to TRUE.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param data the data for the array
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_fromptr2( DType dtype, uint ndim, void *data, ... ){
  va_list ap;
  Array *a; 
  int i;
  uint *size; 
  MALLOC( size, ndim, uint );

  /* get size of dimensions */
  va_start (ap, data ); 
  for( i=0; i<ndim; i++ ){
	 size[i] = (uint)va_arg( ap, uint );
  }
  va_end (ap);                  /* Clean up. */

  a = array_fromptr( dtype, ndim, data, size );

  free( size );

  return a;
}

/** \brief Initialize new array struct. 

	 Memory for the data is not allocated 
	 but set to data. 
	 \warning the  array does NOT take over the responsibility 
	 for the memory. I.e. if you call array_free() it is NOT free'd.
	 If you want this to happen, you need to set array->free_data 
	 to TRUE.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param data the data for the array
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_fromptr( DType dtype, uint ndim, void *data, const uint *size ){
  Array *a; 
  int i;
  MALLOC( a, 1, Array );
  a->dtype=dtype;
  array_SIZEOF_DTYPE( a->dtype_size, dtype );
  a->ndim=ndim;  
  a->free_data=FALSE;
  MALLOC( a->size, ndim, uint );

  /* get size of dimensions */
  a->nbytes=1;
  for( i=0; i<ndim; i++ ){
	 a->size[i] = size[i];
	 a->nbytes *= a->size[i];
  }
  a->nbytes *= a->dtype_size;
  
  a->data = data;

  return a;
}

/** \brief cast memory of type dt in location mem to double.

	 \param out output double
	 \param mem the memory
	 \param dt the DType of mem
 */
void   array_dtype_to_double( double *out, void *mem, DType dt ){
  switch( dt ){
  case CHAR:	
	 *out = *((char*)mem);break;
  case UINT:					
	 *out = *((uint*)mem);break;
  case INT:						
	 *out = *((int*)mem);break;
  case LONG:					
	 *out = *((long*)mem);break;
  case ULONG:					
	 *out = *((ulong*)mem);break;
  case FLOAT:						 
	 *out = *((float*)mem);break;
  case DOUBLE:					
	 *out = *((double*)mem);break;
  default:						
	 warnprintf("Do not know this datatype: %i\n", dt );		
  }
}



/** \brief free all memory associated with the array.
 */
void  array_free( Array *a ){
  if( !a ) return;
  if( a->free_data && a->data ) free( a->data );
  if( a->size ) free( a->size );
  free( a );
}
