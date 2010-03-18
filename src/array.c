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
  char *desc, *tmp, *cptr1, *cptr2, *orgdesc;
  uint *size, **ind;
  int i,j;
  desc = strdup( slicedesc );
  tmp = strdup( slicedesc );

  if( strcount( desc, ',' )!=a->ndim-1 ){
	 errprintf("Slice Description does not contain enough dimensions (need %i)\n",a->ndim );
	 return NULL;
  }
  MALLOC( size, a->ndim, uint );
  MALLOC( ind,  a->ndim, uint*);
  
  /* parse description */
  orgdesc=desc;
  for( i=0; i<a->ndim; i++ ){
	 cptr1=strchr( desc, ',' );
	 if( cptr1 )
		*cptr1='\0';
	 strcpy( tmp, desc );
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

  /* cleaning up */
  free( nbelow );
  free( size );
  free( orgdesc );
  free( tmp );
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
	 for( i=0; i<len; i++ ){
		array_DTYPEPRINT( out, dt, data+i*dtsize );
		fprintf( out, ", " );
	 }
	 if( len<*size )
		fprintf( out, "... " );
	 fprintf( out, "\n" );
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
  char *dt;
  array_DTYPESTRING(dt, a->dtype)
  fprintf( out, "array(%s):", dt );
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
	 if( idx[i]>=a->size[i] ){  /* check array bounds */
		errprintf("Array out of bounds at dim=%i: i=%i, size=%i\n",
					 i, idx[i], a->size[i] );
		return NULL;
	 }
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
Array *array_new ( DType dtype, uint ndim, uint *dims ){
  Array *a; 
  int i;
  MALLOC( a, 1, Array );
  a->dtype=dtype;
  array_SIZEOF_DTYPE( a->dtype_size, dtype );
  a->ndim=ndim;  
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
  int n;
  a = array_new2( dtype, ndim, ap );
  n = a->nbytes/a->dtype_size;
  
  return a;
}

/** \brief Initialize new array struct. 

	 Memory for the data is not allocated 
	 but set to data. \warning the moment you call that function,
	 the array takes over the responsibility for the memory. I.e. if
	 you call array_free() it is free'd.
	 \param dtype the datatype of the array
	 \param ndim number of dimensions
	 \param data the data for the array
	 \param ... the number of elements in each of the dimensions
	 \return array
*/
Array *array_fromptr( DType dtype, uint ndim, void *data, ... ){
  va_list ap;
  Array *a; 
  int i;
  MALLOC( a, 1, Array );
  a->dtype=dtype;
  array_SIZEOF_DTYPE( a->dtype_size, dtype );
  a->ndim=ndim;  
  MALLOC( a->size, ndim, uint );

  /* get size of dimensions */
  a->nbytes=1;
  va_start (ap, data ); 
  for( i=0; i<ndim; i++ ){
	 a->size[i] = (uint)va_arg( ap, uint );
	 a->nbytes *= a->size[i];
  }
  a->nbytes *= a->dtype_size;
  va_end (ap);                  /* Clean up. */
  
  a->data = data;

  return a;
}



/** \brief free all memory associated with the array.
 */
void  array_free( Array *a ){
  if( !a ) return;
  if( a->data ) free( a->data );
  if( a->size ) free( a->size );
  free( a );
}
