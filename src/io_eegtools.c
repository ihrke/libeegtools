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
#include "helper.h"

#include "io_eegtools.h"
#include <stdlib.h>
#include <stdint.h>


/* -------------- READER ---------------- */
/** \brief Read Array from file.

	 The FILE* must be opened in read-binary "rb" mode.
	 
	 The format of the array is described in 
	 \ref arrfileformat .

	 \param out the file pointer to read from (opened read-binary "rb")
	 \return the output array
 */ 
Array* array_from_file( FILE *in ){
  uint16_t dtype;
  uint16_t size_dtype;
  uint16_t ndim;
  uint16_t size;
  uint32_t nbytes;
  uint *dims;
  int i;
  
  ffread( &dtype, 2, 1, in );
  ffread( &size_dtype, 2, 1, in );
  ffread( &ndim, 2, 1, in );
  MALLOC( dims, ndim, uint );
  for( i=0; i<ndim; i++ ){
	 ffread( &size, 2, 1, in );
	 dims[i]=(uint)size;
  }
  Array *a=array_new( dtype, ndim, dims );
  ffread( &nbytes, 4, 1, in );
  if( nbytes!=a->nbytes ){
	 errprintf("Wrong number of bytes, continuing (probably dieing)\n");
  }
  ffread( a->data, 1, nbytes, in );

  free( dims );

  return a;
}

/* -------------- WRITER ---------------- */

/** \brief Write Array to disk.

	 The FILE* must be opened in write-binary "wb" mode.
	 
	 The format of the array is described in 
	 \ref rawfileformat .

	 \param out the file pointer to write to (opened write-binary "wb")
	 \param a the output array
 */ 
void array_to_file( FILE *out, const Array *a ){
  uint16_t ui16;
  uint32_t ui32;
  int i;

  ui16=(uint16_t)a->dtype;
  ffwrite( &ui16, 2, 1, out );
  ui16=(uint16_t)a->dtype_size;
  ffwrite( &ui16, 2, 1, out );
  ui16=(uint16_t)a->ndim;
  ffwrite( &ui16, 2, 1, out );
  for( i=0; i<a->ndim; i++ ){
	 ui16=(uint16_t)a->size[i];
	 ffwrite( &ui16, 2, 1, out );	 
  }
  ui32=(uint32_t)a->nbytes;
  ffwrite( &ui32, 4, 1, out );	 
  /* data */
  ffwrite( a->data, 1, a->nbytes, out );
}
/** \brief Write Array to a given filename.

	 The format of the array is described in
	 \ref rawfileformat .

	 \param fname name of the file
	 \param a the output array
	 \param append if TRUE, open file in 'ab' mode, else in 'wb' mode
 */
void array_to_filename( const char *fname, const Array *a, bool append ){
	FILE *f;
	if( append ){
		if( !(f=fopen( fname, "ab" )) ){
			errprintf("Could not open file '%s'\n",fname);
			return;
		}
	} else {
		if( !(f=fopen( fname, "wb" )) ){
			errprintf("Could not open file '%s'\n",fname);
			return;
		}
	}
	array_to_file( f, a );
	fclose( f );
	return;
}
