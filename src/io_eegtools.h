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

/**\file io_eegtools.h
 \brief \ref status_inprogress io_eegtools Custom file-format for input-output of (EEG)-data.

 
 These IO functions can be used for loslessly saving Array-structs and EEG-structs
 to a file. 

 It is easy to write multiple arrays to the same file:
 \code
 FILE *out;
 out=fopen( NAMEOFFILE, "wb" );
 array_to_file( out, a );
 array_to_file( out, b );
 array_to_file( out, c );
 fclose( f );
 \endcode

 and equally easy to read multiple arrays:
 \code 
 FILE *in;
 in=fopen( NAMEOFFILE, "rb" );
 Array *a=array_from_file( in );
 Array *b=array_from_file( in );
 Array *c=array_from_file( in );
 fclose( f );
 \endcode

 \section arrfileformat File Format for Arrays

 When writing an Array-struct to a file, the following convention is 
 followed:
	 - 2 bytes data-type
	 - 2 bytes sizeof(dtype)
	 - 2 bytes number of dimensions
	 - ndim x 2 bytes size-array (number of elements per dimension)
	 - 4 bytes number of bytes in data-field (remaining)
	 - nbytes data-field

 \section eegtoolsfileformat File Format for EEG-data

 When writing an EEG-struct to a file, the following convention is 
 followed:

 */
#ifndef IO_EEGTOOLS_H
# define IO_EEGTOOLS_H

#include <stdio.h>
#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  
  /* -------------- READER ---------------- */
  Array* array_from_file( FILE *in );
  
  /* -------------- WRITER ---------------- */
  void array_to_file( FILE *out, const Array *a );
  void array_to_filename( const char *fname, const Array *a, bool append );


#ifdef __cplusplus
}
#endif

#endif /* IO_EEGTOOLS_H */
