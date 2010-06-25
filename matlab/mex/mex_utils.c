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

#include <string.h>

#include "mex_utils.h"
#include "linalg.h"



mxArray *create_mex_double( double val )
{
  mxArray *array = mxCreateNumericMatrix( 1, 1, mxDOUBLE_CLASS, mxREAL );
  double *df_ptr = (double*)mxGetData( array );
  *df_ptr = val;
  return array;
}

mxArray *create_mex_string( const char *str )
{
	char *buf = (char*)mxCalloc( strlen( str )+1, sizeof(char) );
	strcpy( buf, str );
	return mxCreateString( buf );
}

char *get_mex_string( const mxArray *arg )
{ 
  /* Input must be a string. */
  if (mxIsChar(arg) != 1)
    mexErrMsgTxt("Input must be a string.");

  /* Input must be a row vector. */
  if (mxGetM(arg) != 1)
    mexErrMsgTxt("Input must be a row vector.");
    
  /* Get the length of the input string. */
  int buflen = (mxGetM(arg) * mxGetN(arg)) + 1;

  /* Allocate memory for input and output strings. */
  char *str = (char*)mxCalloc(buflen, sizeof(char));
  
  /* Copy the string data from prhs[0] into a C string 
   * input_buf. */
  int status = mxGetString(arg, str, buflen);

  if (status != 0) 
	mexWarnMsgTxt("Not enough space. String is truncated.");
 
  return str;
}

double get_mex_double( const mxArray *arg )
{
	double *data = mxGetPr( arg );   
	return data[0];
}

float get_mex_float( const mxArray *arg )
{
	float* data = (float*)mxGetPr( arg );   
	return data[0];
}


int set_mex_field( mxArray *arg, const char *field_name, mxArray *value )
{ 
	int fnum = mxGetFieldNumber( arg, field_name );
	if( fnum == -1 ) {
		fnum = mxAddField( arg, field_name );
	}
	mxSetFieldByNumber( arg, 0, fnum, value );

	return fnum;
}

bool is_mex_scalar( const mxArray *arg )
{
  if( arg == NULL ) return false;
  return mxIsNumeric( arg ) && mxGetN( arg ) == 1 && mxGetM( arg ) == 1;
}

bool is_mex_fid( const mxArray *arg )
{
  if( arg == NULL ) return false;
  return mxIsNumeric( arg ) && (mxGetN( arg ) == 1 || mxGetN( arg ) == 2) && mxGetM( arg ) == 1;
}

bool is_mex_matrix( mxArray *arg ){
  if( arg == NULL ) return false;
  return mxIsDouble( arg ) && mxGetNumberOfDimensions(arg)==2;
}

/** \brief convert MxArray to libeegtools array.
	 
	 Everything must be double-precision.
 */
Array* mex_mxarray_to_array( mxArray * mxa ){
  Array *out;
  int ndim=mxGetNumberOfDimensions(mxa);
  uint *dims=mxGetDimensions(mxa);
  double *ptr=(double*)mxGetPr(mxa);
  out = array_new( DOUBLE, ndim, dims );

  ulong i;
  uint *idx=(uint*)malloc(ndim*sizeof(idx));;
  uint cidx;

  /* col major to row major */
  for( i=0; i<array_NUMEL( out ); i++ ){
	 array_calc_rowindex( i, out->size, out->ndim, idx );
	 cidx=mxCalcSingleSubscript( mxa, ndim, idx);
	 array_INDEX1( out, double, i) = ptr[cidx];
  }
  free( idx );

  array_dimred( out );

  return out;
}

/** \brief convert libeegtools array to MxArray.
	 
	 Everything must be double-precision.
 */
mxArray* mex_array_to_mxarray( Array * a ){
  mxArray *mxa =mxCreateNumericArray(a->ndim, a->size,
												 mxDOUBLE_CLASS, mxREAL );
  
  Array *tmp = array_convert_rowcolmajor( a, TRUE );
  memcpy( mxGetPr(mxa), tmp->data, tmp->nbytes );
  array_free( tmp );
  
  return mxa;
}

/** \brief convert libeegtools array to MxArray.
	 
	 Everything must be int-precision.
*/
mxArray* mex_int_array_to_mxarray( Array * a ){
  mxClassID intclass;
  switch( sizeof(int) ){
  case 1:
	 intclass=mxINT8_CLASS; break;
  case 2:
	 intclass=mxINT16_CLASS; break;
  case 4:
	 intclass=mxINT32_CLASS; break;
  case 8:
	 intclass=mxINT64_CLASS; break;
  default:
	 errprintf("what's that?\n");
	 intclass=mxINT64_CLASS; break;
  }
  mxArray *mxa =mxCreateNumericArray(a->ndim, a->size,
												 intclass, mxREAL );
  
  Array *tmp = array_convert_rowcolmajor( a, TRUE );
  memcpy( mxGetPr(mxa), tmp->data, tmp->nbytes );
  array_free( tmp );
  
  return mxa;
}
