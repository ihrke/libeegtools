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

#ifndef MEX_UTILS_H
#define MEX_UTILS_H

#include "mex.h"
#include "matrix.h"
#include "array.h"

/* convenience stuff */
mxArray *create_mex_double( double val );
mxArray *create_mex_string( const char *str );
char *get_mex_string( const mxArray *arg );
double get_mex_double( const mxArray *arg );
float get_mex_float( const mxArray *arg );
bool is_mex_scalar( const mxArray *arg );
bool is_mex_fid( const mxArray *arg );
int set_mex_field( mxArray *arg, const char *field_name, mxArray *value );
bool is_mex_matrix( mxArray *arg );
bool is_mex_vector( mxArray *arg );
bool mex_have_same_size(mxArray *a, mxArray *b );

char* get_mfile_as_string( const char *fname );

/* array conversion */
Array* mex_mxarray_to_array( mxArray * mxa );
mxArray* mex_array_to_mxarray( Array * a );
mxArray* mex_int_array_to_mxarray( Array * a );
Array* mex_int_mxarray_to_array( mxArray *mxa );

#endif
