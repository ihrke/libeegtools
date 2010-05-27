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

/**\file imageproc.h
 \brief \ref status_inprogress imageproc

 Image Processing Functions.
	
 */
#ifndef IMAGEPROC_H
# define IMAGEPROC_H
#include "definitions.h"
#include "array.h"

#ifdef __cplusplus
extern "C" {
#endif

  /* -------------- FUNCTIONS ---------------- */

  Array*  bresenham_line( int start[2], int end[2] );
  Array*  bresenham_linesegments( const Array *points );


  Array*  disttransform_deadreckoning( const Array *in, Array *dt );


#ifdef __cplusplus
}
#endif

#endif /* IMAGEPROC_H */
