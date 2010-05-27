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

/**\file complex.h
 \brief \ref status_stable complex

 Complex Numbers.
	
 */
#ifndef COMPLEX_H
# define COMPLEX_H
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  typedef struct {
	 double re;
	 double im;
  } Complex;
  

 /* ---------------------------------------------------------------------------- 
	  -- Complex Arithmetic
	  ---------------------------------------------------------------------------- */ 

  Complex complex     ( double re, double im );
  Complex complex_add ( Complex a, Complex b );
  Complex complex_add_dbl ( Complex a, double b );  
  Complex complex_sub ( Complex a, Complex b );
  Complex complex_mul ( Complex a, Complex b );
  Complex complex_mul_double( Complex a, double b );
  double  complex_abs ( Complex a );
  Complex complex_exp ( Complex a );
  Complex complex_conj( Complex a );
  Complex complex_neg ( Complex a );
  Complex complex_div ( Complex a, Complex b);
  Complex complex_sqrt( Complex x );
  Complex complex_bilinear_transform(Complex pz);


#ifdef __cplusplus
}
#endif

#endif /* COMPLEX_H */
