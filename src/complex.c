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

#include "complex.h"
#include <math.h>
#include "mathadd.h"

/** \brief constructor.
	 \param re real part
	 \param im imaginary part
	 \return Complex struct
*/
Complex complex( double re, double im ){
  Complex a;
  a.re = re;
  a.im = im;
  return a;
}

/** \brief complex addition.
 */
Complex complex_add( Complex a, Complex b ){
  a.re += b.re;
  a.im += b.im;
  return a;
}

/** \brief complex subtraction.
 */
Complex complex_sub( Complex a, Complex b ){
  a.re -= b.re;
  a.im -= b.im;
  return a;
}

/** \brief complex multiplication.
 */
Complex complex_mul( Complex a, Complex b ){
  Complex r;
  r.re = (a.re*b.re)-(a.im*b.im);
  r.im = (a.im*b.re)+(a.re*b.im);
  return r;
}

/** \brief complex multiplication with real number.
 */
Complex complex_mul_double( Complex a, double b ){
  a.re *= b;
  a.im *= b;
  return a;
}

/** \brief absolute value of complex number.
 */
double  complex_abs( Complex a ){
  return sqrt( SQR( a.re ) + SQR( a.im ) );
}

/** \brief exponential function of a complex number.

	 this is
	 \f[
	 \exp( a+bi ) = \exp( a )\cdot (\cos b + i\sin b )
	 \f]
 */
Complex complex_exp( Complex a ){
  Complex r;
  r = complex( cos( a.im ), sin( a.im ) );
  r = complex_mul_double( r, exp( a.re ) );
  return r;
}

/** \brief complex conjugate.
 */
Complex complex_conj( Complex a ){
  a.im = -a.im;
  return a;
}

/** \brief complex negate.
 */
Complex complex_neg ( Complex a ){
  a.re = -a.re;
  a.im = -a.im;
  return a;
}

/** \brief complex square root.
 */
Complex complex_sqrt( Complex x ){
  double r = complex_abs(x);
  Complex z = complex(sqrt(0.5 * (r + x.re)),
							 sqrt(0.5 * (r - x.re)));
  if (x.im < 0.0) z.im = -z.im;
  return z;
}
/** \brief complex division.

	 this is
	 \f[ 
	 \frac{(a + bi)}{(c + di)} = \left({ac + bd \over c^2 + d^2}\right) 
	 + \left( {bc - ad \over c^2 + d^2} \right)i
	 \f]
*/
Complex complex_div( Complex a, Complex b){
  Complex r;
  r.re = (a.re*b.re + a.im*b.im)/( SQR( b.re )+SQR( b.im ) );
  r.im = (a.im*b.re - a.re*b.im)/( SQR( b.re )+SQR( b.im ) );
  return r;
}

/** \brief bilinear transform of a complex number.
 */
Complex complex_bilinear_transform(Complex pz){
  return complex_div( complex_add( complex(2.0,0.0), pz ), complex_sub( complex(2.0,0.0), pz ) );
}
