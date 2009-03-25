/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke                                  *
 *   mihrke@uni-goettingen.de                                              *
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

/**\file nonlinear.h
 * Functions using nonlinear systems-theory.

 Phase-Space reconstruction is done as in
 \f[
 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
 \f]
 where \f$ \tau = \nu \Delta t\f$ is a multiple of the sampling step.
	
 * Groups:
 \defgroup nonlinear General
 \{
 \}
 */
#ifndef NONLINEAR_H
# define NONLINEAR_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\addtogroup nonlinear
	*\{
	*/ 

  PhaseSpace* phspace_init ( int m, int tau, double *x, int n );
  void        phspace_free ( PhaseSpace *p );
  void        phspace_print( FILE *out, PhaseSpace *p);

  double      phspace_index_ij( PhaseSpace *p, int i, int j );
  void        phspace_index_i ( PhaseSpace *p, int i, double *x);
  void        phspace_index_j ( PhaseSpace *p, int j, double *x);
  /*\}*/

#ifdef __cplusplus
}
#endif


#endif
