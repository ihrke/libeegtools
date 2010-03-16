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

/**\file regularization.h
	\brief \ref status_unstable Regularization matrices.
 */
#ifndef REGULARIZATION_H
# define REGULARIZATION_H
#include "mathadd.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif


  /** \addtogroup regularization
	*\{
	*/
  double** regularization_linear_points( double **d, int n, 
													  OptArgList *optargs );
  double** regularization_gaussian_line( double **d, int n, 
													  OptArgList *optargs );
  /** \} */

#ifdef __cplusplus
}
#endif


#endif
