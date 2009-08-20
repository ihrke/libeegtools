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

/**\file som.h
 \brief Self-Organizing Maps

 \warning For now this is only experimental but plans are to complete this file.

 Kohonen, T. 1995. Self-Organizing Maps. Springer.

 \todo generalize metric (need to solve Eq. 3.40 (from Kohonen, 1995) for each metric)
 */
#ifndef SOM_H
# define SOM_H

#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif
  
  /**\ingroup som_help 
	\{ */
  Som* som_init( int dimension, int n, int nruns, SOMConnectivityType connectivity );
  void som_free( Som *s );
  void som_print( FILE *f, Som *s );
  /** \} */

  /**\ingroup som_initialize
	\{ */ 
  void som_initialize_random( Som *s, double min, double max );
  void som_initialize_random_samples( Som *s, double **X, int dim, int nsamples );
  /** \} */

  /**\ingroup som_neighbourhood
	\{ */
  double som_neighbourhood_gaussian( int x, int bmu, struct som_struct *s, int t);
  /** \} */

  /**\ingroup som_train
	\{ */
  void som_train_from_data( Som *s, double **X, int dim, int nsamples );
  /** \} */
  
#ifdef __cplusplus
}
#endif

#endif /* SOM_H */
