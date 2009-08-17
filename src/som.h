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

#include "mathadd.h"
#include "definitions.h"


#ifdef __cplusplus
extern "C" {
#endif
  
  /** This is: distance = f( location of node i, location of node j, the model );
	*/
  typedef double(*NeighbourhoodFunction)(int,int);

  double neighbourhood_gaussian( int x, int m );

  /**     giving the dimension/structure of the SOM:
   - 1D_LINEAR - 1-dimensional line, like 1<->2<->3<->...<->n
   - \todo 2D_GRID - grid-like structure such that each node has 8 neighbours
   - \todo 2D_HEXAGONAL
   - \todo CUSTOM - you can define a custom connectivit matrix, where each entry [i,j]
              gives the distance between node i and j
  */
  typedef enum {
    ONED_LINEAR,
    TWOD_GRID,
    TWOD_HEXAGONAL,
    CUSTOM
  } SOMConnectivityType;

  /** Struct for 1-dimensional SOM with trivial topology 
		(neighbouring nodes).
	*/
  typedef struct{
	 double  **m; /**< code-book vectors */
	 int     dimension; /**< dimension of code-book vector (corresponds to data-dimensionality) */
	 int     n; /**< how many code-book vectors? */
	 VectorDistanceFunction *d; /**< distance measure between data and codebook-vectors */
	 
	 int     nruns; /**< number of runs to convergence */
	 int     initial_runs; /**< number of runs with large flexibility (ordering phase) after
									  which it is more restricted; e.g. 0.1*nruns */
	 NeighbourhoodFunction h; /**< neighbourhood function on the line */
	 SOMConnectivityType connectivity; /**< giving the dimension/structure of the SOM */
	 double **custom_connectivity;
  } Som;


  Som* som_init( int dimension, int n );
  void som_free( Som *s );

#ifdef __cplusplus
}
#endif

#endif /* SOM_H */
