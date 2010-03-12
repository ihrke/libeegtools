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

/**\file nnsearch.h
 \brief Fast Nearest-Neighbour searching.

see 

Merkwirth et al. Fast nearest-neighbor searching for nonlinear signal
processing. Physical review E, Statistical physics, plasmas, fluids,
and related interdisciplinary topics (2000) vol. 62 (2 Pt A)
pp. 2089-97
	
 */
#ifndef NNSEARCH_H
# define NNSEARCH_H
#include "mathadd.h"
#include "distances.h"
#include "definitions.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\cond PRIVATE */
  TreeNode*   tnode_init();
  bool tnode_isleaf( TreeNode *C );
  void build_tree_recursive( TreeNode *C, double **D, int N, int *A, int maxel );
  SearchTree* searchtree_init( int n );

  /**\endcond */


  SearchTree* nn_prepare( const double **X, int m, int N, OptArgList *optargs );
  void        nn_search_k( const SearchTree *S, const double *x, int k, 
									int *nn_idx, double *nn_dist );
  void        nn_searck_k_slow( const double **X, int N, int m, 
										  const double *x, int k, int *nn_idx,
										  double *nn_dist, OptArgList *optargs );

#ifdef __cplusplus
}
#endif


#endif
