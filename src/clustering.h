/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke                                  *
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

/**\file clustering.h
 * \brief Clustering functions.
 */

#ifndef CLUSTERING_H
# define CLUSTERING_H
#include <time.h>
#include "mathadd.h"
#include "definitions.h"
#include "averaging.h"
#include "warping.h"
#include "distances.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\addtogroup clustering
	*\{
	*/ 
  Clusters*    kmedoids(const double **dist, int N, int K);
  Clusters*    kmedoids_repeat( const double **dist, int N, int K, int repeat );

  Dendrogram*  agglomerative_clustering(const double **d, int N, LinkageFunction dist);

  /**\addtogroup dendrogram
	*\{
	*/ 

  double dgram_dist_singlelinkage  (const double **d, int N, 
												const Dendrogram *c1, const Dendrogram *c2);
  double dgram_dist_completelinkage(const double **d, int N, 
												const Dendrogram *c1 ,const Dendrogram *c2);
  double dgram_dist_averagelinkage (const double **d, int N, 
												const Dendrogram *c1, const Dendrogram *c2);

  void         dgram_print( Dendrogram *t );
  void         dgram_print_node( Dendrogram *t );
  void         dgram_preorder( const Dendrogram *t, int *vals, int *n );
  Dendrogram*  dgram_init(int val, Dendrogram *left, Dendrogram *right);
  void         dgram_free(Dendrogram *t);
  Dendrogram*  dgram_get_deepest( Dendrogram *c );
  /** \} */


  void      free_cluster(Clusters *c);
  void      print_cluster(FILE *out, const Clusters *c);
  Clusters* init_cluster(int K, int maxN);
  void      copy_cluster(Clusters *dest, const Clusters *src);
  int       compare_clusters(const Clusters *c1, const Clusters *c2);

  /** \} */


  /**\addtogroup gap
	  Reference paper: Tibshirani, 2001
	  \{
  */
  
  GapStatistic* gapstat_init( GapStatistic *g, int K, int B );
  void          gapstat_free( GapStatistic *g );
  void          gapstat_print( FILE *out, GapStatistic *g );
  void          gapstat_calculate( GapStatistic *gap, double **X, int n, int p, 
											  VectorDistanceFunction distfunction, const double** D );
  
  double** gap_get_reference_distribution_simple( const double **X, int n, int p, double **Xr );
  double** gap_get_reference_distribution_svd   ( const double **X, int n, int p, double **Xr );
  double   get_within_scatter (const double **d, int N, const Clusters *c);
  double   get_between_scatter(const double **d, int N, const Clusters *c);
  /** \} */


#ifdef __cplusplus
}
#endif


#endif
