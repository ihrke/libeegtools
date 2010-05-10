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
 * \brief \ref status_unstable Clustering functions.
 */

#ifndef CLUSTERING_H
# define CLUSTERING_H
#include <time.h>
#include "mathadd.h"
#include "definitions.h"
#include "array.h"
#include "averaging.h"
#include "warping.h"
#include "distances.h"

#ifdef __cplusplus
extern "C" {
#endif

  /** \brief Linkage function for hierarchical cluster analysis.
	*/
  typedef double  (*LinkageFunction)   (const Array*,const Dendrogram*,const Dendrogram*);

  /**\addtogroup clustering
	*\{
	*/ 
  Clusters*    kmedoids(const Array *distmat, int K, OptArgList *optargs );
  Clusters*    kmedoids_repeat( const Array *distmat, int K, int repeat );

  Dendrogram*  agglomerative_clustering(const Array *distmat, LinkageFunction dist);

  /**\addtogroup dendrogram
	*\{
	*/ 

  double dgram_dist_singlelinkage  (const Array *d, const Dendrogram *c1, const Dendrogram *c2);
  double dgram_dist_completelinkage(const Array *d, const Dendrogram *c1, const Dendrogram *c2);

  void         dgram_print( Dendrogram *t );
  void         dgram_print_node( Dendrogram *t );
  void         dgram_preorder( const Dendrogram *t, int *vals, int *n );
  Dendrogram*  dgram_init(int val, Dendrogram *left, Dendrogram *right);
  void         dgram_free(Dendrogram *t);
  Dendrogram*  dgram_get_deepest( Dendrogram *c ); 
  int          dgram_num_leaves( const Dendrogram *t );
  Array*       dgram_to_matlab( const Dendrogram *dgram );
  /** \} */


  void      cluster_free(Clusters *c);
  void      cluster_print(FILE *out, const Clusters *c);
  Clusters* cluster_init(int K, int maxN);
  void      cluster_copy(Clusters *dest, const Clusters *src);
  int       cluster_compare(const Clusters *c1, const Clusters *c2);

  double   cluster_within_scatter (const Array *distmat, const Clusters *c);
  double   cluster_between_scatter(const Array *distmat, const Clusters *c);

  /** \} */


  /**\addtogroup gap
	  Reference paper: Tibshirani, 2001
	  
	  \todo rewrite with Array struct
	  \{
  */
  
#if 0
  GapStatistic* gapstat_init( GapStatistic *g, int K, int B );
  void          gapstat_free( GapStatistic *g );
  void          gapstat_print( FILE *out, GapStatistic *g );
  void          gapstat_calculate( GapStatistic *gap, double **X, int n, int p, 
											  VectorDistanceFunction distfunction, const double** D );
  
  double** gap_get_reference_distribution_simple( const double **X, int n, int p, double **Xr );
  double** gap_get_reference_distribution_svd   ( const double **X, int n, int p, double **Xr );

  int      eeg_best_num_clusters_gapstat( const EEG *eeg, VectorDistanceFunction distfunction,
														OptArgList *optargs );
#endif
  /** \} */


#ifdef __cplusplus
}
#endif


#endif
