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
#include "distances.h"

#ifdef __cplusplus
extern "C" {
#endif

  /*-----------------------------------------------------------
	 - CLUSTERING -
	 ---------------------------------------------------------*/

  /** \brief Representation of a complete cluster-analysis.
	*/
  typedef struct{
	 int **clust; /**< indices for the trials in the cluster (Kxn)*/
	 int K;       /**< number of clusters */
	 int *n;      /**< number of trials in each of the K clusters */
  } Clusters;


  /** \brief Representation of a rooted binary tree to hold a dendrogram. 
		
		It is a terminal node if val>0 and left=right=NULL.
		Else it is an intermediate node that has at least one child!=NULL.
	*/
   struct dgram {
	  int val;       /**< content representing object val; if val<0, its an intermediate node */
	  int clustnum;  /**< number of cluster (index) */ 
	  double height; /**< proportional to between sub-cluster distance */
	  struct dgram *left;
	  struct dgram *right;
  };
  typedef struct dgram Dendrogram;

  /** \brief Linkage function for hierarchical cluster analysis.
	*/
  typedef double  (*LinkageFunction)   (const Array*,const Dendrogram*,const Dendrogram*);


  /* ---------------- CLUSTER ANALYSIS ------------------------------- */
  Clusters*    kmedoids(const Array *distmat, int K, OptArgList *optargs );
  Clusters*    kmedoids_repeat( const Array *distmat, int K, int repeat );

  Dendrogram*  agglomerative_clustering(const Array *distmat, LinkageFunction dist);

  /* ---------------- CLUSTER HANDLING ------------------------------- */
  void      cluster_free(Clusters *c);
  void      cluster_print(FILE *out, const Clusters *c);
  Clusters* cluster_init(int K, int maxN);
  void      cluster_copy(Clusters *dest, const Clusters *src);
  int       cluster_compare(const Clusters *c1, const Clusters *c2);

  double   cluster_within_scatter (const Array *distmat, const Clusters *c);
  double   cluster_between_scatter(const Array *distmat, const Clusters *c);

  /* ---------------- DENDROGRAM ------------------------------- */
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

#ifdef __cplusplus
}
#endif


#endif
