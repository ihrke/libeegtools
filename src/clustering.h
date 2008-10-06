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
 * Clustering functions.
 * Groups:
 \defgroup clustering Clustering functions
 \{
    \defgroup distances Distance measures
    \defgroup gap Gap-statistic
	 \defgroup dendrogram Dendrogram (Hierarchical Clustering)
 \}
 */
#ifndef CUSTERING_H
# define CLUSTERING_H
#include <time.h>
#include "mathadd.h"
#include "definitions.h"
#include "averaging.h"

#ifdef __cplusplus
extern "C" {
#endif

  /**\addtogroup clustering
	*\{
	*/ 

  /**\ingroup clustering */
  typedef struct{
	 int **clust; /** indices for the trials in the cluster (Kxn)*/
	 int K;  /** number of clusters */
	 int *n; /** number of trials in each of the K clusters */
  } Clusters;


  /**\ingroup clustering
	  This struct is the representation of a rooted binary tree to
	  hold a dendrogram. 
	  It is a terminal node if val>0 and left=right=NULL.
	  Else it is an intermediate node that has at least one child!=NULL.
	*/
   struct dgram {
	 int val; /** content representing object val; if val<0, its an intermediate node */
	 double height; /** proportional to between sub-cluster distance */
	 struct dgram *left;
	 struct dgram *right;
  };
  typedef struct dgram Dendrogram;

  double** diffmatrix(ModelData *m, double **dm);
  double** eegtrials_diffmatrix_channel(EEGdata_trials *eeg, 
													 double(*dist)(EEGdata*,EEGdata*,int), 
													 int channel);
  void     diffmatrix_standardize(double **d, int N);


  Clusters*    kmedoids(const double **dist, int N, int K);


  /**\addtogroup dendrogram
	*\{
	*/ 

  Dendrogram*  agglomerative_clustering(const double **d, int N, 
													 double(*dist)(double**,int,const Dendrogram*,const Dendrogram*));
  double dgram_dist_singlelinkage  (double **d, int N, const Dendrogram *c1, const Dendrogram *c2);
  double dgram_dist_completelinkage(double **d, int N, const Dendrogram *c1 ,const Dendrogram *c2);
  double dgram_dist_averagelinkage (double **d, int N, const Dendrogram *c1, const Dendrogram *c2);

  void         dgram_print( Dendrogram *t );
  void         dgram_print_node( Dendrogram *t );
  void         dgram_preorder( Dendrogram *t, int *vals, int *n );
  Dendrogram*  dgram_init(int val, Dendrogram *left, Dendrogram *right);
  void         dgram_free(Dendrogram *t);
  /** \} */


  void      free_cluster(Clusters *c);
  void      print_cluster(Clusters *c);
  Clusters* init_cluster(int K, int maxN);
  void      copy_cluster(Clusters *dest, const Clusters *src);


  /** \} */

  /**\addtogroup distances 
	  \ingroup clustering
	  \{*/
  double clusterdist_euclidean_pointwise(EEGdata  *s1, EEGdata *s2, int channel);
  double clusterdist_tw_complete(EEGdata *s1, EEGdata *s2, int channel);
  double clusterdist_tw_markers(EEGdata *s1, EEGdata *s2, int channel);
  /** \} */


  /**\addtogroup gap
	  \ingroup clustering
	  \{
	  Reference paper: Tibshirani, 2001
  */
  int      gap_get_K(const double *gapstat, int k);
  double*  gap_get_gapstat(const double *Wk, const double **Wkref, int B, int k, double *gapstat);
  double** gap_get_reference_distribution(const double **d, int N, int n, double **ref);
  double   gap_get_within_scatter(const double **d, int N, const Clusters *c);
  double*  gap_get_within_scatter_distribution(const double **d, int N, int k, const Clusters **c, double *Wk);
  /** \} */


#ifdef __cplusplus
}
#endif


#endif
