/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
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

#include "clustering.h"


/** standardize matrix d (NxN) by computing
	 \f$
	 \hat{d}_{jk} = \frac{1}{\max d_{jk}} d_{jk}
	 \f$
*/
void diffmatrix_standardize(double **d, int N){
  double maxv;
  int i,j;

  maxv = -DBL_MAX;
  for( i=0; i<N; i++ ){
	 for( j=0; j<N; j++ ){
		if(d[i][j]>maxv)
		  maxv = d[i][j];
	 }
  }
  for( i=0; i<N; i++ ){
	 for( j=0; j<N; j++ ){
		d[i][j] /= maxv;
	 }
  }
    
}

/** return the distance of two ERPs from timewarping the complete ERP 
	 (no time-marker-based warping).
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_tw_complete(EEGdata *s1, EEGdata *s2, int channel){
  double Djk;
  Djk = DTW_get_warpdistance( s1->d[channel], s1->n, s2->d[channel], s2->n, 1.0, 1.0);
  
  return Djk;
}

/** return the distance of two ERPs from timewarping the  ERP 
	 based on the time-markers
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_tw_markers(EEGdata *s1, EEGdata *s2, int channel){
  double Djk;
  Djk = DTW_get_warpdistance( s1->d[channel], s1->n, s2->d[channel], s2->n, 1.0, 1.0);

  return Djk;
}
/** return the distance of two ERPs from euclidean distances.
	 \f$
	 \Delta(s_i, s_j) = ||s_i-s_j||^2 = \sum_t d_t(s_i(t), s_j(t))
	 \f$
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_euclidean_pointwise(EEGdata  *s1, EEGdata *s2, int channel){
  int i;
  double dist;

  massert( s1->n!=s2->n, "n1 != n2, aborting\n" );
  dist = 0.0;
  for(i=0; i<s1->n; i++){
	 dist += SQR( (s1->d[channel][i])-(s2->d[channel][i]) );
  }
  dist = sqrt(dist); ///(double)s1->n;
  return dist;
}

/** Return a matrix of the differences between trials in 
	 the EEGdata_trials struct (NxN) for a given channel.
	 \param eeg 
	 \param dist - the distance function
	 \param channel - which electrode-channel? (index)
	 \param dm - enough memory to hold matrix, or ALLOC_IN_FCT
*/
double** eegtrials_diffmatrix_channel(EEGdata_trials *eeg, 
												  double(*dist)(EEGdata*,EEGdata*,int), 
												  int channel, double **dm){
  int i,j;

  if( dm==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 dm = matrix_init( eeg->ntrials, eeg->ntrials );
  }
  
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=i+1; j<eeg->ntrials; j++ ){
		dm[i][j] = dist( eeg->data[i], eeg->data[j], channel );
		/* dprintf("getting distance (%i x %i)=%f\n", i,j,dm[i][j]); */
		dm[j][i] = dm[i][j];
	 }
  }
  return dm;
}


/** do K-Medoids clustering on the distance-matrix.
	 \param dist
	 \param K 
	 \param clusters - function-allocated memory. contains
	                   K cluster-structs.
*/
Clusters* kmedoids(const double **dist, int N, int K){
  Clusters *C, *Cnew;
  int *medoids;
  int i, j, k, r, minidx, num_not_changed;
  double sum, minsum, mindist, curdist;

  /* memory alloc */
  medoids = (int*) malloc( K*sizeof(int) );
  C    = init_cluster(K, N);
  Cnew = init_cluster(K, N);

  /* initialization step, random init */
  srandom( (unsigned int)time( (time_t *)NULL ) );
  /*random() / (RAND_MAX / N + 1); -- value from 0,...,N-1 */
  for( i=0; i<N; i++ ){
	 r = (random() / (RAND_MAX / K+1));
	 C->clust[r][C->n[r]] = i;
	 (C->n[r])++;
  }

/*   dprintf("initialized cluster\n");  */
/*   print_cluster(C); */

  num_not_changed = 0;
  while(num_not_changed<1){ /* iterate until convergence */
	 for( i=0; i<K; i++ ){ /* minimize distance to one of trials */
		minsum = DBL_MAX;
		minidx = 0;
		for( j=0; j<C->n[i]; j++ ){
		  sum = 0;
		  for( k=0; k<C->n[i]; k++){
			 sum += dist[C->clust[i][j]][C->clust[i][k]];
		  }
		  
		  if( sum<minsum ){
			 /* fprintf(stderr, "sum=%f, minsum=%f\n", sum, minsum);*/
			 minsum = sum;
			 minidx = C->clust[i][j];
		  }
		}
		medoids[i] = minidx;
	 }

	 /* reset cluster */
	 for( i=0; i<Cnew->K; i++ ){
		Cnew->n[i] = 0;
	 }
	 
	 for( i=0; i<N; i++ ){ /* reassign clusters */
		mindist = DBL_MAX;
		minidx = 0;
		for( j=0; j<K; j++ ){ /* look for closest center */
		  curdist = dist[medoids[j]][i];
		  if( curdist<mindist ){
			 mindist=curdist;
			 minidx=j;
		  }
		}
		/*dprintf("Assigning trial %i to cluster %i\n", i, minidx);*/
		Cnew->clust[minidx][Cnew->n[minidx]] = i;
		Cnew->n[minidx] += 1;
	 }
/* 	 if(count>=100000){ */
/* 		fprintf(stderr, "count=%i\n", count); */
/* 		fprintf(stderr, "medoids = (%i, %i)\n", medoids[0], medoids[1]); */
/* 		fprintf(stderr, "C=\n"); */
/* 		print_cluster(C); */
/* 		fprintf(stderr, "Cnew=\n"); */
/* 		print_cluster(Cnew); */
/* 		count = 0; */
/* 	 } else { */
/* 		count++; */
/* 	 } */
	 if(!compare_clusters(C, Cnew))
		num_not_changed++;
	 else
		copy_cluster(C, Cnew);
  }

/*   print_cluster(C); */
/*   fprintf(stderr, "medoids = (%i, %i)\n", medoids[0], medoids[1]); */
  free(medoids);
  free_cluster(Cnew);
  return C;
}

/** Compare whether two cluster structs are "identical", meaning
	 that each cluster contains the same arguments (sequence
	 does not matter). Superficial tests are first carried out, 
	 then intense comparison.
*/
int compare_clusters(const Clusters *c1, const Clusters *c2){
  int i, j, k, l, found, el;
  int c2i, prev_c2i;

  if( c1->K != c2->K ) return 1;
  /* clusters with equal number of trials? */
  for( i=0; i<c1->K; i++ ){
	 found = 0;
	 for( j=0; j<c2->K; j++ ){
		if(c1->n[i] == c2->n[j]){
		  found = 1;
		  break;
		}
	 }
	 if(!found)		
		return 1;
  }

  /* deep comparison */
  for( i=0; i<c1->K; i++ ){ /* each cluster from c1 */ 
	 dprintf("i = %i\n", i);
	 prev_c2i = c1->K+1;
	 c2i = 0;
	 for( j=0; j<c1->n[i]; j++ ){ /* each element from c1 */
		el = c1->clust[i][j];
		/*		dprintf("j = %i\n", j);*/
		for( k=0; k<c2->K; k++ ){
		  for( l=0; l<c2->n[k]; l++ ){
			 if( el == c2->clust[k][l] ){
				c2i = k;
				if(prev_c2i==c1->K+1)
				  prev_c2i = k; /* first run */
			 } 
		  }
		}
		if( c2->n[c2i] != c1->n[i] ){
		  /*		  dprintf("ab: c2->n[%i]=%i, c1->n[%i]=%i\n", c2i, c2->n[c2i], i, c1->n[i]);*/
		  return 1;
		}
		if(prev_c2i!=c2i){
		  /*		  dprintf("ab: prev_c2i=%i, c2i=%i\n", prev_c2i, c2i); */
		  return 1;
		}
	 }
  }
  return 0;
}

/** assumes that dest is already allocated */
void copy_cluster(Clusters *dest, const Clusters *src){
  int i, j;

  dest->K = src->K;
  for( i=0; i<src->K; i++ ){
	 dest->n[i] = src->n[i];
	 for( j=0; j<src->n[i]; j++ ){
		dest->clust[i][j] = src->clust[i][j];
	 }
  }
  
}

void free_cluster(Clusters *c){
  int i;
  free( c->n );
  for( i=0; i<c->K; i++)
	 free(c->clust[i]);
  free(c->clust);
  free(c);
}

Clusters* init_cluster(int K, int maxN){
  Clusters *c;
  int i;

  c = (Clusters*) malloc( sizeof(Clusters) );
  c->K=K;
  c->n = (int*) malloc( K*sizeof(int) );
  c->clust = (int**) malloc( K*sizeof(int*) );
  for( i=0; i<K; i++ ){
	 c->clust[i] = (int*) malloc( maxN*sizeof(int) );
	 c->n[i] = 0;
  }

  return c;
}

void print_cluster(const Clusters *c){
  FILE *o;
  int i, j;
  o = stderr;
  fprintf(o, "Cluster with %i partitions\n", c->K);
  for( i=0; i<c->K; i++ ){
	 fprintf(o, " C[%i]=\{ ", i+1);
	 for( j=0; j<c->n[i]; j++ ){
		fprintf( o, "%i ", c->clust[i][j] );
	 }
	 fprintf( o, "}\n" );
  }
  fprintf( o, "\n" );
}


int      gap_get_K(const double *gapstat, int k);
double*  gap_get_gapstat(const double *Wk, const double **Wkref, int B, int k, double *gapstat);
double** gap_get_reference_distribution(const double **d, int N, int n, double **ref);

/** compute \f$ W_k = \sum_{r=1}^{k}\frac{1}{2n_r}D_r\f$
	 where
	 \f$
	 D_r = \sum_{i,i'\in C_r}d_{ii'}
	 \f$
	 for a given cluster assignment.
	 \param d - difference matrix
	 \param N - number of entries in d (NxN)
	 \param c - cluster-assigment (Clusters - struct)
 */
double   gap_get_within_scatter(const double **d, int N, const Clusters *c){
  double *sums;
  double W;
  int i, j, k;
  print_cluster(c);
  sums = (double*) malloc( c->K*sizeof(double) );
  W = 0;
  for( i=0; i<c->K; i++ ){
	 sums[i] = 0;
	 for( j=0; j<(c->n[i])-1; j++ ){
		for( k=j+1; k<c->n[i]; k++){	 
	/* 	  dprintf("sums[i]=%f\n", sums[i]); */
/* 		  dprintf("i,j,k=(%i,%i,%i), c->clust[i][j]=%i, c->clust[i][k]=%i, d[]=%f\n", */
/* 					 i,j,k, c->clust[i][j], c->clust[i][k], d[c->clust[i][j]][c->clust[i][k]]); */
		  sums[i] += d[c->clust[i][j]][c->clust[i][k]];
		}
	 }
	 sums[i] /= (double)2*c->n[i];
	 dprintf("sums[i]=%f\n", sums[i]);
	 W += sums[i];
  }

  free(sums);

  return W;
}
double*  gap_get_within_scatter_distribution(const double **d, int N, int k, const Clusters **c, double *Wk);


/** Agglomerative clustering.
	 \param d distance matrix for the N objects
	 \param N number of objects
	 \param dist function giving the between sub-cluster distance; 
	             defined are dgram_dist_singlelinkage(), dgram_dist_completelinkage(),
					 dgram_dist_averagelinkage()
 */
Dendrogram* agglomerative_clustering(const double **d, int N, 
												  double(*dist)(const double**,int,const Dendrogram*,const Dendrogram*)){
  Dendrogram **nodes, *tmp; /* at the lowest level, we have N nodes */
  double min_d, cur_d;
  int min_i=0, min_j=0;
  int num_nodes;
  int i, j;
 
  nodes = (Dendrogram**) malloc( N*sizeof( Dendrogram* ) );

  for( i=0; i<N; i++ ){ /* initialize terminal nodes */
	 nodes[i] = dgram_init( i, NULL, NULL );
  }
  num_nodes = N;

  while( num_nodes > 1 ){
	 dprintf("nmodes = %i\n", num_nodes );
	 min_d = DBL_MAX;
	 for( i=0; i<num_nodes-1; i++ ){
		for( j=i+1; j<num_nodes; j++ ){ 
		  cur_d = dist( d, N, nodes[i], nodes[j] );
		  if( cur_d<=min_d ){
			 min_d=cur_d;
			 min_i=i;
			 min_j=j;
		  }
		}
	 }
	 
	 /* we know min_i and min_j, now */
	 tmp = dgram_init( -1, nodes[min_i], nodes[min_j] ); /* link the two nodes */
	 nodes[min_i] = tmp; /* remove min_i */
	 tmp = nodes[num_nodes-1]; /* swap with last entry in current list */
	 nodes[min_j] = tmp;
	 num_nodes--;
	 dprintf("Linking (%i,%i)\n", min_i, min_j );
  }
  tmp = nodes[0]; /* this is the root-node now */
  free( nodes );
  return tmp;
}

/** Allocate memory for a single Dendrogram Node and set its value to val.
	 left and right can be set to NULL.
*/
Dendrogram* dgram_init(int val, Dendrogram *left, Dendrogram *right ){
  Dendrogram *c;
  c = (Dendrogram*) malloc( sizeof(Dendrogram) );
  c->val = val;
  c->height=0.0;
  c->left = left;
  c->right= right;
  return c;
}


/** frees the complete Dendrogram referred to by t (all children)
 */
void dgram_free(Dendrogram *t){
  if( t->right==NULL && t->left==NULL ){ /* base case */
	 free( t );
	 return;
  } else { /* recurse into sub-trees */
	 if( t->right!=NULL ){
		dgram_free( t->right );
	 } 
	 if( t->left!=NULL ){
		dgram_free( t->left );
	 }
  }
}
/** \cond PRIVATE */
void _dgram_preorder_recursive( const Dendrogram *t, int *vals, int *n ){
  if( t->left==NULL && t->right==NULL ){
	 if( t->val<0 ){
		errprintf("t->val<0 (%i) even though it's a terminal node\n", t->val);
		return;
	 }
	 vals[(*n)++] = t->val;
	 return;
  } 
  if( t->left!=NULL ){
	 //	 dprintf("walking  left subtree\n");
	 _dgram_preorder_recursive( t->left, vals, n );
  }
  if( t->right!=NULL ){
	 _dgram_preorder_recursive( t->right, vals, n );
  }
  return;
}
/** \endcond */

/** preorder traversal of t. Store all non-negative elements in val and return 
	 the number of these elements in n[0];
	 \param t the tree
	 \param vals output
	 \param n output (number of elements in vals)
*/
void         dgram_preorder( const Dendrogram *t, int *vals, int *n ){
  int m;

  m = 0;
  //  dprintf(" starting recursive walk, m=%i\n", m);
  //  dgram_print_node(t);
  _dgram_preorder_recursive( t, vals, &m );

  /* fprintf(stderr, "vals=["); */
  /* for( i=0; i<m; i++ ){ */
  /* 	 fprintf(stderr, "%i ", vals[i]); */
  /* } */
  /* fprintf(stderr, "]\n"); */
  //  dprintf(" done with recursive walk, m=%i\n", m);
  *n = m;

  return;
}
/** single linkage clustering:
	 \f[
	 d_{SL}(G, H) = \min_{i\in G, j\in H} d_{ij}
	 \f]
 */
double dgram_dist_singlelinkage  (const double **d, int N, const Dendrogram *c1, const Dendrogram *c2){
  double dist;
  int *el1, *el2;
  int n1, n2;
  int i, j;

  dist = DBL_MAX;
  
  el1 = (int*) malloc( N*sizeof( int ) ); /* N suffices */
  el2 = (int*) malloc( N*sizeof( int ) ); 

  dgram_preorder( c1, el1, &n1 );
  dgram_preorder( c2, el2, &n2 );

  for( i=0; i<n1; i++ ){
	 for( j=0; j<n2; j++ ){
		if( el1[i]<0 || el1[i]>=N || el2[j]<0 || el2[j]>=N ){
		  errprintf("something's wrong, el[i] not in range 0-%i\n", N);
		  return -1;
		}
		if( d[el1[i]][el2[j]]<dist ){
		  dist = d[el1[i]][el2[j]];
		}
	 }
  }
  
  //dprintf("SL-distance = %f\n", dist);

  free(el1);
  free(el2);
  return dist;
}

/** complete linkage clustering:
	 \f[
	 d_{CL}(G, H) = \max_{i\in G, j\in H} d_{ij}
	 \f]
 */
double dgram_dist_completelinkage(const double **d, int N, const Dendrogram *c1, const Dendrogram *c2){
  double dist;
  int *el1, *el2;
  int n1, n2;
  int i, j;

  dist = DBL_MIN;
  
  el1 = (int*) malloc( N*sizeof( int ) ); /* N suffices */
  el2 = (int*) malloc( N*sizeof( int ) ); 

  dgram_preorder( c1, el1, &n1 );
  dgram_preorder( c2, el2, &n2 );

  for( i=0; i<n1; i++ ){
	 for( j=0; j<n2; j++ ){
		if( el1[i]<0 || el1[i]>=N || el2[j]<0 || el2[j]>=N ){
		  errprintf("something's wrong, el[i] not in range 0-%i\n", N);
		  return -1;
		}
		if( d[el1[i]][el2[j]]>dist ){
		  dist = d[el1[i]][el2[j]];
		}
	 }
  }
  
  dprintf("SL-distance = %f\n", dist);

  free(el1);
  free(el2);
  return dist;
}


/** average linkage clustering:
	 \f[
	 d_{AL}(G, H) = \frac{1}{N_G N_H} \sum_{i\in G} \sum_{j\in H} d_{ij}
	 \f]
*/
double dgram_dist_averagelinkage (const double **d, int N, const Dendrogram *c1, const Dendrogram *c2){
  errprintf(" NOT IMPLEMENTED YET!!\n");
  return -1;
}

/** \cond PRIVATE */
#define END_NODE 0
#define INTERMEDIATE_NODE 1
#define NULL_NODE 2
int _dgram_get_deepest_recursive( Dendrogram *c, Dendrogram **candidate, int *depth, int curdepth ){
  int status_left=NULL_NODE;
  int status_right=NULL_NODE;

  if( c->left==NULL && c->right==NULL ){
	 return END_NODE;
  }

  if( c->left!=NULL ){
	 status_left  = _dgram_get_deepest_recursive( c->left,  candidate, depth, curdepth+1 );  
  }
  if( c->right!=NULL ){
	 status_right = _dgram_get_deepest_recursive( c->right, candidate, depth, curdepth+1 );  
  }

  if( status_left==END_NODE && status_right==END_NODE ){

	 if( curdepth>=*depth ){
		*candidate = c;
		*depth = curdepth;
	 }
	 dprintf("END_NODE visited, cd=%i, d=%i, c=%p, candidate=%p\n", curdepth, *depth, c, candidate);
  }

  return INTERMEDIATE_NODE;
}
/** \endcond */

/** \return a pointer to the deepest non-leaf node in the Dgram c 
	 (left and right tree are non-NULL)
 */
Dendrogram* dgram_get_deepest( Dendrogram *c ){
  Dendrogram **d, *tmp;
  int depth;
  d = (Dendrogram **) malloc( sizeof( Dendrogram* ) );

  depth = 0;
  _dgram_get_deepest_recursive( c, d, &depth, 0 );
  dprintf("deepest=%p, depth=%i\n", *d, depth);
  tmp = *d;
  free(d);
  return tmp;
}


/** \cond PRIVATE */
void _dgram_print_recursive( Dendrogram *t, int level ){
  int i;

  fprintf( stdout, "'%p(%4i)' - ", t, t->val );
  if( t->left ){ /* there is a left tree, print it */
	 _dgram_print_recursive( t->left, level+1 );
  }  
  /* we need to add spaces */
  fprintf( stdout, "\n" );
  for( i=0; i<level; i++)
	 fprintf( stdout, "           " );

  if( t->right ){
	 _dgram_print_recursive( t->right, level+1 );
  }
}
/** \endcond */

void         dgram_print( Dendrogram *t ){
  _dgram_print_recursive( t, 0 );
  fprintf( stdout, "\n");
}


void dgram_print_node( Dendrogram *t ){
  fprintf( stdout, "t=%p, val=%i, height=%f, left=%p, right=%p\n", t, t->val, t->height, t->left, t->right );
}
#undef END_NODE 
#undef INTERMEDIATE_NODE 
#undef NULL_NODE
