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



/** run the kmedoids function a couple of times and pick the best
	 in terms of within scatter.
	 \param dist ance matrix
	 \param N number of observations
	 \param K number of medoids to compute
	 \param repeat number of repetitions
	 \return clusters - function-allocated memory. contains
	                    K cluster-structs.
 */
Clusters* kmedoids_repeat( const double **dist, int N, int K, int repeat ){
  Clusters *C;
  Clusters *Cnew;
  int i;
  double ratio, tmp;

  C = init_cluster( K, N );
  ratio=DBL_MAX;
  for( i=0; i<repeat; i++ ){
	 Cnew = kmedoids( dist, N, K );
	 tmp =  get_within_scatter ( dist, N, Cnew );
	 tmp /= get_between_scatter( dist, N, Cnew );

	 if( tmp<ratio ){
		dprintf("Run %i: %f\n", i, tmp);
		ratio=tmp; 
		copy_cluster( C, Cnew );
	 }
	 free_cluster( Cnew );
  }

  return C;
}

/** do K-Medoids clustering on the distance-matrix.
	 \param dist ance matrix
	 \param N number of observations
	 \param K number of medoids to compute
	 \return clusters - function-allocated memory. contains
	                   K cluster-structs.
*/
Clusters* kmedoids(const double **dist, int N, int K){
  Clusters *C, *Cnew;
  int *medoids;
  int i, j, k, r, minidx, num_not_changed;
  double sum, minsum, mindist, curdist;
  const gsl_rng_type * T;
  gsl_rng * randgen;

  if( K>N ){
	 errprintf("K>=N: %i>=%i\n", K, N );
	 return NULL;
  }

  /* memory alloc */
  medoids = (int*) malloc( K*sizeof(int) );
  C    = init_cluster(K, N);
  Cnew = init_cluster(K, N);

  /* initialization step, random init */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  randgen = gsl_rng_alloc (T);

  int *permut;
  permut = (int*)malloc( K*sizeof(int) );
  for( i=0; i<K; i++ )
	 permut[i]=i;
  vector_shuffle_int( permut, K ); /* random permut */
  for( i=0; i<N; i++ ){
	 if( i<K ){ /* first each partition gets a guy */
		r = permut[i];
	 } else {
		r = (random() / (RAND_MAX / K+1));
	 }
	 C->clust[r][C->n[r]] = i;
	 (C->n[r])++;
  }
  free( permut );

  dprintf("initialized cluster\n");  
  print_cluster(stderr, C); 

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

	 if(!compare_clusters(C, Cnew))
		num_not_changed++;
	 else
		copy_cluster(C, Cnew);
  }

  free(medoids);
  free_cluster(Cnew);
  gsl_rng_free (randgen);

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
	 /* dprintf("i = %i\n", i); */
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

void print_cluster(FILE *o,const Clusters *c){
  int i, j;

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


/** initialize GapStatistic struct.
	 \param g - a pointer to allocated mem, or NULL -> own memory is allocated
	 \param K - max num clusters
	 \param B - num rep ref dist
 */
GapStatistic* gapstat_init( GapStatistic *g, int K, int B ){
  int i;

  if( g==ALLOC_IN_FCT ){
	 g = (GapStatistic*)malloc( sizeof(GapStatistic) );
  }
  g->K = K;
  g->B = B;
  g->khat = -1;
  dprintf("allocating memory\n");
  g->gapdistr = (double*)malloc( K*sizeof(double) );
  g->sk = (double*)malloc( K*sizeof(double) );
  g->Wk= (double*)malloc( K*sizeof(double) );
  g->Wkref= (double**) malloc( B*sizeof(double*) );
  for( i=0; i<B; i++ )
	 g->Wkref[i] = (double*) malloc( K*sizeof(double) );
  g->progress=NULL;

  return g;
}

void gapstat_free( GapStatistic *g ){
  int i;
  free( g->Wk );
  for( i=0; i<g->B; i++ )
	 free( g->Wkref[i] );
  free( g->Wkref );
  free( g->gapdistr );
  free( g->sk );
  free( g );
}

/** calculate the gap-statistic. the struct contains all important information.
	 \todo fix this, there appears to be an NaN error
	 \param gap
	 \param X is nxp where n is number of observations and p num features (e.g. 
	          trials x voltage
	 \param D is the distance matrix for X (nxn)
 */
void gapstat_calculate( GapStatistic *gap, double **X, int n, int p, 
								VectorDistanceFunction distfunction, const double** D ){
  int i, k, b;
  Clusters *C;
  double **Xr; /* reference distribution */
  double **Dr; /* distance in ref dist */
  double l_bar; /* (1/B sum_b( log(Wkref[b]) )) */

  /* init */
  Xr = (double**) malloc( n*sizeof( double* ) );
  Dr = (double**) malloc( n*sizeof( double* ) );
  for( i=0; i<n; i++ ){
	 Xr[i] = (double*) malloc( p*sizeof( double ) );
	 Dr[i] = (double*) malloc( n*sizeof( double ) );
  }

  if( gap->progress ){
	 gap->progress( PROGRESSBAR_INIT, (gap->B+1)*gap->K );
  }

  /* calculation */
  for( k=1; k<=gap->K; k++ ){ /* data within-scatter */
	 if( gap->progress )
		gap->progress( PROGRESSBAR_CONTINUE_LONG, k );
	 C = kmedoids_repeat( D, n, k, 50 );
	 //	 print_cluster( stderr, C );
	 gap->Wk[k-1] = get_within_scatter( D, n, C );
	 free_cluster( C );
  }

  for( b=0; b<gap->B; b++ ){ /* monte Carlo for ref-data*/
	 for( k=1; k<=gap->K; k++ ){
		if( gap->progress ){
		  gap->progress( PROGRESSBAR_CONTINUE_LONG, (b+1)*gap->K+k );
		}
		
		Xr = gap_get_reference_distribution_simple( (const double **) X, n, p, Xr );
		Dr = vectordist_distmatrix( distfunction, (const double**)Xr, n, p, Dr, NULL, NULL );
		C = kmedoids_repeat( (const double**)Dr, n, k, 50 );
		//		print_cluster(stderr, C );
		gap->Wkref[b][k-1] = get_within_scatter( (const double**)Dr, n, C );
		free_cluster( C );
	 }
  }
  
  for( k=1; k<=gap->K; k++ ){ /* gap distr */
	 gap->gapdistr[k-1]=0;
	 for( b=0; b<gap->B; b++ ){
		gap->gapdistr[k-1] += log( gap->Wkref[b][k-1] );
	 }
	 l_bar = ((gap->gapdistr[k-1])/(double)gap->B);
	 gap->gapdistr[k-1] = l_bar - log( gap->Wk[k-1] );

	 /* need std and sk */
	 gap->sk[k-1] = 0;
	 for( b=0; b<gap->B; b++ ){
		gap->sk[k-1] += SQR( log( gap->Wkref[b][k-1] ) - l_bar );
	 }
	 gap->sk[k-1] = sqrt( gap->sk[k-1]/(double)gap->B );
	 gap->sk[k-1] = gap->sk[k-1] * sqrt( 1.0 + 1.0/(double)gap->B );
  }
  
  /* find best k */
  for( k=0; k<=gap->K-1; k++ ){
	 if( gap->gapdistr[k] >= ((gap->gapdistr[k+1])-(gap->sk[k+1])) ){
		gap->khat = k+1;
		break;
	 }
  }
  if( gap->progress )
	 gap->progress( PROGRESSBAR_FINISH, 0 );

  /* free */
  for( i=0; i<n; i++ ){
	 free(Xr[i]);
	 free(Dr[i]);
  }
  free( Xr );
  free( Dr );
}

void gapstat_print( FILE *out, GapStatistic *g ){
  int i, j;
  double mean; 

  fprintf(stderr, "GapStatistic:\n"
			 " K=%i\n"
			 " B=%i\n"
			 " khat=%i\n", g->K, g->B, g->khat );
  fprintf(stderr, " Gap       = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.2f) ", i+1, g->gapdistr[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " sk        = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.4f) ", i+1, g->sk[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " Wk        = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.4f) ", i+1, g->Wk[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " <Wkref>_B = [ ");
  for(i=0; i<g->K; i++ ){
	 mean = 0;
	 for(j=0; j<g->B; j++){
		mean += g->Wkref[j][i];
	 }
	 mean /= (double)g->B;
	 fprintf(stderr, "(%i, %.4f) ", i+1, mean);
  }
  fprintf(stderr, "]\n");

}

/** Compute a reference distribution for data X.
	 \param X original data
	 \param n,p dimensions of X
	 \param Xr n x p matrix with reference data
 */
double** gap_get_reference_distribution_simple( const double **X, int n, int p, double **Xr ){
  int i, j;
  double minf, maxf;

  for( i=0; i<p; i++ ){ /* features */
	 minf=DBL_MAX;
	 maxf=DBL_MIN;
	 for( j=0; j<n; j++ ){ 		  /* find min/max for feature i */
		if( X[j][i]<minf )
		  minf = X[j][i];
		if( X[j][i]>maxf )
		  maxf = X[j][i];
	 }
	 
	 for( j=0; j<n; j++ ){ /* uniform random values in [minf,maxf] */
		Xr[j][i] = (drand48()*maxf)+minf;
		if( isnan(Xr[j][i]) ){
		  errprintf(" Xr[%i][%i] is nan\n", j, i);
		}
	 }
  }

  return Xr;
}

/** Compute a reference distribution for data X. 
	 CAUTION, not tested.
	 \param X original data
	 \param n,p dimensions of X
	 \param Xr n x p matrix with reference data
	 \todo does not yet work for n<p, since libgsl does not support this for SVD
 */
double** gap_get_reference_distribution_svd( const double **X, int n, int p, double **Xr ){
  int i, j;
  double minf, maxf;
  gsl_matrix *gX, *gXd, *V;
  gsl_vector *work, *S;

  /* convert to gsl */
  gX = gsl_matrix_alloc( n, p );
  gXd= gsl_matrix_alloc( n, p );
  V  = gsl_matrix_alloc( n, p );
  S  = gsl_vector_alloc( p );
  work=gsl_vector_alloc( p );

  for( i=0; i<n; i++ ){
	 for( j=0; j<p; j++ ){
		gsl_matrix_set( gX, i, j, X[i][j] );
	 }
  }
  gsl_matrix_memcpy( gXd, gX );

  gsl_linalg_SV_decomp( gXd, V, S, work);

  gsl_blas_dgemm(CblasNoTrans, /* matrix mult */
					  CblasNoTrans,
					  1.0, gX, V, 0.0, gXd);

  for( i=0; i<p; i++ ){ /* features */
	 minf=DBL_MAX;
	 maxf=DBL_MIN;
	 for( j=0; j<n; j++ ){ 		  /* find min/max for feature i */
		if( gsl_matrix_get( gXd, j, i )<minf )
		  minf = gsl_matrix_get( gXd, j, i );
		if( gsl_matrix_get( gXd, j, i )>maxf )
		  maxf = gsl_matrix_get( gXd, j, i );
	 }
	 
	 for( j=0; j<n; j++ ){ /* uniform random values in [minf,maxf] */
		gsl_matrix_set( gX, j, i, ((drand48()*maxf)+minf) );
	 }
  }

  gsl_matrix_transpose_memcpy( gXd, V );
  gsl_blas_dgemm(CblasNoTrans,
					  CblasNoTrans,
					  1.0, gX, gXd, 0.0, V);

  for( i=0; i<n; i++ ){
	 for( j=0; j<p; j++ ){
		Xr[j][i] = gsl_matrix_get( V, i, j );
	 }
  }

  gsl_matrix_free( gX );
  gsl_matrix_free( gXd );
  gsl_matrix_free( V );
  gsl_vector_free( work );
  gsl_vector_free( S );

  return Xr;
}


/** compute \f$ W_k = \sum_{r=1}^{k}\frac{1}{2n_r}D_r\f$
	 where
	 \f$
	 D_r = \sum_{i,i'\in C_r}d_{ii'}
	 \f$
	 for a given cluster assignment.
	 \param d - difference matrix (NxN)
	 \param N - number of entries in d (NxN)
	 \param c - cluster-assigment (Clusters - struct)
 */
double   get_within_scatter(const double **d, int N, const Clusters *c){
  double *sums;
  double W;
  int i, j, k;
  //  print_cluster( stderr, c);
  sums = (double*) malloc( c->K*sizeof(double) );
  W = 0;
  for( i=0; i<c->K; i++ ){ /* cluster loop */
	 sums[i] = 0;
	 for( j=0; j<(c->n[i])-1; j++ ){
		for( k=j+1; k<c->n[i]; k++){	 
		  sums[i] += d[c->clust[i][j]][c->clust[i][k]];
		}
	 }
	 sums[i] /= (double)2*c->n[i];
	 /* dprintf("sums[i]=%f\n", sums[i]); */
	 W += sums[i];
  }

  free(sums);

  return W;
}

/** compute the between-scatter between clusters 
	 (sum of distances of means of clusters)
	 \param d distance matrix (NxN)
	 \param c the partioning of d
 */
double get_between_scatter( const double **d, int N, const Clusters *c ){
  double W;
  int k1, k2, i, j;
  int num_comp; 

  W = 0.0; 
  num_comp=0;
  for( k1=0; k1<c->K-1; k1++ ){
	 for( k2=k1+1; k2<c->K; k2++ ){ /* compare cluster k1 and k2 */
		for( i=0; i<c->n[k1]; i++ ){
		  for( j=0; j<c->n[k2]; j++ ){
			 W += d[c->clust[k1][i]][c->clust[k2][j]];
			 num_comp++;
		  }
		}
	 }
  }
  W /= (double)2*num_comp;

  return W;
}

/** Agglomerative clustering.
	 \param d distance matrix for the N objects
	 \param N number of objects
	 \param dist function giving the between sub-cluster distance; 
	             defined are dgram_dist_singlelinkage(), dgram_dist_completelinkage(),
					 dgram_dist_averagelinkage()
 */
Dendrogram* agglomerative_clustering(const double **d, int N, LinkageFunction dist ){
  Dendrogram **nodes, *tmp; /* at the lowest level, we have N nodes */
  double min_d, cur_d=0.0;
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
  
  //  dprintf("SL-distance = %f\n", dist);

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


/** return the best number of clusters as returned by the
	 Gap-Statistic.
	 To get the data-distribution, all channels in eeg are averaged.
	 \param optargs may contain
	 - "max_num_clusters=int", maximum number of clusters, default=10
	 - "num_ref_dist=int", number of reference distribution iterations, default=ntrials
	 - <b>optargs is passed to:</b> distfunction
	 \return best number of clusters or 0 in case of error
 */

int      eeg_best_num_clusters_gapstat( const EEG *eeg, VectorDistanceFunction distfunction, OptArgList *optargs ){
  GapStatistic *gap;
  double x;
  int max_num_clusters;
  int num_ref_dist;
  EEG *meaneeg;
  double **D;
  int bestclust=0;

  max_num_clusters = 10;
  num_ref_dist = eeg->ntrials;
  if( optarglist_has_key( optargs, "max_num_clusters" ) ){
	 x = optarglist_scalar_by_key( optargs, "max_num_clusters" );
	 if( !isnan( x ) ) max_num_clusters=(int)x;
  }
  if( optarglist_has_key( optargs, "num_ref_dist" ) ){
	 x = optarglist_scalar_by_key( optargs, "num_ref_dist" );
	 if( !isnan( x ) ) num_ref_dist=(int)x;
  }

  meaneeg = eeg_average_channels( eeg );
  eeg_print( stderr, meaneeg, 3 );
  D = eeg_distmatrix( meaneeg, distfunction, ALLOC_IN_FCT, optargs );

  gap = gapstat_init( NULL, max_num_clusters, num_ref_dist );
  gapstat_calculate( gap, meaneeg->data[0], meaneeg->ntrials, meaneeg->n, 
							distfunction, D );

  bestclust = gap->khat;
#ifdef DEBUG
  gapstat_print( stderr, gap );
#endif 

  matrix_free( D, meaneeg->ntrials );
  eeg_free( meaneeg );
  gapstat_free( gap );

  return bestclust;
}
