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
#include "linalg.h"
#include <time.h>
#include "eeg.h"


/** run the kmedoids function a couple of times and pick the best
	 in terms of within scatter.
	 \param dist ance matrix
	 \param K number of medoids to compute
	 \param repeat number of repetitions
	 \return clusters - function-allocated memory. contains
	                    K cluster-structs.
 */
Clusters*    kmedoids_repeat( const Array *distmat, int K, int repeat ){
  Clusters *C;
  Clusters *Cnew;
  int i;
  double ratio, tmp;
  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix) return NULL;
  int N = distmat->size[0];

  long initial_seed=(long)time(NULL);
  OptArgList *opts=optarglist("seed=long", initial_seed );

  C = cluster_init( K, N );
  ratio=DBL_MAX;
  for( i=0; i<repeat; i++ ){
	 optarglist_optarg_by_key( opts, "seed" )->data_scalar=initial_seed+i*10;
	 Cnew = kmedoids( distmat, K, opts );
	 tmp = 1.0;
	 
	 tmp =  cluster_within_scatter ( distmat, Cnew ); 
	 tmp /= cluster_between_scatter( distmat, Cnew ); 

	 if( tmp<ratio ){
		dprintf("Run %i: %f\n", i, tmp);
		ratio=tmp; 
		cluster_copy( C, Cnew );
	 }
	 cluster_free( Cnew );
  }
  optarglist_free( opts );

  return C;
}

/** do K-Medoids clustering on the distance-matrix.
	 \param dist ance matrix
	 \param K number of medoids to compute
    \param optargs may contain:
	 - <tt>seed=int</tt> seed for initializing the cluster-configuration 
	 (if 0, use time(NULL))
	 \return clusters - function-allocated memory. contains
	                   K cluster-structs.	
*/
Clusters*    kmedoids(const Array *distmat, int K, OptArgList *optargs ){
  Clusters *C, *Cnew;
  int *medoids;
  int i, j, k, r, minidx, num_not_changed;
  double sum, minsum, mindist, curdist;
  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix) return NULL;
  int N = distmat->size[0];

  /* optional arguments */
  unsigned long seed=0;
  double x;
  optarg_PARSE_SCALAR( optargs, "seed", seed, unsigned long, x );

  if( K>N ){
	 errprintf("K>=N: %i>=%i\n", K, N );
	 return NULL;
  }

  /* memory alloc */
  medoids = (int*) malloc( K*sizeof(int) );
  C    = cluster_init(K, N);
  Cnew = cluster_init(K, N);

  /* initialization step, random init */ 
  const gsl_rng_type * T;
  gsl_rng * randgen;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  randgen = gsl_rng_alloc (T);
  if( seed==0 )
	 seed=(unsigned long)time(NULL);
  gsl_rng_set (randgen, seed );

  Array *permut=array_new2( UINT, 1, K );
  for( i=0; i<K; i++ )
	 array_INDEX1( permut, uint, i )=i;
  array_shuffle( permut, seed );
  for( i=0; i<N; i++ ){
	 if( i<K ){ /* first each partition gets a guy */
		r = array_INDEX1( permut, uint, i );
	 } else {
		r = gsl_rng_uniform_int( randgen, K );
	 }
	 C->clust[r][C->n[r]] = i;
	 (C->n[r])++;
  }
  array_free( permut );

  dprintf("initialized cluster\n");  
#ifdef DEBUG
  cluster_print(stderr, C); 
#endif

  /*------------------ computation ----------------------------*/
  num_not_changed = 0;
  while(num_not_changed<1){ /* iterate until convergence */
	 for( i=0; i<K; i++ ){ /* minimize distance to one of trials */
		minsum = DBL_MAX;
		minidx = 0;
		for( j=0; j<C->n[i]; j++ ){
		  sum = 0;
		  for( k=0; k<C->n[i]; k++){
			 sum += mat_IDX( distmat, C->clust[i][j], C->clust[i][k] );
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
		  curdist = mat_IDX( distmat, medoids[j], i );
		  if( curdist<mindist ){
			 mindist=curdist;
			 minidx=j;
		  }
		}
		/*dprintf("Assigning trial %i to cluster %i\n", i, minidx);*/
		Cnew->clust[minidx][Cnew->n[minidx]] = i;
		Cnew->n[minidx] += 1;
	 }

	 if(!cluster_compare(C, Cnew))
		num_not_changed++;
	 else
		cluster_copy(C, Cnew);
  }
  /*------------------ /computation ----------------------------*/

  free(medoids);
  cluster_free(Cnew);
  gsl_rng_free (randgen);

  return C;
}

/** Compare whether two cluster structs are "identical", meaning
	 that each cluster contains the same arguments (sequence
	 does not matter). Superficial tests are first carried out, 
	 then intense comparison.
*/
int cluster_compare(const Clusters *c1, const Clusters *c2){
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
void cluster_copy(Clusters *dest, const Clusters *src){
  int i, j;

  dest->K = src->K;
  for( i=0; i<src->K; i++ ){
	 dest->n[i] = src->n[i];
	 for( j=0; j<src->n[i]; j++ ){
		dest->clust[i][j] = src->clust[i][j];
	 }
  }
  
}

void cluster_free(Clusters *c){
  int i;
  free( c->n );
  for( i=0; i<c->K; i++)
	 free(c->clust[i]);
  free(c->clust);
  free(c);
}

Clusters* cluster_init(int K, int maxN){
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

void cluster_print(FILE *o,const Clusters *c){
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

#if 0
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
	 //	 cluster_print( stderr, C );
	 gap->Wk[k-1] = get_within_scatter( D, n, C );
	 cluster_free( C );
  }

  for( b=0; b<gap->B; b++ ){ /* monte Carlo for ref-data*/
	 for( k=1; k<=gap->K; k++ ){
		if( gap->progress ){
		  gap->progress( PROGRESSBAR_CONTINUE_LONG, (b+1)*gap->K+k );
		}
		
		Xr = gap_get_reference_distribution_simple( (const double **) X, n, p, Xr );
		Dr = vectordist_distmatrix( distfunction, (const double**)Xr, n, p, Dr, NULL, NULL );
		C = kmedoids_repeat( (const double**)Dr, n, k, 50 );
		//		cluster_print(stderr, C );
		gap->Wkref[b][k-1] = get_within_scatter( (const double**)Dr, n, C );
		cluster_free( C );
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

#endif

/** \brief get within-scatter within clusters.
	 compute \f$ W_k = \sum_{r=1}^{k}\frac{1}{2n_r}D_r\f$
	 where
	 \f$
	 D_r = \sum_{i,i'\in C_r}d_{ii'}
	 \f$
	 for a given cluster assignment.
	 \param d - distance matrix (NxN)
	 \param c - cluster-assigment (Clusters - struct)
	 \return the within-scatter
 */
double   cluster_within_scatter (const Array *distmat, const Clusters *c){
  double *sums;
  double W;
  int i, j, k; 
  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix) return -1;

  //  cluster_print( stderr, c);
  sums = (double*) malloc( c->K*sizeof(double) );
  W = 0;
  int numcomp=0;
  for( i=0; i<c->K; i++ ){ /* cluster loop */
	 sums[i] = 0;
	 for( j=0; j<(c->n[i])-1; j++ ){
		for( k=j+1; k<c->n[i]; k++){	 
		  sums[i] += mat_IDX( distmat, c->clust[i][j], c->clust[i][k]);
		  numcomp++;
		}
	 }
	 //	 sums[i] /= (double)2*c->n[i];
	 /* dprintf("sums[i]=%f\n", sums[i]); */
	 W += sums[i];
  }
  W /= (double)2*numcomp;

  free(sums);

  return W;
}

/** \brief get within-scatter between clusters.
	 (sum of distances of means of clusters)
	 \param d distance matrix (NxN)
	 \param c the partioning of d
	 \return the between-scatter
 */ 
double   cluster_between_scatter(const Array *distmat, const Clusters *c){
  double W;
  int k1, k2, i, j;
  int num_comp; 
  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix) return -1;

  W = 0.0; 
  num_comp=0;
  for( k1=0; k1<c->K-1; k1++ ){
	 for( k2=k1+1; k2<c->K; k2++ ){ /* compare cluster k1 and k2 */
		for( i=0; i<c->n[k1]; i++ ){
		  for( j=0; j<c->n[k2]; j++ ){
			 W += mat_IDX( distmat, c->clust[k1][i], c->clust[k2][j] );
			 num_comp++;
		  }
		}
	 }
  }
  W /= (double)2*num_comp;

  return W;
}


/** \brief Agglomerative clustering.
	 
	 The function was tested against MATLAB's hierarchical cluster-analysis.

	 \param distmat distance matrix for the N objects
	 \param distfct function giving the between sub-cluster distance; 
	             defined are dgram_dist_singlelinkage(), dgram_dist_completelinkage(),
					 dgram_dist_averagelinkage()
    \return Dendrogram 
 */
Dendrogram*  agglomerative_clustering(const Array *distmat, LinkageFunction distfct){
  Dendrogram **nodes, *tmp; /* at the lowest level, we have N nodes */
  double min_d, cur_d=0.0;
  int min_i=0, min_j=0;
  int num_nodes;
  int i, j;
  bool ismatrix;
  matrix_CHECKSQR( ismatrix, distmat );
  if( !ismatrix) return NULL;
  int N=distmat->size[1];

  nodes = (Dendrogram**) malloc( N*sizeof( Dendrogram* ) );

  int clustidx=0;
  for( i=0; i<N; i++ ){ /* initialize terminal nodes */
	 nodes[i] = dgram_init( i, NULL, NULL );
	 nodes[i]->clustnum=clustidx++;
  }
  num_nodes = N;

  while( num_nodes > 1 ){
	 dprintf("nmodes = %i\n", num_nodes );
	 min_d = DBL_MAX;
	 for( i=0; i<num_nodes-1; i++ ){
		for( j=i+1; j<num_nodes; j++ ){ 
		  cur_d = distfct( distmat, nodes[i], nodes[j] );
		  if( cur_d<=min_d ){
			 min_d=cur_d;
			 min_i=i;
			 min_j=j;
		  }
		}
	 }
	 
	 /* we know min_i and min_j, now */
	 tmp = dgram_init( -1, nodes[min_i], nodes[min_j] ); /* link the two nodes */
	 tmp->height=min_d;
	 tmp->clustnum=clustidx++;
	 nodes[min_i] = tmp; /* remove min_i */
	 tmp = nodes[num_nodes-1]; /* swap with last entry in current list */
	 nodes[min_j] = tmp;
	 num_nodes--;
	 dprintf("Linking (%i,%i), height=%f\n", min_i, min_j, min_d );
  }
  tmp = nodes[0]; /* this is the root-node now */
  free( nodes );
  return tmp;
}

/** \cond PRIVATE */
void _dgram_walk_matlab( const Dendrogram *t, Array *a ){
  if( t->val>=0 ){ /* leave */
	 return;
  }
  mat_IDX( a, t->clustnum, 0 ) = t->left->clustnum+1;
  mat_IDX( a, t->clustnum, 1 ) = t->right->clustnum+1;
  mat_IDX( a, t->clustnum, 2 ) = t->height; 
  _dgram_walk_matlab( t->left, a );
  _dgram_walk_matlab( t->right, a );
}
/** \endcond */


/** \brief convert a Dendrogram to a N-1 x 3 array as used by 
	 MATLAB.

	 This can be used for conveniently plotting a dendrogram in 
	 MATLAB.
	 
	 From the MATLAB-manual about the output format:
	 (http://www.mathworks.com/access/helpdesk/help/toolbox/stats/linkage.html)
	 
	 \verbatim
	 The output, Z, is an (m-1)-by-3 matrix containing cluster tree
	 information. The leaf nodes in the cluster hierarchy are the objects
	 in the original data set, numbered from 1 to m. They are the singleton
	 clusters from which all higher clusters are built. Each newly formed
	 cluster, corresponding to row i in Z, is assigned the index m+i, where
	 m is the total number of initial leaves.
	 
	 Columns 1 and 2, Z(i,1:2), contain the indices of the objects that
	 were linked in pairs to form a new cluster. This new cluster is
	 assigned the index value m+i. There are m-1 higher clusters that
	 correspond to the interior nodes of the hierarchical cluster tree.
	 
	 Column 3, Z(i,3), contains the corresponding linkage distances between
	 the objects paired in the clusters at each row i.
	 
	 For example, consider a case with 30 initial nodes. If the tenth
	 cluster formed by the linkage function combines object 5 and object 7
	 and their distance is 1.5, then row 10 of Z will contain the values
	 (5, 7, 1.5). This newly formed cluster will have the index
	 10+30=40. If cluster 40 shows up in a later row, that means this newly
	 formed cluster is being combined again into some bigger cluster.
	 \endverbatim
	 
	 \param dgram a Dendrogram
	 \param N number of terminal 
	 \return a N-1 x 3 Array in the format described above
*/
Array* dgram_to_matlab( const Dendrogram *dgram ){
  int N=dgram_num_leaves( dgram );
  Array *mat=array_new2( DOUBLE, 2, N+N-1, 3 );
  _dgram_walk_matlab( dgram, mat );
  char buf[200];
  sprintf( buf, "%i-%i,:", N, N+N-2 );
  Array *rmat=array_slice( mat, buf );
  array_free(mat);
  return rmat;
}

/** Allocate memory for a single Dendrogram Node and set its value to val.
	 left and right can be set to NULL.
*/
Dendrogram* dgram_init(int val, Dendrogram *left, Dendrogram *right ){
  Dendrogram *c;
  c = (Dendrogram*) malloc( sizeof(Dendrogram) );
  c->val = val;
  c->height=-1.0;
  c->clustnum=0;
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

/** \brief preorder traversal of Dendrogram. 

	 Store all non-negative elements in val and return 
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
/** \brief single linkage clustering.
	 \f[
	 d_{SL}(G, H) = \min_{i\in G, j\in H} d_{ij}
	 \f]
 */
double dgram_dist_singlelinkage  (const Array *d, const Dendrogram *c1, const Dendrogram *c2){
  double dist;
  int *el1, *el2;
  int n1, n2;
  int i, j;
  bool ismatrix;
  matrix_CHECK( ismatrix,d );
  if( !ismatrix ) return -1;
  if( d->size[0]!=d->size[1] ){	 
	 errprintf("distmat must be N x N, got (%i,%i)\n",d->size[0],d->size[1]);
	 return -1; 
  }
  int N = d->size[0];

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
		if( mat_IDX( d, el1[i], el2[j]) < dist ){
		  dist = mat_IDX( d, el1[i], el2[j]);
		}
	 }
  }
  
  //dprintf("SL-distance = %f\n", dist);

  free(el1);
  free(el2);
  return dist;
}

/** \brief complete linkage clustering.
	 \f[
	 d_{CL}(G, H) = \max_{i\in G, j\in H} d_{ij}
	 \f]
 */
double dgram_dist_completelinkage(const Array *d, const Dendrogram *c1, const Dendrogram *c2){
  double dist;
  int *el1, *el2;
  int n1, n2;
  int i, j;
  bool ismatrix;
  matrix_CHECK( ismatrix,d );
  if( !ismatrix ) return -1; 
  if( d->size[0]!=d->size[1] ){	 
	 errprintf("distmat must be N x N, got (%i,%i)\n",d->size[0],d->size[1]);
	 return -1; 
  }
  int N = d->size[0];

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
		if( mat_IDX( d, el1[i], el2[j] )>dist ){
		  dist = mat_IDX( d, el1[i], el2[j]);
		}
	 }
  }
  
  //  dprintf("SL-distance = %f\n", dist);

  free(el1);
  free(el2);
  return dist;
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

  fprintf( stdout, "'%p(%4i,%.2f,%4i)' - ", t, t->clustnum, t->height, t->val );
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

/** \brief return the number of leaves in the Dendrogram.
 */ 
int dgram_num_leaves( const Dendrogram *t ){
  int r=0;
  if( t==NULL ){
	 return 0;
  } else if( t->val!=-1 ){
	 r=1;
  }
  r += dgram_num_leaves( t->left );
  r += dgram_num_leaves( t->right );
  return r;
}

void dgram_print_node( Dendrogram *t ){
  fprintf( stdout, "t=%p, val=%i, height=%f, left=%p, right=%p\n", t, t->val, t->height, t->left, t->right );
}
#undef END_NODE 
#undef INTERMEDIATE_NODE 
#undef NULL_NODE

#if 0
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
							distfunction, (const double**)D );

  bestclust = gap->khat;
#ifdef DEBUG
  gapstat_print( stderr, gap );
#endif 

  dblpp_free( D, meaneeg->ntrials );
  eeg_free( meaneeg );
  gapstat_free( gap );

  return bestclust;
}
#endif
