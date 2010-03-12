/***************************************************************************
 *   Copyright (C) 2010 by Matthias Ihrke   *
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
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "nnsearch.h"
#include <gsl/gsl_sort.h>
#include "pqueue.h"
#include "mathadd.h"

/** Initialize the search tree for fast-nearest-neighbour searching.

	 \param data matrix of m variables with N measurements each
	 \param m dimensionality
	 \param N number of measurements
	 \param optargs may contain:
	  - 'metric=void*' the distance metric used (default=euclidean)
	  - 'max_numel_terminal_node=int' the maximum number of elements in a 
	     terminal node (default=50)
	  - optional arguments for the called function (e.g. metric)
	 \return the searchtree-struct used to look for nearest neighbours
 */
SearchTree* nn_prepare( const double **X, int N, int m, OptArgList *optargs ){
  SearchTree *S;
  double **D;
  int *A;
  VectorDistanceFunction distfct;
  int max_numel;
  void *tmp;
  double x;
  
  /* --------------------------------------*/
  distfct = vectordist_euclidean;
  max_numel = 50;
  if( optarglist_has_key( optargs, "metric" ) ){
	 tmp = optarglist_ptr_by_key( optargs, "metric" );
	 if( tmp )
		distfct = (VectorDistanceFunction) tmp;
  }
  if( optarglist_has_key( optargs, "max_numel_terminal_node" ) ){
	 x = optarglist_scalar_by_key( optargs, "max_numel_terminal_node" );
	 if( !isnan( x ) )
		max_numel=(int)x;
  }
  dprintf("maxnumel=%i\n", max_numel);
  /* -------------------------------------- */

  dprintf("preparing search tree, %i points of dim %i\n", N, m);
  S = searchtree_init(N);
  A=S->A;
  D = vectordist_distmatrix( distfct, X,
									  N, m, NULL, NULL, optargs );
  S->d=X;
  S->distfct = distfct;
  S->optargs = optargs;
  S->m=m;
  matrix_print( X, N, m );
  matrix_print( (const double**)D, N, N );
  S->root = tnode_init(); 
  S->root->c = (int)(random() / (RAND_MAX / N+1)); 
  S->root->start=0;
  S->root->end=N-1;
  S->root->R = vector_max( D[S->root->c], N, NULL );
  dprintf("root center: %i, radius=%f\n", S->root->c, S->root->R );
  
  build_tree_recursive( S->root, D, N, A, max_numel );
  dprintf(" S->N=%i\n", S->N);
  return S;
}

/** Searching the k nearest neighbour of vector x in X.
	 This is a slow implementation:
	 - calculate all distances from x to any point in X
	 - sort the distances
	 - return the indices with the k smallest distances

	 \param X dataset to query
	 \param N,m dimensions of X
	 \param x query point; must be of dimension m
	 \param k number of neighbours to look for
	 \param nn_idx memory for holding the indices for the k NN's 
	 \param nn_dist memory for holding distances for k NN's
	 \param optargs may contain:
	  - 'metric=void*' the distance metric used (default=euclidean)
	  - optional arguments for the called function (e.g. metric)
 */
void nn_search_k_slow( const double **X, int N, int m, const double *x, int k, 
							  int *nn_idx, double *nn_dist, OptArgList *optargs ){
  VectorDistanceFunction distfct;
  void *tmp;
  int i;
  size_t *permut;  
  double *ds;

  permut=(size_t*)malloc(N*sizeof(size_t));
  ds = (double*)malloc( N*sizeof(double) );

  /* --------------------------------------*/
  distfct = vectordist_euclidean;
  if( optarglist_has_key( optargs, "metric" ) ){
	 tmp = optarglist_ptr_by_key( optargs, "metric" );
	 if( tmp )
		distfct = (VectorDistanceFunction) tmp;
  }
  /* --------------------------------------*/
 
  for( i=0; i<N; i++ ){
	 ds[i] = distfct( x, X[i], m, optargs);
  }
  gsl_sort_index( permut, ds, 1, N );    

  for( i=0; i<k; i++ ){
	 nn_idx[i]=permut[i];
	 nn_dist[i] = ds[permut[i]];
  }

  dprintf("Done\n");
  /* freeing */
  free( ds );
  free( permut );
}

/** Searching the k nearest neighbour of vector x in searchtree.
	 This implementation uses a searchtree-based method.
	 
	 \param S the searchtree, previously initialized by nn_prepare
	 \param x query point; must be of dimension m in searchtree S
	 \param k number of neighbours to look for
	 \param nn_idx memory for holding the indices for the k NN's 
	 \param nn_dist memory for holding distances for k NN's
 */
void nn_search_k( const SearchTree *S, const double *x, int k, int *nn_idx,
						double *nn_dist ){
  int i;
  size_t *permut; 
  PriorityQueue *pq=pq_init();
  double dmin, dtmp;

  permut=(size_t*)malloc(k*sizeof(size_t));

  dprintf("searching for %i NN's\n", k );
  /* initialize nearest neighbours */
  for( i=0; i<k; i++ ){ /* random init */
	 nn_idx[i] = (int)(random() / (RAND_MAX / S->N+1)); 
	 nn_dist[i]= S->distfct( x, S->d[nn_idx[i]], S->m, S->optargs );
  }
#ifdef DEBUG
  dprintf(" initial m_k\n");
  vector_print_int( nn_idx, k );
  dprintf(" initial d(m_k,q)\n");
  vector_print( nn_dist, k );
#endif

  /* insert root node */
  dmin = S->distfct( x, S->d[S->root->c], S->m, S->optargs )-S->root->R;
  pq_insert( pq, (void*)S->root, dmin );

  TreeNode *C;
  PQnode *pqn;
  int count, ki;
  double dcq, dclq, dcrq;
  double pdmin; 
  while( 1 ){
	 /* pop next high node */
	 pqn = pq_pop_node(pq);
	 C = (TreeNode*)pqn->content;
	 dmin = pqn->priority;
	 pqnode_free( pqn );
	 dcq = S->distfct( x, S->d[C->c], S->m, S->optargs );

	 /* check for termination */
	 for( count=0,i=0; i<k; i++ ){
		if( dmin<nn_dist[i] ){
		  count++;
		}
	 }
	 if( count==k ){ 
		break;
	 }
	 
	 /* terminal node? */
	 if( tnode_isleaf( C ) ){ /* scan within the cluster */
		for( ki=0; ki<k; ki++ ){
		  for( i=0; i<C->end-C->start+1; i++ ){
			 if( nn_dist[ki] < ABS( dcq - C->cdist[i] ) ){
				dtmp = S->distfct( x, S->d[S->A[C->start+i]], S->m, S->optargs );
				if( dtmp<nn_dist[ki] ){
				  nn_dist[ki]=dtmp;
				  nn_idx[ki] =S->A[C->start+i];
				}
			 }
		  }
		}
	 } else { /* insert child clusters into PQ */
		pdmin = dmin;
		dclq=S->distfct( x, S->d[C->left->c], S->m, S->optargs );
		dcrq=S->distfct( x, S->d[C->right->c], S->m, S->optargs );
		if( (dtmp=(dclq-dcrq-C->left->g)/2.0) > pdmin ){
		  dmin = dtmp;
		}
		if( (dtmp=(dclq - C->left->R)) > pdmin ){
		  dmin = dtmp;
		}
		pq_insert( pq, C->left, dmin );
		dmin=pdmin;
		if( (dtmp=(dcrq-dclq-C->right->g)/2.0) > pdmin ){
		  dmin = dtmp;
		}
		if( (dtmp=(dcrq - C->right->R)) > pdmin ){
		  dmin = dtmp;
		}
		pq_insert( pq, C->right, dmin);
	 }
  }

  /* sort nn arrays */
  gsl_sort_index( permut, nn_dist, 1, k );    
  int *nnidx_tmp=(int*)malloc( k*sizeof(int));
  double *nndist_tmp=(double*)malloc(k*sizeof(double));
  memcpy( nnidx_tmp, nn_idx, k*sizeof(int));
  memcpy( nndist_tmp, nn_dist, k*sizeof(double));
  for( i=0; i<k; i++ ){
	 nn_idx[i]  = nnidx_tmp[permut[i]];
	 nn_dist[i] = nndist_tmp[permut[i]];
  }

  /* cleaning */
  free( nnidx_tmp );
  free( nndist_tmp );
  pq_free( pq );
  free( permut );
}

/**\cond PRIVATE
 */

void build_tree_recursive( TreeNode *C, double **D, int N, int *A, int maxel ){
  TreeNode *L, *R;
  int il,ir;
  int i;

  C->left=tnode_init();
  C->right=tnode_init();
  L = C->left; R = C->right;

  /* left child far away from C */
  dprintf("calc left child c\n");
  L->c = A[C->start];
  for( i=C->start; i<=C->end; i++ ){
	 if( D[L->c][C->c] < D[A[i]][C->c] ){
		L->c = A[i];
	 }
  }
  /* right child far away from sister */
  dprintf("calc right child c\n");
  R->c = A[C->start];
  for( i=C->start; i<C->end; i++ ){
	 if( D[R->c][L->c] < D[A[i]][L->c] ){
		R->c = A[i];
	 }
  }
  dprintf("New Child centers: %i, %i\n", L->c, R->c );

  /* give children data */	 
  L->start  = C->start;
  L->end    = C->start;
  R->start  = C->end;
  R->end    = C->end;
	 
  il=L->end;
  ir=R->start;  
  double tmp;
  while(il <= ir){
	 while(D[L->c][A[il]] < D[R->c][A[il]]){
		if( (tmp=D[R->c][A[il]]-D[L->c][A[il]]) < L->g ){
		  L->g=tmp;
		}
		il++;
	 }
	 while(D[L->c][A[ir]] > D[R->c][A[ir]]){
		if( (tmp=D[L->c][A[ir]]-D[R->c][A[ir]]) < R->g ){
		  R->g=tmp;
		}
		ir--;
	 }
	 vector_print_int( A, N );
	 dprintf("ptrs: il,ir=(%i,%i)\n", il, ir);
		
	 if( il<ir ){
		SWAPT( int, A[il], A[ir] );
		L->end=il;
		R->start=ir;
		il++; ir--;
	 } else {
		L->end=il-1;
		R->start=ir+1;
		break;
	 }
  }

  dprintf("new cranges: (%i-%i), (%i-%i)\n", 
			 L->start, L->end, R->start, R->end );

  /* calculate enclosing radii */
  L->R = 0;
  for( i=L->start; i<=L->end; i++ ){
	 if( D[A[i]][L->c]>L->R )
		L->R = D[A[i]][L->c];
  }
  R->R = 0;
  for( i=R->start; i<=R->end; i++ ){
	 if( D[A[i]][R->c]>R->R )
		R->R = D[A[i]][R->c];
  }

  vector_print_int( A, N );
  dprintf("L->R=%f, R->R=%f\n", L->R, R->R);

  /* terminal node? */
  if( L->end-L->start > maxel ){
	 dprintf("left recursion\n");
	 build_tree_recursive( L, D, N, A, maxel );
  } else {
	 L->cdist = (double*)malloc( (L->end - L->start + 1)*sizeof(double));
	 for( i=0; i<L->end-L->start+1; i++ ){
		L->cdist[i] = D[L->c][A[L->start+i]];
	 }
  }
  if( L->end-L->start > maxel ){
	 dprintf("right recursion\n");
	 build_tree_recursive( R, D, N, A, maxel );
  } else {
	 R->cdist = (double*)malloc( (R->end - R->start + 1)*sizeof(double));
	 for( i=0; i<R->end-R->start+1; i++ ){
		R->cdist[i] = D[R->c][A[R->start+i]];
	 }
  }
}

SearchTree* searchtree_init( int n ){
  int i;
  SearchTree *S=(SearchTree*)malloc(sizeof(SearchTree));
  S->A = (int*)malloc(n*sizeof(int));
  for( i=0; i<n; i++ ){
	 S->A[i] = i;
  }
  S->m = 0;
  S->N = n;
  S->root = NULL;
  return S;
}

TreeNode* tnode_init(){
  TreeNode *t=(TreeNode*)malloc(sizeof(TreeNode));
  t->c=-1;
  t->R=-1;
  t->g=DBL_MAX;
  t->start=-1;
  t->end=-1;
  t->cdist=NULL;
  t->left=NULL;
  t->right=NULL;
  return t;
}

bool tnode_isleaf( TreeNode *C ){
  if( C->cdist )
	 return TRUE;
  else
	 return FALSE;
}


/**\endcond */
