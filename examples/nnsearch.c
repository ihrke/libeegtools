#include "optarg.h"
#include "mathadd.h"
#include "nnsearch.h"
#include <time.h>

int main(){
  double **X;
  double *q;
  int N, p;
  int i,j;
  N=11;
  p=3;
  X = (double**)malloc( N*sizeof(double*) );
  q = (double*) malloc( p*sizeof(double ) );
  for( i=0; i<N; i++){
	 X[i] = (double*)malloc( p*sizeof(double*) );
	 for( j=0; j<p; j++ ){
		X[i][j] = drand48();
	 }
  }
  for( i=0; i<p; i++ ){
	 q[i] = drand48();
  }
  srandom(time(NULL));
  OptArgList *o;
  o=optarglist("max_numel_terminal_node=int", 2);
  optarglist_print( o, stderr );

  SearchTree *S;
  S = nn_prepare( (const double**)X, N, p, o );

  dprintf("S->N=%i\n", S->N );
  int k=N/5;
  int* nn_idx = (int*)malloc( k*sizeof(int));
  double *nnd = (double*)malloc( k*sizeof(double));
  dprintf("looking for k=%i neighbours of vector\n",k);
  dblp_print( q, p );

  nn_search_k( S, q, k, nn_idx, nnd );
  dprintf("Finished searching, result: \n");
  dblp_print( nnd, k );
  dblp_print_int(nn_idx, k);

  dprintf("slow searching:\n");
  nn_search_k_slow( (const double**)X, N, p, (const double*)q, k, nn_idx, nnd, NULL );
  dprintf("Finished searching, result: \n");
  dblp_print( nnd, k );
  dblp_print_int(nn_idx, k);

  double *ds = (double*)malloc( N*sizeof(double) );
  dblpp_print( (const double**)X, N, p );
  for( i=0; i<N; i++ ){
	 ds[i] = vectordist_euclidean( q, X[i], p, NULL);
  }
  dblp_print( ds, N );

  for( i=0; i<N; i++){
	 free( X[i] );
  }
  free( ds );
  free( X );
  free( q );
  free( nn_idx );
  free( nnd );
}
