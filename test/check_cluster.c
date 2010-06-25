/*
** check_cluster.c
** 
** Unit-Test suite using libcheck
**
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
*/

#include<stdlib.h>
#include<check.h>
#include<math.h>
#include "check_all.h"
#include "array.h"
#include "distances.h"
#include "helper.h"
#include "clustering.h"
#include "io.h"

START_TEST (test_kmedoids)
{
  int n=100, K=4;
  int i;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  Clusters *C=kmedoids( dmat, K, NULL );
  //  cluster_print( stderr, C );

  fail_unless( C->K==K );
  for( i=0; i<C->K; i++ ){
	 fail_unless( C->n>0 );
  }

  cluster_free( C );
  array_free( dat );
  array_free( dmat );
}
END_TEST

START_TEST (test_kmedoids_repeat)
{
  int n=100, K=4;
  int rep=10;
  int i;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  Clusters *C =kmedoids_repeat( dmat, K, rep );
  Clusters *C2=kmedoids( dmat, K, NULL );
  //  cluster_print( stderr, C );

  fail_unless( C->K==K );
  fail_unless( C2->K==K );
  for( i=0; i<C->K; i++ ){
	 fail_unless( C->n>0 );
	 fail_unless( C2->n>0 );
  }

  fail_if( cluster_within_scatter( dmat, C )>
			  cluster_within_scatter( dmat, C2 ) );
  fail_if( cluster_between_scatter( dmat, C )<
			  cluster_between_scatter( dmat, C2 ) );
  dprintf("Crep=(%f,%f) vs. C=(%f,%f)\n",
			 cluster_within_scatter( dmat, C ),
			 cluster_between_scatter( dmat, C ),
			 cluster_within_scatter( dmat, C2 ),
			 cluster_between_scatter( dmat, C2 ) );

  cluster_free( C );
  cluster_free( C2 );

  array_free( dat );
  array_free( dmat );
}
END_TEST

START_TEST (test_scatter)
{
  int n=100, K=4;
  int i;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  Clusters *C=kmedoids( dmat, K, NULL );

  double wit=cluster_within_scatter( dmat, C );
  double bet=cluster_between_scatter( dmat, C );

  dprintf("within=%f, between=%f\n", wit, bet );
  fail_if( wit>bet );

  fail_unless( C->K==K );
  for( i=0; i<C->K; i++ ){
	 fail_unless( C->n>0 );
  }

  cluster_free( C );
  array_free( dat );
  array_free( dmat );
}
END_TEST

START_TEST (test_agglomerative)
{
  int n=100;
  int i;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  
  Dendrogram *dgram=agglomerative_clustering( dmat, dgram_dist_singlelinkage );
  //  dgram_print( dgram );

  dgram_free( dgram );
  array_free( dat );
  array_free( dmat );
}
END_TEST

START_TEST (test_numleaves)
{
  int n=100;
  int i;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  
  Dendrogram *dgram=agglomerative_clustering( dmat, dgram_dist_singlelinkage );
  int nl=dgram_num_leaves( dgram );
  dprintf("numleaves=%i\n",nl);
  fail_if( nl!=n );
  //  dgram_print( dgram );

  dgram_free( dgram );
  array_free( dat );
  array_free( dmat );
}
END_TEST

#ifdef MATIO
START_TEST (test_dgram_to_matlab)
{
  /* Test in MATLAB as follows:

	  load test.mat
	  Y=pdist( dat, 'euclidean');
	  Z=linkage(Y, 'single');
	  figure;
	  subplot( 2, 1, 1 );
	  dendrogram( Z, 'colorthreshold','default');
	  subplot( 2, 1, 2 );
	  dendrogram( dend, 'colorthreshold','default');
	*/

  int n=100;
  int i, succ;
  Array *dat = array_randunif( 0, 2, n, 2 );
  Array *dmat= matrix_distmatrix( vectordist_euclidean,
											 dat, NULL, NULL );
  
  Dendrogram *dgram=agglomerative_clustering( dmat, dgram_dist_singlelinkage );
  Array *ml=dgram_to_matlab( dgram );
  array_print( ml, -1, stderr );
  fail_if( ml->size[0]!=n-1 );
  fail_if( ml->size[1]!=3 );

  succ=write_array_matlab( dat, "dat", "test.mat", FALSE);
  fail_if( succ );
  succ=write_array_matlab( dmat, "dmat", "test.mat", TRUE);
  fail_if( succ );
  succ=write_array_matlab( ml, "dend", "test.mat", TRUE);
  fail_if( succ );

  dgram_free( dgram );
  array_free( dat );
  array_free( ml );
  array_free( dmat );
}
END_TEST
#endif

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_cluster_suite (void){
  Suite *s = suite_create ("Clustering-Functions");

  TCase *tc_core = tcase_create ("ClusterCore");
  tcase_add_test (tc_core, test_kmedoids);
  tcase_add_test (tc_core, test_kmedoids_repeat);
  tcase_add_test (tc_core, test_scatter);
  tcase_add_test (tc_core, test_agglomerative);
  tcase_add_test (tc_core, test_numleaves);

#ifdef MATIO
  tcase_add_test (tc_core, test_dgram_to_matlab);
#endif

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

