/*
** check_average.c
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
#include "averaging.h"


START_TEST (test_hierarchical)
{
  int N=10,
	 C=2,
	 n=12;
  Array *a=array_randunif( SEED_RANDOM, 3, N, C, n );
  Array *b=array_randunif( SEED_RANDOM, 2, N, n );

  Array *asmall=array_slice( a, ":,1,:");
  Array *da=matrix_distmatrix( vectordist_euclidean, asmall, NULL, NULL ) ;  
  Array *db=matrix_distmatrix( vectordist_euclidean, b, NULL, NULL ) ;  


  Array *avga=hierarchical_average( a, da, average_example, NULL );
  Array *avgac=array_new2( DOUBLE, 2, C, n );
  int i,j, k;

  for( j=0; j<C; j++ ){
	 for( k=0; k<n; k++ ){
		for( i=0; i<N; i++ ){
		  array_INDEX2(avgac,double,j,k) +=
			 array_INDEX3( a,double, i,j,k );
		}
		array_INDEX2(avgac,double,j,k) /= (double)N;
	 }
  }

  /* array_print( avgac, -1, stderr ); */
  /* array_print( avga, -1, stderr ); */

  /* hierarchical average with average_example is the same as
	  a pointwise avg */
  for( i=0; i<array_NUMEL( avga ); i++ ){
	 fail_unless( cmpdouble(array_INDEX1( avga, double, i ),
									array_INDEX1( avgac, double, i ),3 )==0 );
  }

  Array *avgb=hierarchical_average( b, db, average_example, NULL );
  Array *avgbc=matrix_mean(b, 0 );

  /* array_print( avgbc, -1, stderr ); */
  /* array_print( avgb, -1, stderr ); */

  for( i=0; i<array_NUMEL( avgb); i++ ){
	 fail_unless( cmpdouble(array_INDEX1( avgb, double, i ),
									array_INDEX1( avgbc, double, i ),3 )==0 );
  }

  array_free( a );
  array_free( asmall );
  array_free( b );
  array_free( da );
  array_free( avga );
  array_free( avgac );
  array_free( db );
  array_free( avgb );
  array_free( avgbc );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_average_suite (void){
  Suite *s = suite_create ("Average-Functions");

  TCase *tc_core = tcase_create ("AvgCore");
  tcase_add_test (tc_core, test_hierarchical);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

