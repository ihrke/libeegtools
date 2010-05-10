/*
** check_warping.c
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

START_TEST (test_cumulate)
{
  int n=10;
  Array *a=array_new_dummy( DOUBLE, 2, n, 2 );
  Array *b=array_new_dummy( DOUBLE, 2, n, 2 );
  Array *c=array_randunif( 0, 2, n, 10 );
  Array *d=array_randunif( 0, 2, n, 10 );
  array_reverse( b );

  Array *m=distmatrix_signaldist( vectordist_euclidean,
											 a, b, NULL, NULL );
  fail_if( m->ndim<2 || m->size[0]!=n || m->size[1]!=n );
  matrix_dtw_cumulate( m, FALSE, NULL );  
  int i, j;
  for( i=0; i<n; i++ ){ /* monotonicity test */
	 fail_unless( mat_IDX( m, n-1, n-1 )>=mat_IDX(m,i,i),
					  "last element (%i,%i)=%f not maximum (%i,%i)=%f",
					  n-1,n-1,mat_IDX( m, n-1, n-1 ),i,i, mat_IDX(m,i,i));
  }

  Array *m2=distmatrix_signaldist( vectordist_euclidean,
											  c, d, NULL, NULL );
  fail_if( m2->ndim<2 || m2->size[0]!=n || m2->size[1]!=n );
  matrix_dtw_cumulate( m2, FALSE, NULL );  
  for( i=0; i<n; i++ ){ /* monotonicity test */
	 fail_unless( mat_IDX( m2, n-1, n-1 )>=mat_IDX(m2,i,i),
					  "last element (%i,%i)=%f not maximum (%i,%i)=%f",
					  n-1,n-1,mat_IDX( m2, n-1, n-1 ),i,i, mat_IDX(m2,i,i));
  }

  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
  array_free( m );
  array_free( m2 );
}
END_TEST

START_TEST (test_backtrack)
{
  int n=10;
  int i;
  Array *a=array_randunif( 3423, 2, n, 2 );
  Array *b=array_randunif( 250939, 2, n, 2 );
  Array *m=distmatrix_signaldist( vectordist_euclidean,
											 a, b, NULL, NULL );
  fail_if( m->ndim<2 || m->size[0]!=n || m->size[1]!=n );
  matrix_dtw_cumulate( m, FALSE, NULL );  

  Array *P = matrix_dtw_backtrack( m );

  /* start/end condition */
  fail_unless( array_INDEX2( P, uint, 0,0)==0 );
  fail_unless( array_INDEX2( P, uint, 1,0)==0 );
  fail_unless( array_INDEX2( P, uint, 0,P->size[1]-1)==n-1 );
  fail_unless( array_INDEX2( P, uint, 1,P->size[1]-1)==n-1 );
  
  for( i=1; i<P->size[1]; i++ ){ /* monotonicity test */
	 fail_unless( array_INDEX2( P, uint, 0, i )>=
					  array_INDEX2( P, uint, 0, i-1 ) );
	 fail_unless( array_INDEX2( P, uint, 1, i )>=
					  array_INDEX2( P, uint, 1, i-1 ) );
  }  
  /* check for no gaps */
  for( i=1; i<P->size[1]; i++ ){
	 fail_if( array_INDEX2( P, uint, 0, i )-
				 array_INDEX2( P, uint, 0, i-1 )>1 );
	 fail_if( array_INDEX2( P, uint, 1, i )-
				 array_INDEX2( P, uint, 1, i-1 )>1 );
	 fail_if( array_INDEX2( P, uint, 0, i )-
				 array_INDEX2( P, uint, 0, i-1 )==0 &&
				 array_INDEX2( P, uint, 1, i )-
				 array_INDEX2( P, uint, 1, i-1 )==0 );
  }
  array_free( a );
  array_free( b );
  array_free( m );
  array_free( P );
}
END_TEST

START_TEST (test_slopeconstraint)
{
  int n=1000;
  int i;
  Array *a=array_randunif( 3423, 2, n, 2 );
  Array *b=array_randunif( 250939, 2, n, 2 );

  Array *m=distmatrix_signaldist( vectordist_euclidean,
											  a, b, NULL, NULL );
  fail_if( m->ndim<2 || m->size[0]!=n || m->size[1]!=n );

  OptArgList *opts=optarglist("slope_constraint=int", SLOPE_CONSTRAINT_NONE );
  Array *m1=matrix_dtw_cumulate( m, TRUE, opts );  
  optarglist_optarg_by_key(opts, "slope_constraint")->data_scalar=SLOPE_CONSTRAINT_LAX;
  Array *m2=matrix_dtw_cumulate( m, TRUE, opts );
  optarglist_optarg_by_key(opts, "slope_constraint")->data_scalar=SLOPE_CONSTRAINT_MEDIUM;
  Array *m3=matrix_dtw_cumulate( m, TRUE, opts );  
  optarglist_optarg_by_key(opts, "slope_constraint")->data_scalar=SLOPE_CONSTRAINT_SEVERE;
  Array *m4=matrix_dtw_cumulate( m, TRUE, opts );  

  Array *P1 = matrix_dtw_backtrack( m1 );
  Array *P2 = matrix_dtw_backtrack( m2 );
  Array *P3 = matrix_dtw_backtrack( m3 );
  Array *P4 = matrix_dtw_backtrack( m4 );

  /* check for no gaps */
  for( i=1; i<P1->size[1]; i++ ){
	 fail_if( array_INDEX2( P1, uint, 0, i )-
				 array_INDEX2( P1, uint, 0, i-1 )>1 );
	 fail_if( array_INDEX2( P1, uint, 1, i )-
				 array_INDEX2( P1, uint, 1, i-1 )>1 );
	 fail_if( array_INDEX2( P1, uint, 0, i )-
				 array_INDEX2( P1, uint, 0, i-1 )==0 &&
				 array_INDEX2( P1, uint, 1, i )-
				 array_INDEX2( P1, uint, 1, i-1 )==0 );
  }
  for( i=1; i<P2->size[1]; i++ ){
	 fail_if( array_INDEX2( P2, uint, 0, i )-
				 array_INDEX2( P2, uint, 0, i-1 )>1 );
	 fail_if( array_INDEX2( P2, uint, 1, i )-
				 array_INDEX2( P2, uint, 1, i-1 )>1 );
	 fail_if( array_INDEX2( P2, uint, 0, i )-
				 array_INDEX2( P2, uint, 0, i-1 )==0 &&
				 array_INDEX2( P2, uint, 1, i )-
				 array_INDEX2( P2, uint, 1, i-1 )==0 );
  }
  for( i=1; i<P3->size[1]; i++ ){
	 fail_if( array_INDEX2( P3, uint, 0, i )-
				 array_INDEX2( P3, uint, 0, i-1 )>1 );
	 fail_if( array_INDEX2( P3, uint, 1, i )-
				 array_INDEX2( P3, uint, 1, i-1 )>1 );
	 fail_if( array_INDEX2( P3, uint, 0, i )-
				 array_INDEX2( P3, uint, 0, i-1 )==0 &&
				 array_INDEX2( P3, uint, 1, i )-
				 array_INDEX2( P3, uint, 1, i-1 )==0 );
  }
  for( i=1; i<P4->size[1]; i++ ){
	 fail_if( array_INDEX2( P4, uint, 0, i )-
				 array_INDEX2( P4, uint, 0, i-1 )>1 );
	 fail_if( array_INDEX2( P4, uint, 1, i )-
				 array_INDEX2( P4, uint, 1, i-1 )>1 );
	 fail_if( array_INDEX2( P4, uint, 0, i )-
				 array_INDEX2( P4, uint, 0, i-1 )==0 &&
				 array_INDEX2( P4, uint, 1, i )-
				 array_INDEX2( P4, uint, 1, i-1 )==0 );
  }
  /* start/end condition */
  fail_unless( array_INDEX2( P1, uint, 0,0)==0 );
  fail_unless( array_INDEX2( P1, uint, 1,0)==0 );
  fail_unless( array_INDEX2( P1, uint, 0,P1->size[1]-1)==n-1 );
  fail_unless( array_INDEX2( P1, uint, 1,P1->size[1]-1)==n-1 );

  /* slope constraints check */
  dprintf("sizes: %i, %i, %i, %i\n", P1->size[1], P2->size[1], P3->size[1], P4->size[1] );
  fail_unless( P1->size[1]>=P2->size[1] );
  fail_unless( P2->size[1]>=P3->size[1] );
  fail_unless( P3->size[1]>=P4->size[1] );

  array_free( a );
  array_free( b );
  array_free( m );
  array_free( m1 );  
  array_free( m2 );
  array_free( m3 );  
  array_free( m4 );
  array_free( P1 );
  array_free( P2 );
  array_free( P3 );
  array_free( P4 );
}
END_TEST


/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_warping_suite (void){
  Suite *s = suite_create ("Warping-Functions");

  TCase *tc_core = tcase_create ("WarpCore");
  tcase_add_test (tc_core, test_cumulate);
  tcase_add_test (tc_core, test_backtrack);
  tcase_add_test (tc_core, test_slopeconstraint);
 

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

