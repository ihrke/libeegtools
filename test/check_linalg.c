/*
** check_list.c
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
#include "linalg.h"
#include "helper.h"

START_TEST (test_gsl)
{
  gsl_matrix *g;
  Array *a = array_new_dummy( DOUBLE, 2, 10, 11 );

  g = matrix_to_gsl( a, TRUE );
  fail_unless( memcmp( g->data, a->data, a->nbytes )==0 );
  gsl_matrix_free( g );

  g = matrix_to_gsl( a, FALSE );
  fail_unless( g->data==a->data );
  gsl_matrix_free( g );

  Array *b = array_new_dummy( INT, 2, 11, 12 );
  g = matrix_to_gsl( b, TRUE );
  fail_unless( cmpdouble( g->data[1], 1.0, 3 )==0 );
  gsl_matrix_free( g );

  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_mean)
{
  Array *a = array_new_dummy( DOUBLE, 2, 12, 3 );
  Array *v  = matrix_mean( a, 0 );
  Array *v2 = matrix_mean( a, 1 );

  /* array_print( a, -1, stdout ); */
  /* array_print( v, -1, stdout ); */
  /* array_print( v2, -1, stdout ); */

  bool isvec;
  vector_CHECK( isvec, v );
  fail_unless( isvec );
  fail_unless( v->size[0]==3 );
  fail_unless( cmpdouble( array_INDEX1( v, double, 0 ), 16.5, 3 )==0 );
  vector_CHECK( isvec, v2 );
  fail_unless( isvec );
  fail_unless( v2->size[0]==12 );
  fail_unless( cmpdouble( array_INDEX1( v2, double, 0 ), 1.0, 3 )==0 );
}
END_TEST

START_TEST (test_getrow)
{
  Array *a = array_new_dummy( DOUBLE, 2, 12, 13 );
  Array *v  = matrix_get_row( a, 0, TRUE );
  Array *v2 = matrix_get_row( a, 10, FALSE );

  /* array_print( a, -1, stdout ); */
  /* array_print( v, -1, stdout ); */
  /* array_print( v2, -1, stdout ); */

  bool isvec;
  vector_CHECK( isvec, v );
  fail_unless( isvec );
  fail_unless( v->free_data==TRUE );

  vector_CHECK( isvec, v2 );
  fail_unless( isvec );
  fail_unless( v2->free_data==FALSE );

  array_free( a );
  array_free( v );
  array_free( v2 );
}
END_TEST

START_TEST (test_getcol)
{
  Array *a = array_new_dummy( DOUBLE, 2, 12, 13 );
  Array *v  = matrix_get_col( a, 0 );
  Array *v2 = matrix_get_col( a, 10 );

  /* array_print( a, -1, stdout ); */
  /* array_print( v, -1, stdout ); */
  /* array_print( v2, -1, stdout ); */

  bool isvec;
  vector_CHECK( isvec, v );
  fail_unless( isvec );
  fail_unless( v->free_data==TRUE );

  vector_CHECK( isvec, v2 );
  fail_unless( isvec );
  fail_unless( v2->free_data==TRUE );

  array_free( a );
  array_free( v );
  array_free( v2 );
}
END_TEST

START_TEST (test_mmult)
{
  /*
	 compare with MATLAB code:
	 d1=3,d2=4,d3=5;
	 A=reshape(0:d1*d2-1, d2, d1)';
	 B=reshape(0:d2*d3-1, d3, d2)';
	 Z=(A*B)';
	 num2str(Z(:)', '%.2f, ')
	*/
  int d1=3,
	 d2=4,
	 d3=5;
  double correct[]={70.00, 76.00, 82.00, 88.00, 94.00, 190.00, 212.00, 
						  234.00, 256.00, 278.00, 310.00, 
						  348.00, 386.00, 424.00, 462.00};
  Array *a = array_new_dummy( DOUBLE, 2, d1, d2 );
  Array *b = array_new_dummy( DOUBLE, 2, d2, d3 );
  Array *c = matrix_mult( a, b );
  bool ismat;
  matrix_CHECK( ismat, c );

  int i;
  for( i=0; i<array_NUMEL( c ); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( c,double,i),
									 correct[i], 3 )==0 );
  } 

  /* array_print( a, -1, stdout ); */
  /* array_print( b, -1, stdout ); */
  /* array_print( c, -1, stdout ); */

  fail_unless( c->size[0]==d1 );
  fail_unless( c->size[1]==d3 );

  array_free( a );
  array_free( b );
  array_free( c );
}
END_TEST

START_TEST (test_transpose)
{
  int i,j;
  Array *a = array_new_dummy( DOUBLE, 2, 12, 3 );
  Array *b = matrix_transpose( a, TRUE );

  fail_unless( b->size[0]==a->size[1] );
  fail_unless( b->size[1]==a->size[0] );
  for( i=0; i<a->size[0]; i++ ){
	 for( j=0; j<a->size[1]; j++ ){

		fail_unless( cmpdouble( array_INDEX2( a,double,i,j),
										array_INDEX2( b,double,j,i), 3 )==0 );
	 }
  }

  /* array_print(a,-1,stdout); */
  /* array_print(b,-1,stdout); */

  matrix_transpose( a, FALSE );

  for( i=0; i<array_NUMEL(a); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( a,double,i),
									 array_INDEX1( b,double,i), 3 )==0 );
  }
  /* array_print(a,-1,stdout); */

  array_free( a );
  array_free( b );
}
END_TEST

START_TEST (test_pca)
{ 
  /*
	 compare with MATLAB code:
	 d1=12,d2=3;
	 X=reshape(0:d1*d2-1, d2, d1)';
	 [Y var P]=pca_svd( X' );
	 num2str(Y(:)', '%.2f, ')
	 num2str(var(:)', '%.2f, ')
	*/
  Array *a = array_new_dummy( DOUBLE, 2, 12, 3 );
  Array *vars, *pa;
  double correctX[]={-8.62, -0.00, -0.00, -7.05, -0.00, -0.00, -5.48, -0.00, -0.00, -3.92, -0.00, -0.00, -2.35, 0.00, -0.00, -0.78, -0.00, -0.00, 0.78, 0.00, 0.00, 2.35, 0.00, 0.00, 3.92, 0.00, 0.00, 5.48, 0.00, 0.00, 7.05, 0.00, 0.00, 8.62, 0.00, 0.00};
  double correctv[]={351.00, 0.00, 0.00};

  array_print( a, -1, stdout );

  pa = matrix_pca( a, &vars, TRUE );

  int i;
  for( i=0; i<array_NUMEL( pa ); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( pa,double,i),
									 correctX[i], 2 )==0,
					  "Y(%i): %f!=%f", i,array_INDEX1( pa,double,i),
					  correctX[i]);
  } 
  for( i=0; i<array_NUMEL( vars ); i++ ){
	 fail_unless( cmpdouble( array_INDEX1( vars,double,i),
									 correctv[i], 2 )==0 );
  } 

  array_print( pa, -1, stdout );
  array_print( vars, -1, stdout );

  array_free( a );
  array_free( vars );
  array_free( pa );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_linalg_suite (void){
  Suite *s = suite_create ("Linear Algebra-Functions");

  TCase *tc_core = tcase_create ("LinalgCore");
  tcase_add_test (tc_core, test_gsl );
  tcase_add_test (tc_core, test_mean );  
  tcase_add_test (tc_core, test_getrow );  
  tcase_add_test (tc_core, test_getcol );  
  tcase_add_test (tc_core, test_mmult );  
  tcase_add_test (tc_core, test_transpose );  

  tcase_add_test (tc_core, test_pca );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

