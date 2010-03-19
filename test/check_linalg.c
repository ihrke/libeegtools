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

START_TEST (test_pca)
{
  Array *a = array_new_dummy( DOUBLE, 2, 12, 3 );
  Array *vars, *p;
  p = matrix_pca( a, &vars, TRUE );
  array_print( p, -1, stdout );
  array_print( vars, -1, stdout ); 
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

  tcase_add_test (tc_core, test_pca );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

