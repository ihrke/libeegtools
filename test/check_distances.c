/*
** check_distances.c
** 
** Unit-Test suite using libcheck
**
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Wed Oct 10 20:28:00 2007 Matthias Ihrke
** Last update Thu Oct 11 13:30:19 2007 Matthias Ihrke
*/

#include<stdlib.h>
#include<check.h>
#include<math.h>
#include "array.h"
#include "distances.h"
#include "check_all.h"
#include "helper.h"
#include "mathadd.h"

     
START_TEST (test_dist_point_line)
{
  double p[2]={ -0.40, 10.32};
  double x[2]={ -12.1, 3.10 };
  double y[2]={ 100.1, -10.0};
 
  fail_if( cmpdouble( 8.5281, dist_point_line( p,x,y ), 4 ),
			  "8.5281 (MATLAB) != '%f'", dist_point_line( p,x,y ));
  x[0]=10.32;
  fail_if( cmpdouble( 5.5966, dist_point_line( p,x,y ), 4 ),
			  "5.5966 (MATLAB) != '%f'", dist_point_line( p,x,y ));
  y[1] = 1000.0;
  fail_if( cmpdouble( 11.3244, dist_point_line( p,x,y ), 4 ),
			  "11.3244 (MATLAB) != '%f'", dist_point_line( p,x,y ));
  
}
END_TEST


START_TEST (test_sdistmat_euclidean)
{
  int n=10;
  Array *a=array_new_dummy( DOUBLE, 1, n );
  Array *b=array_new_dummy( DOUBLE, 1, n );
  array_reverse(b);
  Array *c=distmatrix_signaldist( vectordist_euclidean,
											 a, b, NULL, NULL );
  Array *d=array_new2( DOUBLE, 2, n,n );
  d = distmatrix_signaldist( vectordist_euclidean,
									  a, b, d, NULL );
  int i;
  for( i=0; i<array_NUMEL( c ); i++ ){
	 fail_if( array_INDEX1(c, double, i) < 0 );
	 fail_if( array_INDEX1(d, double, i) < 0 );
	 fail_unless( array_INDEX1(c, double, i)==array_INDEX1(d, double, i) );
  }
  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
}
END_TEST

START_TEST (test_vectordistmatrix)
{ 
  Array *a=array_new_dummy( DOUBLE, 2, 100, 20 );
  Array *b=matrix_distmatrix( vectordist_euclidean, a, NULL, NULL ) ;  
  Array *c=array_new2( DOUBLE, 2, 100, 100 );
  c=matrix_distmatrix( vectordist_euclidean,
							a, c, NULL );
  int i;
  for( i=0; i<array_NUMEL( b ); i++ ){
	 fail_if( array_INDEX1(b, double, i) < 0 );
	 fail_if( array_INDEX1(c, double, i) < 0 );
	 fail_unless( array_INDEX1(b, double, i)==array_INDEX1(c, double, i) );
  }
  array_free( a );
  array_free( b );
  array_free( c );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/


Suite * init_distance_suite (void){
  Suite *s = suite_create ("Distance-Functions");

  TCase *tc_core = tcase_create ("DistanceCore");
  tcase_add_test (tc_core, test_dist_point_line);
  tcase_add_test (tc_core, test_sdistmat_euclidean);
  tcase_add_test (tc_core, test_vectordistmatrix);


  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

