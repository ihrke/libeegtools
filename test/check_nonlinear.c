/*
** check_nonlinear.c
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
#include "helper.h"
#include "mathadd.h"
#include "nonlinear.h"
#include "recurrence_plot.h"



   
#define TESTN 100
START_TEST (test_predict_simple)
{
  double x[TESTN], y[TESTN];
  int i;
  for( i=0; i<TESTN; i++){
	 x[i] = PI*i/(double)TESTN;
	 y[i] = sin(x[i]);
  }
  TimeDelayReconstruction *p = tdelay_init( 3, 1, y, TESTN );
  double sample[3]={ 0.5, 0.6, 0.7 };  
  dprintf("predict=%f\n", tdelay_predict_simple( p, sample, 2, 0.1 ) );

  tdelay_free( p );
}
END_TEST

START_TEST (test_prediction_error_simple)
{
  double x[TESTN], y[TESTN], z[TESTN];
  int i;
  for( i=0; i<TESTN; i++){
	 x[i] = 2*PI*i/(double)TESTN;
	 y[i] = sin(x[i]);
	 z[i] = cos(x[i]);
  }

  TimeDelayReconstruction *p = tdelay_init( 3, 1, y, TESTN );
  dprintf("crmse=%f\n", 
			 tdelay_simple_nonlinear_prediction_error( p, z, TESTN, 
																	  1, 0.1 ) );

  tdelay_free( p );
}
END_TEST

START_TEST (test_recplot)
{
  int i;
  Array *a=array_randunif( 0, 2, 100, 10 );
  Array *b=array_randunif( 0, 2, 100, 10 );
  Array *R=recplot( a, b, NULL, 0.5, NULL );
  bool ismatrix;
  matrix_CHECK(ismatrix,R);
  fail_unless( ismatrix );
  for( i=0; i<array_NUMEL(R); i++ ){
	 fail_unless( cmpdouble( array_INDEX1(R,double,i), 0.0, 3 )==0 ||
					  cmpdouble( array_INDEX1(R,double,i), 1.0, 3 )==0 );
  }

  Array *c=array_randunif( 0, 1, 100 );
  Array *d=array_randunif( 0, 1, 100 );
  R=recplot( c, d, R, 0.5, NULL );
  matrix_CHECK(ismatrix,R);
  fail_unless( ismatrix );
  for( i=0; i<array_NUMEL(R); i++ ){
	 fail_unless( cmpdouble( array_INDEX1(R,double,i), 0.0, 3 )==0 ||
					  cmpdouble( array_INDEX1(R,double,i), 1.0, 3 )==0 );
  }

  OptArgList *o = optarglist( "fan=int", 10 );
  R=recplot( a, b, R, 0.5, o );
  matrix_CHECK(ismatrix,R);
  fail_unless( ismatrix );
  for( i=0; i<array_NUMEL(R); i++ ){
	 fail_unless( cmpdouble( array_INDEX1(R,double,i), 0.0, 3 )==0 ||
					  cmpdouble( array_INDEX1(R,double,i), 1.0, 3 )==0 );
  }

  optarglist_free( o );
  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
  array_free( R );
}
END_TEST

START_TEST (test_loscalc_dtw)
{
	 int i, n=100;
	 Array *a=array_randunif( 0, 2, n, 10 );
	 Array *b=array_randunif( 0, 2, n, 10 );
	 Array *R=recplot( a, b, NULL, 0.5, NULL );
	 bool ismatrix;
	 matrix_CHECK(ismatrix,R);
	 fail_unless( ismatrix );

	 Array *los=recplot_los_dtw_noise( R );
	 bool ispath;
	 warppath_CHECK(ispath,los);
	 fail_unless( ispath );

	 /* start/end condition */
	 fail_unless( array_INDEX2( los, uint, 0,0)==0 );
	 fail_unless( array_INDEX2( los, uint, 1,0)==0 );
	 fail_unless( array_INDEX2( los, uint, 0,los->size[1]-1)==n-1 );
	 fail_unless( array_INDEX2( los, uint, 1,los->size[1]-1)==n-1 );

	 for( i=1; i<los->size[1]; i++ ){ /* monotonicity test */
		fail_unless( array_INDEX2( los, uint, 0, i )>=
						 array_INDEX2( los, uint, 0, i-1 ) );
		fail_unless( array_INDEX2( los, uint, 1, i )>=
						 array_INDEX2( los, uint, 1, i-1 ) );
	 }
	 /* check for no gaps */
	 for( i=1; i<los->size[1]; i++ ){
		fail_if( array_INDEX2( los, uint, 0, i )-
					array_INDEX2( los, uint, 0, i-1 )>1 );
		fail_if( array_INDEX2( los, uint, 1, i )-
					array_INDEX2( los, uint, 1, i-1 )>1 );
		fail_if( array_INDEX2( los, uint, 0, i )-
					array_INDEX2( los, uint, 0, i-1 )==0 &&
					array_INDEX2( los, uint, 1, i )-
					array_INDEX2( los, uint, 1, i-1 )==0 );
	 }

	 array_free( a );
	 array_free( b );
	 array_free( R );
	 array_free( los );
}
END_TEST

START_TEST (test_loscalc_marwan)
{
	 int i, n=100;
	 Array *a=array_randunif( 0, 2, n, 10 );
	 Array *b=array_randunif( 0, 2, n, 10 );
	 Array *R=recplot( a, b, NULL, 0.5, NULL );
	 bool ismatrix;
	 matrix_CHECK(ismatrix,R);
	 fail_unless( ismatrix );

	 Array *los=recplot_los_marwan( R, 2, 2 );
	 bool ispath;
	 warppath_CHECK(ispath,los);
	 fail_unless( ispath );

	 /* start/end condition */
	 fail_unless( array_INDEX2( los, uint, 0,0)==0 );
	 fail_unless( array_INDEX2( los, uint, 1,0)==0 );
	 fail_unless( array_INDEX2( los, uint, 0,los->size[1]-1)==n-1 );
	 fail_unless( array_INDEX2( los, uint, 1,los->size[1]-1)==n-1 );

	 for( i=1; i<los->size[1]; i++ ){ /* monotonicity test */
		  fail_unless( array_INDEX2( los, uint, 0, i )>=
							array_INDEX2( los, uint, 0, i-1 ) );
		  fail_unless( array_INDEX2( los, uint, 1, i )>=
							array_INDEX2( los, uint, 1, i-1 ) );
	 }
	 /* check for no gaps */
	 for( i=1; i<los->size[1]; i++ ){
		  fail_if( array_INDEX2( los, uint, 0, i )-
					  array_INDEX2( los, uint, 0, i-1 )>1 );
		  fail_if( array_INDEX2( los, uint, 1, i )-
					  array_INDEX2( los, uint, 1, i-1 )>1 );
		  fail_if( array_INDEX2( los, uint, 0, i )-
					  array_INDEX2( los, uint, 0, i-1 )==0 &&
					  array_INDEX2( los, uint, 1, i )-
					  array_INDEX2( los, uint, 1, i-1 )==0 );
	 }

	 array_free( a );
	 array_free( b );
	 array_free( R );
	 array_free( los );
}
END_TEST
/* template
START_TEST (test_)
{
}
END_TEST
*/
Suite * init_nonlinear_suite (void){
  Suite *s = suite_create ("Nonlinear-Functions");

  TCase *tc_core = tcase_create ("NonlinearCore");
  tcase_add_test (tc_core, test_predict_simple);
  tcase_add_test (tc_core, test_prediction_error_simple);
  tcase_add_test (tc_core, test_recplot);
  tcase_add_test (tc_core, test_loscalc_dtw);
  tcase_add_test (tc_core, test_loscalc_marwan);


  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

