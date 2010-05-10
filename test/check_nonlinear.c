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
  PhaseSpace *p = phspace_init( 3, 1, y, TESTN );
  double sample[3]={ 0.5, 0.6, 0.7 };  
  dprintf("predict=%f\n", phspace_predict_simple( p, sample, 2, 0.1 ) );

  phspace_free( p );
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

  PhaseSpace *p = phspace_init( 3, 1, y, TESTN );
  dprintf("crmse=%f\n", 
			 phspace_simple_nonlinear_prediction_error( p, z, TESTN, 
																	  1, 0.1 ) );

  phspace_free( p );
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

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

