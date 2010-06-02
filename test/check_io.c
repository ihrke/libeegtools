/*
** check_io.c
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
#include "io.h"

START_TEST (test_arrtomatlab)
{
  Array *a=array_new_dummy( DOUBLE, 3, 3, 4, 5 );
  Array *b=array_new_dummy( DOUBLE, 2, 3, 4 );

  bool succ;
  succ=write_array_matlab( a, "a", "test.mat", FALSE);
  fail_if( succ );
  succ=write_array_matlab( b, "b", "test.mat", TRUE);
  fail_if( succ );


  Array *c=read_array_matlab( "test.mat", "a");
  Array *d=read_array_matlab( "test.mat", "b");

  ulong i;
  for( i=0; i<array_NUMEL(a); i++ ){
	 fail_unless( cmpdouble( array_INDEX1(a,double,i), 
									 array_INDEX1(c,double,i), 3)==0);
  }
  for( i=0; i<array_NUMEL(b); i++ ){
	 fail_unless( cmpdouble( array_INDEX1(b,double,i), 
									 array_INDEX1(d,double,i), 3)==0);
  }

  /* array_print( a, -1, stderr ); */
  /* array_print( c, -1, stderr ); */

  array_free( a );
  array_free( b );
  array_free( c );
  array_free( d );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_io_suite (void){
  Suite *s = suite_create ("IO-Functions");

  TCase *tc_core = tcase_create ("IOCore");
  tcase_add_test (tc_core, test_arrtomatlab);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

