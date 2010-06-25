/*
** check_optarg.c
** 
** Unit-Test suite using libcheck
**
**
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
*/

#include <stdlib.h>
#include <check.h>
#include "check_all.h"
#include "optarg.h"
#include <math.h>
     
START_TEST (test_optarg_init)
{
  double *test=(double*)malloc(10*sizeof(double));
  OptArgList *opts = optarglist( " hi =double , test = double *, du=int ", 0.5, test, 10);
  optarglist_print( opts, stderr );

  fail_if( strcmp( opts->args[0].key, "hi" ) );
  fail_if( strcmp( opts->args[0].type, "double" ), "'%s' not '%s'", opts->args[0].type, "double");
  fail_if( cmpdouble( opts->args[0].data_scalar, 0.5, 2 ) );
  fail_if( strcmp( opts->args[2].key, "du" ) );
  fail_if( strcmp( opts->args[2].type, "int" ) );
  fail_if( ((int) opts->args[2].data_scalar )!=10 );
  fail_if( strcmp( opts->args[1].key, "test" ) );
  fail_if( strcmp( opts->args[1].type, "double*" ) );
  fail_if( opts->args[1].data_ptr != test );

  optarglist_free( opts );

  opts = optarglist( "hi=int", 5 );
  optarglist_print( opts, stderr );
  fail_if( strcmp( opts->args[0].key, "hi" ) );


  free(test);
}
END_TEST

     
START_TEST (test_optarg_get_by_key)
{
  double *test=(double*)malloc(10*sizeof(double));
  OptArgList *opts = optarglist( "hi =double,du=int , test = double *  ",0.5,10,test);

  OptArg *arg=optarglist_optarg_by_key( opts, "test" );

  fail_if( strcmp( arg->type, "double*" ) );
  fail_if( arg->scalar );
  fail_if( arg->data_ptr!=test );

  double *test2;
  test2 = optarglist_ptr_by_key( opts, "test" );
  fail_if( test2!=test );

  double x = optarglist_scalar_by_key( opts, "gibtsnich");
  fail_if( !isnan( x ) );

  x = optarglist_scalar_by_key( opts, "hi");
  fail_if( isnan( x ) );
  fail_if( cmpdouble( x, 0.5, 2 ) );
  optarglist_free( opts );
  free(test);
}
END_TEST

START_TEST (test_optarg_append)
{
  double *test=(double*)malloc(10*sizeof(double));
  OptArgList *opts = optarglist( "hi =double,du=int , test = double *  ",0.5,10,test);
  optarglist_print( opts, stderr );

  OptArg *new1 = optarg_scalar( "newone", 10.0 );
  OptArg *new2 = optarg_ptr   ( "newtwo", test );
  OptArgList *newo  = optarglist_append_arg( opts, new1 );
  OptArgList *newo2 = optarglist_append_arg( newo, new2 );

  optarglist_print( newo, stderr );
  optarglist_print( newo2, stderr );


  optarglist_free( opts );
  optarglist_free( newo );
  optarglist_free( newo2 );
  free( new1 );
  free( new2 );
  free(test);
}
END_TEST

START_TEST (test_parse_scalar)
{
  OptArgList *opts=optarglist( "test=long,"
										 "test2=void*", (long)6, NULL );
  long test=0;
  double x;
  optarg_PARSE_SCALAR( opts, "test", test, long, x );
  fail_if( test!=6 );
  test = 0;
  optarg_PARSE_SCALAR( opts, "test3", test, long, x );
  fail_if( test!=0 );  
  optarg_PARSE_SCALAR( opts, "test2", test, long, x );
  fail_if( test!=0 );
}
END_TEST

START_TEST (test_parse_ptr)
{
  int test[10];
  double *test2=(double*)malloc( 10*sizeof(double));
  OptArgList *opts=optarglist( "test=void*,"
										 "test2=void*", (void*)test,
										 (void*) test2);

  void *ptr;
  int *t1;
  double *t2;
  optarg_PARSE_PTR( opts, "test", t1, int*, ptr );
  fail_if( t1!=test );
  optarg_PARSE_PTR( opts, "test3", t1, int*, ptr );
  fail_if( test!=t1 );  

  optarg_PARSE_PTR( opts, "test2", t2, double*, ptr );
  fail_if( t2!=test2 );

  free( test2 );
}
END_TEST

START_TEST (test_delete_arg)
{
  int test[10];
  double *test2=(double*)malloc( 10*sizeof(double));
  OptArgList *opts=optarglist( "test=void*,"
										 "test2=void*,"
										 "dritter=int", (void*)test,
										 (void*) test2, 3);
  OptArg *tmp=optarglist_optarg_by_key( opts, "dritter" );
  OptArgList *o=optarglist_delete_arg( opts, tmp );
  
  fail_unless( o->nargs == opts->nargs-1 );
  fail_if(  optarglist_has_key( o, tmp->key ) );

  optarglist_free( opts );
  optarglist_free( o );
}
END_TEST

START_TEST (test_inittmp)
{  
  double *test=(double*)malloc(10*sizeof(double));
  OptArgList *opts = optarglisttmp( " hi =double , test = double *, du=int ", 0.5, test, 10);

  fail_if( strcmp( opts->args[0].key, "hi" ) );
  fail_if( strcmp( opts->args[0].type, "double" ), "'%s' not '%s'", opts->args[0].type, "double");
  fail_if( cmpdouble( opts->args[0].data_scalar, 0.5, 2 ) );
  fail_if( strcmp( opts->args[2].key, "du" ) );
  fail_if( strcmp( opts->args[2].type, "int" ) );
  fail_if( ((int) opts->args[2].data_scalar )!=10 );
  fail_if( strcmp( opts->args[1].key, "test" ) );
  fail_if( strcmp( opts->args[1].type, "double*" ) );
  fail_if( opts->args[1].data_ptr != test );

  optarglist_free( opts );

  opts = optarglisttmp( "hi=int", 5 );
  fail_if( strcmp( opts->args[0].key, "hi" ) );
  optarglist_free( opts );

  free(test);
}
END_TEST

START_TEST (test_removeflag)
{
  double *test=(double*)malloc(10*sizeof(double));
  OptArgList *opts = optarglisttmp( " hi =double , test = double *, du=int ", 0.5, test, 10);

  //  optarglist_print( opts, stderr );  
  bool flag;
  OptArgList *opts2=optarglist_remove_freeflag( opts, &flag );
  fail_unless( flag );
  fail_unless( opts2!=opts );
  optarglist_free( opts2 );

  opts = optarglist( " hi =double , test = double *, du=int ", 0.5, test, 10);
  opts2=optarglist_remove_freeflag( opts, &flag );
  fail_if( flag );
  fail_if( opts2!=opts );
  optarglist_free( opts2 );

  free( test );
}
END_TEST

/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_optarg_suite (void){
  Suite *s = suite_create ("Optarg-Functions");

  TCase *tc_core = tcase_create ("OptargCore");
  tcase_add_test (tc_core, test_optarg_init );
  tcase_add_test (tc_core, test_optarg_get_by_key );
  tcase_add_test (tc_core, test_optarg_append );
  tcase_add_test (tc_core, test_parse_scalar);
  tcase_add_test (tc_core, test_parse_ptr);
  tcase_add_test (tc_core, test_delete_arg);
  tcase_add_test (tc_core, test_inittmp);
  tcase_add_test (tc_core, test_removeflag);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

