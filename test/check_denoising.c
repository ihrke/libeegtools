/*
** check_denoising.c
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
#include "check_all.h"
#include "helper.h"
#include "mathadd.h"
#include "filter.h"

     
START_TEST (test_sampled_line)
{
  double d[10], diff;
  int i;
  sampled_line(d, 10, -10.0, 12.0);
  diff=round((double)ABS(d[1]-d[0]));

  for(i=1; i<10; i++){
    fail_if(round(fabs(d[i]-d[i-1]))!=diff);
  }
  fail_if(d[0]!=-10);
  fail_if(d[9]!=12);
}
END_TEST

START_TEST (test_filter)
{
	char *script="load filter.mat\n"
				 "plot( a, 'b' )\n"
				 "hold on\n"
				 "plot( b, 'r' )\n";
	int n=10000;
	Array *a=array_randunif( 3423, 1, n );
	Array *b=filter( a, "LpBu4/4", 1000, true );

	fail_unless( b );

	write_array_matlab( a, "a", CHECKDATADIR"filter.mat", false);
	write_array_matlab( b, "b", CHECKDATADIR"filter.mat", true);
	write_script( CHECKDATADIR"test_filter.m", script);

	array_free( a ) ;
	array_free( b );
}
END_TEST



Suite * init_denoising_suite (void){
  Suite *s = suite_create ("Denoising Functions");

  TCase *tc_core = tcase_create ("DenoiseCore");
  tcase_add_test (tc_core, test_sampled_line);
  tcase_add_test (tc_core, test_filter);
  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

