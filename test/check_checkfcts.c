/*
** check_checkfcts.c
** 
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Thu Oct 11 13:18:26 2007 Matthias Ihrke
** Last update Thu Oct 11 13:41:15 2007 Matthias Ihrke
*/

#include<check.h>
#include "checkfcts.h"
#include "check_all.h"

START_TEST (test_cmpdouble)
{
  fail_if(cmpdouble(1.0, 1.0, 1));
  fail_if(cmpdouble(1.0, 1.02, 1));
  fail_unless(cmpdouble(1.0, 1.02, 2)==-1);
  fail_unless(cmpdouble(-1.000001, -2.0, 3)==1); 
  fail_if(cmpdouble(-1.000001, -1.000000002, 3));
  fail_unless(cmpdouble(-1.000001, -1.000000002, 10));
}
END_TEST

START_TEST (test_isequal_doublearray_binfile)
{
  double d[20] = {-1.1905, -1.2245, -1.0235, 0.5122, -0.7338,
		  -0.4516, -0.4515, -0.1916, -0.1333, 1.3741,
		  -1.2610, -1.4282, -1.6060, -1.3679, 0.2498,
		  -0.2019, -0.7183, -0.5378, -0.6745, 0.0662};
  fail_unless(isequal_doublearray_binfile(d, 20, "testdata/bindoubles1.dat", 4));
  d[0] = -1.1904;
  fail_if    (isequal_doublearray_binfile(d, 20, "testdata/bindoubles1.dat", 5));
}
END_TEST

Suite * init_checkfct_suite (void){
  Suite *s = suite_create ("Check Functions");

  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_cmpdouble);
  tcase_add_test (tc_core, test_isequal_doublearray_binfile);
  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

