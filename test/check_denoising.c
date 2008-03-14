/*
** check_denoising.c
** 
** Unit-Test suite using libcheck
**
**  use DB() to generate output
**  use ML() to plot with matlab
**  use DB( ML ( code ) ) to conditionally plot
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
#include "check_all.h"
#include "helper.h"

#define DB(code) code
     
START_TEST (test_sampled_line)
{
  double d[10], diff;
  int i;
  sampled_line(d, 10, -10.0, 12.0);
  diff=round(fabs(d[1]-d[0]));
  DB( vprint_vector("d", d, 10) );

  for(i=1; i<10; i++){
    DB( printf("%f ", d[i]-d[i-1]) );
    fail_if(round(fabs(d[i]-d[i-1]))!=diff);
  }
  fail_if(d[0]!=-10);
  fail_if(d[9]!=12);
}
END_TEST


START_TEST (test_lininterp)
{
  double s[8] = {3.0, 7.0, 3.0, 7.0, 9.0, 3.0, 2.0, 1.0};
  double r[8];
  double snew[20];
  double rnew[20];

  sampled_line(r, 8, 0, 10);
  sampled_line(rnew, 20, 1, 8);
  lininterp(r, s, 8, rnew, snew, 20);
  DB(
     ML( Engine *m;
	 m = ml_init();
	 ml_plot(m, r, s, 8, "r", 1); 
	 engEvalString(m, "hold on;");
	 ml_plot(m, rnew, snew, 20, "bo", 0);
	 ml_wait(m);
	 ml_close(m);
	 );
     );
}
END_TEST


Suite * init_denoising_suite (void){
  Suite *s = suite_create ("Denoising Functions");

  TCase *tc_core = tcase_create ("Core");
  tcase_add_test (tc_core, test_sampled_line);
  tcase_add_test (tc_core, test_lininterp);
  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

