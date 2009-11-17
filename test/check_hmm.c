/*
** check_hmm.c
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
#include "hmm.h"

     
START_TEST (test_hmm_alloc_cphmm)
{
  CPHiddenMarkovModel *m =cphmm_alloc( 10, 1000, 2000, 10, 5 );
  fail_if( m->K!=10 );
  fail_if( m->n!=1000 );
  fail_if( m->M!=2000 );
  fail_if( m->Q!=10 );
  fail_if( m->J!=5 );
  m->q[0] = 1;  m->q[m->Q-1] = 10;
  m->u[0] = 1;  m->u[m->K-1] = 10;
  m->d[0] = 1;  m->d[m->J-1] = 10;
  m->z[0] = 1;  m->z[m->M-1] = 10;

  m->tau[m->K-1][m->n-1]=1;
  m->phi[m->K-1][m->n-1]=1;
}
END_TEST



Suite * init_hmm_suite (void){
  Suite *s = suite_create ("HMM-Functions");

  TCase *tc_core = tcase_create ("HMMcore");
  tcase_add_test (tc_core, test_hmm_alloc_cphmm);

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}
