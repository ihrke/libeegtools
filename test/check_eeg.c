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

#include <stdlib.h>
#include <check.h>
#include "check_all.h"
#include <math.h>
#include "definitions.h"
#include "eeg.h"


/* template
START_TEST (test_)
{
}
END_TEST
*/

Suite * init_eeg_suite (void){
  Suite *s = suite_create ("EEG-Functions");

  TCase *tc_core = tcase_create ("EEGCore");
  //  tcase_add_test (tc_core, test_ );

  tcase_set_timeout(tc_core, 20);
  suite_add_tcase (s, tc_core);
  return s;
}

