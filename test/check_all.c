/*
** check_all.c
** 
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Thu Oct 11 13:26:29 2007 Matthias Ihrke
** Last update Thu Oct 11 14:04:09 2007 Matthias Ihrke
*/

#include "check_all.h"

int main (void){
  int number_failed;
  Suite *s = init_denoising_suite ();
  Suite *s2= init_checkfct_suite ();

  SRunner *sr = srunner_create (s);
  srunner_add_suite(sr, s2);
  srunner_set_log(sr, "test.log");
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);
  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
