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
  Suite *cfct = init_checkfct_suite ();
  Suite *den  = init_denoising_suite();
  Suite *dist = init_distance_suite ();
  Suite *hmm  = init_hmm_suite      ();
  Suite *list = init_list_suite     ();
  Suite *other= init_other_suite    ();
  Suite *nlin = init_nonlinear_suite();
  Suite *arr  = init_array_suite    ();


  SRunner *sr = srunner_create (cfct);
  srunner_add_suite( sr, den  );
  srunner_add_suite( sr, dist );
  srunner_add_suite( sr, hmm  );
  srunner_add_suite( sr, list );
  srunner_add_suite( sr, other);
  srunner_add_suite( sr, nlin );
  srunner_add_suite( sr, arr  );


  srunner_set_log(sr, "test.log");
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);

  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
