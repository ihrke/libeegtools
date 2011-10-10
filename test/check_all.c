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
  Suite *list = init_list_suite     ();
  Suite *other= init_other_suite    ();
  Suite *nlin = init_nonlinear_suite();
  Suite *arr  = init_array_suite    ();
  Suite *lalg = init_linalg_suite   ();
  Suite *warp = init_warping_suite  ();
  Suite *clust= init_cluster_suite  ();
  Suite *avg  = init_average_suite  ();
  Suite *io   = init_io_suite  ();
  Suite *optarg=init_optarg_suite  ();
  Suite *eeg  = init_eeg_suite  ();


  SRunner *sr = srunner_create (cfct);
  srunner_add_suite( sr, dist );
/*
  srunner_add_suite( sr, lalg  );
  srunner_add_suite( sr, warp  );
  srunner_add_suite( sr, clust  );
  srunner_add_suite( sr, avg  );
  srunner_add_suite( sr, list );
  srunner_add_suite( sr, eeg  );
  srunner_add_suite( sr, optarg );
  srunner_add_suite( sr, io  ); 
  srunner_add_suite( sr, nlin );
  srunner_add_suite( sr, other);
*/
  srunner_add_suite( sr, den  );
  srunner_add_suite( sr, arr  );


#ifdef EXPERIMENTAL
  Suite *hmm  = init_hmm_suite      ();
  /* srunner_add_suite( sr, hmm  ); */
#endif


  srunner_set_log(sr, "test.log");
  srunner_run_all (sr, CK_VERBOSE);
  number_failed = srunner_ntests_failed (sr);

  srunner_free (sr);
  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
