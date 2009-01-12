/** \example t_test.c
 *
 * \brief Testing ground
 *
 * Compilation:
 *\code
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "reader.h"
#include "writer.h"
#include "definitions.h"
#include "helper.h"
#include "clustering.h"
#include "warping.h"
#include <libplotter/cplotter.h>
#include <gsl/gsl_spline.h>
#define PL(code)

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int padtw_test(int argc, char **argv);
int interpolation_test(int argc, char **argv);

int main(int argc, char **argv){ 
  dprintf("in t_test\n");
  /* return padtw_test(argc, argv); */
  return interpolation_test(argc, argv);
}

int interpolation_test(int argc, char **argv){
  EEGdata_trials *eeg;
  EEGdata *new;
  double **Delta;
  int N;
  double *s1, *s2;
  int n, newn;

  /* get data */				  
  oprintf("Reading from '%s'\n", argv[1]);
  eeg=read_eegtrials_from_raw( argv[1] );
  print_eegdata_trials(stderr, eeg);
  n = eeg->nsamples;
  newn = atoi(argv[2]);
  s1 = eeg->data[0]->d[0];

  s2 = resample_gsl( s1, n, newn, ALLOC_IN_FCT, gsl_interp_cspline );

  plot_format( NULL, s1, n, "r-");
  plot_format( NULL, s2, newn, "b.");
  plot_show();

  return 0;
}

int padtw_test(int argc, char **argv){
  EEGdata_trials *eeg;
  EEGdata *new;
  double **Delta;
  int N;

 
  /* get data */				  
  oprintf("Reading from '%s'\n", argv[1]);
  eeg=read_eegtrials_from_raw( argv[1] );
  print_eegdata_trials(stderr, eeg);
  N = eeg->ntrials;

  Delta = eegtrials_diffmatrix_channel( eeg, clusterdist_euclidean_pointwise, 0, ALLOC_IN_FCT);

  /* DO IT */
  new = init_eegdata( eeg->data[0]->nbchan, eeg->nsamples, eeg->nmarkers_per_trial );
  eegtrials_PADTW( eeg, Delta, N, 1, new );
  write_eegdata_ascii_file( "test.out", new );


  /* cleaning up */
  free_eegdata_trials( eeg );
  free_eegdata( new );
  matrix_free( Delta, N );

  return 0;
}
