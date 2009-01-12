/** \example t_adtw_markers.c
 *
 * \brief Computing the ADTW of the first two trials within a raw-file.
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
#include "helper.h"
#include "definitions.h"
#include "clustering.h"
#include "averaging.h"
#include "mathadd.h"
#include "warping.h"

#include <libplotter/cplotter.h>

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  FILE *f;
  char buffer[255];
  EEGdata_trials *eeg;
  EEGdata *new, *s1, *s2;
  double **G, **d;
  double maxdist;
  WarpPath *P;
  int n, i, channel;
  
  /* get data */				  
  oprintf("Reading from '%s'\n", argv[1]);
  eeg=read_eegtrials_from_raw( argv[1] );
  print_eegdata_trials(stderr, eeg);
 
  n = eeg->data[0]->n;
  channel = 0;
  s1 = eeg->data[0];
  s2 = eeg->data[1];

  /* computation */  
  G = matrix_init(n,n);
  d = matrix_init(n,n);

  dprintf(" init d,G\n");
  d = eeg_distmatrix_euclidean_derivative_channel(  s1, s2, 0, d, 1,1 );
  dprintf(" d\n");
  G = eeg_gaussian_regularization_bresenham_channel( s1, s2, 0, atof(argv[2]), G);

  dprintf(" computed d_ij and G_f\n");
  maxdist = matrix_max( (const double**)d, n, n, NULL, NULL );
  scalar_minus_matrix( maxdist, d, n, n );
  matrix_dottimes_matrix( d, (const double**)G, n, n );
  scalar_minus_matrix( maxdist, d, n, n );
  
  P = DTW_path_from_square_distmatrix( d, n, ALLOC_IN_FCT );
  new = eeg_ADTW_from_path( s1, s2, ALLOC_IN_FCT, 0, P );

  plot_image( d, n, n, "hot" ); 
  plot_format_int( P->upath, P->spath, P->J+P->K, "b-");

  plot_switch( plot_add( ) ); 
  plot_format( eeg->times, s1->d[0], n, "g");
  plot_format( eeg->times, s2->d[0], n, "b");
  plot_format( eeg->times, new->d[0], n, "r");


  plot_switch( plot_add( ) );
  plot_image( G, n, n, "hot" ); 
  plot_show(); 
  
  write_double_matrix_ascii_file( "test.dst", (const double**) G, n, n );

  /* cleaning up */
  free_eegdata_trials( eeg );
  free_eegdata( new );
  free_warppath( P );

  for( i=0; i<n; i++ ){
	 free( d[i] );
	 free( G[i] );
  }
  free( d );
  free( G );
  
  return 0; 
}
