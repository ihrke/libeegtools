/** \file t_dtw.c
 *
 * \brief Testing DTW with restrictions
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

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
#define PL(code) (code)
#else
#define PL(code)
#endif

void usage(const char *);

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  FILE *f;
  int c, i, j, k, chan1, chan2;
  int n1, n2;
  EEGdata_trials *eeg;
  WarpPath *P;
  double **dist;
  double K;
  double pos[2], size[2];
  const char *cmap = "hot";
  chan1 = 0;
  chan2 = 0;

  if(argc<3){
	 usage(argv[0]);
  }
  fprintf( stderr, "reading file          : %s\n", argv[1] );
  fprintf( stderr, "restriction parameter : %f\n", K=atof( argv[2] ) );
  fprintf( stderr, "using electrode       : (%i, %i)\n", chan1, chan2 );

  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);

  dist = matrix_init(eeg->data[0]->n, eeg->data[1]->n);

  for( i=0; i<2; i++){
	 dprintf("TRIAL: i=%i\n", i);
	 n2 = eeg->data[i  ]->n;
	 n1 = (eeg->data[i+1]->n)-100;
	 
	 if(i==0){
		P = init_warppath(n1, n2);
	 }

	 dist = DTW_build_distmatrix( eeg->data[i  ]->d[chan1], n1, 
											eeg->data[i+1]->d[chan2], n2, 1.0, 1.0, dist);
	 pos[0] = -n1-150; pos[1] = 0;
	 size[0] = n1;  size[1] = n2;
	 //	 dist[100][10]=100000;
	 PL( plot_image_position( dist, n1, n2, cmap, pos, size) );

	 dist = DTW_build_restricted_cumdistmatrix( eeg->data[i  ]->d[chan1], n1, 
															  eeg->data[i+1]->d[chan2], n2, 
															  K, dist );  
	 PL( plot_image( dist, n1, n2, cmap) );
	 P = DTW_path_from_cumdistmatrix(dist, n1, n2, P);

	 PL( plot_format_int( P->upath, P->spath, (P->J)+(P->K), "b." ) );  

	 double pos[4] = {0, -100, 600, -1};
	 double size[4]= {0, 600, -50, +50 };
	 int sb;
	 PL( sb = plot_subplot_create(pos, size) );
	 PL( plot_subplot_select(sb) );
	 PL( plot_format(NULL, eeg->data[i], n1, "r.") );
	 PL( plot_format(NULL, eeg->data[i+1], n2, "b.") );
	 PL( plot_subplot_select(-1) );

	 PL( plot_add() );
	 dprintf("(n1,n2) = (%i,%i)\n", n1, n2);
  }
  
  PL( plot_show() );

  free_eegdata_trials( eeg );
  matrix_free(dist, eeg->data[0]->n);
  free_warppath(P);

  return 0;
}


 void usage(const char *p){
	fprintf(stdout, "Usage: %s <eeg raw-file> K\n", p);
	exit(-1);
 }
