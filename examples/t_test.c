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

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
#define PL(code) (code)
#else
#define PL(code)
#endif

void usage(const char *);

void dtwtest1(int argc, const char **argv);

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  dtwtest1(argc, argv);

  return 0;
}


 void usage(const char *p){
	fprintf(stdout, "Usage: %s <eeg raw-file> K\n", p);
	exit(-1);
 }


void dtwtest1(int argc, const char **argv){
 FILE *f;
  double *t;
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
  t=NULL;


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

  for( i=0; i<1; i++){
	 n1 = eeg->markers[i  ][1];
	 n2 = eeg->markers[i+1][1];
	 dprintf("TRIAL: i=%i, (n1,n2)=(%i,%i)\n", i, n1,n2);

	 if(i==0){
		P = init_warppath(NULL, n1, n2);
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


	 t = ADTW_from_path( eeg->data[i]->d[chan1], n1, eeg->data[i+1]->d[chan2], n2, P, t );

	 double pos[4] = {0, -100, 500, -1};
	 double size[4]= {0, 500, -50, +50 };
	 int sb;
	 PL( sb = plot_subplot_create(pos, size) );
	 PL( plot_subplot_select(sb) );
	 PL( plot_format(NULL, eeg->data[i]->d[chan1], n1, "r-") );
	 PL( plot_format(NULL, eeg->data[i+1]->d[chan2], n2, "b-") );
	 PL( plot_format(NULL, t, (n1+n2)/2, "g-") );
	 PL( plot_subplot_select(-1) );

	 PL( plot_add() );
	 dprintf("(n1,n2) = (%i,%i)\n", n1, n2);
	 free(t);
	 t=NULL;
  }
  
  PL( plot_show() );

  free_eegdata_trials( eeg );
  matrix_free(dist, eeg->data[0]->n);
  free_warppath(P);

}
