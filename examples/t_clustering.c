/** \example t_clustering.c
 * Testing clustering
 *
 * Compilation:
 *\code
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "reader.h"
#include "definitions.h"
#include "clustering.h"
#include "helper.h"
#include <libgen.h> /* basename */

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
#define PL(code) (code)
#else
#define PL(code)
#endif
/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  char filename[255];
  int i, j, chan;
  EEGdata_trials *eeg;
  double **dist;
  double(*distmetric)(EEGdata*,EEGdata*,int);
  chan = 36;

  fprintf(stderr, "reading file   : %s\n", argv[1]);
  fprintf(stderr, "using distance : %s\n", argv[2]);
  fprintf(stderr, "using electrode: %i\n", chan);

  if(!strcmp(argv[2], "tw")){
	 fprintf(stderr, "... recognized TW-distance\n");
	 distmetric=clusterdist_tw_complete;
  } else if(!strcmp(argv[2], "euclid")){
	 fprintf(stderr, "... recognized euclidean-distance\n");
	 distmetric=clusterdist_euclidean_pointwise;
  } else {
	 fprintf(stderr, "... metric not recognized, aborting\n");
	 return -1;
  }
  sprintf( filename, "clustdat_%s_%s.tsv", basename( argv[1] ), argv[2] );
  fprintf( stderr, " writing file: %s\n", filename);
  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);

  /* get distances */
  dist = eegtrials_diffmatrix_channel( eeg, distmetric, chan );
  f = fopen(filename, "w");
  for(i=0; i<eeg->ntrials; i++){
	 for(j=0; j<eeg->ntrials; j++){
		fprintf(f, "%.2f\t", dist[i][j]);
	 }
	 fprintf(f, "\n");
  }
  fclose(f);

  for(i=0; i<eeg->ntrials; i++){
	 PL( plot_format(eeg->times, (eeg->data[i])->d[chan], (eeg->data[i])->n, "r") );
  } 								 

  PL(	plot_show() );
  free_eegdata_trials(eeg);
  return 0;
}
