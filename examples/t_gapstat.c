/** \example t_gapstat.c
 *
 * \brief Testing gap statistic
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

void test_cmpclust();
void test_dummycluster();

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  FILE *f;
  int i, j, k, chan;
  EEGdata_trials *eeg;
  double **dist;
  double(*distmetric)(EEGdata*,EEGdata*,int);
  Clusters *C;
  int maxk=10;
  double xWk[200], Wk[200], logWk[200];
  chan = 0;

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


  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);

  /* get distances */
  dist = eegtrials_diffmatrix_channel( eeg, distmetric, chan );

  diffmatrix_standardize(dist, eeg->ntrials);
  write_double_matrix_ascii( "diffeeg.txt", dist, eeg->ntrials, eeg->ntrials );

  /* test_cmpclust();
	  test_dummycluster();*/

  for( k=1; k<maxk; k++ ){
	 C = kmedoids( dist, eeg->ntrials, k );
	 print_cluster( C );
	 xWk[k] = (double)k;
	 Wk[k] = gap_get_within_scatter( dist, eeg->ntrials, C );
	 logWk[k] = log( Wk[k] );
	 free_cluster( C );
	 dprintf( "W[%i]=%f\n", k, Wk[k] );
  }
  PL( plot_format( &(xWk[1]), &(Wk[1]), maxk-1, "rO-") ); 
  PL( plot_format( &(xWk[1]), &(logWk[1]), maxk-1, "bO-") );
  PL( plot_show() );
  free_eegdata_trials( eeg );

  return 0;
}


void test_dummycluster(){
  double **d;
  int i, x, y;
  Clusters *C;

  int l = 25;
  d = read_double_matrix_ascii( "pam.txt", l, l, NULL );
  
  for(y=0; y<l; y++){
	 for( x= 0; x<l; x++){
		fprintf(stderr, "%f\t", d[y][x]);
	 }
	 fprintf(stderr, "\n");
  }

  C = kmedoids(d, l, 2);

  for( i=0; i<l; i++ )
	 free( d[i] );
  free( d );
  free_cluster(C);
}



void test_cmpclust(){
  Clusters *c1, *c2;
  int i;

  c1 = init_cluster(3, 10);
  c2 = init_cluster(3, 10);

  c1->n[0]=2;
  c1->clust[0][0]=5;
  c1->clust[0][1]=1;

  c1->n[1]=4;
  c1->clust[1][0]=9;
  c1->clust[1][1]=4;
  c1->clust[1][2]=8;
  c1->clust[1][3]=2;

  c1->n[2]=3;
  c1->clust[2][0]=3;
  c1->clust[2][1]=7;
  c1->clust[2][2]=6;

  c2->n[1]=2;
  c2->clust[1][0]=5;
  c2->clust[1][1]=1;

  c2->n[0]=4;
  c2->clust[0][0]=2;
  c2->clust[0][1]=8;
  c2->clust[0][2]=4;
  c2->clust[0][3]=9;

  c2->n[2]=3;
  c2->clust[2][0]=3;
  c2->clust[2][1]=7;
  c2->clust[2][2]=6;


  fprintf(stderr, "compare_clusters(c1,c2)=%i\n", compare_clusters(c1, c2));

  fprintf(stderr, "C1=\n");
  print_cluster(c1);
  fprintf(stderr, "C2=\n");
  print_cluster(c2);

  free_cluster(c1);
  free_cluster(c2);
}
