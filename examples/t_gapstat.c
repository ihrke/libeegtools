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

#ifdef PLOTTER
#include <libplotter/cplotter.h>
#endif

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  //  oldgapstat(argc, argv);
  int i,j;
  double mean[100];
  double logWk[100];
  EEGdata_trials *eeg;
  GapStatistic *gapstat_handle;
  double **X, **D;
  int K, B;
  int n,p;
  unsigned long int myseed;

  myseed = (unsigned long int)time( (time_t *)NULL );
  srandom( myseed );
  

  if( argc<4 ){
	 fprintf( stderr, "Usage: t_gapstat test.raw max_num_clusters num_rep\n");
	 return -1;
  }
  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);
  
  /* compose nxp-matrix (observations x features) */
  X = matrix_init( eeg->ntrials, eeg->nsamples );
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=0; j<eeg->nsamples; j++ ){
		/* just a single channel for now */
		X[i][j] = eeg->data[i]->d[0][j];
	 }
  }
  n = eeg->ntrials; p = eeg->nsamples;
  
/*   n=100; */
/*   p = 2; */
/*   X = read_double_matrix_ascii( argv[1], p, n, X); */
  
  fprintf( stderr, "Distance: ");  
  D = vectordist_distmatrix( vectordist_euclidean, X, n, p, 
									  ALLOC_IN_FCT, progressbar_rotating, 
									  (void*)signaldist_euclidean  );
  write_double_matrix_ascii_file( "out.txt", D, n, n );


  K = atoi(argv[2]);
  B = atoi(argv[3]);
  gapstat_handle = gapstat_init( NULL, K, B );
  gapstat_handle->progress = progressbar_rotating;

  fprintf( stderr, "Gap-Statistic: ");
  gapstat_calculate( gapstat_handle, X, n, p, vectordist_euclidean, D );
  printf(" done gapstat\nBest K=%i\n", gapstat_handle->khat);
  gapstat_print(stderr, gapstat_handle);


#ifdef PLOTTER
  plot_format( NULL, gapstat_handle->gapdistr+1, gapstat_handle->K-1, "r-o" );

  for( i=0; i<gapstat_handle->K; i++ ){
	 mean[i] = 0;
	 for(j=0; j<gapstat_handle->B; j++){
		mean[i] += log(gapstat_handle->Wkref[j][i]);
	 }
	 mean[i] /= (double)gapstat_handle->B;

	 logWk[i] = log(gapstat_handle->Wk[i]);
  }

  /* plot_format( NULL, mean, gapstat_handle->K, "g" );  */
  /* plot_format( NULL, logWk, gapstat_handle->K, "b" );  */

  plot_show();
#endif
}

