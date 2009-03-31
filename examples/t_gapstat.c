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

//#include "config.h"
//#define HAVE_LIBPLOTTER
#ifdef HAVE_LIBPLOTTER
#include <libplotter/cplotter.h>
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
  //  oldgapstat(argc, argv);
  int i,j;
  double mean;
  EEGdata_trials *eeg;
  GapStatistic *gapstat_handle;
  double **X, **D;
  int K, B;
  int n,p;
  unsigned long int myseed;

  myseed = (unsigned long int)time( (time_t *)NULL );
  srandom( myseed );
  
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
  
  D = vectordist_distmatrix( vectordist_euclidean, X, n, p, ALLOC_IN_FCT, NULL );

  printf("starting gapstat\n");
  K = atoi(argv[2]);
  B = atoi(argv[3]);
  gapstat_handle = gapstat_init( NULL, K, B );

  printf(" now really!\n");
  gapstat_calculate( gapstat_handle, X, n, p, vectordist_euclidean, D );
  printf(" done gapstat\nBest K=%i\n", gapstat_handle->khat);
  gapstat_print(stderr, gapstat_handle);

  FILE *f;
  f = fopen( "gapstat.txt", "w");
  for( i=0; i<gapstat_handle->K; i++ ){
	 fprintf(f,  "%10f\t", gapstat_handle->gapdistr[i]);
  }

  fprintf(f, "\n");
  for( i=0; i<gapstat_handle->K; i++ ){
	 fprintf(f,  "%10f\t", gapstat_handle->Wk[i]);
  }

  fprintf(f, "\n");
  for( i=0; i<gapstat_handle->K; i++ ){
	 mean = 0;
	 for(j=0; j<gapstat_handle->B; j++){
		mean += log(gapstat_handle->Wkref[j][i]);
	 }
	 mean /= (double)gapstat_handle->B;

	 fprintf(f,  "%10f\t", mean);
  }

  fclose(f);

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
