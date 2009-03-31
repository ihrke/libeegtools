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
#include <time.h>
#include <string.h> /* memcpy */
#include "reader.h"
#include "writer.h"
#include "definitions.h"
#include "helper.h"
#include "nonlinear.h"
#include "recurrence_plot.h"

#include <libplotter/cplotter.h>

void test_embedding( int argc, char **argv );
void test_recplot  ( int argc, char **argv );


/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  test_recplot( argc, argv );
  //test_embedding( argc, argv );
  return 0;
}
void test_embedding( int argc, char **argv ){
  int i,j;
  EEGdata_trials *eeg;
  PhaseSpace *p1;
  int tau;
  double *x1, *x2;

  if( argc<3 ) {
	 fprintf( stderr, "prog test.raw tau\n" );
	 return;
  } 
  
  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);
  tau = atoi( argv[2] );
  x1 = (double*)malloc( eeg->nsamples*sizeof(double) );
  x2 = (double*)malloc( eeg->nsamples*sizeof(double) );

  p1 = phspace_init( 2, tau, eeg->data[2]->d[0], eeg->nsamples );
  phspace_index_j( p1, 0, x1 );
  phspace_index_j( p1, 1, x2 );
  plot_format( x1, x2, eeg->nsamples, "r-" );
  plot_show( );
  
  free( x1 );
  free( x2 );
  phspace_free( p1 );
  free_eegdata_trials( eeg );
}
void test_recplot( int argc, char **argv ){
  int i,j;
  EEGdata_trials *eeg;
  PhaseSpace *p1, *p2;
  RecurrencePlot *R;
  double epsilon;
  int m, tau, fan;
  double *s1, *s2;
  WarpPath *P1, *P2;

  long seed = (long) time( NULL ); // intitializing the random number generator
  srand48(seed);

  
  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);
  if( argc<5 ) {
	 fprintf( stderr, "prog test.raw m tau fan\n" );
	 return;
  }

  m = atoi( argv[2] );
  tau=atoi( argv[3] );
  fan  = atoi( argv[4] );
  

  int trial1 = randint(0,eeg->ntrials-1);
  int trial2 = randint(0,eeg->ntrials-1);
  trial1=0;
  trial2=0;
  oprintf("Comparing Trials %i <--> %i\n", trial1, trial2 );
  s1 = eeg->data[0]->d[trial1];
  s2 = eeg->data[0]->d[trial2];
  p1 = phspace_init( m, tau, s1, eeg->nsamples );
  p2 = phspace_init( m, tau, s2, eeg->nsamples );
  phspace_print( stderr, p1 );
  phspace_print( stderr, p2 );

  R = recplot_init( p1->xn, p2->xn, fan, RPLOT_FAN );
  recplot_calculate( R, p1, p2 );
  recplot_print( stderr, R );

  P1 = recplot_los_marwan( R, 100,100 );


  double **d, **D;
  d = matrix_init( R->m, R->n );
  D = matrix_init( R->m, R->n );
  matrix_copy( R->R, d, R->m, R->n ); 
  scalar_minus_matrix( 1.0, d, R->m, R->n );
  matrix_copy( d, D, R->m, R->n ); 
  dtw_cumulate_matrix( D, R->m, R->n );

  P2 = dtw_backtrack( (const double**) D, R->m, R->n, NULL );
  matrix_normalize_by_max( D, R->m, R->n );
  //P2 = recplot_los_dtw   ( R );

  /* plot_format( NULL, s1, eeg->nsamples, "r"); */
  /* plot_format( NULL, s2, eeg->nsamples, "b"); */
  plot_image( d, R->m, R->n, "hot");
  plot_format_int( P1->t1, P1->t2, P1->n, "r");
  plot_format_int( P2->t1, P2->t2, P2->n, "g");
  plot_switch(plot_add());
   plot_image( D, R->m, R->n, "hot");
  double position[2] = {0,0};
  double size[2] = {R->m, R->n};
  double scaling[2] = {0.0, 0.01};
  write_double_matrix_ascii_file( "test.txt", D, R->m, R->n );
  plot_image_position_scaling_nocopy(D, R->m, R->n, "hot",
												 position, size, scaling);

  plot_format_int( P1->t1, P1->t2, P1->n, "r");
  plot_format_int( P2->t1, P2->t2, P2->n, "g");
  plot_switch(plot_add());
  plot_image( R->R, R->m, R->n, "hot");
  plot_format_int( P1->t1, P1->t2, P1->n, "r.o");
  plot_format_int( P2->t1, P2->t2, P2->n, "g.o");
  plot_show();


  phspace_free( p1 );
  phspace_free( p2 );
  recplot_free( R );
  free_eegdata_trials( eeg );
}
