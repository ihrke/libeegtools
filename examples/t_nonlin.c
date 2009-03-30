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
  
  s1 = eeg->data[0]->d[0];
  s2 = eeg->data[0]->d[1];
  p1 = phspace_init( m, tau, s1, eeg->nsamples );
  p2 = phspace_init( m, tau, s2, eeg->nsamples );
  phspace_print( stderr, p1 );
  phspace_print( stderr, p2 );

  R = recplot_init( p1->xn, p2->xn, fan, RPLOT_FAN );
  recplot_calculate( R, p1, p2 );


  plot_format( NULL, s1, eeg->nsamples, "r");
  plot_format( NULL, s2, eeg->nsamples, "b");
  plot_switch(plot_add());
  plot_image( R->R, R->m, R->n, "hot");
  plot_show();


  phspace_free( p1 );
  phspace_free( p2 );
  recplot_free( R );
  free_eegdata_trials( eeg );
}
