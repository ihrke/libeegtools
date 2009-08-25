/** \example t_som.c
 *
 * \brief Testing SOM
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
#include "som.h"

#ifdef PLOTTER
#include <libplotter/cplotter.h>
#endif 

void som_train_from_data2( Som *s, double **X, int dim, int nsamples ){
  int i, j, t;
  int bmu; /* best matching unit */
  double h, tmp; 
  double bmu_score; /* best matching unit score */
  double *input;

  if( dim != s->dimension ){
	 errprintf( "Data dimension and SOM-dimension do not match... exit\n");
	 return;
  }

  if( s->progress ){
	 (*(s->progress))( PROGRESSBAR_INIT, s->nruns );
  }
  double *x, *y;
  x = vector_init( NULL, s->n, 0.0 );
  y = vector_init( NULL, s->n, 0.0 );
  plot_format_nocopy( x, y, s->n, "r.-" );
  plot_set_limits( -500, 2000, -20, 20 );
  plot_show_in_thread();
  getc(stdin);

  for( t=1; t<=s->nruns; t++ ){
	 if( s->progress ){
		s->progress( PROGRESSBAR_CONTINUE_LONG, t-1 );
	 }

	 /* choose a sample */
	 input = X[gsl_rng_uniform_int( s->rng, nsamples )];//(t-1) % nsamples];

	 /* find BMU */
	 bmu_score = DBL_MAX;
	 for( i=0; i<s->n; i++ ){	 
		if( (tmp=s->distancefct( input, s->m[i], dim, NULL ))<bmu_score ){
		  bmu_score = tmp;
		  bmu = i;
		}
	 }
	 
	 if( t%10==0 ){
		for( i=0; i<s->n; i++ ){
		  x[i] = s->m[i][0];
		  y[i] = s->m[i][1];
		}
		
		plot_update();
		//usleep(1000);
	 }
  
  
	 /* adapt weights (learning) */
	 for( i=0; i<s->n; i++ ){
		if( s->progress ){
		  s->progress( PROGRESSBAR_CONTINUE_SHORT, i );
		}

		h = s->neighbourhoodfct( i, bmu, s, t );
		for( j=0; j<dim; j++ ){
		  s->m[i][j] += h*(input[j]-s->m[i][j]);
		}
	 }
  }

  getc(stdin);
}

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  EEGdata_trials *eeg;
  Som *s;
  double **X;
  int i,j;
  int ntrials;
  double *avg;


  if( argc<2 ){
	 errprintf("Need a raw-file\n");
	 return -1;
  }
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);
  ntrials = 10;// eeg->ntrials;

  avg = vector_init( NULL, eeg->nsamples, 0.0 );
  for( i=0; i<ntrials; i++ ){
	 plot_format_nocopy( eeg->times, eeg->data[i]->d[0], eeg->nsamples, "y.-" ); 
	 for( j=0; j<eeg->nsamples; j++ ){
		avg[j] += eeg->data[i]->d[0][j];
	 }
  }
  for( i=0; i<eeg->nsamples; i++ ){
	 avg[i] /= (double)ntrials;
  }
  plot_format_nocopy( eeg->times, avg, eeg->nsamples, "b.-" ); 

  oprintf("Using %i trials starting from 0\n", ntrials);
  
  X = matrix_init( eeg->nsamples*ntrials, 2 );
  for( i=0; i<ntrials; i++ ){
	 for( j=0; j<eeg->nsamples; j++ ){
		X[i*ntrials+j][0] = eeg->times[j];
		X[i*ntrials+j][1] = eeg->data[i]->d[0][j];
	 }
  }
  
  s = som_init( 2, 2*eeg->nsamples, 100000, ONED_LINEAR );
  s->progress = progressbar_rotating;

  som_print( stdout, s );

  som_initialize_random_samples( s, X, 2, ntrials*eeg->nsamples );
  for( i=0; i<eeg->nsamples; i++ ){
	 s->m[2*i  ][0] = eeg->times[i];
	 s->m[2*i+1][0] = eeg->times[i];
	 s->m[2*i  ][1] = X[i][1];
	 s->m[2*i+1][1] = X[i][1];
  }
  
  som_train_from_data( s, X, 2, ntrials*eeg->nsamples );
  
  
  return 0;
}
