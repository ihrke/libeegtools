/** \example t_timefreq.c
 *
 * \brief Testing Time-Frequency representations (Spectrogram).
 * 
 * Compilation:
 *\code
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include <argp.h>
#include "reader.h"
#include "writer.h"
#include "helper.h"
#include "definitions.h"
#include "clustering.h"
#include "averaging.h"
#include "mathadd.h"
#include "time_frequency.h"
#include "warping.h"

#include <libplotter/cplotter.h>


/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  FILE *f;
  char buffer[255];
  EEGdata_trials *eeg;
  Spectrogram *spectgram;
  double *win;
  double *data;
  int N,i;
  int win_length=257;


  double pb[2] = {20.0, 50.0};
  double A = 40;
  double srate = 500.0;
  double alpha;
  int M;
  estimate_kaiseralpha_bandpass( pb, A, srate, &alpha, &M );
  dprintf("pb=(%f,%f), A=%f, srate=%f, alpha=%f, M=%i\n",
			 pb[0],pb[1], A, srate, alpha, M );

  /* get data */
  dprintf("%s\n", argv[1]);
  eeg=read_eegtrials_from_raw( argv[1] );
  print_eegdata_trials(stderr, eeg);
  
  data = eeg->data[0]->d[0]; 
  N = eeg->data[0]->n;
 
 /* N = 2000; */
 /*  data = read_double_vector_ascii( argv[1], N, NULL ); */

  /* win = window_gaussian( NULL, win_length, 0.3 ); */
  win = window_hanning( NULL, win_length );
  oprintf( "w=%i\n", win_length);
  /* spectgram = spectrogram_stft( eeg->data[0]->d[0], eeg->data[0]->n, 500, */
  /* 										  win, win_length,  */
  /* 										  100, 100, NULL ); */
  double corner_freqs[2] = {0.0,20.0};
  spectgram = spectrogram_stft( data, N, srate,
  										  win, win_length,
  										  1000, N, corner_freqs, NULL );

  spectrogram_compute_powerspectrum( spectgram );

  /* plot_format( NULL, win, win_length, "r" ); */
  /* plot_switch( plot_add() ); */
  /* /\* plot_format( eeg->times, eeg->data[0]->d[0], eeg->data[0]->n, "r" ); *\/ */
  /* plot_format( NULL, data, N, "r" ); */
  /* plot_switch( plot_add() ); */
  plot_image( spectgram->powerspect, spectgram->N_time, spectgram->N_freq, "jet" );
  plot_show();


  /* cleaning up */
  free_eegdata_trials( eeg );
  spectrogram_free( spectgram );
  free( win );
  /*  free( data );*/

  return 0;
}
