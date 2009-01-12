/** \example t_rawreader.c
 *
 * \brief Testing reader
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
#include "helper.h"

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <libplotter/cplotter.h>
#define PL(code) (code)
#else
#define PL(code)
#endif
/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  int i, j, chan, t1, t2;
  EEGdata_trials *eeg;
  EEGdata *target;
  chan = 1;
  t1 = 1;
  t2 = 2;

  dprintf("rawreader\n\n");
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);

  target = init_eegdata( eeg->data[t1]->nbchan, eeg->data[t1]->n, eeg->nmarkers_per_trial );
  //eeg_ADTW_markers_channel( eeg->data[t1], eeg->data[t2], target, chan, 1 /*theta*/ );

  PL( plot_format( eeg->times, eeg->data[t1]->d[chan], eeg->data[t1]->n, "r" ) );
  PL( plot_format( eeg->times, eeg->data[t2]->d[chan], eeg->data[t2]->n, "r" ) );  
  PL( plot_format( eeg->times, target->d[chan], target->n, "b" ) );  

/*   for(i=0; i<eeg->ntrials; i++){ */
/* 	 PL( plot_format(eeg->times, (eeg->data[i])->d[chan], (eeg->data[i])->n, "r") ); */
/*   } */
  PL( plot_show() );

  free_eegdata(target);
  free_eegdata_trials(eeg);
  return 0;
}
