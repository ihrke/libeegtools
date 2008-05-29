/** \file t_rawreader.c
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
  int i, j, chan;
  EEGdata_trials *eeg;
  chan = 1;

  dprintf("rawreader\n\n");
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);


  
  for(i=0; i<eeg->ntrials; i++){
	 PL( plot_format(eeg->times, (eeg->data[i])->d[chan], (eeg->data[i])->n, "r") );
  }
  PL( plot_show() );

  free_eegdata_trials(eeg);
  return 0;
}
