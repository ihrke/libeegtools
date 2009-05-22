/** \example t_testfilter.c
 *
 * \brief Testing filter functions
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
#include "denoising.h"
#include "writer.h"

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  int i, j;
  EEGdata_trials *eeg;
  
  /* get data */
  eeg=read_eegtrials_from_raw(argv[1]);
  print_eegdata_trials(stderr, eeg);
  
  for( i=0; i<eeg->ntrials; i++ ){
	 eeg_filter_fidlib( eeg->data[i], eeg->sampling_rate, "BpBu4/10-20" );
  }

  write_eegtrials_to_raw_file( eeg, "out.raw");

  return 0;
}
