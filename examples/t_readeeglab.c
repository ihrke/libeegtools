/** \example t_readeeglab.c
 *
 * \brief Testing EEGlab-reader
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

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  int i, j;
  EEGdata_trials *eeg;

  eeg=read_segmented_eeg_from_eeglabset(argv[1]);
  print_eegdata_trials(stderr, eeg);

  free_eegdata_trials(eeg);
  return 0;
}
