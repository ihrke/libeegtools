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
  EEG *eeg;

  eeg=read_eeglab_file( "/scratch.local/ihrke/data/eegtast_19_epoched.set" );

  return 0;
}
