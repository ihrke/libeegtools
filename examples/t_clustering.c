/** \example t_clustering.c
 * Testing clustering
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "reader.h"
#include "definitions.h"
#include "clustering.h"
#include "helper.h"
#include "eeg.h"

#include <libgen.h> /* basename */

#ifdef PLOTTER
#include <libplotter/cplotter.h>
#endif

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  EEG *eeg;
  eeg_simple_average( EEG *eeg );
  return 0;
}
