/** \example t_hmm.c
 *
 * \brief Testing Hidden Markov Models
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
#include "hmm.h"

#ifdef PLOTTER
#include <libplotter/cplotter.h>
#endif 

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  CPHiddenMarkovModel *m;
  m = cphmm_alloc( 10, 10, 10, 10, 10 );
  

  return 0;
}
