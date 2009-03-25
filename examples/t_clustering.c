/** \example t_clustering.c
 * Testing clustering
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
#include "clustering.h"
#include "helper.h"
#include <libgen.h> /* basename */

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
  return 0;
}
