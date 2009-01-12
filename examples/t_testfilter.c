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

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <libplotter/cplotter.h>
#define PL(code) (code)
#else
#define PL(code)
#endif
double fct(double x){
  return x*x;
}

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  FILE *f;
  int i, j;

  double *x, *y;
  int n;
  n=100;
  x = (double*)malloc(n*sizeof(double));
  y = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++){
	 x[i] = (double)i/(double)n;
	 y[i] = fct(x[i]);
  }
  

  PL( plot_format(x, y, n, "b.") );
  moving_average( y, n, 20 );
  /*running_median( y, n, 20 );*/
  PL( plot_format(x, y, n, "r.") );

  /*  plot_fct(-10, 10, 0.1, fct, "--g");*/
  /*      plot_newplot();*/
  
  
  PL( plot_show() );
  printf("Freeing memory\n");
  free(x);
  free(y);

  return 0;
}
