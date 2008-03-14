/** \file t_denoise.c
 *
 * \brief Single Trial estimation of ERP's using wavelet denoising.
 *
 *  - using GSL
 *  - MATLAB wrapper: wrap_wavelet.m
 *
 * Format of input file:
 *   - binary double format;
 *   - for all trials: n points s_{c,i} 
 * 
 * Format of output file:
 *   - denoised signals, one after the other of length n
 *
 * Compilation:
 *\code
 *   make waveest
 *   (gcc -g -Wall -pedantic -o waveest waveest.c
 *-I/home/ihrke/local/include -L /home/ihrke/local/lib -lgsl
 *-lgslcblas -lm)
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "helper.h"
#include "denoising.h"
#include "mathadd.h" 
#include "argspec.h"

const  char * argp_program_version = "waveest 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Waveest -- single Trial estimation of ERP's";
static char   args_doc[] = "<input-file>";

/* the cmd-line parser */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct cmdargs     args; /* the globally available cmd-line args */


/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){
  extern int *verbosity;
  FILE *f, *g;
  double *data, *dptr; 
  int i, j, J;

  arg_set_default_values(&args);
  verbosity = &(args.verbose);
  
  /* Parse the arguments */
  argp_parse (&argp, argc, argv, 0, 0, &args);

#ifdef DEBUG
  arg_print_argset(stderr, &args);
#endif

  J = (int)round(glog((double)args.n, 2));
  if(((int)pow(2, J))<args.n)
    v_printf(0, " Warning: Signal length %i is not a power of 2"
			  "-> extending signal to length %i\n",
			  args.n, (int)pow(2,++J));
  
    
  v_printf(10, "Using file '%s', reading %i doubles, output to '%s'\n", 
	   args.infile, args.n, args.output);
  v_printf(10, "L=%i, e=%i, N=%i, n=%i\n", args.L, args.extension, args.N, args.n);
  data = (double*)malloc(((int)pow(2,J)) * sizeof(double));

  if((f = fopen(args.infile, "rb"))==NULL) errormsg(ERR_IO, 1);
  if((g = fopen(args.output, "wb"))==NULL) errormsg(ERR_IO, 1);

  for(j=0; j<args.N; j++){
    if(fread(data, sizeof(double), args.n, f)<args.n) 
      errormsg(ERR_IO, 1);

    /* signal extension if necessary */
    dptr = args.sigextfct(data, args.n, (int)pow(2,J));
    
    generic_denoising(data, (int)pow(2,J), args.L, args.cleanfct, args.eta);
    
    v_printf(10, " Successfully finished denoising function\n");
    if(fwrite(dptr, sizeof(double), args.n, g)<args.n)
      errormsg(ERR_IO, 1);
  }

  fclose(f);
  fclose(g);
  
  free(data);
  return 0;
}
