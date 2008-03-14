/**\file t_warpavg.c
 * \brief Test the complete Warpavg Procedure.
 *
 * Format of input file:
 * 		- N,n, zero
 *       - following N doubles are the sRi sampline points
 *       - then N*n doubles -- ui
 *
 * Output:
 *
 * Compilation:
 *   make t_warpavg
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include <argp.h>
#include "helper.h"
#include "denoising.h"
#include "argspec.h"

/* #define DEBUG_LOCAL  */

const  char * argp_program_version = "t_warpavg 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "";
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
  double **ui, *wavg, *sRid, tmp;
  /*int *sRi;*/
  int **markers;
  int i,j,k, zero;
 

  arg_set_default_values(&args);
  verbosity = &(args.verbose);
  
  /* Parse the arguments */
  argp_parse (&argp, argc, argv, 0, 0, &args);

#ifdef DEBUG
  arg_print_argset(stderr, &args);
#endif
   
  v_printf(10, "Using file '%s', reading %i doubles, output to '%s'\n", 
	   args.infile, args.n, args.output);

  if((f = fopen(args.infile, "rb"))==NULL) errormsg(ERR_IO, 1);
  if(args.output){
    if((g = fopen(args.output, "wb"))==NULL) 
      errormsg(ERR_IO, 1);
  }
  if(fread(&tmp, sizeof(double), 1, f)<1)              errormsg(ERR_IO, 1);
  args.N = (int)tmp;
  if(fread(&tmp, sizeof(double), 1, f)<1)              errormsg(ERR_IO, 1);
  args.n = (int)tmp;

#ifdef DEBUG
  arg_print_argset(stderr, &args);
#endif

  
  /* allocating mem */
  sRid    = (double*) malloc(args.N*sizeof(double));
  markers = (int **)   malloc(args.N*sizeof(int*));
  for(i=0; i<args.N; i++) markers[i]=(int*)malloc(2*sizeof(int));
  ui      = (double**)malloc(args.N*sizeof(double*));
  for(i=0; i<args.N; i++)
    ui[i]=(double*)malloc(args.n*sizeof(double));
  wavg  = (double*)calloc(args.n, sizeof(double));

  /* get Ri */
  if(fread(&tmp, sizeof(double), 1, f)<1)              errormsg(ERR_IO, 1);
  zero = (int)tmp;
  if(fread(sRid, sizeof(double), args.N, f)<args.N)    errormsg(ERR_IO, 1);
  for(i=0; i<args.N; i++){
	  markers[i][0]=zero;
	  markers[i][1]=(int)sRid[i];
  }
  

/* get the ui */
  for(i=0; i<args.N; i++)
    if(fread(ui[i], sizeof(double), args.n, f)<args.n) 
      errormsg(ERR_IO, 1);

  
  wavg = warpaverage(ui, args.N, args.n, markers, wavg);
  
  /* cleaning up */
  fclose(f); 
  if(args.output) fclose(g);
  for(i=0; i<args.N; i++){
	  free(markers[i]);
	  free(ui[i]);
  }
  free(ui);
  free(markers);  free(sRid);
  free(wavg);
  return 0;
}

