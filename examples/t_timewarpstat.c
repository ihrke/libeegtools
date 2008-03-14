/** \file t_timewarpstat.c
 * \brief Testing Timewarping and report statistics.
 *
 * Format of input file:
 *   - binary double format;
 *   - first double is the length of "real" signal u (Nu)
 *   - next come Nu doubles for the real signal u
 *   - next N points (n1, ..., nN) give the no. 
 *     of sampling points for the following N samples (u_i) 
 *   - then for all trials: n_i points which give the signal u_i
 * 
 * Output:
 *   if output is specified, a double file with the N signals,
 *   all of length Nu 
 *
 * Compilation:
 *   make t_timewarpstat
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include <argp.h>
#include "helper.h"
#include "denoising.h"
#include "argspec.h"

/*#define DEBUG_LOCAL*/

const  char * argp_program_version = "t_timewarpstat 0.1";
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
  double *u, *ui;
  double *lens;
  double lenu;
  int maxn;
  int i, j;
  int phandle;
  double *stat_rmse, *stat_snr; /* arrays of length N, filled with rmse's and snr's */

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
  /* get u */
  fread(&lenu, sizeof(double), 1, f);
  u = (double*)malloc(lenu*sizeof(double));
  if(fread(u, sizeof(double), lenu, f)<lenu)
    errormsg(ERR_IO, 1);

  /* get the n_1, ..., n_N */
  lens  = (double*)malloc(args.N*sizeof(double));
  if(fread(lens, sizeof(double), args.N, f)<args.N) 
    errormsg(ERR_IO, 1);

#ifdef DEBUG_LOCAL
  printf("lens=\n"); for(j=0; j<args.N; j++) printf("%2.2f ", lens[j]); printf("\n");
  printf("u=\n"); for(j=0; j<lenu; j++) printf("%2.2f ", u[j]); printf("\n");
#endif

  maxn = (int)ceil(maxel(lens, args.N));
  ui = (double*)malloc(maxn*sizeof(double));


  /* initialize statistics */
  if((stat_rmse = (double*)malloc(args.N * sizeof(double)))==NULL) errormsg(ERR_MEM, 1);
  if((stat_snr  = (double*)malloc(args.N * sizeof(double)))==NULL) errormsg(ERR_MEM, 1);

  dprintf("maxn=%i, lenu=%i\n", (int)maxn, (int)lenu);
  for(i=1; i<=args.N; i++){ /* trial counter */
    if(fread(ui, sizeof(double), lens[i-1], f)<lens[i-1]) 
      errormsg(ERR_IO, 1);

#ifdef DEBUG_LOCAL
    printf("ui(%i)=\n", i); for(j=0; j<lens[i-1]; j++) printf("%2.2f ", ui[j]); printf("\n");
#endif


    /* -------- Warping ------------- */
    timewarp(u, lenu, ui, lens[i-1], args.theta1, args.theta2);
    stat_rmse[i-1]=rmse(u, ui, lenu);
    stat_snr[i-1] =snr (u, ui, lenu);

    if(args.output)
      fwrite(ui, sizeof(double), lenu, g);
  }


 /* TODO: statistics for the SNR/RMSE */
  printf("+-----------------+\n| Statistics:     |\n+-----------------+\n");
  printf(" avg(RMSE) = %.4f\n", gsl_stats_mean(stat_rmse, 1, args.N));
  printf(" avg(SNR)  = %.4f\n", gsl_stats_mean(stat_snr , 1, args.N));
  printf(" sd(RMSE)  = %.4f\n", gsl_stats_sd(stat_rmse, 1, args.N));
  printf(" sd(SNR)   = %.4f\n", gsl_stats_sd(stat_snr, 1, args.N));

  fclose(f); 
  if(args.output) fclose(g);
  free(stat_rmse); 
  free(stat_snr);
  free(lens);
  free(u);
  free(ui);

  return 0;
}

