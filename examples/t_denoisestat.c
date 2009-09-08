/** \example t_denoisestat.c
 *
 *  Single Trial estimation of ERP's using wavelet
 * denoising. Returns statistics of multiple runs.
 *
 *  - using GSL
 *
 * Format of input file:
 *   - binary double format;
 *   - first n points are "times" array (range)
 *   - for all trials: n points u_{c,i}, then n points s_{c,i}
 * 
 * Format of output file:
 *   - denoised signals, one after the other of length n
 *   - outputs RMSE and SNR of the signals
 *   - use with paramscan.py
 *
 * Compilation:
 *   make evalfct
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include <argp.h>
#include "helper.h"
#include "denoising.h"
#include "argspec.h"

const  char * argp_program_version = "evalfct 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "evalfct -- Parameter scanning";
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
  double *data, *orgdata, *real, *times, *dptr; 
  double *stat_rmse, *stat_snr; /* arrays of length N, filled with rmse's and snr's */
  int i, j, J;
  int phandle;

  arg_set_default_values(&args);
  verbosity = &(args.verbose);
  
  /* Parse the arguments */
  argp_parse (&argp, argc, argv, 0, 0, &args);

#ifdef DEBUG
  arg_print_argset(stderr, &args);
#endif
   
  J = (int)round(glog((double)args.n, 2));
  if(((int)pow(2, J))<args.n)
    v_printf(0, " Warning: Signal length %i is not a power of 2 -> extending signal to length %i\n",
	     args.n, (int)pow(2,++J));
  
    
  v_printf(10, "Using file '%s', reading %i doubles, output to '%s'\n", 
	   args.infile, args.n, args.output);
  v_printf(10, "L=%i, e=%i, N=%i, n=%i\n", args.L, args.extension, args.N, args.n);
  data  = (double*)malloc(((int)pow(2,J)) * sizeof(double));
  orgdata  = (double*)malloc(((int)pow(2,J)) * sizeof(double));
  real  = (double*)malloc(((int)pow(2,J)) * sizeof(double));
  times = (double*)malloc(((int)pow(2,J)) * sizeof(double));
  stat_rmse = (double*)malloc(args.N * sizeof(double));
  stat_snr  = (double*)malloc(args.N * sizeof(double));

  if((f = fopen(args.infile, "rb"))==NULL) errormsg(ERR_IO, 1);
  if(args.output)
    if((g = fopen(args.output, "wb"))==NULL) errormsg(ERR_IO, 1);

  if(fread(times, sizeof(double), args.n, f)<args.n) 
    errormsg(ERR_IO, 1);

  for(j=0; j<args.N; j++){
    if(fread(real, sizeof(double), args.n, f)<args.n) 
      errormsg(ERR_IO, 1);
    if(fread(data, sizeof(double), args.n, f)<args.n) 
      errormsg(ERR_IO, 1);
    
    orgdata = memcpy(orgdata, data, args.n*sizeof(double));
    /* signal extension if necessary */
    dptr = args.sigextfct(data, args.n, (int)pow(2,J));

    /*   for(i=0; i<args.n; i++) */
    /*       printf("%.2f ", real[i]); */
    /*     printf("\n"); */
	 wavelet_denoising(data, (int)pow(2,J), args.L, args.cleanfct, args.eta);

    stat_rmse[j]=rmse(real, dptr, args.n);
    stat_snr[j] =snr (real, dptr, args.n);
    v_printf(10, " Successfully finished denoising function\n");
    if(args.output)
      if(fwrite(dptr, sizeof(double), args.n, g)<args.n)
	errormsg(ERR_IO, 1);
  }

  /* TODO: statistics for the SNR/RMSE */
  printf("+-----------------+\n| Statistics:     |\n+-----------------+\n");
  printf(" avg(RMSE) = %.2f\n", gsl_stats_mean(stat_rmse, 1, args.N));
  printf(" avg(SNR)  = %.2f\n", gsl_stats_mean(stat_snr , 1, args.N));
  printf(" sd(RMSE)  = %.2f\n", gsl_stats_sd(stat_rmse, 1, args.N));
  printf(" sd(SNR)   = %.2f\n", gsl_stats_sd(stat_snr, 1, args.N));
  fclose(f);
  if(args.output) fclose(g);

  free(times);
  free(orgdata);
  free(data); free(real);
  free(stat_rmse); free(stat_snr);
  return 0;
}

