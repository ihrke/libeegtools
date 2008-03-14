/*
*  C Implementation: loocv
*
* Description: 
*
*
* Author: Matthias Ihrke <mihrke@uni-goettingen.de>, (C) 2008
*
* Copyright: See COPYING file that comes with this distribution
*
*/
/** \file t_loocv.c
 *
 * \brief Leave One Out Cross-Validion.
 *
 *  - using GSL
 *  - MATLAB wrapper: 
 *
 *  * Format of input file:
 * 		- 1st double is N
 * 		- then double n
 *       - 1st double is R_c (in real time)
 *       - following N doubles are the R_{c,i} (in real time)
 *       - then n doubles -- times array 
 *       - then n doubles -- u_c
 *       - then N*n doubles -- s_{c,i}
 * 
 * Format of output file:
 *   - denoised signals, one after the other of length n
 *
 * Compilation:
 *\code
 *
 *\endcode
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include "helper.h"
#include "denoising.h"
#include "averaging.h"
#include "mathadd.h" 
#include "argspec.h"

const  char * argp_program_version = "loocv 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Leave One Out Cross-Validion";
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
	double N,n;
	double *cv;
	ModelData d;
	int i;
	
	arg_set_default_values(&args);
	verbosity = &(args.verbose);
	/* Parse the arguments */
	argp_parse (&argp, argc, argv, 0, 0, &args);
	#ifdef DEBUG
	arg_print_argset(stderr, &args);
	#endif
	d.den_params = (DenoisingParameters*)malloc(sizeof(DenoisingParameters));
	d.tw_params = (TimewarpParameters*)malloc(sizeof(TimewarpParameters));
	d.den_params->L = args.L;
	d.den_params->cleanfct = args.cleanfct;
	d.den_params->eta = args.eta;
	d.den_params->sigextfct= args.sigextfct;
	d.tw_params->theta1=args.theta1;
	d.tw_params->theta2=args.theta2;
	
	if((f = fopen(args.infile, "rb"))==NULL) errormsg(ERR_IO, 1);
	if(args.output){
		if((g = fopen(args.output, "wb"))==NULL) 
			errormsg(ERR_IO, 1);
	}
	if(fread(&N, sizeof(double), 1, f)<1)       errormsg(ERR_IO, 1);
	if(fread(&n, sizeof(double), 1, f)<1)       errormsg(ERR_IO, 1);
	dprintf("Got N=%.2f=%i from file\n", N,(int)N);
	dprintf("Got n=%.2f=%i from file\n", n,(int)n);
	d.N = (int)N;
	d.n = (int)n;
	d.Ri = (double*)malloc(n*sizeof(double));
	d.times=(double*)malloc(n*sizeof(double));
	d.u = (double*)malloc(n*sizeof(double));
	cv = (double*)malloc(n*sizeof(double));
	d.si    = (double**)malloc(d.N*sizeof(double*));
	for(i=0; i<d.N; i++)
		d.si[i]=(double*)malloc(d.n*sizeof(double));
	
	if(fread(&(d.R), sizeof(double), 1, f)<1)       errormsg(ERR_IO, 1);
	if(fread(d.Ri, sizeof(double), d.N, f)<d.N) errormsg(ERR_IO, 1);
	if(fread(d.times, sizeof(double), d.n, f)<d.n) errormsg(ERR_IO, 1);
	if(fread(d.u, sizeof(double), d.n, f)<d.n)     errormsg(ERR_IO, 1);
	/* get the si */
	for(i=0; i<d.N; i++)
		if(fread(d.si[i], sizeof(double), d.n, f)<d.n) 
			errormsg(ERR_IO, 1);
	 
	/* --- starting up -- */
	cv = loocv(&d, cv, iterative_warpavg);
	
	if(args.output)
		if(fwrite(cv, sizeof(double), d.n, g)<d.n)
			errormsg(ERR_IO, 1);
	
	if(f)	 fclose(f);
	if(args.output) 
		fclose(g);	
	free(d.Ri); free(d.times); 
	free(d.u);
	for(i=0; i<d.N; i++)
		free(d.si[i]);
	free(d.si);
	return 0;
}
