/*
** argspec.h
** 
** Made by Matthias Ihrke
** Login   <mihrke@localhost>
** 
** Started on  Fri Jul  6 00:37:26 2007 Matthias Ihrke
** Last update Sun Oct 14 19:12:17 2007 Matthias Ihrke
*/

#ifndef   	ARGSPEC_H_
# define   	ARGSPEC_H_
#include <argp.h>
#include "mathadd.h"
#include "denoising.h"
#include "helper.h"

/* ---------------------------------------------------------------------------- 
   -- globals for argp                                                       -- 
   ---------------------------------------------------------------------------- */
static struct argp_option options[] = {
  {"verbose",  'v', "level",      0,  "Produce verbose output according to level (0-100)" },
  {"output",   'o', "FILE", 0,
   "Output denoised curves to FILE (default: no output)" },
  {"n", 'n', "int", 0, "number of doubles to read from file"}, 
  {"extension", 'e', "fctname", 0, "signal extension ('zero', 'zeror', 'smooth', 'sym' (default))"},
  {"function", 'f', "fctname",0, "estimation function (valid: 'ti', 'sureshrink', 'heursure', 'conventional')"},
  {"eta", 't', "fctname", 0, "hard or soft thresholding ('s', 'h')"},
  {"reslevel", 'L', "level",0, "first level to use thresholding (default: 5)"},
  {"N", 'N', "trials",0, "number of trials to read and write *default=1)"},
/*   {"plot", 'p', "device",0, "any libplot device ('X', 'ps', 'gif'...) or 0 (default)"}, */
  {"theta1", 'g', "double", 0, "weight for the amplitude in d (timewarp); default=1"},
  {"theta2", 'h', "double", 0, "weight for the gradient in d (timewarp); default=1"},
  { 0 }
};

/** available cmd-line args */
struct cmdargs {
  char *infile;
  int   n;        /** number of doubles to read from file */
  int   N;        /** number of trials */
  int   L;        /** first resolution level to use thresholding on */
  int   verbose;  /** global verbosity level (1-100), used by
		     v_printf() */
  double theta1;
  double theta2;
/*   char  *plotdev;  /\* libplot device *\/ */
  int   extension;/** use zero-padding (0), symmetrization (1) */
  char  *output;   /** file to print the output to */
  double (*cleanfct)(const double*, int); /** cleaning function to use */
  double (*eta)(double, double); /** hard or soft thresholding */
  double* (*sigextfct)(double*, int, int); /** extension fct: 'sym', 'zeros' */
};

/* ---------------------------------------------------------------------------- 
   -- Functions                                                              -- 
   ---------------------------------------------------------------------------- */
error_t parse_opt (int key, char *arg, struct argp_state *state);
int arg_set_default_values(struct cmdargs *arg);
int arg_print_argset(FILE* out, const struct cmdargs *arg);

#endif 	    /* !ARGSPEC_H_ */
