/*
** argspec.c
** 
** Made by (Matthias Ihrke)
** Login   <mihrke@localhost>
** 
** Started on  Fri Jul  6 00:38:09 2007 Matthias Ihrke
** Last update Wed Oct 10 20:49:30 2007 Matthias Ihrke
*/

#include "argspec.h"

/* ----------------------------------------------------------------- */

int arg_set_default_values(struct cmdargs *arg){
  arg->verbose=100;
  arg->output=NULL;
  arg->N = 2;
/*   arg->plotdev = NULL; */
  arg->cleanfct = sureshrink;
  arg->sigextfct= sigext_sym;
  arg->eta = eta_s;
  arg->L = 5;
  arg->extension = 1;
  arg->theta1=1.0;
  arg->theta2=1.0;
  return 0;
}

int arg_print_argset(FILE* out, const struct cmdargs *arg){
  fprintf(out, "+Argset content:\n");
  fprintf(out, "  infile =%s\n", arg->infile);
  fprintf(out, "  verbose=%i\n", arg->verbose);
  fprintf(out, "  output =%s\n", arg->output);
/*   fprintf(out, "  plotdev=%s\n", arg->plotdev); */
  fprintf(out, "  N      =%i\n", arg->N);
  fprintf(out, "  L      =%i\n", arg->L);
  fprintf(out, "  sigext =%i\n", arg->extension);
  fprintf(out, "  theta1  =%.2f\n", arg->theta1);
  fprintf(out, "  theta2  =%.2f\n", arg->theta2);
  return 0;
}

/* Parse a single cmd-line option. */
error_t parse_opt (int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct cmdargs *arguments = state->input;

  switch (key){
  case 'v':
    arguments->verbose = atoi(arg);
    break;
  case 'o':
    v_printf(10, " Output to file '%s'\n", arg);
    arguments->output = (char*)malloc((strlen(arg)+1)*sizeof(char));
    arguments->output = strcpy(arguments->output, arg);
    break;
  case 'n':
    arguments->n = atoi(arg);
    break;
  case 'N':
    arguments->N = atoi(arg);
    break;
  case 'L':
    arguments->L = atoi(arg);
    break;
  case 'g':
    arguments->theta1 = strtod(arg, (char**)NULL);
    break;
  case 'h':
    arguments->theta2 = strtod(arg, (char**)NULL);
    break;
/*   case 'p': */
/*     v_printf(10, " Plotdev = '%s'\n", arg); */
/*     if(strcmp(arg, "0")==0) */
/*       arguments->plotdev = NULL; */
/*     else { */
/*       arguments->plotdev = (char*)malloc((strlen(arg)+1)*sizeof(char)); */
/*       arguments->plotdev = strcpy(arguments->plotdev, arg); */
/*     } */
/*     break; */
  case 'e':
    v_printf(10, " Using extension function '%s'\n", arg);
    if(strcmp(arg, "zero")==0)
      arguments->sigextfct = sigext_zeros;
    else if(strcmp(arg, "zeror")==0)
      arguments->sigextfct = sigext_zerosr;
    else if(strcmp(arg, "sym")==0)
      arguments->sigextfct = sigext_sym;
    else if(strcmp(arg, "smooth")==0)
      arguments->sigextfct = sigext_smooth;
    else
      v_printf(0, "Warning: Function '%s' not available, using default\n", arg);
    break;
  case 'f':
    v_printf(10, " Using estimation function '%s'\n", arg);
    if(strcmp(arg, "ti")==0)
      arguments->cleanfct = translation_invariant_thresholding;
    else if(strcmp(arg, "sureshrink")==0)
      arguments->cleanfct = sureshrink;
    else if(strcmp(arg, "heursure")==0)
      arguments->cleanfct = heuristic_sure;
    else if(strcmp(arg, "conventional")==0)
      arguments->cleanfct = conventional_thresholding;
    else
      v_printf(0, "Warning: Function '%s' not available, using default\n", arg);
    break;
  case 't':
    v_printf(10, " Using thresholding function (eta) '%s'\n", arg);
    if(strcmp(arg, "s")==0)
      arguments->eta = eta_s;
    else if(strcmp(arg, "h")==0)
      arguments->eta = eta_h;
    else
      v_printf(0, "Warning: Eta '%s' not available, using default\n", arg);
    break;

  case ARGP_KEY_ARG:
    if (state->arg_num >= 1) /* Too many arguments. */
      argp_usage (state);
    arguments->infile = arg;
    break;

  case ARGP_KEY_END:
    if (state->arg_num < 1) /* Not enough arguments. */
      argp_usage (state);
    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}
