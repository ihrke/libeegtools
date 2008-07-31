/** \file t_dtw.c
 *
 * \brief Testing DTW with restriction parameter and multiple electrodes.
 *
 * Compilation:
 *\code
 *\endcode
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* memcpy */
#include <argp.h>
#include "reader.h"
#include "writer.h"
#include "definitions.h"
#include "helper.h"
#include "clustering.h"

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
static int plotit=0;
#define PL(...) if(plotit){ __VA_ARGS__ }
#else
#define PL(...)
#endif

/* ------------- Argument parsing ------------------------------ */
error_t parse_opt (int key, char *arg, struct argp_state *state);
static struct argp_option options[] = {
  {"output",   'o', "FILE", 0,
   "Output file (default: no output)" }, 
  {"theta",   't', "double [0,1]", 0,
   "restriction parameter (default=1)" },
#ifdef HAVE_LIBPLOTTER
  {"plot",  'p', 0, OPTION_ARG_OPTIONAL,  
	"Use Flag for Plotting [0|1]" },
#endif
  { 0 }
};

/* available cmd-line args */
struct cmdargs {
  char *input;
  double theta;    /** restriction parameter */
  char  *output;   /** file to print the output to */
};

const  char * argp_program_version = "t_dtw 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Testing DTW with restriction parameter and multiple electrodes.";
static char   args_doc[] = "<input-file>";

/* the cmd-line parser */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct cmdargs     args; /* the globally available cmd-line args */

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  EEGdata_trials *eeg;
  WarpPath *P;
  double R;
  double **dist, /* distance for two trials */
	 **wdist; /* distance of two warppaths */
  int trial, chan, n1, n2, num_chan;
  int num_trials;

#ifdef HAVE_LIBPLOTTER
  extern int plotit;
  double position[2], size[2]={10,10};
  const char *cmap="cool";
  int plotmatrix=8;
#endif

  /* Parse the arguments */
  args.theta=1.0;
  args.output="stdout";
  argp_parse (&argp, argc, argv, 0, 0, &args);
  
  fprintf( stderr, "reading file          : %s\n", args.input );
  fprintf( stderr, "restriction parameter : %f\n", args.theta );
  fprintf( stderr, "output to             : %p\n", args.output);
  fprintf( stderr, "Plotting enabled      : %i\n", plotit);

  /* get data */
  eeg=read_eegtrials_from_raw(args.input);
  print_eegdata_trials(stderr, eeg);

  num_trials = 2;
  num_chan   = eeg->data[0]->nbchan;
  dprintf("num_chan,num_trial = (%i,%i)\n", num_chan, num_trials);

  /* nxn matrix (maximum necessary) */
  dist  = matrix_init( eeg->data[0]->n, eeg->data[0]->n );
  wdist = matrix_init( num_chan, num_chan );
  P     = init_warppath( eeg->data[0]->n, eeg->data[0]->n );

  /* main loop */
  for( trial=0; trial<num_trials; trial++ ){ 
	 dprintf("Matching trial %i <-> %i\n", trial, trial+1);

	 for( chan=0; chan<num_chan; chan++){
		dprintf("  Channel: %i\n", chan);
		n1 = eeg->markers[trial  ][1];
		n2 = eeg->markers[trial+1][1];
		dprintf("(n1,n2)=(%i,%i)\n", n1,n2);

		dist = DTW_build_restricted_cumdistmatrix( eeg->data[trial  ]->d[chan], n1, 
																 eeg->data[trial+1]->d[chan], n2, 
																 args.theta, dist );  

		P = DTW_path_from_cumdistmatrix(dist, n1, n2, P);

		PL( 
			int xpos = chan/plotmatrix;
			int ypos = chan - (plotmatrix*xpos);
			int sb;
			double sbpos[4], sbsize[4];
			double psize = 10;
			position[0] = channelcoords_64[chan].x; //(double)xpos;
			position[1] = channelcoords_64[chan].y; //(double)ypos;
			sbpos[0] = position[0]; sbpos[1]=position[1];
			sbpos[2] = sbpos[0]+psize;  sbpos[3]=sbpos[1]+psize;
			sbsize[0]=0; sbsize[1]=n1;
			sbsize[2]=0; sbsize[3]=n2;
			dprintf("position=(%f,%f)\n", position[0], position[1]);
			dprintf("sbpos, size = (%f,%f,%f,%f), (%f,%f,%f,%f)\n", 
					  sbpos[0], sbpos[1], sbpos[2], sbpos[3],
					  sbsize[0],  sbsize[1], sbsize[2], sbsize[3]);
			plot_image_position( dist, n1, n2, cmap, position, size); 
			sb = plot_subplot_create(sbpos, sbsize);
			plot_subplot_select(sb);
			plot_format_int( P->upath, P->spath, (P->J)+(P->K), "b." );
			plot_subplot_select(-1);
			 );

	 } /* chan loop */
	 PL( plot_add(); );
  } /* trial loop */
  
  /* cleaning up */
  PL( plot_show(); );
 
  dprintf("Freeing Memory\n");
  matrix_free(dist, eeg->data[0]->n);
  matrix_free(wdist, num_chan);
  free_eegdata_trials( eeg );
  free_warppath(P);

  return 0;
}

/* ----------------------------------------------------------------- */
/* ----- Parse a single cmd-line option. --------------------------- */

error_t parse_opt (int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct cmdargs *arguments = state->input;

#ifdef HAVE_LIBPLOTTER
  extern int plotit;
#endif

  switch (key){
  case 'o':
    dprintf(" Output to file '%s'\n", arg);
    arguments->output = (char*)malloc((strlen(arg)+1)*sizeof(char));
    arguments->output = strcpy(arguments->output, arg);
    break;
  case 't':
    dprintf(" theta=%s\n", arg);
    arguments->theta = atof(arg);
    break;
#ifdef HAVE_LIBPLOTTER
  case 'p':
    dprintf(" Plotting enabled\n");
	 plotit=1;
    break;
#endif

  case ARGP_KEY_ARG:
    if (state->arg_num >= 1){ /* Too many arguments. */
		errprintf("Too many args\n");
      argp_usage (state);
	 }
    arguments->input = arg;
    break;

  case ARGP_KEY_END:
    if (state->arg_num < 1){ /* Not enough arguments. */
		errprintf("Too few args\n");
      argp_usage (state);
	 }
    break;

  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}
