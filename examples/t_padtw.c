/** \file t_padtw.c
 *
 * \brief Computing the PADTW of a raw-file.
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
#include "helper.h"
#include "definitions.h"
#include "clustering.h"

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
static int plotit=0;
#define PL(...) if(plotit){ __VA_ARGS__ }
void plot_eegdata( double *times, EEGdata *s, ChannelInfo *channelcoords, const char *format );
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
  EEGdata *new;
  Dendrogram *T, *Tsub;
  double **d, /* distance matrix, NxN */
	 **cumd; /* cumulated distance matrix */
  int num_chan,
	 nsamples, nmarkers,
	 N; /* number of trials */
  int i, j, idx1, idx2,
	 trials_left,
	 chan;

#ifdef HAVE_LIBPLOTTER
  extern int plotit;
  double position[2], size[2]={10,10};
  const char *cmap="cool";
  int plotmatrix=8;
#endif

  /* Parse the arguments */
  args.theta=1.0;
  args.output=NULL;
  argp_parse (&argp, argc, argv, 0, 0, &args);
  
  fprintf( stderr, "reading file          : %s\n", args.input );
  fprintf( stderr, "restriction parameter : %f\n", args.theta );
  fprintf( stderr, "output to             : %p\n", args.output);
  fprintf( stderr, "Plotting enabled      : %i\n", plotit);

  /* get data */
  eeg=read_eegtrials_from_raw(args.input);
  print_eegdata_trials(stderr, eeg);
  PL( 
	  plot_eegdata( eeg->times, eeg->data[0], channelcoords_64, "r-" ); 
	  plot_eegdata( eeg->times, eeg->data[1], channelcoords_64, "b-" ); 
	  plot_show(); 
		);
 
  /* prepare settings */
  num_chan = eeg->data[0]->nbchan;
  nsamples = eeg->data[0]->n;
  nmarkers = eeg->nmarkers_per_trial;
  N = 10;//eeg->ntrials;


  /* build a huge distance matrix over the average of the channels
	  for all trials 
  */
  d    = matrix_init( N, N );
  cumd = matrix_init( N, N );
  for( chan=0; chan<num_chan; chan++ ){
	 for( i=0; i<N-1; i++ ){
		for( j=i+1; j<N; j++ ){
		  d[i][j] = clusterdist_euclidean_pointwise( eeg->data[i], eeg->data[j], 0 /*channel*/ );
		  d[j][i] = d[i][j];
		}
	 }
	 matrix_add_matrix( cumd,  (const double**)d, N, N );
  }
  matrix_divide_scalar( cumd, N, N, (double)num_chan );
  matrix_print(d, N, N);

  /* now build hierarchical dendrogram */
  T = agglomerative_clustering( (const double**)d, N, dgram_dist_singlelinkage );
  dgram_print( T );

  /* now walk the tree to find pairs of trials to match */
  trials_left = N;
  while( trials_left >= 2 ){
	 Tsub = dgram_get_deepest( T );
	 oprintf("trials_left = %i, Tsub=%p\n", trials_left, Tsub);
	 oprintf("Tsub=%p, val=%i, lval=%i, rval=%i\n", Tsub, Tsub->val, Tsub->left->val, Tsub->right->val );
	 idx1 = Tsub->left->val;
	 idx2 = Tsub->right->val;
	
	 if( eeg->data[idx1]==NULL || eeg->data[idx2]==NULL ){
		errprintf( "try to touch a NULL-node: eeg->data[idx1]=%p, eeg->data[idx2]=%p\n", 
					  eeg->data[idx1], eeg->data[idx2] );
		return;
	 }

	 new = init_eegdata( num_chan, nsamples+1000, nmarkers );

	 eeg_ADTW_markers_channel(eeg->data[idx1], eeg->data[idx2], new, 0, args.theta);

	 /* remove the previous trials */
	 free_eegdata(eeg->data[idx1]);
	 free_eegdata(eeg->data[idx2]);
	 eeg->data[idx1] = new; /* ADTW goes to idx1 */
	 eeg->data[idx2] = NULL; /* do not touch this again */ 
	 
	 /* replace node by leaf representing ADTW(idx1, idx2) */
	 Tsub->val = idx1;
	 Tsub->left = NULL;
	 Tsub->right = NULL;

	 trials_left--;
  }

  /* cleaning up */
  dprintf("Freeing Memory\n");
  dgram_free( T );
  matrix_free( d, N );
  matrix_free( cumd, N );

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

/* ----------------------------------------------------------------------
	LIBPLOTTER
   ----------------------------------------------------------------------*/
#ifdef HAVE_LIBPLOTTER
void plot_eegdata( double *times, EEGdata *s, ChannelInfo *channelcoords, const char *format ){
  int chan, sb;
  double sbpos[4], sbsize[4];
  double psize = 10;
  double position[2], size[2]={10,10};

  for( chan=0; chan<s->nbchan; chan++ ){
	 position[0] = channelcoords[chan].y; 
	 position[1] = channelcoords[chan].x; 

	 sbpos[0] = position[0]; sbpos[1]=position[1];
	 sbpos[2] = sbpos[0]+psize;  sbpos[3]=sbpos[1]+psize;
	 sbsize[0]=times[0]; sbsize[1]=times[s->n-1];//s->n;
	 sbsize[2]=-50; sbsize[3]=50;

	 dprintf("position=(%f,%f)\n", position[0], position[1]);
	 dprintf("sbpos, size = (%f,%f,%f,%f), (%f,%f,%f,%f)\n", 
				sbpos[0], sbpos[1], sbpos[2], sbpos[3],
				sbsize[0],  sbsize[1], sbsize[2], sbsize[3]);

	 sb = plot_subplot_create(sbpos, sbsize);
	 plot_subplot_select(sb);
	 plot_format( times, s->d[chan], s->n, format );
	 plot_subplot_select(-1);
  }
}
#endif

