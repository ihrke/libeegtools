/** \example t_dtw.c
 *
 * \brief Testing DTW with restriction parameter and multiple electrodes.
 * 
 * The program computes WarpPathes for all channels, for trial n vs. n+1.
 * Then the (mean) difference matrix of all these warppathes is given as output.
 *
 
 There is an R-script

     misc/cluster_channels_with_pathdist.R

 which plots clusters given the output of that program.


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
#include "chaninfo.h"

#include "config.h"
#ifdef HAVE_LIBPLOTTER
#include <libplotter/cplotter.h>
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
  {"clusters", 'c', "int", 0, "number of clusters (default=2)"}, 
#ifdef HAVE_LIBPLOTTER
  {"plot",  'p', 0, OPTION_ARG_OPTIONAL,  
	"Use Flag for Plotting [0|1]" },
#endif
  { 0 }
};

/* available cmd-line args */
struct cmdargs {
  char *input;
  int num_clust; /** number of clusters K to compute */
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
  WarpPath **P;
  double R;
  double **dist, /* distance for two trials */
	 **wdist, /* distance of two warppaths */
	 **cwdist; /* cumulated (averaged) distance over trials */
  int trial, chan, /* counter main loop */
	 n1, n2, n,
	 num_chan, 
	 c1, c2, i; /* counter diff */
  int num_trials;
  Clusters *C;

#ifdef HAVE_LIBPLOTTER
  extern int plotit;
  double position[2], size[2]={10,10};
  const char *cmap="cool";
  int plotmatrix=8;
#endif

  /* Parse the arguments */
  args.theta=1.0;
  args.num_clust = 2;
  args.output=NULL;
  argp_parse (&argp, argc, argv, 0, 0, &args);
  
  fprintf( stderr, "reading file          : %s\n", args.input );
  fprintf( stderr, "restriction parameter : %f\n", args.theta );
  fprintf( stderr, "output to             : %p\n", args.output);
  fprintf( stderr, "Plotting enabled      : %i\n", plotit);
  fprintf( stderr, "num_clust             : %i\n", args.num_clust);

  /* get data */
  eeg=read_eegtrials_from_raw(args.input);
  print_eegdata_trials(stderr, eeg);

  num_trials = eeg->ntrials-1;
  num_chan   = eeg->data[0]->nbchan;
  n = eeg->data[0]->n;
  dprintf("num_chan,num_trial = (%i,%i)\n", num_chan, num_trials);

  /* allocating stuff */
  dist  = matrix_init( eeg->data[0]->n, eeg->data[0]->n );  /* nxn matrix (maximum necessary) */
  wdist = matrix_init( num_chan, num_chan );
  cwdist = matrix_init( num_chan, num_chan );
  P     = (WarpPath **) malloc( num_chan * sizeof(WarpPath*));
  for( chan=0; chan<num_chan; chan++ ){
	 P[chan] = init_warppath( P[chan], eeg->data[chan]->n, eeg->data[chan]->n );
  }

  FILE *out;
  out = fopen( "path.txt", "w");

  /* main loop */
  for( trial=0; trial<num_trials; trial++ ){ 
	 oprintf("Matching trial %i <-> %i\n", trial, trial+1);


	 for( chan=0; chan<num_chan; chan++){
		dprintf("  Channel: %i\n", chan);
		n1 = eeg->markers[trial  ][1];
		n2 = eeg->markers[trial+1][1];
		dprintf("(n1,n2)=(%i,%i)\n", n1,n2);

		dist = DTW_build_restricted_cumdistmatrix( eeg->data[trial  ]->d[chan], n1, 
																 eeg->data[trial+1]->d[chan], n2, 
																 args.theta, dist );  

		P[chan] = DTW_path_from_cumdistmatrix( (const double**) dist, n1, n2, P[chan]);
		if(trial==4){
		  write_int_vector_ascii( out, P[chan]->upath, n+n );
		  write_int_vector_ascii( out, P[chan]->spath, n+n );
		}
	 } /* chan loop */
	 if(trial==4) 		  
		fclose(out);

	 /* channels are known, now compare them */
	 for( c1=0; c1<num_chan-1; c1++ ){
		for( c2=c1+1; c2<num_chan; c2++){
		  wdist[c1][c2] = DTW_distance_between_paths(P[c1], P[c2]);
		  wdist[c2][c1] = wdist[c1][c2];
		}
	 }
	 
	 /* cumulate for averaging over trials */
	 matrix_add_matrix( cwdist, (const double**)wdist, num_chan, num_chan ); 
  } /* trial loop */

  matrix_divide_scalar( cwdist, num_chan, num_chan, (double) num_trials );

  /* clustering */
  C = kmedoids( (const double**) cwdist, num_chan, args.num_clust );
  print_cluster( C );
	 
#ifdef HAVE_LIBPLOTTER
  if( plotit ){
	 int sb;
	 int k;
	 //plot_image( wdist, num_chan, num_chan, "jet" ); 
	 double sbpos[4], sbsize[4];
	 double sbx[64], sby[64];
	 const char clustformat[][20] = {"r.0", "g.0", "b.0", "m.0", "y.0"};
	 sbpos[0] = -10; sbpos[1]=-10;
	 sbpos[2] = 0;   sbpos[3]=0;
	 sbsize[0]= -100; sbsize[1]=100;
	 sbsize[2]= -100; sbsize[3]=100;
	 dprintf("sbpos, size = (%f,%f,%f,%f), (%f,%f,%f,%f)\n", 
				sbpos[0], sbpos[1], sbpos[2], sbpos[3],
				sbsize[0],  sbsize[1], sbsize[2], sbsize[3]);
	 /* sb = plot_subplot_create(sbpos, sbsize); */
	 /* plot_subplot_select(sb); */
		 
	 for( k=0; k<C->K; k++ ){
		dprintf("Prepare cluster %i/%i\n", k, C->K);
		for( i=0; i<C->n[k]; i++ ){
		  dprintf(" C->clust[k][i]=%i,  channelcoords_64[ C->clust[k][i] ].x=%f\n", 
					 C->clust[k][i],channelcoords_64[ C->clust[k][i] ].x ); 
		  sbx[i] = channelcoords_64[ C->clust[k][i] ].x;
		  sby[i] = channelcoords_64[ C->clust[k][i] ].y;
		}
		dprintf("Plot cluster %i/%i\n", k, C->K);
		plot_format( sbx, sby, C->n[k], clustformat[k] );
	 }

	 /* plot_subplot_select(-1); */
  }

#endif
	 
  free_cluster( C );

  /* cleaning up */
  PL( plot_show(); );
 

  /* output */
  if(!args.output){
	 fprintf(stdout, "# ");
	 for( i=0; i<argc; i++ ){
		fprintf( stdout, "%s ", argv[i]);
	 }
	 fprintf( stdout, "\n");
	 write_double_matrix_ascii( stdout, (const double**)cwdist, num_chan, num_chan );
  } else {
	 FILE *f;
	 if((f= fopen(args.output, "w"))==NULL) 	 
		errormsg(ERR_IO, 1);
	 fprintf(f, "# ");
	 for( i=0; i<argc; i++ ){
		fprintf( f, "%s ", argv[i]);
	 }
	 fprintf( f, "\n");
	 write_double_matrix_ascii( f, (const double**)cwdist, num_chan, num_chan );
	 fclose(f);
  }

  dprintf("Freeing Memory\n");
  matrix_free(dist, eeg->data[0]->n);
  matrix_free(wdist, num_chan); 
  matrix_free(cwdist, num_chan); 

  for( chan=0; chan<num_chan; chan++ ){
	 free_warppath(P[chan]);
  }
  free(P);
  free_eegdata_trials( eeg );

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

  case 'c':
	 dprintf(" num_clusters=%s\n", arg);
	 arguments->num_clust=atoi(arg);
	 break;
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
/*
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
*/
