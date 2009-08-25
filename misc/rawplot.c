/** \cond plotter */

#include "definitions.h"


/* usual includes */
#include <stdio.h>

#ifdef PLOTTER

#include <argp.h>
/* my includes */
#include "helper.h"
#include "chaninfo.h"
#include "tools.h"
#include <libplotter/cplotter.h>

/* function definitions */
void plot_eegdata_nocopy( double *times, EEGdata *s, ChannelInfo *channelcoords, 
								  int *chan_idx, int num_chans, const char *format );
void plot_eegdata_notopology_nocopy( double *times, EEGdata *s, int *chan_idx, 
												 int num_chans, const char *format );

/* ------------- Argument parsing ------------------------------ */
error_t parse_opt (int key, char *arg, struct argp_state *state);

static struct argp_option options[] = {
  {"channel", 'c', "vector", 0, 
	"plot channels in vector 'idx, idx2,...', default 'all'"},
  {"trials", 't', "vector", 0, "plot trials in vector, default 'all'"},
  {"separate-plots", 's',  0, OPTION_ARG_OPTIONAL, "Plot trials in separate plots. " 
  "This is much slower but allows to scroll back and forth between trials."
  "Else, only forward scrolling is possible."},
  {"notopology", 'n', 0, OPTION_ARG_OPTIONAL,  "plot all channels on top of each other (default)"},
  {"topology_file", 'f', "filename", 0, "plot using this topology (overrides -n)"},
  { 0 }
};

/* available cmd-line args */
struct cmdargs {
  char *input;
  int *channels;
  int num_channels;
  int num_trials;
  int *trials;
  int sep_plots;
  int notopology;
  char *topology_file;
};
void args_free ( struct cmdargs );
void args_print( struct cmdargs );
void args_setdefault( struct cmdargs* );
const  char * argp_program_version = "rawplot 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Plotting Contents of RAW-file.";
static char   args_doc[] = "<input-file>";

/* the cmd-line parser */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct cmdargs     args; /* the globally available cmd-line args */

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main( int argc, char **argv ){
  EEGdata_trials *eeg;
  EEGdata *trial;
  int org_ntrials, org_nchans, org_nsamples, 
	 org_nmarkers; /* as read from file */
  int i, t;
  ChannelInfo *chaninfo;

  chaninfo = channelcoords_64;

  /* Parse the arguments */
  args_setdefault( &args );
  argp_parse (&argp, argc, argv, 0, 0, &args);
  args_print( args );

  /* get data */
  eeg=read_eegtrials_from_raw( (const char*) args.input );
  print_eegdata_trials(stderr, eeg);

  /* eegtrials_remove_baseline( eeg, -500, 0 ); */

  org_ntrials  = eeg->ntrials;
  org_nchans   = eeg->data[0]->nbchan;
  org_nsamples = eeg->data[0]->n;
  org_nmarkers = eeg->nmarkers_per_trial;
  if( args.topology_file ){
	 chaninfo = read_chaninfo_ced( args.topology_file, ALLOC_IN_FCT );
  } 

  trial = init_eegdata( org_nchans, org_nsamples, org_nmarkers );
  oprintf(" This is trial 0\n");
  copy_similar_eegdata( trial, eeg->data[0] );

  /* if num_channels==0, use all chan-indices */
  if( args.num_channels==0 ){
	 args.channels=(int*)malloc( org_nchans*sizeof(int) );
	 args.num_channels = org_nchans;
	 for( i=0; i<org_nchans; i++ )
		args.channels[i] = i;
  }

  /* if num_trials==0, use all trial-indices */
  if( args.num_trials==0 ){
	 args.trials=(int*)malloc( org_ntrials*sizeof(int) );
	 args.num_trials = org_ntrials;
	 for( i=0; i<org_ntrials; i++ )
		args.trials[i] = i;
  }


  if( args.notopology ){
	 plot_eegdata_notopology_nocopy( eeg->times, trial, args.channels, 
												args.num_channels, "r.-" );
  } else {
	 plot_eegdata_nocopy( eeg->times, trial, chaninfo, 
								 args.channels, args.num_channels, "r-" );
  }
  plot_show_in_thread();

  for( i=1; i<args.num_trials; i++ ){	 
	 plot_update();
	 if( fgetc( stdin ) == 'q' ) break;
	 t = args.trials[i]; /* actual trial index */
	 oprintf(" This is trial %i\n", t );
	 copy_similar_eegdata( trial, eeg->data[i] );
  }
  
  plot_update();
  fgetc( stdin );
  
  args_free( args );
  free_eegdata_trials( eeg );
  free_eegdata( trial );
  if( args.topology_file ){
	 free( chaninfo );
  }
  return 0;
}

/* ----------------------------------------------------------------- */
/* ----- Parse a single cmd-line option. --------------------------- */

error_t parse_opt (int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct cmdargs *arguments = state->input;
  int count,i,tmp; 
  char *comma;

  switch (key){
  case 'c': /* which channels? */
	 dprintf(" using channels '%s'\n", arg);
	 if( !strcmp( arg, "all" ) ) break;
	 count = strcount( arg, ',' );
	 arguments->num_channels = count+1;
	 arguments->channels = (int*) malloc( arguments->num_channels*sizeof(int) );
	 for( i=0; i<arguments->num_channels; i++ ){
		if( (comma = strchr( arg, ',' ))!=NULL )
		  comma[0] = '\0';
		arguments->channels[i] = atoi( arg );
		dprintf( "a->c[%i] = %i\n", i, arguments->channels[i]);
		arg = comma+1;
	 }
    break;

  case 't': /* which trials? */
	 dprintf(" plotting trials '%s'\n", arg);
	 if( !strcmp( arg, "all" ) ) break;
	 count = strcount( arg, ',' );
	 arguments->num_trials = count+1;
	 arguments->trials = (int*) malloc( arguments->num_trials*sizeof(int) );
	 for( i=0; i<arguments->num_trials; i++ ){
		if( (comma = strchr( arg, ',' ))!=NULL )
		  comma[0] = '\0';
		arguments->trials[i] = atoi( arg );
		dprintf( "a->c[%i] = %i\n", i, arguments->trials[i]);
		arg = comma+1;
	 }
    break;

  case 'n': /* notopology */
	 arguments->notopology=1;
	 break;
  case 's': /* separate-plots */
	 arguments->sep_plots=1;
	 errprintf("sep-plots NOT IMPLEMENTED\n");
	 break;
	 
  case 'f': /* topology-file */
	 arguments->notopology=0;
	 arguments->topology_file=arg;
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

void args_setdefault(  struct cmdargs *arg){
  arg->input = NULL;
  arg->channels = NULL;
  arg->num_channels=0;
  arg->num_trials=0;
  arg->trials=NULL;
  arg->notopology=1;
  arg->sep_plots=0;
  arg->topology_file=NULL;
}

void args_print(  struct cmdargs arg){
  int i;

  oprintf( "file=%s\n", arg.input );
  oprintf( "  notopology=%i\n", arg.notopology );
  oprintf( "  top-file=%s\n", arg.topology_file );
  oprintf( "  sep_plots =%i\n", arg.sep_plots );
  if( arg.num_channels==0 ){
	 oprintf( " all channels\n");
  } else {
	 oprintf( "channels=");
	 for( i=0; i<arg.num_channels; i++ ){
		printf( "%i,", arg.channels[i] );
	 }
	 printf( "\b\  (%i)\n", arg.num_channels );
  }
  if( arg.num_trials==0 ){
	 oprintf( " all trials\n");
  } else {
	 oprintf( "trials=");
	 for( i=0; i<arg.num_trials; i++ ){
		printf( "%i,", arg.trials[i] );
	 }
	 printf( "\b\  (%i)\n", arg.num_trials );
  }
  
}
void args_free(  struct cmdargs arg){
  if(arg.channels){
	 free(arg.channels);
  }
  if(arg.trials){
	 free(arg.trials);
  }
}


void plot_eegdata_nocopy( double *times, EEGdata *s, 
								  ChannelInfo *channelcoords, 
								  int *chan_idx, int num_chans,
								  const char *format ){
  int chan, c, sb;
  double sbpos[4], sbsize[4];
  double psize = 10;
  double position[2], size[2]={10,10};
  int nm, i;  
  double xlinex[2], yliney[2], xliney[2]={0};
	 
  xlinex[0] = times[0]; xlinex[1] = times[s->n-1];
  yliney[0] = -100; yliney[1] = 100;


  for( c=0; c<num_chans; c++ ){
	 chan = chan_idx[c];
	 position[0] = channelcoords[chan].y; 
	 position[1] = channelcoords[chan].x; 

	 sbpos[0] = position[0]; sbpos[1]=position[1];
	 sbpos[2] = sbpos[0]+psize;  sbpos[3]=sbpos[1]+psize;
	 sbsize[0]=times[0]; sbsize[1]=times[s->n-1];//s->n;
	 sbsize[2]=-100; sbsize[3]=100;

	 /*	 dprintf("position=(%f,%f)\n", position[0], position[1]);
	 dprintf("sbpos, size = (%f,%f,%f,%f), (%f,%f,%f,%f)\n", 
				sbpos[0], sbpos[1], sbpos[2], sbpos[3],
				sbsize[0],  sbsize[1], sbsize[2], sbsize[3]);
	 */
	 sb = plot_subplot_create(sbpos, sbsize);
	 plot_subplot_select(sb); 
	 plot_format( xlinex, xliney, 2, "k-" );
	 plot_format( xliney, yliney, 2, "k-" );
	 plot_format_nocopy( times, s->d[chan], s->n, format );

	 nm = s->nmarkers;
	 for( i=0; i<nm; i++ ){
		/* dprintf(" i=%i, (x,y)=(%*/
		plot_format_nocopy( &(times[s->markers[i]]), &(s->d[chan][s->markers[i]]), 1, "kO" );
	 }
	 plot_subplot_select(-1);
  }
}

void plot_eegdata_notopology_nocopy( double *times, EEGdata *s, 
												 int *chan_idx, int num_chans,
												 const char *format ){
  int chan, c, sb, i;
  double xlinex[2], yliney[2], xliney[2]={0};
	 
  xlinex[0] = times[0]; xlinex[1] = times[s->n-1];
  yliney[0] = -100; yliney[1] = 100;
  
  plot_format( xlinex, xliney, 2, "k-" );
  plot_format( xliney, yliney, 2, "k-" );
  for( c=0; c<num_chans; c++ ){
	 chan = chan_idx[c];
	 plot_format_nocopy( times, s->d[chan], s->n, format );
	 
	 for( i=0; i<s->nmarkers; i++ ){
		plot_format_nocopy( &(times[s->markers[i]]), &(s->d[chan][s->markers[i]]), 1, "kO" );
	 }
  }
}
#else

int main(){
  fprintf( stderr, "Sorry, no Plotter installed\n");
  return -1;
}


#endif
/** \endcond */
