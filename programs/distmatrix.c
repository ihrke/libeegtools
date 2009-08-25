/** \file distmatrix.c
 *
 * \brief Takes a .raw input file and computes the distance matrix between
          the trials. Storage as ascii-file.
*/

#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include "reader.h"
#include "definitions.h"
#include "writer.h"
#include "warping.h"
#include "averaging.h"
#include "mathadd.h"

#define BUFSIZE 255

/* ------------- Argument parsing ------------------------------ */
error_t parse_opt (int key, char *arg, struct argp_state *state);

static struct argp_option options[] = {
  {"output", 'o', "filename", 0, "specifies output file (default '<input>.dst')"}, 
  {"metric", 'm', "string", 0, "{'euclid', 'normalize_euclid', 'dtw'}"},
  {"channels", 'c', "chan1,chan2,...", 0, "channels to use for calculation, default all"},
  { 0 }
};

/* available cmd-line args */
struct cmdargs {
  char output_filename[BUFSIZE];
  char eeg_filename[BUFSIZE];
  VectorDistanceFunction  metric;
  PointwiseDistanceFunction dtwdist;
  int num_channels;
  int *channels;
};
void args_setdefault(  struct cmdargs *arg );
void args_free ( struct cmdargs );

const  char * argp_program_version = "distmatrix 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Takes a .raw input file and computes the distance matrix between"
  "the trials. Storage as ascii-file.";

static char   args_doc[] = "<input-file>";

/* the cmd-line parser */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct cmdargs     args; /* the globally available cmd-line args */

int main( int argc, char **argv ){ 
  EEGdata_trials *eeg;
  double **Delta, **tmp;
  int N;
  int i,j;

  /* Parse the arguments */
  args_setdefault( &args );
  argp_parse (&argp, argc, argv, 0, 0, &args);
  if( strlen( args.output_filename )==0 ){
	 strcpy( args.output_filename, args.eeg_filename );
	 sprintf( args.output_filename, "%s.dst", args.output_filename );
  }

  /* get data */
  eeg=read_eegtrials_from_raw( (const char*) args.eeg_filename );
  print_eegdata_trials(stderr, eeg);

  N = eeg->ntrials;
  
  /* if num_channels==0, use all chan-indices */
  if( args.num_channels==0 ){
	 args.channels=(int*)malloc( eeg->data[0]->nbchan*sizeof(int) );
	 args.num_channels = eeg->data[0]->nbchan;
	 for( i=0; i<args.num_channels; i++ )
		args.channels[i] = i;
  }

  Delta = matrix_init( N, N );
  tmp   = matrix_init( N, N );
  
  for( j=0; j<args.num_channels; j++ ){ /* average over channels */
	 fprintf( stderr, "  Channel %i\n", args.channels[j] );
	 if( args.channels[j] < 0 || args.channels[j]>=eeg->data[0]->nbchan ){
		warnprintf("  Channel %i not in EEG-file, skipping...\n", args.channels[j] );
		continue;
	 }
	 if( args.metric == vectordist_dtw ){
		tmp = eegtrials_distmatrix_channel( eeg, args.metric, args.channels[j], tmp, signaldist_euclidean_derivative );
	 } else {
		tmp = eegtrials_distmatrix_channel( eeg, args.metric, args.channels[j], tmp, NULL );
	 }
	 matrix_add_matrix( Delta, (const double**)tmp, N, N );
  }
  matrix_divide_scalar( Delta, N, N, (double)args.num_channels );

  fprintf( stderr, "Writing Distance matrix to '%s'\n", args.output_filename );
  write_double_matrix_ascii_file( args.output_filename, (const double**)Delta, N, N );
	
  /* cleaning up */
  matrix_free( Delta, eeg->ntrials );
  matrix_free( tmp  , eeg->ntrials );
  free_eegdata_trials( eeg );
  args_free( args );

  return 0;
}
/* ------------- Argument parsing ------------------------------ */

error_t parse_opt (int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct cmdargs *arguments = state->input;
  int count, i,tmp;
  char *comma;

  switch (key){
  case 'o': /* output */
	 fprintf(stderr, "Output to '%s'\n", arg );
	 strncpy( arguments->output_filename, arg, BUFSIZE-1 );
	 arguments->output_filename[BUFSIZE-1]='\0';
	 break; 
  case 'm': /* metric */
	 fprintf(stderr, "Metric: '%s'\n", arg );
	 if( !strcmp(arg, "euclid") ){
		fprintf(stderr, "Recognized euclidean-metric\n" );
		arguments->metric=vectordist_euclidean;
	 } else if( !strcmp( arg, "euclid_normalized" ) ){
		arguments->metric=vectordist_euclidean_normalized;
	 } else if( !strcmp( arg, "dtw" ) ){
		arguments->metric=vectordist_dtw;
	 } else {
		fprintf( stderr, "Did not recognize this metric, choose default\n");
	 }
	 break;
  case 'c': /* channels */
	 fprintf(stderr," using channels '%s'\n", arg);
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
  case ARGP_KEY_ARG:	
	 fprintf(stderr, "EEG file: %s\n", arg );
	 strncpy( arguments->eeg_filename, arg, BUFSIZE-1 );
	 arguments->eeg_filename[BUFSIZE-1] = '\0';
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

void args_setdefault(  struct cmdargs *arg ){
  arg->output_filename[0]='\0';
  arg->metric=vectordist_euclidean; 
  arg->channels = NULL;
  arg->num_channels=0;
}
  
void args_free ( struct cmdargs args ){
  int i;
  if( args.channels ){
	 free( args.channels );
  }
}
