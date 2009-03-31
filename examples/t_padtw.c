/** \example t_padtw.c
 *
 * \brief Computing the Hierarchical DTW-average of many raw-files.
 Takes a list of "similar" .raw-files (same condition across
 subjects), computes the PADTW for each of them and merges them
 into a new .raw file with one trace for each subject.
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
/* #include <libplotter/cplotter.h> */

/* ------------- Argument parsing ------------------------------ */
error_t parse_opt (int key, char *arg, struct argp_state *state);

static struct argp_option options[] = {
  {"output", 'o', "filename", 0, "specifies output file (default 'padtw_output.raw')"}, 
  {"sigma", 's', "double in [0,1]", 0, "maximal sigma for regularization [0,1]; default 0.2"},
  {"metric", 'm', "string", 0, "{'locfreq', 'derivative_euclid'}"},
  {"writedist", 'w', 0, 0, "write distance matrix Delta to <filename>.dst (default 0)"},
  {"readdist", 'r', 0, 0, "try to read distance matrix Delta from <filename>.dst\n"
  "If it is not found, it is computed (warning message). (default 0)"},
  { 0 }
};

/* available cmd-line args */
struct cmdargs {
  char *output_filename;
  char **eeg_filenames;
  int  num_eeg_files;
  double sigma;
  int  metric; 					  /* 0-locfreq, 1-derivative_euclid */
  int readdist;					  /* flag 0,1 */
  int writedist;					  /* flag 0,1 */
};
void args_setdefault(  struct cmdargs *arg, int max_number_of_eeg_files );
void args_free ( struct cmdargs );

const  char * argp_program_version = "padtw 0.1";
const  char * argp_program_bug_address = "<mihrke@uni-goettingen.de>";
static char   doc[] = "Takes a list of 'similar' .raw-files (same condition across subjects), "
	 "computes the PADTW for each of them and merges them into a new "
	 ".raw file with one trace for each subject.  No time-markers are "
	 "written in this file.";
static char   args_doc[] = "<input-files>";

/* the cmd-line parser */
static struct argp argp = { options, parse_opt, args_doc, doc };
struct cmdargs     args; /* the globally available cmd-line args */


/* ---------------------------------------------------------------------------- 
   -- MAIN                                                           -- 
   ---------------------------------------------------------------------------- */
int main( int argc, char **argv ){
  EEGdata_trials *eeg, *outeeg;
  double **Delta=NULL, **tmp;
  char dstfname[255];
  char *rawfile;
  int i, j;
  SettingsHierarchicalDTW settings;
  int num_channels, num_markers, num_samples;
  double corner_freqs[2]={0.0, 25.0}; /* we have filtered data */

  /* Parse the cmd-line arguments */
  args_setdefault( &args, argc );
  argp_parse (&argp, argc, argv, 0, 0, &args);


  for( i=0; i<args.num_eeg_files; i++ ){ /* for each .raw file from the list */
	 rawfile = args.eeg_filenames[i];
	 fprintf( stderr, "Processing file '%s'\n", rawfile );
	 eeg=read_eegtrials_from_raw( rawfile );
	 print_eegdata_trials( stderr, eeg );
	 if( i==0 ){					  /* first file */
		num_channels = eeg->data[0]->nbchan;
		num_markers = eeg->nmarkers_per_trial;
		num_samples = eeg->nsamples;  
		/* prepare output struct */
		outeeg = init_eegdata_trials( args.num_eeg_files, /* trials */
												num_markers,					  /* markers */
												num_channels,
												num_samples,
												eeg->times );
		fprintf(stderr, "OutEEG:\n");
		print_eegdata_trials( stderr, outeeg );
	 } else {						  /* check whether other files are similar to first file */
		if( (num_channels != eeg->data[0]->nbchan) ||
			 (num_markers != eeg->nmarkers_per_trial) ||
			 (num_samples != eeg->nsamples) ){
		  fprintf( stderr, "ERROR: Raw-file '%s' does not match the first file '%s'\n Aborting...\n",
					  rawfile, args.eeg_filenames[0] ); 
		}
	 }

	 /* compute or read distance matrix between trials */	
	 strncpy( dstfname, rawfile, strlen(rawfile)-3 );
	 dstfname[strlen(rawfile)-3+0]='d';
	 dstfname[strlen(rawfile)-3+1]='s';
	 dstfname[strlen(rawfile)-3+2]='t';
	 if( args.readdist ){
		Delta=read_double_matrix_ascii( dstfname, eeg->ntrials, eeg->ntrials, ALLOC_IN_FCT );
	 } 
	 if( !Delta ) {				  /* failed to read, or need to calculate */
		fprintf(stderr, "compute distance matrix between trials\n" );
		Delta = matrix_init( eeg->ntrials, eeg->ntrials );
		tmp   = matrix_init( eeg->ntrials, eeg->ntrials );
		for( j=0; j<num_channels; j++ ){ /* average over channels */
		  fprintf( stderr, "  Channel %i\n", j );
		  tmp = eegtrials_distmatrix_channel( eeg, vectordist_euclidean, j, tmp );
		  matrix_add_matrix(Delta, (const double**)tmp, eeg->ntrials, eeg->ntrials );
		}
		matrix_divide_scalar( Delta, eeg->ntrials, eeg->ntrials, (double)eeg->ntrials ); 
	 }
	 if( args.writedist ){
		fprintf( stderr, "Writing Distance matrix to '%s'\n", dstfname );
		write_double_matrix_ascii_file( dstfname, (const double**)Delta, eeg->ntrials, eeg->ntrials );
	 }

	 /* compute PADTW */
	 fprintf(stderr, " compute PADTW No. %i for file %s\n", i, rawfile );
	 
	 settings = init_dtw_hierarchical( eeg );
	 settings.sigma=args.sigma;
	 settings.corner_freqs[0] = corner_freqs[0];
	 settings.corner_freqs[1] = corner_freqs[1];
	 if( args.metric==0 ){		  /* locfreq */
		settings.pointdistance=eeg_distmatrix_stft_channel;
	 } else {
		settings.pointdistance=eeg_distmatrix_euclidean_derivative_channel;
	 }
	 settings.progress = progressbar_rotating;
	 print_settings_hierarchicaldtw( stderr, settings );
	 fprintf( stderr, "--> Progress: " );
	 eegtrials_dtw_hierarchical( eeg, (const double**)Delta, eeg->ntrials, outeeg->data[i], settings );

	 /* cleaning up */
	 matrix_free( Delta, eeg->ntrials );
	 matrix_free( tmp  , eeg->ntrials );
	 free_eegdata_trials( eeg );
  }
  
  /* write output to file */
  fprintf(stderr, "Writing Results to '%s'\n", args.output_filename );  
  for( i=0; i<outeeg->ntrials; i++ ){
	 for( j=0; j<outeeg->nmarkers_per_trial; j++ ){
		outeeg->markers[i][j]=outeeg->data[i]->markers[j];
	 }
  }
  print_eegdata_trials( stderr, outeeg );
  write_eegtrials_to_raw_file( outeeg, args.output_filename );
  fprintf( stderr, " Done\n");

  free_eegdata_trials( outeeg );
  args_free( args );
  return 0;
}




/* ------------- Argument parsing ------------------------------ */

error_t parse_opt (int key, char *arg, struct argp_state *state){
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct cmdargs *arguments = state->input;

  switch (key){
  case 'o': /* output */
	 fprintf(stderr, "Output to '%s'\n", arg );
	 arguments->output_filename=arg;
	 break; 
  case 's': /* output */
	 fprintf(stderr, "Maximal sigma '%s'\n", arg );
	 arguments->sigma=atof(arg);
	 break;
  case 'r':
	 fprintf(stderr, "Read distance file\n");
	 arguments->readdist=1;
	 break;
  case 'w':
	 fprintf(stderr, "Write distance file\n");
	 arguments->writedist=1;
	 break;
  case 'm': /* output */
	 fprintf(stderr, "Metric: '%s'\n", arg );
	 if( !strcmp(arg, "locfreq") ){
		fprintf(stderr, "Recognized local-frequency-metric\n" );
		arguments->metric=0;
	 } else if( !strcmp( arg, "derivative_euclid" ) ){
		arguments->metric=1;
	 } else {
		fprintf( stderr, "DId not recognize this metric, choose default\n");
	 }
	 break;
  case ARGP_KEY_ARG:	
	 fprintf(stderr, "EEG file No. %i: %s\n", arguments->num_eeg_files+1,  arg );
	 arguments->eeg_filenames[arguments->num_eeg_files] = (char*) malloc( (strlen(arg)+1)*sizeof( arg ) );
	 strcpy( arguments->eeg_filenames[arguments->num_eeg_files], arg );
	 arguments->num_eeg_files += 1;
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

void args_setdefault(  struct cmdargs *arg, int max_number_of_eeg_files ){
  arg->output_filename="padtw_output.raw";
  arg->num_eeg_files=0;
  arg->eeg_filenames=(char**) malloc( max_number_of_eeg_files*sizeof( char* ) );
  arg->sigma=0.2;
  arg->metric=0;
  arg->readdist=0;
  arg->writedist=0;
}
  
void args_free ( struct cmdargs args ){
  int i;
  for( i=0; i<args.num_eeg_files; i++ ){
	 free( args.eeg_filenames[i] );
  }
  if( args.eeg_filenames ){
	 free( args.eeg_filenames );
  }
}
