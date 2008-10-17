/** \example t_padtw.c
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
#define noPL(...)
#ifdef HAVE_LIBPLOTTER
#include <cplotter.h>
static int plotit=0;
#define PL(...) if(plotit){ __VA_ARGS__ }
void plot_eegdata( double *times, EEGdata *s, ChannelInfo *channelcoords, const char *format );
#else
#define PL(...)
#endif


void eeg_ADTW_markers_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta);
WarpPath* eeg_DTW_get_paths_by_markers2( const EEGdata *s1, const EEGdata *s2, int channel, double theta, WarpPath *P );
void eeg_ADTW_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta );
double* ADTW_from_path2(const double *u, int J, const double *s, int K, const WarpPath *P, double *avg);

/* ------------- Argument parsing ------------------------------ */
error_t parse_opt (int key, char *arg, struct argp_state *state);

static struct argp_option options[] = {
  {"output",   'o', "FILE", 0,
   "Output file (default: no output)" }, 
  {"distmatrix", 'm', "FILE", 0, "file containing distance matrix"},
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
  char *distfile; /** filename for distmatrix */
  double theta;    /** restriction parameter */
  char  *output;   /** file to print the output to */
};
void args_free( struct cmdargs );
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
  FILE *out;
  char buffer[255];
  EEGdata_trials *eeg;
  EEGdata *new, *avg, *savg;
  Dendrogram *T, *Tsub;
  double **d, /* distance matrix, NxN */
	 **cumd; /* cumulated distance matrix */
  int num_chan, use_channel,
	 nsamples, nmarkers,
	 N; /* number of trials */
  int i, j, idx1, idx2, c,
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
  args.distfile=NULL;
  argp_parse (&argp, argc, argv, 0, 0, &args);
  
  fprintf( stderr, "reading file          : %s\n", args.input );
  fprintf( stderr, "restriction parameter : %f\n", args.theta );
  fprintf( stderr, "output to             : %p\n", args.output);
  fprintf( stderr, "distfile              : %p\n", args.distfile);
  fprintf( stderr, "Plotting enabled      : %i\n", plotit);

  /* get data */
  eeg=read_eegtrials_from_raw(args.input);
  print_eegdata_trials(stderr, eeg);

  /* prepare settings */
  num_chan = eeg->data[0]->nbchan;
  nsamples = eeg->data[0]->n;
  nmarkers = eeg->nmarkers_per_trial;
  N = eeg->ntrials;
  use_channel=0;

  /* build a huge distance matrix over the average of the channels
	  for all trials; build only if no distance-file is provided
  */
  if( !args.distfile ){ /* compute directly */
	 d    = matrix_init( N, N );
	 cumd = matrix_init( N, N );
	 for( chan=0; chan<num_chan; chan++ ){
		for( i=0; i<N-1; i++ ){
		  for( j=i+1; j<N; j++ ){
			 //oprintf("Channel=%i, clusterdist( %i, %i ) of %i trials\n", chan, i, j, N);
			 d[i][j] = clusterdist_euclidean_pointwise( eeg->data[i], eeg->data[j], chan /*channel*/ );
			 //d[i][j] = clusterdist_tw_complete( eeg->data[i], eeg->data[j], chan /*channel*/ );
			 d[j][i] = d[i][j];
		  }
		}
		matrix_add_matrix( cumd,  (const double**)d, N, N );
	 }
	 
	 matrix_divide_scalar( cumd, N, N, (double)num_chan );
	 //matrix_print(cumd, N, N);
	 strcpy( buffer, (const char*)args.output );
	 strcat( buffer, "_dist.txt" );
	 write_double_matrix_ascii_file( buffer, cumd, N, N );
  } else { /* read from file */
	 cumd = read_double_matrix_ascii( args.distfile, N, N, NULL );
  }


  /* now build hierarchical dendrogram */
  T = agglomerative_clustering( (const double**)cumd, N, dgram_dist_completelinkage );
  //dgram_print( T );

  /* for comparison, get the simple avg */
  savg = eegtrials_simple_average( eeg, NULL );
  
  
  /* now walk the tree to find pairs of trials to match */
  trials_left = N;
  while( trials_left >= 2 ){
	 Tsub = dgram_get_deepest( T );
	 idx1 = Tsub->left->val;
	 idx2 = Tsub->right->val;

	 oprintf("trials_left = %i, Tsub=%p\n", trials_left, Tsub);
	 oprintf("Tsub=%p, val=%i, lval=%i, rval=%i, d[%i,%i]=%f\n", Tsub, Tsub->val, Tsub->left->val, 
				Tsub->right->val,  idx1,  idx2, d[idx1][idx2]);

	 if( eeg->data[idx1]==NULL || eeg->data[idx2]==NULL ){
		errprintf( "try to touch a NULL-node: eeg->data[idx1]=%p, eeg->data[idx2]=%p\n", 
					  eeg->data[idx1], eeg->data[idx2] );
		return;
	 }

	 new = init_eegdata( num_chan, nsamples, nmarkers );

	 /* loop this for all channels ! */
	 for( c=0; c<num_chan; c++){
		dprintf("Channel=%i\n", c);
		eeg_ADTW_channel2( eeg->data[idx1], eeg->data[idx2], new, c, args.theta );
		//eeg_ADTW_markers_channel2( eeg->data[idx1], eeg->data[idx2], new, c, args.theta );
	 }

	 /* remove the previous trials */
	 //	 free_eegdata(eeg->data[idx1]);
	 //	 free_eegdata(eeg->data[idx2]);
	 eeg->data[idx1] = new; /* ADTW goes to idx1 */
	 eeg->data[idx2] = NULL; /* do not touch this again */ 

	 /* replace node by leaf representing ADTW(idx1, idx2) */
	 Tsub->val = idx1;
	 Tsub->left = NULL;
	 Tsub->right = NULL;

	 trials_left--;
  }

  avg = eeg->data[idx1];

  if(args.output){
	 out = fopen( args.output, "w" );
	 write_eegdata_ascii( out, savg );
	 write_eegdata_ascii( out, avg );
	 fclose( out );
  }

  PL( 
	  //plot_eegdata( eeg->times, avg, channelcoords_64, "r-" ); 
	  //plot_eegdata( eeg->times, savg, channelcoords_64, "b-" ); 
	  plot_format(NULL, avg->d[use_channel], avg->n, "r");
	  plot_format(NULL, savg->d[use_channel], savg->n, "b");
	  plot_show(); 
		);
 



  /* cleaning up */
  dprintf("Freeing Memory\n");
  free_eegdata(savg);
  dgram_free( T );
  matrix_free( d, N );
  matrix_free( cumd, N );
  args_free( args );

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
  case 'm':
	 dprintf(" using distmatrix-file '%s'\n", arg);
	 arguments->distfile = (char*)malloc((strlen(arg)+1)*sizeof(char));
	 arguments->distfile = strcpy(arguments->distfile, arg);
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

void args_free(  struct cmdargs arg){
  if(arg.output){
	 free(arg.output);
  }
  if( arg.distfile ){
	 free(arg.distfile);
  }
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
  double *tmpx, *tmpy;
  int nm, i;
  tmpx = (double*) malloc( (s->nmarkers+2)*sizeof( double ) );
  tmpy = (double*) malloc( (s->nmarkers+2)*sizeof( double ) );

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

	 nm = s->nmarkers;
	 tmpx[0] = times[0];
	 tmpy[0] = s->d[chan][0];
	 tmpx[nm-1] = times[s->n-1];
	 tmpy[nm-1] = s->d[chan][s->n-1];
	 for( i=1; i<nm+2; i++){
		tmpx[i] = times[s->markers[i-1]];
		tmpy[i] = s->d[chan][s->markers[i-1]];
	 }
	 plot_format( tmpx, tmpy, nm+2, "ko.");
	 plot_subplot_select(-1);
  }
  free(tmpx);
  free(tmpy);
}
#endif



/** compute the marker-based ADTW for one channel in s1,s2 and put it into target.
	 \param s1,s2 sources
	 \param target output
	 \param theta restriction parameter (see \ref timewarping)
 */
void eeg_ADTW_markers_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta){
  int nparts, part, i, N, J, K, n;
  int pos_1, pos_2, pos_t;
  WarpPath *Pptr, *Pptr2;
  unsigned long *mark1, *mark2, *markt; /* convenience, to include 0 and n */

  if( eegdata_cmp_settings( s1, s2 ) ){
	 errprintf("cannot ADTW the to EEGdata structs, because they are too different\n");
	 return;
  }

  target->nmarkers = s1->nmarkers;
  n = s1->n;
  N = target->nmarkers;
  nparts = target->nmarkers+1;

  Pptr = (WarpPath*) malloc( nparts*sizeof( WarpPath ) );
  for( part=0; part<nparts; part++ ){
	 Pptr2 = Pptr+part;
	 Pptr2 = init_warppath( Pptr2, n, n );
  }

  /* add trivial markers for convenience in the loop */
  mark1 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark2 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  markt = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark1[0]   = 0;   mark2[0]   = 0;   markt[0]   = 0;
  mark1[N+1] = n-1; mark2[N+1] = n-1; markt[N+1] = n-1;
  memcpy( mark1+1, s1->markers, N*sizeof(unsigned long) );
  memcpy( mark2+1, s2->markers, N*sizeof(unsigned long) );

  /* set new markers as (J+K)/2+1 */
  target->n = n;
  for( i=1; i<N+1; i++ ){
		target->markers[i-1] = ( ( (mark1[i])+
											(mark2[i]) ) )/2 + 1;
		markt[i] = target->markers[i-1];
  }

  for( i=0; i<N+2; i++){
	 //	 dprintf( "s1->markers[%i] = %i \t s2->markers[%i] = %i\n", i, s1->markers[i], i, s2->markers[i]);
	 dprintf("nm1[%i] = %i | nm2[%i] = %i | nmt[%i] = %i\n", i, mark1[i], i, mark2[i], i, markt[i]);
  }


  /* DEBUG */
  print_eegdata( stdout, target );

  	/* we get all paths simultaneously */
  Pptr = eeg_DTW_get_paths_by_markers2( s1, s2, channel, theta, Pptr );


  for( part=0; part<nparts; part++ ){
	 dprintf("Part No.=%i/%i\n",  part, nparts);
	 J = mark1[part+1]-mark1[part];
	 K = mark2[part+1]-mark1[part];
	 pos_1 = mark1[part];
	 pos_2 = mark2[part];
	 pos_t = markt[part];
	 dprintf("(J, K, p1, p2, pt)=(%i,%i,%i,%i,%i)\n", J, K,pos_1,pos_2,pos_t);
	 ADTW_from_path( s1->d[channel]+pos_1, J, /* source 1 */
						  s2->d[channel]+pos_2, K, /* source 2 */
						  Pptr+part,               /* which path */
						  target->d[channel]+pos_t );  /* write it to */
  }
  
  free( mark1 ); 
  free( mark2 );
  free( markt );
  for( i=0; i<nparts; i++ ){
	 dprintf("i=%i, Pptr+i=%p\n",i, (Pptr+i));
	 free( (Pptr+i)->upath );
	 free( (Pptr+i)->spath );
  }
  free(Pptr);
}


/** compute multiple warp-Pathes by computing restricted warppathes 
	 with restriction parameter \f$\theta\f$ for each of the segments 
	 from marker 0...i, i...i+1, ..., N-1...N (see \ref timewarping).
	 \param s1,s2 data
	 \param channel channel to use in the dataset
	 \param theta restriction parameter (Chiba-Band)
	 \param P pointer to Warppath-structs (nmarker many); if NULL, own memory is allocated
 */
WarpPath* eeg_DTW_get_paths_by_markers2( const EEGdata *s1, const EEGdata *s2, int channel, double theta, WarpPath *P ){
  int i,j, N,n, offset1, offset2;
  int J,K;
  WarpPath *Pptr;
  unsigned long *mark1, *mark2; /* convenience, to include 0 and n */
  double **dist;

  if( s1->nmarkers != s2->nmarkers || s1->n != s2->n || s1->nbchan < channel || s2->nbchan<channel){
	 errprintf("not the same number of markers (%i, %i) or something else\n", s1->nmarkers, s2->nmarkers );
	 return P;
  }
  if( theta>1 ){
	 errprintf( "Theta > 1, choose theta=1\n");
	 theta=1;
  } else if(theta<0){
	 errprintf( "Theta < 0, choose theta=0\n");
	 theta=0;
  }
  N = s1->nmarkers;
  n = s1->n;

  /* add trivial markers for convenience in the loop */
  mark1 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark2 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark1[0]   = 0;   mark2[0]   = 0;
  mark1[N+1] = n-1; mark2[N+1] = n-1;
  memcpy( mark1+1, s1->markers, N*sizeof(unsigned long) );
  memcpy( mark2+1, s2->markers, N*sizeof(unsigned long) );
  dist = matrix_init( n, n );

  if( P==NULL ){
	 warnprintf("Warning: memory allocated within function\n");
	 P = (WarpPath*) malloc( (N+1)*sizeof( WarpPath ) );
	 for( i=0; i<N+1; i++){	
		J = mark1[i+1] - mark1[i];
		K = mark2[i+1] - mark2[i];
		Pptr = P+i;
		Pptr = init_warppath( Pptr, J, K );
	 }
  }

  // plot_format( NULL, s1->d[channel], s1->n, "r");
  // plot_format( NULL, s2->d[channel], s1->n, "b");
  //plot_show();

  for( i=0; i<N+2; i++){
	 //	 dprintf( "s1->markers[%i] = %i \t s2->markers[%i] = %i\n", i, s1->markers[i], i, s2->markers[i]);
	 dprintf("nm1[%i] = %i | nm2[%i] = %i\n", i, mark1[i], i, mark2[i]);
  }
  N += 2;
  Pptr = P;
  const char *cols[] = {"r.", "g.", "b.", "k."};
  for( i=1; i<N; i++ ){
	 J = mark1[i] - mark1[i-1];
	 K = mark2[i] - mark2[i-1];
	 offset1 = mark1[i-1];
	 offset2 = mark2[i-1];
	 dprintf("i=%i: mark1[%i]=%i, mark1[%i]=%i,  J,K = (%i,%i), Pptr=%p\n", i, i-1, mark1[i-1], i, mark1[i], J, K, Pptr);

	 /* prepare warppath */
	 dprintf("resttingin Pptr=%p to J=%i, K=%i\n", Pptr, J, K);
	 reset_warppath( Pptr, J, K );
	 

	 dist = DTW_build_restricted_cumdistmatrix( s1->d[channel]+offset1, J, 
															  s2->d[channel]+offset2, K, 
															  theta, dist );  
	 //	 plot_image( dist, J, K, "bone");
	 Pptr = DTW_path_from_cumdistmatrix( (const double**) dist, J, K, Pptr);

	 dprintf("Pptr=%p 1st=(%i,%i), 2nd=(%i,%i)\n", Pptr, Pptr->upath[0], Pptr->spath[0], Pptr->upath[1], Pptr->spath[1]);
	 noPL(	 plot_format_int( Pptr->upath, Pptr->spath, Pptr->K+Pptr->J, cols[i-1]);
	 plot_format( NULL, s1->d[channel]+offset1, J, "k" );
	 plot_format( NULL, s2->d[channel]+offset2, K, "k" );
	 plot_switch(plot_add());
				 );
	 Pptr = P+i;
  }

  //  plot_show();
  free(mark1); free(mark2);
  matrix_free(dist, n);

  return P;
}

/** compute the  ADTW for one channel in s1,s2 and put it into target. Ignore time-markers.
	 \param s1,s2 sources
	 \param target output
	 \param channel index to use
	 \param theta restriction parameter (see \ref timewarping)
 */
void eeg_ADTW_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta ){
  double **dist;
  WarpPath *P;

  dist = matrix_init( s1->n, s2->n );

  dist = DTW_build_restricted_cumdistmatrix( s1->d[channel], s1->n,  s2->d[channel], s2->n, theta, dist );
  P = DTW_path_from_cumdistmatrix( dist, s1->n, s2->n, NULL );
  target->d[channel] = ADTW_from_path2( s1->d[channel], s1->n, s2->d[channel], s2->n, P, target->d[channel] );

  noPL(  
		 plot_format( NULL, s1->d[channel], s1->n, "b" );
		 plot_format( NULL, s2->d[channel], s2->n, "b" );
		 plot_format_int( P->upath, P->spath, P->J+P->K, "r.");
		 plot_format( NULL, target->d[channel], target->n, "g" );

		 plot_show();  
			);

  matrix_free( dist, s1->n );
  free_warppath( P );
}

double* ADTW_from_path2(const double *u, int J, const double *s, int K, const WarpPath *P, double *avg){
	double *tmp;
	int i, idx;

	tmp = (double*)calloc(K+J, sizeof(double));
	if(avg==NULL){
	  dprintf("Allocating own memory\n");
	  avg = (double*)calloc((K+J)/2+1, sizeof(double));
	}
	
	idx = 0;
	for(i=0; i<J+K; i++){
	  if( P->upath[i]==0 && P->spath[i]==0 ){
		 continue;
	  }
	  tmp[idx++] = (u[P->upath[i]] + s[P->spath[i]])/2.0;
	}
	//	dprintf("J=%i,K=%i, idx=%i\n", J, K, idx);

	tmp = flip_array( tmp, idx );
	avg = resample_nearest_neighbour( tmp, idx, (J+K)/2+1, avg );
	//	avg = resample_linear( tmp, idx, (J+K)/2+1, avg );

	noPL(
		  plot_format( NULL, u, J, "k" );	
		  plot_format( NULL, s, K, "k" );	
		  plot_format( NULL, tmp, idx, "b" );	
		  plot_format( NULL, avg, (J+K)/2+1, "g" );	
		  plot_show();
		  );
	free(tmp);
	return avg;
}
