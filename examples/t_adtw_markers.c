/** \example t_adtw_markers.c
 *
 * \brief Computing the ADTW of the first two trials within a raw-file.
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
#include "averaging.h"
#include "mathadd.h"

#include <cplotter.h>

void eeg_ADTW_markers_channel3(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta);

void eeg_ADTW_markers_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta);
WarpPath* eeg_DTW_get_paths_by_markers2( const EEGdata *s1, const EEGdata *s2, int channel, double theta, WarpPath *P );
void eeg_ADTW_channel2(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta );
double* ADTW_from_path2(const double *u, int J, const double *s, int K, const WarpPath *P, double *avg);


double** test( const double *s1, int n1, const double *s2, int n2, double theta, double **d){
  
}

/* ---------------------------------------------------------------------------- 
   -- main routine                                                           -- 
   ---------------------------------------------------------------------------- */
int main(int argc, char **argv){ 
  FILE *f;
  char buffer[255];
  EEGdata_trials *eeg;
  EEGdata *new;

  /* get data */				  
  eeg=read_eegtrials_from_raw( argv[1] );
  print_eegdata_trials(stderr, eeg);
  /* new = init_eegdata( eeg->data[0]->nbchan, eeg->data[0]->n, eeg->nmarkers_per_trial ); */

  /* computation */  
  /*  test(*/
  /* eeg_ADTW_markers_channel3( eeg->data[0], eeg->data[1], new, 0, 0.5 ); */


  /* write */
  f = fopen( "test.raw", "w");
  write_raw_header( f, 1, 1, new->n, new->nmarkers );
  ffwrite( eeg->times, sizeof(double), eeg->data[0]->n, f );
  write_raw_markers( f, new );
  ffwrite( new->d[0], sizeof(double), new->n, f );
  fclose(f);
  fprintf( stderr, "WRITTEN\n");

  /* cleaning up */
  free_eegdata_trials( eeg );
  free_eegdata( new );

  return 0;
}

/** compute the marker-based ADTW for one channel in s1,s2 and put it into target.
	 \param s1,s2 sources
	 \param target output
	 \param theta restriction parameter (see \ref timewarping)
 */
void eeg_ADTW_markers_channel3(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double R){
  int i,j,k;
  double avgu, avgs, rmsu=0.0, rmss=0.0;
  double snorm, snormp; /* variable for a point of the normalized
									signal and the preceding point */
  double unorm, unormp;
  int theta;
  double left,down,downleft;
  const double *u, *s;
  int J,K;
  double **d;
  int *bham;
  char c;
  WarpPath *P;

  u = s1->d[channel];
  s = s2->d[channel];
  J = s1->n;
  K = s2->n;

  d = (double**) malloc(J*sizeof(double*));
  for(j=0; j<J; j++)
	 d[j]=(double*)malloc(K*sizeof(double));

  bham = (int*)malloc( MAX( J,K )*2*sizeof( int ) );
  bresenham(0,0, J-1, K-1, bham);
    
  if( R>1 ){
	 dprintf("Restriction R=%f too large, using 1.0\n");
	 R = 1.0;
  } else if( R<0 ){
	 dprintf("Restriction R=%f < 0, aborting\n");
	 return;
  }

  theta = (int)floor( R*MIN( K, J) );
  /* theta = (int)floor( ( R*sqrt( SQR( J )+ SQR( K ) ) )/2.0 ); */
  dprintf("theta=%i pixels\n", theta);


  avgu = gsl_stats_mean(u, 1, J);
  avgs = gsl_stats_mean(s, 1, K);
  for(j=0; j<J; j++) rmsu += SQR( u[j] );
  rmsu = sqrt(rmsu/(double)J);
  for(k=0; k<K; k++) rmss += SQR( s[k] );
  rmss = sqrt(rmss/(double)K);

  for(j=0; j<J; j++){ // set everything to NAN for restrictions
    for(k=0; k<K; k++){
  		d[j][k] = NAN;
  	 }
  }

  plot_image( d, K, J, "hot" );
  
  int b = 1;
  if( K>J ) b=0;

  int lower_corridor, upper_corridor;


  /* computing d_jk */
  for( j=0; j<MAX( J, K ); j++ ){ /* J>K */
	 lower_corridor = MAX( 0, bham[(2*j)+b]-theta );
	 upper_corridor = MIN( bham[(2*j)+b]+theta, K );
	 /* dprintf("b=%i, bham=(%i,%i), j=%i, corridor: (%i, %i)\n", */
	 /* 			b, bham[(2*j)+0], bham[(2*j)+1], j, lower_corridor, upper_corridor); */
    for( k= lower_corridor; k<upper_corridor; k++ ){ /* corridor */
      snorm = (s[k]-avgs)/rmss;
      unorm = (u[j]-avgu)/rmsu;
      (k==0) ? (snormp = 0) : (snormp = ((s[k-1]-avgs)/rmss));
      (j==0) ? (unormp = 0) : (unormp = ((u[j-1]-avgu)/rmsu));

		/* swapping necessary for K>J (reswapped after cumulation) */
		if(K>J) swap2i(&j, &k);

		d[j][k] = fabs(unorm - snorm) + fabs( (unorm-unormp) - (snorm-snormp) ); /* (1) */

		/* cumulate matrix, NAN is treated as inf */
		if(k==0 && j==0) ;
      else if(k==0 && j>0) 
		  d[j][k] += d[j-1][k];
      else if(k>0 && j==0) 
		  d[j][k] += d[j][k-1];
      else { /* j,k > 0 */
		  (isnan( d[j-1][k  ] ) )?(left     = DBL_MAX):(left     = d[j-1][k  ]);
		  (isnan( d[j  ][k-1] ) )?(down     = DBL_MAX):(down     = d[j  ][k-1]);
		  (isnan( d[j-1][k-1] ) )?(downleft = DBL_MAX):(downleft = d[j-1][k-1]);

		  d[j][k] += MIN(MIN(left, down), downleft);
		}

		/* reswap */
		if(K>J) swap2i(&j, &k);
    }
  }
  plot_switch( plot_add() );
  plot_image( d, K, J, "hot" );


  


  plot_show( );

  /* cleaning up */
  free( bham );
  for( j=0; j<J; j++)
	 free( d[j] );
  free( d );
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
	//avg = resample_nearest_neighbour( tmp, idx, (J+K)/2+1, avg );
		avg = resample_linear( tmp, idx, (J+K)/2+1, avg );
	free(tmp);
	return avg;
}
