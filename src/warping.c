/***************************************************************************
 *   Copyright (C) 2008 by Matthias Ihrke   *
 *   mihrke@uni-goettingen.de   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "warping.h"

/** This function applies regularization matrix R to distance matrix
	 d using
	 \f[
	 \hat{d}_{jk} = \mbox{max}(d_{jk}) - (  [\mbox{max}(d_{jk})-d_jk] \cdot R_{jk} )
	 \f]
	 \param d distance matrix 
	 \param R regularization matrix
	 \param M,N dimensions of d and R
 */
void      dtw_regularize_matrix( double **d, const double **R, int M, int N ){
  double maxdist;

  maxdist = matrix_max( (const double**)d, M, N, NULL, NULL );
  scalar_minus_matrix( maxdist, d, M, N );
  matrix_dottimes_matrix( d, (const double**)R, M, N );
  scalar_minus_matrix( maxdist, d, M, N ); 
}

/** cumulate a distance matrix d to give
	 \f[
	 \D_{jk} = \mathbf{d}_{jk}+\min{\{\D_{j,k-1}, \D_{j-1,k}, \D_{j-1, k-1}\}}
	 \f]
	 \param d input/output matrix
	 \param M,N dimensions of d
 */
void      dtw_cumulate_matrix  ( double **d, int M, int N ){
  /* j = 0,...,M-1
	  k = 0,...,N-1
  */
  int j, k;
  
   /* computing D_jk */
  for(j=0; j<M; j++){
    for(k=0; k<N; k++){
      if(k==0 && j==0) ;
      else if(k==0 && j>0)
		  d[j][k] += d[j-1][k];
      else if(k>0 && j==0)
		  d[j][k] += d[j][k-1];
      else /* j,k > 0 */
		  d[j][k] += MIN(MIN(d[j][k-1], d[j-1][k]), d[j-1][k-1]);
    }
  }
}
/**

 */
WarpPath* dtw_backtrack        ( const double **d, int M, int N, WarpPath *P ){ 
  int j,k; /* j = 0,...,M-1
				  k = 0,...,N-1  */
  int idx;
  double left, down, downleft;

  if( P==NULL ){
	 P = init_warppath( NULL, M, N );
  }

  /* Backtracking */
  j=M-1; k=N-1;

  idx = 1;
  P->t1[0] = j;
  P->t2[0] = k;
  while( j>0 || k>0 ){
	 if( k==0 ){ /* base cases */
		j--;
	 } else if( j==0 ){
		k--;
	 } else { /* min( d[j-1][k], d[j-1][k-1], d[j][k-1] ) */
	  
		left     = d[j-1][k  ];
		down     = d[j  ][k-1];
		downleft = d[j-1][k-1];

		(isnan( d[j-1][k  ] ) )?(left     = DBL_MAX):(left     = d[j-1][k  ]);
		(isnan( d[j  ][k-1] ) )?(down     = DBL_MAX):(down     = d[j  ][k-1]);
		(isnan( d[j-1][k-1] ) )?(downleft = DBL_MAX):(downleft = d[j-1][k-1]);


		/* this order of comparisons was chosen to ensure, that diagonal
			steps are preferred in case of equal distances
		*/
		if( downleft<=left ){	  /* downleft or down */
		  if( downleft<=down ){  
			 k--; j--;			 	  /* downleft */
		  } else {
			 k--;						  /* down */
		  } 
		} else {						  /* left or down */
		  if( down<=left ){
			 k--;						  /* down */
		  } else {
			 j--; 					  /* left */
		  }
		}

	 /* 	if( left<=downleft ){	  /\* left or down *\/ */
	 /* 	  if( left <= down ) */
	 /* 		 j--;						  /\* left *\/ */
	 /* 	  else */
	 /* 		 k--;						  /\* down *\/ */
	 /* 	} else {						  /\* downleft or down *\/ */
	 /* 	  if( downleft <= down ){ */
	 /* 		 k--; j--;				  /\* downleft *\/ */
	 /* 	  } else  */
	 /* 		 k--;						  /\* down *\/ */
	 /* 	} */


	 }	/* if */

	 P->t1[idx] = j;
	 P->t2[idx] = k;
	 dprintf("(j, k), idx = (%i,%i), %i\n", j, k, idx);
	 if(idx>0 && P->t1[idx]>=P->t1[idx-1] && P->t2[idx]>=P->t2[idx-1]  ){
		warnprintf("P=%p: idx=%i\n", P, idx);
	 }
	 idx++;
  }
  P->n=idx-1;
  //  print_warppath( stderr, P );

  /* now the warppath has wrong order, reverse */
  for( j=0; j<(P->n)/2; j++ ){
	 swap2i( &(P->t1[j]), &(P->t1[P->n-1-j]) );
	 swap2i( &(P->t2[j]), &(P->t2[P->n-1-j]) );
  }
  //  print_warppath( stderr, P );

  return P;
}
  

/** Add magnitude of signals according to warppath P.
	 Data is resampled to have length (J+K)/2+1 using cspline interpolation.
 * \param path - contains warppath
 * \param weights - weights in average, for using it with hierarchical averaging
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* warp_add_signals_by_path(const double *s1, int n1, 
											const double *s2, int n2, 
											const WarpPath *P, double *avg, 
											const double weights[2]){
  double *x, *y, *xp;
  int newn, i;
  
  newn = (n1+n2)/2;			  /* length of average signal */
  x = (double*)malloc( (n1+n2)*sizeof(double) );
  y = (double*)malloc( (n1+n2)*sizeof(double) );
  xp= (double*)malloc( newn*sizeof(double) );
  if(avg==NULL){
	 dprintf("Allocating own memory\n");
	 avg = (double*)calloc(newn, sizeof(double));
  }

  for( i=0; i<newn; i++ ){
	 xp[i] = (double)i;
  }


  for(i=0; i<P->n; i++){
	 x[i] = (double)((P->t1[i])+(P->t2[i]))/2.0; /* average latencies */
	 y[i] = weights[0]*(s1[P->t1[i]]) + weights[1]*(s2[P->t2[i]]); /* average magnitude */
	 dprintf( "x[%i]=(%f,%f), (%i, %i)\n", i, x[i], y[i], (P->t1[i]), (P->t2[i]) );
	 if( i>0 && x[i]<=x[i-1] ){
		warnprintf("x not monotonic at P=%p x[%i]=(%f,%f), (%i, %i)\n", P, i, x[i], y[i], (P->t1[i]), (P->t2[i]) );
	 }
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, P->n);
  gsl_spline_init (spline, x, y, P->n);

  for( i=0; i<newn; i++ ){
	 dprintf("xp[%i]=%f\n", i, xp[i] ); 
	 avg[i] = gsl_spline_eval( spline, xp[i], acc );
	 dprintf(" avg[%i] = (%f, %f)\n", i,xp[i], avg[i] );
  }

  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free( x );
  free( y );
  free( xp );

  return avg;
}


/** \todo this function is pretty huge. Modularize?

	 \param eeg_in 
	 \param distmatrix Distance matrix between the trials in eeg_in
	 \param N number of trials (columns/rows) in distmatrix
	 \param max_sigma param for regularization
	 \param corner_freqs for local-frequency metric; maximal {0, eeg_in->sampling_rate/2.0}
	 \param dont_touch_eeg flag; if >0, duplicate complete eeg-struct within the function. 
	                       Else, eeg is modified and cleared within the function.
	 \param out memory for the final average or NULL (own memory is allocated)
	 \return pointer to final average
 */  
EEGdata* eegtrials_dtw_hierarchical( EEGdata_trials *eeg_in, const double **distmatrix, int N, 
												 EEGdata *out, SettingsHierarchicalDTW settings ){
  Dendrogram *T, *Tsub; 
  EEGdata_trials *eeg;
  EEGdata *new, *s1, *s2;
  WarpPath *P;
  int num_chan,
	 nsamples, nmarkers, n;
  int i, idx1, idx2, c,
	 trials_left;
  double **G, **d;
  double maxdist;
  double srate;
  double weights[2]={1.0,1.0}; 			  /* this is for recursive averaging */
  int    *indices;				  /* n-array containing number averagings for each trial */
  int *channels;
  int nchan,chan;

  if( N>eeg_in->ntrials ){
	 errprintf(" distance matrix contains more trials than eegdata. this is fatal\n");
	 return NULL;
  }

  if(settings.channels==NULL){
	 channels = linspace(0,num_chan-1);
	 nchan = num_chan;
  } else {
	 channels = settings.channels;
	 nchan = settings.num_channels;
  }

  /* build hierarchical dendrogram */
  T = agglomerative_clustering( (const double**)distmatrix, N, settings.linkage );

  if( settings.dont_touch_eeg ){
	 dprintf("Cloning EEGData\n");
	 eeg = clone_eegdata_trials( eeg_in );
  } else {
	 dprintf("using eeg_in and deleting it!\n");
	 eeg = eeg_in;
  }
  
  /* prepare settings */
  num_chan = eeg->data[0]->nbchan;
  nsamples = eeg->data[0]->n;
  nmarkers = eeg->nmarkers_per_trial;
  srate    = eeg->sampling_rate;
  n = nsamples;
  G = matrix_init( n, n );
  d = matrix_init( n, n );  
  indices = (int*)malloc( N*sizeof(int) );
  for( i=0; i<N; i++ ){
	 indices[i] = 1;
  }
  if( !out ){
	 warnprintf( "allocating eegdata within function \n");
	 out = init_eegdata( num_chan, nsamples, nmarkers );
  }
  P = init_warppath( ALLOC_IN_FCT, n, n );

  idx1=0; 
  idx2=0;
  /* now walk the tree to find pairs of trials to match */
  trials_left = N;
  if( settings.progress ){
	 settings.progress( PROGRESSBAR_INIT, N );
  }
  while( trials_left >= 2 ){
	 if( settings.progress ){
		settings.progress( PROGRESSBAR_CONTINUE_LONG, N-trials_left );
	 }
	 
	 Tsub = dgram_get_deepest( T );
	 idx1 = Tsub->left->val;
	 idx2 = Tsub->right->val;
	 dprintf("trials_left = %i, Tsub=%p\n", trials_left, Tsub);
	 dprintf("Tsub=%p, val=%i, lval=%i, rval=%i, d[%i,%i]=%f\n", Tsub, Tsub->val, Tsub->left->val, 
				Tsub->right->val,  idx1,  idx2, d[idx1][idx2]);


	 weights[0] = indices[idx1]/(double)(indices[idx1]+indices[idx2]);
	 weights[1] = indices[idx2]/(double)(indices[idx1]+indices[idx2]);
	 dprintf("indices[%i,%i]=(%i,%i)\n", idx1, idx2, indices[idx1], indices[idx2]);
	 dprintf("weights=(%f,%f)\n", weights[0], weights[1] );

	 if( eeg->data[idx1]==NULL || eeg->data[idx2]==NULL ){
		errprintf( "try to touch a NULL-node: eeg->data[idx1]=%p, eeg->data[idx2]=%p\n", 
					  eeg->data[idx1], eeg->data[idx2] );
		return NULL;
	 }
	 
	 /* prepare average */
	 new = init_eegdata( num_chan, nsamples, nmarkers );
	 s1  = eeg->data[idx1];
	 s2  = eeg->data[idx2];

	 dprintf(" compute G\n");
	 if( settings.regularize ){
		G = settings.regularize( s1, s2, settings.sigma, G); /* we need G only 
																				  once for all channels */
		if(!G){
		  errprintf("Regularization did not work\n");
		  return NULL;
		}
	 }


  
	 /* loop this for all channels */
	 for( c=0; c<nchan; c++){	
		chan = channels[c];

		if( settings.progress ){
		  //		  oprintf("Trials (%i,%i): Channel=%i\n", idx1, idx2, chan); 
		  settings.progress( PROGRESSBAR_CONTINUE_SHORT, 0 );
		}

		dprintf(" compute d\n");
		if(!(d = settings.pointdistance( s1, s2, chan, d, (void*)(&settings) )) ){
		  errprintf("pointdistance faulty\n");
		  return NULL;
		}
		dprintf(" ...done\n");  
		if( settings.regularize ){
		  dprintf(" Regularize d\n");
		  dtw_regularize_matrix( d, G, n, n );
		  dprintf(" ...done\n");  
		}
		dprintf(" Compute path\n");  
		dtw_cumulate_matrix( d, n, n );
		P = dtw_backtrack( (const double**) d, n, n, P );
		/* P = DTW_path_from_square_distmatrix( (const double**)d, n, P ); */
		dprintf(" ...done\n");  
		dprintf(" Warpavg\n");  
		new = eeg_warp_add_signals_by_path( s1, s2, new, chan, P, weights );
		dprintf(" ...done\n");  
	 }

	 /* remove the previous trials */
	 free_eegdata(s1);
	 free_eegdata(s2);
	 eeg->data[idx1] = new; /* ADTW goes to idx1 */
	 eeg->data[idx2] = NULL; /* do not touch this again */ 

	 indices[idx1] += indices[idx2];
	 indices[idx2]=-1; 			  /* never to be used again */
	 dprintf("new indices[%i,%i]=(%i,%i)\n", idx1, idx2, indices[idx1], indices[idx2]);

	 /* replace node by leaf representing ADTW(idx1, idx2) */
	 Tsub->val = idx1;
	 Tsub->left = NULL;
	 Tsub->right = NULL;

	 trials_left--;
  }
  for( i=0; i<N; i++ ){
	 dprintf("indices[%i]=%i\n", i, indices[i]);
  }

  copy_similar_eegdata( out, eeg->data[idx1] );
  
  /* scale by number of trials */
  if( settings.progress ){
	 settings.progress( PROGRESSBAR_FINISH, 0 );
  }

  /* cleaning up */
  dprintf("Freeing Memory\n");
  matrix_free( d, n );
  matrix_free( G, n );
  dgram_free( T );
  free( indices );
  free_warppath( P );
  if( !settings.channels ){
	 free( channels );
  }


  return out;
}


EEGdata* eeg_warp_add_signals_by_path(const EEGdata *s1, const EEGdata *s2, 
												 EEGdata *target, int channel, 
												 const WarpPath *P, const double weights[2]){
  int i;

  if(target==NULL){
	 target = init_eegdata(s1->nbchan, (s1->n+s2->n)/2, 0);
	 warnprintf("ALLOC: allocated memory in function!\n");
  }
  warp_add_signals_by_path( s1->d[channel], s1->n,
									s2->d[channel], s2->n,
									P, target->d[channel], weights );
  
  /* set time-markers */
  for( i=0; i<s1->nmarkers; i++ ){
	 target->markers[i] = (s1->markers[i]+s2->markers[i])/2;
	 dprintf( "New marker for target[%i]=%i\n", i, (int)target->markers[i] );
  }

  return target;
}


/** Warp-average according to Gibbons+Stahl 2007.
	 They proposed to stretch or
	 compress the single signals in order to match the average reaction
	 time, by simply moving the sampling points in time according to a
	 linear, quadratic, cubic or to-the-power-of-four function. 
	 Formally, they adjusted the time-axis by letting
	 \f[
	 \phi_t^{-1}(t) = t + \frac{t^k}{R^k_t}(E(R) - R_t)
	 \f]
	 where \f$ R_t \f$ denotes the reaction time of the current trial and E is
	 the expected value (the mean reaction time across trials). In their
	 work, Gibbons et al. studied this approach for \f$ k \in \{1,2,3,4\} \f$.

	 Individual trials are warped according to \f$ \phi \f$ and also the
	 averages obtained from different individuals.
	 Warping takes place between stimulus-onset-marker and response-marker
	 \param eeg_in input
	 \param target or ALLOC_IN_FCT
	 \param stmulus_marker gives the index indicating which of the markers within eeg_in  
	                       is the stimulus-onset
	 \param response_marker gives the index indicating which of the markers within eeg_in  
	                       is the response-onset						  
	 \param k parameter for gibbon's method
	 \return pointer to target or newly allocated memory
*/  
EEGdata* eegtrials_gibbons( const EEGdata_trials *eeg, EEGdata *target, 
									 int stimulus_marker, int response_marker, double k ){
  double *times, *new_times;
  double meanRT=0, curRT;
  double step;
  int n, N, i, j, l;
  int numchan, c;
  int stim_onset, resp_onset;

  if( stimulus_marker>=eeg->nmarkers_per_trial ||
		response_marker>=eeg->nmarkers_per_trial ||
		stimulus_marker>=response_marker ){
	 errprintf( "stimulus or response marker not correct (%i, %i), abort fct\n", 
					 stimulus_marker, response_marker );
	 return NULL;
  }

  /* for convenience */
  times = eeg->times;  
  N = eeg->ntrials;
  n = eeg->nsamples;
  numchan = eeg->data[0]->nbchan;


  if( target==ALLOC_IN_FCT ){
	 warnprintf(" allocating memory in fct\n");
	 target = init_eegdata( numchan, n, eeg->nmarkers_per_trial );
  } else {
	 if( eegdata_cmp_settings( target, eeg->data[0] ) ){
		errprintf( " target and input not similar enough\n" );
		return NULL;
	 }
	 reset_eegdata( target );	  /* set to zero */
  }

  /* gsl interpolation allocation */
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n); /* gibbons uses linear */

  /* computing mean reaction time */
  for( i=0; i<N; i++ ){
	 meanRT += times[eeg->markers[i][response_marker]] 
		- times[eeg->markers[i][stimulus_marker]];
  }
  meanRT /= (double)N; 

  /* compute warping for each trial */
  new_times = (double*) malloc( n*sizeof(double) );
  memcpy( new_times, times, n*sizeof(double) ); /* copy times */
  for( i=0; i<N; i++ ){
	 stim_onset = eeg->markers[i][stimulus_marker];
	 resp_onset = eeg->markers[i][response_marker];
	 curRT = times[resp_onset]-times[stim_onset];

	 for( j=stim_onset; j<resp_onset; j++ ){
		new_times[j] = times[j] + pow( times[j], k )
		  /pow( curRT, k )*( meanRT-curRT ); /* gibbons eq. */
	 }
	 /* warp the rest of the segments linearly */
	 step =(times[n-1]-curRT)/(double)(n-1-resp_onset);
	 l = 0;
	 for( j=resp_onset; j<n; j++ ){
		new_times[j]=meanRT + l*step;
	 }

	 for( c=0; c<numchan; c++ ){ /* loop channels */
		gsl_spline_init (spline, new_times, eeg->data[i]->d[c], n);
		for( j=0; j<n; j++ ){	  /* and samples */
		  target->d[c][j] += gsl_spline_eval ( spline, times[j], acc );
		}
	 }
  }
  for( c=0; c<numchan; c++ ){
	 for( j=0; j<n; j++ ){
		target->d[c][j] /= (double)N; /* average */
	 }
  }
  
  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free( new_times );

  return target;
}

/** \cond PRIVATE */
void _settings_guess_from_eegtrials( const EEGdata_trials *eeg, SettingsHierarchicalDTW *settings ){
  int n = eeg->nsamples;
  settings->winlength = MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 );
  settings->N_time=n;
  settings->N_freq=(settings->winlength)*4;
  settings->sampling_rate=eeg->sampling_rate;
  settings->corner_freqs[0]=0.0;
  settings->corner_freqs[1]=settings->sampling_rate/2.0;
}
/** \endcond */

SettingsHierarchicalDTW init_dtw_hierarchical( const EEGdata_trials *eeg ){	 
  SettingsHierarchicalDTW settings;
  settings.regularize=NULL; /* eeg_regularization_gaussian_line; */
  settings.linkage=dgram_dist_completelinkage;
  settings.dont_touch_eeg=1;
  settings.pointdistance=eeg_distmatrix_stft_channel;
  settings.theta1=1.0;
  settings.theta2=1.0;
  settings.corner_freqs[0]=-1;
  settings.corner_freqs[1]=-1;
  settings.winfct=window_hanning;
  settings.winlength=-1;
  settings.sampling_rate=-1;
  settings.N_freq=-1;
  settings.N_time=-1;
  settings.sigma=0.2;
  settings.channels=NULL;
  settings.num_channels = -1;
  settings.m = 10;
  settings.tau = 15;
  settings.FAN = (int)0.05*eeg->nsamples;
  settings.progress = NULL;
  _settings_guess_from_eegtrials( eeg, &settings );
  return settings;
}

void  print_settings_hierarchicaldtw( FILE *out, SettingsHierarchicalDTW s ){
  int i;

  fprintf( out, "SettingsHierarchicalDTW '%p':\n"
			  "  regularize     = %p\n"
			  "  sigma          = %f\n"
			  "  linkage        = %p\n"
			  "  dont_touch_eeg = %i\n"
			  "  pointdistance  = %p\n"
			  "  theta          = (%f,%f)\n"
			  "  corner_freqs   = (%f,%f)\n"
			  "  sampling_rate  = %f\n"
			  "  winfct         = %p\n"
			  "  winlength      = %i\n"
			  "  N_freq         = %i\n"
			  "  N_time         = %i\n"
			  "  m              = %i\n"
			  "  tau            = %i\n"
			  "  FAN            = %i\n"
			  "  channels       = {", &s, s.regularize, s.sigma, s.linkage,
			  s.dont_touch_eeg, s.pointdistance, s.theta1, s.theta2,
			  s.corner_freqs[0], s.corner_freqs[1], s.sampling_rate, 
			  s.winfct, s.winlength, s.N_freq, s.N_time, s.m, s.tau, s.FAN );
  for( i=0; i<s.num_channels; i++ ){
	 fprintf( out, "%i, ", s.channels[i] );
  }
  fprintf( out, "}\n"
			  "  num_channels  = %i\n"
			  "  progress      = %p\n", s.num_channels, s.progress);
}



WarpPath* init_warppath(WarpPath *path, int n1, int n2){
  if(path==NULL){
	 path = (WarpPath*)malloc(sizeof(WarpPath));
	 path->t1=NULL;
	 path->t2=NULL;
  }
  path->n = (n1+n2);
  path->n1= n1;
  path->n2= n2;

  if( path->t1 ){
	 free( path->t1 );
  } 
  if( path->t2 ){
	 free( path->t2 );
  }
  path->t1 = (int*)calloc(n1+n2, sizeof(int));
  path->t2 = (int*)calloc(n1+n2, sizeof(int));

  return path;
}

void      reset_warppath(WarpPath *P, int n1, int n2){
  memset( P->t1, 0, P->n*sizeof(int) );
  memset( P->t2, 0, P->n*sizeof(int) );
  P->n1 = n1; P->n2 = n2;
  P->n = 0;
  memset( P->t1, 0, n1*sizeof(int) );
  memset( P->t2, 0, n2*sizeof(int) );
}

void free_warppath(WarpPath *p){
	free(p->t1);
	free(p->t2);
	free(p);
}

void print_warppath( FILE *out, WarpPath *P ){
  int howmany, i;
  howmany=3;

  fprintf( out, "WarpPath '%p':\n"
			  " n1 = %i\n"
			  " n2 = %i\n"
			  " n  = %i\n"
			  " P  = ", P, P->n1, P->n2, P->n );
  if( P->n>=howmany ){
	 for( i=0; i<howmany; i++ ){
		fprintf( out, " (%i -> %i),", P->t1[i], P->t2[i] );
	 }
  } else if( P->n>0 ){
	 for( i=0; i<P->n; i++ ){
		fprintf( out, " (%i -> %i),", P->t1[i], P->t2[i] );
	 }
  } else {
	 fprintf( out, "<empty>");
  }
  fprintf( out, "\n" );
}
  /* ------------------------------------------------------------------------------
	  ------------------------------------------------------------------------------
	  -----------------------OBSOLETE-----------------------------------------------
	  ------------------------------------------------------------------------------
	  ------------------------------------------------------------------------------ */

/** \cond OBSOLETE */
/** Construct the warppath from a square distance matrix.

	 Algorithm:
	 - cumulate the pointwise distance-matrix d by choosing
	 \f[
	 \mathbf{D}_{jk} = \mathbf{d}_{jk}+\min{\{\mathbf{D}_{j,k-1}, \mathbf{D}_{j-1,k}, \mathbf{D}_{j-1, k-1}\}}
	 \f]
	 \param d -  distance matrix 
	 \param n - nxn matrix d
	 \param path - pointer to WarpPath-struct of length (J+K) or NULL (allocated in function)
	 \return WarpPath struct
*/
WarpPath* DTW_path_from_square_distmatrix(const double **d, int n, WarpPath *path){
  int j,k;
  int idx;
  double left=0, down=0, downleft=0;
  double **D;

  if( path==NULL ){
	 path = init_warppath(NULL, n, n);
	 warnprintf( "allocating warppath in function\n");
  }
  dprintf("path=%p,path->n1=%i, path->n2=%i\n",path, path->n1, path->n2);
  reset_warppath(path, n, n);
  
  D = matrix_init( n, n );
  matrix_copy( d, D, n, n );

  /* computing D_jk */
  for(j=0; j<n; j++){
    for(k=0; k<n; k++){
      if(k==0 && j==0) ;
      else if(k==0 && j>0)
		  D[j][k] += D[j-1][k];
      else if(k>0 && j==0)
		  D[j][k] += D[j][k-1];
      else /* j,k > 0 */
		  D[j][k] += MIN(MIN(D[j][k-1], D[j-1][k]), D[j-1][k-1]);
    }
  }
  dprintf("Djk = %f\n", D[n-1][n-1] );
  /* Backtracking */
  j=n-1; k=n-1;

  idx = (j+k);
  path->t1[idx] = j;
  path->t2[idx] = k;
  while( j>0 || k>0 ){
	 if( k==0 ){ /* base cases */
		j--;
	 } else if( j==0 ){
		k--;
	 } else { /* min( d[j-1][k], d[j-1][k-1], d[j][k-1] ) */
	  
		left     = D[j-1][k  ];
		down     = D[j  ][k-1];
		downleft = D[j-1][k-1];
		
		(isnan( D[j-1][k  ] ) )?(left     = DBL_MAX):(left     = D[j-1][k  ]);
		(isnan( D[j  ][k-1] ) )?(down     = DBL_MAX):(down     = D[j  ][k-1]);
		(isnan( D[j-1][k-1] ) )?(downleft = DBL_MAX):(downleft = D[j-1][k-1]);

		if( left<=downleft ){
		  if( left <= down )
			 j--;
		  else
			 k--;
		} else {
		  if( downleft <= down ){
			 k--; j--;
		  } else 
			 k--;
		}
	 }
	 
	 /* if( left>=downleft && downleft<=down ){ /\* diagonal step *\/ */
	 /* 	path->t1[idx] = j+1; */
	 /* 	path->t2[idx] = k; */
	 /* 	idx--; */
	 /* 	path->t1[idx] = j; */
	 /* 	path->t2[idx] = k; */
	 /* } else { */
		path->t1[idx] = j;
		path->t2[idx] = k;
	 /* } */
	 idx--;
	 /* dprintf("(j, k), idx = (%i,%i), %i\n", j, k, idx); */
  }
  dprintf("idx=%i, length(path)=%i\n", idx, 2*(n-1)-idx );

  /* free */
  matrix_free( D, n );

  dprintf("leaving\n");
  return path;
}



/** build the cumulated dissimilarity matrix D
	 \param u,J 1st signal
	 \param s,K 2nd signal
	 \param R - restriction band ( in [0,1] )
	 \param d -- if NULL, function allocates the memory
	 \return pointer to JxK matrix D
 */
double** DTW_build_restricted_cumdistmatrix(const double *u, int J, 
														  const double *s, int K, 
														  double R, double **d){ 
  int j,k;
  double avgu, avgs, rmsu=0.0, rmss=0.0;
  double snorm, snormp; /* variable for a point of the normalized
									signal and the preceding point */
  double unorm, unormp;
  int theta;
  double left,down,downleft;

  int *bham;

  if( d==NULL ){
	 d = (double**) malloc(J*sizeof(double*));
	 for(j=0; j<J; j++)
		d[j]=(double*)malloc(K*sizeof(double));
  }

  bham = (int*)malloc( MAX( J,K )*2*sizeof( int ) );
  bresenham(0,0, J-1, K-1, bham);
    
 
  if( R>1 ){
	 dprintf("Restriction R=%f too large, using 1.0\n", R);
	 R = 1.0;
  } else if( R<0 ){
	 dprintf("Restriction R=%f < 0, aborting\n", R);
	 return d;
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
 
  /* for( i=0; i<MAX( J,K ); i++ ){ */
  /* 	 d[bham[(2*i)+0]][bham[(2*i)+1]] = NAN; */
  /* } */
  
  free(bham);
  return d;
}

/** build the cumulated dissimilarity matrix D
	 \param u,J 1st signal
	 \param s,K 2nd signal
	 \param theta1/theta2 - weights for metric
	 \param d -- if NULL, function allocates the memory
	 \return pointer to JxK matrix D
 */
double** DTW_build_cumdistmatrix(const double *u, int J, const double *s, int K, 
											double theta1, double theta2, double **d){  
  int j,k;
  double avgu, avgs, rmsu=0.0, rmss=0.0;
  double snorm, snormp; /* variable for a point of the normalized
									signal and the preceding point */
  double unorm, unormp;

  if( d==NULL ){
	 d = (double**) malloc(J*sizeof(double*));
	 for(j=0; j<J; j++)
		d[j]=(double*)malloc(K*sizeof(double));
  }
    
  avgu = gsl_stats_mean(u, 1, J);
  avgs = gsl_stats_mean(s, 1, K);
  for(j=0; j<J; j++) rmsu += pow(u[j], 2);
  rmsu = sqrt(rmsu/(double)J);
  for(k=0; k<K; k++) rmss += pow(s[k], 2);
  rmss = sqrt(rmss/(double)K);
  
  /* computing D_jk */
  for(j=0; j<J; j++){
    for(k=0; k<K; k++){
      snorm = (s[k]-avgs)/rmss;
      unorm = (u[j]-avgu)/rmsu;
      (k==0) ? (snormp = 0) : (snormp = ((s[k-1]-avgs)/rmss));
      (j==0) ? (unormp = 0) : (unormp = ((u[j-1]-avgu)/rmsu));
      d[j][k] = theta1 * fabs(unorm - snorm) + theta2 * fabs( (unorm-unormp) - (snorm-snormp) ); /* (1) */
      if(k==0 && j==0) ;
      else if(k==0 && j>0)
		  d[j][k] += d[j-1][k];
      else if(k>0 && j==0)
		  d[j][k] += d[j][k-1];
      else /* j,k > 0 */
		  d[j][k] += MIN(MIN(d[j][k-1], d[j-1][k]), d[j-1][k-1]);
    }
  }

  return d;
}  





/** construct the warppath from distance matrix. Use DTW_build_cumdistmatrix() for that.
	 \param d - cumulated distmatrix
	 \param J,K - dims for d
	 \param path - pointer to WarpPath-struct of length (J+K) or NULL (allocated in function)
	 \return WarpPath struct
*/
WarpPath* DTW_path_from_cumdistmatrix(const double **d, int J, int K, WarpPath *path){
  int j, k;
  int idx;
  double left, down, downleft;

  if( path==NULL ){
	 path = init_warppath(NULL, J, K);
  }
  dprintf("path=%p,path->n1=%i, path->n2=%i\n",path, path->n1, path->n2);
  reset_warppath(path, J, K);


  /* Backtracking */
  j=J-1; k=K-1;

  idx = 1;
  path->t1[0] = j;
  path->t2[0] = k;
  while( j>0 || k>0 ){
	 if( k==0 ){ /* base cases */
		j--;
	 } else if( j==0 ){
		k--;
	 } else { /* min( d[j-1][k], d[j-1][k-1], d[j][k-1] ) */
	  
		left     = d[j-1][k  ];
		down     = d[j  ][k-1];
		downleft = d[j-1][k-1];

		(isnan( d[j-1][k  ] ) )?(left     = DBL_MAX):(left     = d[j-1][k  ]);
		(isnan( d[j  ][k-1] ) )?(down     = DBL_MAX):(down     = d[j  ][k-1]);
		(isnan( d[j-1][k-1] ) )?(downleft = DBL_MAX):(downleft = d[j-1][k-1]);

		if( left<=downleft ){
		  if( left <= down )
			 j--;
		  else
			 k--;
		} else {
		  if( downleft <= down ){
			 k--; j--;
		  } else 
			 k--;
		}
	 }

	 
	 path->t1[idx] = j;
	 path->t2[idx] = k;
	 idx++;
	 /*dprintf("(j, k), idx = (%i,%i), %i\n", j, k, idx);*/
  }
  dprintf("leaving\n");
  return path;
}

/** \endcond */
