/***************************************************************************
 *   Copyright (C) 2008/2009 by Matthias Ihrke   *
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

#include "distances.h"


/** build a distance matrix D from data X which consists of n observations with
	 p features each. 
	 Observations are compared using f.
	 \param f distance function
	 \param X data (nxp)
	 \param D output matrix or NULL -> own memory allocation
	 \param progress NULL or a progressbar
	 \return the distance matrix between all columns in X
 */
double** vectordist_distmatrix( VectorDistanceFunction f, 
										  const double **X, int n, int p, double **D, 
										  ProgressBarFunction progress,
										  void *userdata ){
  int i,j;
  int idx;
  
  if( D==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 D = matrix_init( n, n );
  }
  if( progress ){
	 progress( PROGRESSBAR_INIT, n*(n-1)/2 );
  }
  
  idx=1;
  for( i=0; i<n; i++ ){
	 for( j=i+1; j<n; j++ ){
		if( progress ){
		  progress( PROGRESSBAR_CONTINUE_LONG, idx );
		}
		D[i][j] = f( (double*)X[i], (double*)X[j], p, userdata );
		D[j][i] = D[i][j];
		if( isnan( D[j][i] ) ){
		  errprintf("D[%i][%i] is nan\n", j, i);
		}
		idx++;
	 }
  }
  if( progress ){
	 progress( PROGRESSBAR_FINISH, 0 );
  }

  return D;
}

/** Euclidean distance between vector x1 and x2.
	 \param x1, x2 vectors of size p
	 \param userdata is ignored
*/
double vectordist_euclidean( double *x1, double *x2, int p, void *userdata ){
  double d;
  int i;

  d=0.0;
  for( i=0; i<p; i++ ){
	 d += SQR( x1[i]-x2[i] );
  }
  d = sqrt( d );

  return d;
}
/** Euclidean distance between normalized vector x1 and x2.
	 Normalization is defined by
	 \f[
	 \hat{x}_i = \frac{x_i-\langle x_i \rangle_i}{\sqrt{\langle x_i^2\rangle_i}}
	 \f]

	 \param x1, x2 vectors of size p
	 \param userdata is ignored
 */
double   vectordist_euclidean_normalized( double *x1, double *x2, int p, void *userdata ){
  int i;
  double d;
  double avg1, avg2,
	 rms1=0, rms2=0;
  double norm1,norm2; /* normalized signal at time t */
	 
  avg1 = gsl_stats_mean(x1, 1, p);
  avg2 = gsl_stats_mean(x2, 1, p);
  
  for(i=0; i<p; i++){
	 rms1 += SQR( x1[i] );
	 rms2 += SQR( x2[i] );
  }
  
  rms1 = sqrt(rms1/(double)p);
  rms2 = sqrt(rms2/(double)p);

  d=0.0;
  for( i=0; i<p; i++ ){
	 norm1 = (x1[i]-avg1)/rms1;
	 norm2 = (x2[i]-avg2)/rms2;

	 d += SQR( norm1-norm2 );
  }
  d = sqrt( d );

  return d;
}
/** Distance which sums the absolute deviation from the main 
	 diagonal of the dtw-warp-path computed on x1 and x2.
	 \param x1, x2 vectors of size p
	 \param userdata is (PointwiseDistanceFunction, userdata for distfunction)
 */
double   vectordist_dtw( double *x1, double *x2, int p, void *userdata ){
  double d;
  double **D;
  PointwiseDistanceFunction f;
  WarpPath *P;
  double diag1[2], diag2[2], path[2];
  int i;

  f = *((PointwiseDistanceFunction)userdata);
  userdata += sizeof(PointwiseDistanceFunction);

  D = f( x1, p, x2, p, NULL, userdata );
  dtw_cumulate_matrix( D, p, p );
  P = dtw_backtrack( D, p, p, NULL );

  d = 0;
  diag1[0] = 0;   
  diag1[1] = 0;
  diag2[0] = p;
  diag2[0] = p; /* line through (0,0), (p,p) */
  for( i=0; i<P->n; i++ ){
	 path[0] = (double)P->t1[i];
	 path[1] = (double)P->t2[i];
	 d += dist_point_line( path, diag1, diag2 );
  }

  matrix_free( D, p );
  free_warppath( P );

  return d;
}

/** \todo
	 compute the cumulated sum along the regularized warping path.
	 \param userdata (int)userdata[0] is nmarkers; 
	       ((int)userdata)[1] and following is an 2 x nmarkers matrix 
			 given the corresponding markers in the signals.
 */
double   vectordist_regularized_dtw( double *x1, double *x2, int p, void *userdata ){
  /* int nmarkers; */
  /* int **markers; */
  /* double **G, **D; */
  /* double maxsigma; */
  /* double Djk; */

  /* nmarkers = *((int*)userdata); */
  /* userdata = ((int*)userdata)+1; */
  /* maxsigma = *((double*)userdata); */
  /* userdata = ((double*)userdata)+1; */
  /* markers = (int**)userdata; */
  
  /* G = regularization_gaussian_line( markers[0], markers[1], nmarkers, p, maxsigma, NULL ); */
  /** TODO **/
  errprintf("Not implemented yet, HUA, HUA!\n");
  return -1;
}


/** pointwise distance matrix for euclidean metric.
	 \param userdata is ignored
 */
double** signaldist_euclidean( const double *s1, int n1, const double *s2, int n2,
										 double **d, void *userdata ){
  int i,j;
  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( n1, n2 );
  }

  for( i=0; i<n1; i++){
	 for( j=0; j<n2; j++ ){
		d[i][j] = SQR( fabs( s1[i] - s2[j] ) );
	 }
  }
  return d;
}


/** calculate
	 \f[
	 d(s_1(t_1), s_2(t_2))_{{derivative}} := \theta_1|\tilde{s_1}(t_1) - \tilde{s_2}(t_2)| + \theta_2|\tilde{s_1}'(t_1) - \tilde{s_2}'(t_2)|
	 \f]
	 where
	 \f[ 
	 \tilde{s}(t) := \frac{s(t) - \langle s(t)
	 \rangle_t}{\sqrt{\langle s(t)^2 \rangle_t}}
	 \f]
	 Derivatives are approximated with s'(t) = s(t)-s(t-1)
	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param userdata is (theta1,theta2) both as double values
	 \return d or NULL (error)
 */
double** signaldist_euclidean_derivative( const double *s1, int n1, const double *s2, int n2, double **d, void *userdata ){
  int i, j, k;
  double avg1, avg2,
	 rms1=0, rms2=0;
  double norm1,norm2, /* normalized signal at time t */
	 norm1p,norm2p;    /* normalized signal at time t-1 */
  double theta1,theta2;

  if( !userdata ){
	 theta1=1.0;
	 theta2=1.0;
  } else {
	 theta1 = *((double*)userdata);
	 userdata = ((double*)userdata)+1;
	 theta2 = *((double*)userdata);
  }

  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( n1, n2 );
  }

  avg1 = gsl_stats_mean(s1, 1, n1);
  avg2 = gsl_stats_mean(s2, 1, n2);
  for(i=0; i<MAX(n1, n2); i++){
	 if( i<n1 )
		rms1 += SQR( s1[i] );
	 if( i<n2 )
		rms2 += SQR( s2[i] );
  }

  rms1 = sqrt(rms1/(double)n1);
  rms2 = sqrt(rms2/(double)n2);

  for( j=0; j<n1; j++ ){
	 for( k=0; k<n2; k++ ){
      norm1 = (s1[j]-avg1)/rms1;
      norm2 = (s2[k]-avg2)/rms2;
      (j==0) ? (norm1p = 0) : (norm1p = ((s1[j-1]-avg1)/rms1));
      (k==0) ? (norm2p = 0) : (norm2p = ((s2[k-1]-avg2)/rms2));
		d[j][k] = theta1*fabs(norm1 - norm2) + theta2*fabs( (norm1-norm1p) - (norm2-norm2p) ); /* (1) */
	 }
  }
  return d;
}


/** 
	 \todo at the moment, N_time is required to be n, the number of
	       sampling points in s1 and s2. This needs to be fixed by interpolation

	 calculate
	 \f[  
	 d_{{STFT}}(s_1(t_1), s_2(t_2)) := || STFT\{s_1(t_1)\} -
	 \STFT\{s_2(t_2)\}||_{\circ}.
	 \f]
	 where 
	 \f[
	 STFT\{s(t)\}(\omega) = \int_{-\infty}^{+\infty} s(t)w(t-\tau)e^{-i\omega t} dt
	 \f]
	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param params is a
	 \return d or NULL (error)
 */
double** signaldist_stft( const double *s1, int n1, const double *s2, int n2, double **d, void *userdata ){
  int i, j;
  Spectrogram *sp1, *sp2;
  double *window;
  double *d1, *d2,
	 mean1, mean2;
  
  /* extract from userdata */
  double sample_frequency;
  WindowFunction winfct; 
  int winlength;
  int N_freq;
  int N_time;
  double corner_freqs[2];

  /* get params from void pointer */
  sample_frequency = *((double*)(userdata));
  userdata        += sizeof(double);
  winfct           = *((WindowFunction*)(userdata));
  userdata        += sizeof(WindowFunction);
  winlength        = *((int*)(userdata));
  userdata        += sizeof(int);
  N_freq           = *((int*)(userdata));
  userdata        += sizeof(int);
  N_time           = *((int*)(userdata));
  userdata        += sizeof(int);
  corner_freqs[0]  = *((double*)(userdata));
  userdata        += sizeof(double);
  corner_freqs[1]  = *((double*)(userdata));
  userdata        += sizeof(double);

  dprintf("got sfreq=%f, winfct=%p, winlength=%i, N_f=%i, N_t=%i, cf=(%f,%f)\n",
			 sample_frequency, winfct, winlength, N_freq, N_time, corner_freqs[0], corner_freqs[1] );

  if( d==ALLOC_IN_FCT ){
	 warnprintf("allocating matrix in fct\n");
	 d=matrix_init( n1, n2 );
  }

  /* remove mean from the signals */
  dprintf("Normalization\n");
  d1 = vector_init( NULL, n1, 0.0 );
  d2 = vector_init( NULL, n2, 0.0 );
  mean1 = gsl_stats_mean( s1, 1, n1 );
  mean2 = gsl_stats_mean( s2, 1, n2 );
  for( i=0; i<MAX(n1,n2); i++ ){	  
	 if( n1<i )
		d1[i] = s1[i]-mean1;
	 if( n2<i )
		d2[i] = s2[i]-mean2;
  }

  window = winfct(ALLOC_IN_FCT, winlength);
  sp1 = spectrogram_stft( d1, n1, sample_frequency,
								  window, winlength, N_freq, N_time, 
								  corner_freqs, ALLOC_IN_FCT );
  sp2 = spectrogram_stft( d2, n2, sample_frequency,
								  window, winlength, N_freq, N_time,
								  corner_freqs, ALLOC_IN_FCT );
  if( !sp1 || !sp2 ){
	 errprintf("Spectrogram broken\n");
	 return NULL;
  }
  spectrogram_compute_powerspectrum( sp1 );
  spectrogram_compute_powerspectrum( sp2 );

  dprintf(" Compute distances (%i vectors)\n", sp1->N_freq);
  for( i=0; i<N_time; i++ ){
	 for( j=0; j<N_time; j++ ){
		d[i][j] = vector_euclidean_distance( sp1->powerspect[i], sp2->powerspect[j], sp1->N_freq );
	 }
  }
  dprintf("\nDone\n");

  /* free */
  free_spectrogram( sp1 );
  free_spectrogram( sp2 );
  free( d1 );
  free( d2 );
  free( window );
  
  return d;
}

/** calculate
	 \f$ 
	 d_{euclid}(s_1(t_1),s_2(t_2)) = |s_1(t_1)-s_2(t_2)|^2
	 \f$
	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param params is ignored
	 \return d or NULL (error)
 */
double** eeg_distmatrix_euclidean_channel( const EEGdata *s1, const EEGdata *s2, 
														 int channel, double **d, void *userdata ){

  d = signaldist_euclidean( s1->d[channel], s1->n, s2->d[channel], s2->n, d, userdata );

  return d;
}

/** calculate
	 \f[
	 d(s_1(t_1), s_2(t_2))_{{derivative}} := \theta_1|\tilde{s_1}(t_1) - \tilde{s_2}(t_2)| + \theta_2|\tilde{s_1}'(t_1) - \tilde{s_2}'(t_2)|
	 \f]
	 where
	 \f[ 
	 \tilde{s}(t) := \frac{s(t) - \langle s(t)
	 \rangle_t}{\sqrt{\langle s(t)^2 \rangle_t}}
	 \f]
	 Derivatives are approximated with s'(t) = s(t)-s(t-1)
	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param params is a SettingsHierarchicalDTW struct that should contain valid values for theta1,theta2 weights for distance
	 \return d or NULL (error)
 */
double** eeg_distmatrix_euclidean_derivative_channel( const EEGdata *s1, const EEGdata *s2, 
																		int channel, double **d, 
																		void *params ){ 
  /* get params from void pointer */
  SettingsHierarchicalDTW *settings;
  double theta[2];
  settings = (SettingsHierarchicalDTW*)params;
  if( !settings ){
	 theta[0]=1.0;
	 theta[1]=1.0;
  } else {
	 theta[0]=settings->theta1;
	 theta[1]=settings->theta2;
  }

  if( eegdata_cmp_settings( s1, s2 ) ){
	 errprintf(" ERROR: two EEGdata sets are not similar enough\n");
	 return NULL;
  }

  d = signaldist_euclidean_derivative( s1->d[channel], s1->n, s2->d[channel], s2->n, d, &theta );
													
  return d;
}

/** 
	 \todo at the moment, N_time is required to be n, the number of
	       sampling points in s1 and s2. This needs to be fixed by interpolation

	 calculate
	 \f[  
	 d_{{STFT}}(s_1(t_1), s_2(t_2)) := || STFT\{s_1(t_1)\} -
	 \STFT\{s_2(t_2)\}||_{\circ}.
	 \f]
	 where 
	 \f[
	 STFT\{s(t)\}(\omega) = \int_{-\infty}^{+\infty} s(t)w(t-\tau)e^{-i\omega t} dt
	 \f]
	 This function calls signaldist_stft() with the appropriate variables packed 
	 together.

	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param params is a SettingsHierarchicalDTW struct that should contain valid:
	   - (sample_frequency of the signal),
	   - (winfct which windowing function? one of window_*()),
	   - (winlength lenght of window in samples),
	   - (N_freq how many frequencies),
	   - (N_time how many time-points),
	   - (corner_freqs - corner frequencies (lower, upper) for the returned
	     spectrum at each time-sample in Hz; maximal would be
	     {0, srate/2}) 
	 fields
	 \return d or NULL (error)
 */
double** eeg_distmatrix_stft_channel( const EEGdata *s1, const EEGdata *s2, 
												  int channel, double **d, void *params ){
  void *pack, *orgpack;
  SettingsHierarchicalDTW *settings;
  settings = (SettingsHierarchicalDTW*)params;

  /* packing stuff for signaldist-fct */
  pack = (void*)malloc( 4*sizeof(double) + 3*sizeof(int) + sizeof(WindowFunction) );
  orgpack = pack;
  memcpy( pack, &(settings->sampling_rate   ), sizeof(double) );
  pack += sizeof(double);
  memcpy( pack, &(settings->winfct          ), sizeof(WindowFunction) );
  pack += sizeof(WindowFunction);
  memcpy( pack, &(settings->winlength       ), sizeof(int) );
  pack += sizeof(int);
  memcpy( pack, &(settings->N_freq          ), sizeof(int) );
  pack += sizeof(int);
  memcpy( pack, &(settings->N_time          ), sizeof(int) );
  pack += sizeof(int);
  memcpy( pack, &(settings->corner_freqs    ), 2*sizeof(double) );

  d = signaldist_stft( s1->d[channel], s1->n, s2->d[channel], s2->n, d, orgpack );
  free( orgpack );
  return d;
}

double** eeg_distmatrix_recplot_losdtwnoise_channel( const EEGdata *s1,const  EEGdata *s2, 
																	  int channel, double **d, void *params ){
  int i,j; 
  double noise;
  double noiseamp = 0.01;
  RecurrencePlot *R;
  PhaseSpace *p1, *p2;
  SettingsHierarchicalDTW *settings;
  settings = (SettingsHierarchicalDTW*)params;

  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( s1->n, s2->n );
  }
  srand((long)time(NULL));
  p1 = phspace_init( settings->m, settings->tau, s1->d[channel], s1->n );
  p2 = phspace_init( settings->m, settings->tau, s2->d[channel], s2->n );

  R = recplot_init( s1->n, s2->n, settings->FAN, RPLOT_FAN );
  recplot_calculate( R, p1, p2 );

  matrix_copy( (const double**)R->R, d, R->m, R->n ); 
  
  /* flip binary matrix */
  scalar_minus_matrix( 2.0, d, R->m, R->n );
  
  /* add small amount of noise */
  for( i=0; i<R->m; i++ ){
	 for( j=0; j<R->n; j++ ){
		noise = noiseamp*(((double)rand()) / RAND_MAX);		
		if( d[i][j]<2 )
		  d[i][j] += noise;
	 }
  }

  return d;
}


/** Compute trial-to-trial distance matrix for all trials in eeg-struct.
	 Average over all channels in the EEG-struct.

	 \todo implement that the called distance function incorporates time-markers or whatever is needed
	 via the userdata
	 \param eeg the EEG-data
	 \param f the function used to compare two ERPs
	 \param d user allocated memory, or NULL -> own memory is alloc'ed
	 \param userdata userdata is passed to the vectordist_*() function
 */
double** eeg_distmatrix( EEG *eeg, VectorDistanceFunction f, double **d, void *userdata ){
  int i,j;
  int c;

  if( d==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 d = matrix_init( eeg->ntrials, eeg->ntrials );
  }
  


  for( i=0; i<eeg->ntrials; i++ ){
	 d[i][i] = 0.0;
	 for( j=i+1; j<eeg->ntrials; j++ ){
		for( c=0; c<eeg->nbchan; c++ ){
		  d[i][j] += f( eeg->data[c][i], 
							 eeg->data[c][i],
							 eeg->n, userdata );
		  if( isnan( d[i][j] ) ){
			 errprintf("D[%i][%i] is nan\n", i, j);
		  }
		}
		d[i][j] /= (double)eeg->nbchan;
		d[j][i] = d[i][j];
	 }
  }
  
  return d;
}

double pointdist_euclidean(double x, double y){
  return fabs(x-y);
}

/** Compute distance from point p to the line through points x and y. 
	 \code
	 d = abs(det([y-x;p-x]))/norm(y-x);
	 \endcode

	 \param p point in question (x,y)-plane
	 \param x,y points specify the line x=(x_1,x_2), y=(y_1, y_2)
 */
double   dist_point_line(double *p, double *x, double *y){
  double det;
  double d;
  det = ((y[0]-x[0])*(p[1]-x[1])) - ((y[1]-x[1])*(p[0]-x[0]));
  //dprintf("det=%f\n", det);
  d = (double)ABS( det ) / (double)sqrt( SQR(y[0]-x[0]) + SQR(y[1]-x[1]) );

  return d;
}

/** computes the distance between two warping-functions. 
	 It is assumed that both pathes have the same starting and endpoint.

	 Algorithm:
	 \code
	 1. compute distance transform for p1
	 2. line-integral of p2 through DT(p1) 
	 \endcode
 */
double   pathdist_euclidean_dt(WarpPath *p1, WarpPath *p2){
  int **I;
  int n,m,
	 i,j;
  double dist;
  double **d;
  
  n = p1->t1[p1->n-1]+1;
  m = p1->t2[p1->n-1]+1;

  dprintf("calculate distance with (%i x %i) points\n", n, m );

  I = (int**)malloc( n*sizeof(int*) );
  for( i=0; i<n; i++){
	 I[i] = (int*)malloc( m*sizeof(int) );
	 for( j=0; j<m; j++ ){
		I[i][j] = 0;
	 }
  }

  for( i=0; i<p1->n; i++ ){
	 //dprintf(" (%i , %i)\n", p1->t1[i], p1->t2[i] );
	 I[p1->t1[i]][p1->t2[i]]=1;
  }
  
  dprintf(" init done \n");
	 
  d = disttransform_deadreckoning( I, n, m, ALLOC_IN_FCT );
  matrix_divide_scalar( d, n, m, matrix_max( d, n, m, NULL, NULL ) );
  dprintf(" dt done \n");

  dist = 0.0;
  for( i=0; i<p2->n; i++ ){
	 //	 dprintf(" (%i , %i)\n", p2->t1[i], p2->t2[i] );
	 if(  p2->t1[i]>=n ||  p2->t2[i]>=n ){
		warnprintf("points in P2 do not fall in P1,  (%i , %i) -> ignored\n", p2->t1[i], p2->t2[i] );
	 } else {
		dist += d[p2->t1[i]][p2->t2[i]];
	 }
  }
  
  //  dist /= (double)n*m;
  dprintf(" dist=%f compute done \n", dist);
 
  for( i=0; i<n; i++ ){
	 free( I[i] );
  }
  free(I); 
  //matrix_free( d, m );

  return dist;
}
