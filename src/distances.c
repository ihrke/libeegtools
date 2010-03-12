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
#include "optarg.h"
#include <gsl/gsl_statistics.h>

/** build a distance matrix D from data X which consists of n observations with
	 p features each. 
	 Observations are compared using f.
	 \param f distance function
	 \param X data (nxp)
	 \param D output matrix or NULL -> own memory allocation
	 \param progress NULL or a progressbar
	 \param optargs may contain arguments for the vectordistancefunction f
	 \return the distance matrix between all columns in X
 */
double** vectordist_distmatrix( VectorDistanceFunction f, 
										  const double **X, int n, int p, double **D, 
										  ProgressBarFunction progress,
										  OptArgList *optargs ){
  int i,j;
  int idx;
  
  dprintf("data dim is (%i,%i)\n", n, p);
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
		D[i][j] = f( (double*)X[i], (double*)X[j], p, optargs );
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
	 \param optargs is ignored
*/
double vectordist_euclidean( const double *x1, const double *x2, int p, OptArgList *optargs ){
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
	 \param optargs is ignored
 */
double   vectordist_euclidean_normalized( const double *x1, const double *x2, int p, OptArgList *optargs ){
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
	 \param optargs may contain:
	 <ul>
	 <li> <tt>pointdistance=void*</tt> a PointwiseDistanceFunction, default is \c signaldist_euclidean
	 <li> optargs for the distfunction (optargs is passed 'as is' to this function)
	 </ul>
 */
double   vectordist_dtw( const double *x1, const double *x2, int p, OptArgList *optargs ){
  double d;
  double **D;
  PointwiseDistanceFunction f;
  void *tmp;
  WarpPath *P;
  double diag1[2], diag2[2], path[2];
  int i;

  f = signaldist_euclidean;
  if( optarglist_has_key( optargs, "distfct" ) ){
	 tmp = optarglist_ptr_by_key( optargs, "distfct" );
	 if( tmp )
		f = (PointwiseDistanceFunction) tmp;
  }

  D = f( x1, p, x2, p, NULL, optargs );
  dtw_cumulate_matrix( D, p, p, NULL );
  P = dtw_backtrack( (const double**)D, p, p, NULL );

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
	 \param optargs should contain the markers
 */
double   vectordist_regularized_dtw( const double *x1, const double *x2, int p, OptArgList *optargs ){
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
	 \param optargs is ignored
 */
double** signaldist_euclidean( const double *s1, int n1, const double *s2, int n2,
										 double **d, OptArgList *optargs ){
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
	 \param optargs can contain:
	 - <tt>theta1=double</tt>, default is \c 1.0
	 - <tt>theta2=double</tt>, default is \c 1.0
	 \return d or NULL (error)
 */
double** signaldist_euclidean_derivative( const double *s1, int n1, const double *s2, int n2, double **d, OptArgList *optargs ){
  int i, j, k;
  double avg1, avg2,
	 rms1=0, rms2=0;
  double norm1,norm2, /* normalized signal at time t */
	 norm1p,norm2p;    /* normalized signal at time t-1 */
  double theta1,theta2;
  double x;

  theta1=1.0;
  theta2=1.0;

  if( optarglist_has_key( optargs, "theta1" ) ){
	 x = optarglist_scalar_by_key( optargs, "theta1" );
	 if( !isnan( x ) )
		theta1=x;
  }
  if( optarglist_has_key( optargs, "theta2" ) ){
	 x = optarglist_scalar_by_key( optargs, "theta2" );
	 if( !isnan( x ) )
		theta2=x;
  }
  dprintf("Using theta=(%f,%f)\n", theta1, theta2 );


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
	 \param optargs may contain:
	 - <tt>sample_frequency=double</tt> of the signal, default is \c 500 (should really be provided!)
	 - <tt>winfct=void*</tt> windowing function, default is \c window_hanning
	 - <tt>winlength=int</tt> size of the window, default is <tt> MAX( SQR( sqrt(next_pow2( n ))-3 ), 5 )</tt>
	 - <tt>N_freq=int</tt> number of frequency bins, default is \c winlength*4
	 - <tt>N_time=int</tt> number of time bins, default is \c n
	 - <tt>corner_freqs=double*</tt> (array with two double entries), default is \c (0.0,250.0)
	 \return d or NULL (error)
 */
double** signaldist_stft( const double *s1, int n1, const double *s2, int n2, double **d, OptArgList *optargs ){
  int i, j;
  Spectrogram *sp1, *sp2;
  double *window;
  double *d1, *d2,
	 mean1, mean2;
  double x;
  void *ptr;

  double sample_frequency;
  WindowFunction winfct; 
  int winlength;
  int N_freq;
  int N_time;
  double corner_freqs[2];

  /* defaults */
  sample_frequency = 500.0;
  winfct = window_hanning;
  winlength =  MAX( SQR( sqrt(next_pow2( n1 ))-3 ), 5 );
  N_freq = winlength*4;
  N_time = n1;
  corner_freqs[0] = 0.0;
  corner_freqs[1] = 250.0;
  
  /* get params */
  if( optarglist_has_key( optargs, "sample_frequency" ) ){
	 x = optarglist_scalar_by_key( optargs, "sample_frequency" );
	 if( !isnan( x ) ) sample_frequency=x;
  }
  if( optarglist_has_key( optargs, "winfct" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "winfct" );
	 if( ptr ) winfct = (WindowFunction)ptr;
  }
  if( optarglist_has_key( optargs, "winlength" ) ){
	 x = optarglist_scalar_by_key( optargs, "winlength" );
	 if( !isnan( x ) ) winlength=(int)x;
  }
  if( optarglist_has_key( optargs, "N_freq" ) ){
	 x = optarglist_scalar_by_key( optargs, "N_freq" );
	 if( !isnan( x ) ) N_freq=(int)x;
  }
  if( optarglist_has_key( optargs, "N_time" ) ){
	 x = optarglist_scalar_by_key( optargs, "N_time" );
	 if( !isnan( x ) ) N_time=(int)x;
  }
  if( optarglist_has_key( optargs, "corner_freqs" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "corner_freqs" );
	 if( ptr ){
		corner_freqs[0] = ((double*)ptr)[0];
		corner_freqs[1] = ((double*)ptr)[1];
	 }
  }

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
								  corner_freqs, NULL, ALLOC_IN_FCT );
  sp2 = spectrogram_stft( d2, n2, sample_frequency,
								  window, winlength, N_freq, N_time,
								  corner_freqs, NULL, ALLOC_IN_FCT );
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
  spectrogram_free( sp1 );
  spectrogram_free( sp2 );
  free( d1 );
  free( d2 );
  free( window );
  
  return d;
}

/** builds the pointwise distance matrix based on the
	 line-of-synchrony as found in cross-recurrence plots. For details,
	 see 
	 
	 Matthias Ihrke, Hecke Schrobsdorff and J. Michael Herrmann:
	 Recurrence-Based Synchronization of Single Trials for EEG-Data
	 Analysis. Lecture Notes on Computer Science, Proceedings of IDEAL
	 2009. 

	 \param s1,n1,s2,n2 the two signals
	 \param d memory or ALLOC_IN_FCT
	 \param optargs may contain:
	 - <tt>m=int</tt> the embedding dimension for the phase-space reconstruction of the signals, default=10
	 - <tt>tau=int</tt> time-lag for phase-space rec., default=15
	 - <tt>FAN=int</tt> fixed amount of neighbours for each point in the recplot, default=(int)0.05*n1
	 - <tt>noiseamp=double</tt> factor to use for adding noise to the recplot, default=\c 0.01
 */
double** signaldist_recplot_los( const double *s1, int n1, const double *s2, int n2, 
											double **d, 
											OptArgList *optargs ){
  int i,j; 
  double noise;
  double noiseamp;
  RecurrencePlot *R;
  PhaseSpace *p1, *p2;
  int m, tau, FAN;
  double x;

  /* defaults */
  m = 10;
  tau =15 ;
  FAN = (int)0.05*n1;
  noiseamp = 0.01;

  /* override */
  if( optarglist_has_key( optargs, "m" ) ){
	 x = optarglist_scalar_by_key( optargs, "m" );
	 if( !isnan( x ) ) m=(int)x;
  }
  if( optarglist_has_key( optargs, "tau" ) ){
	 x = optarglist_scalar_by_key( optargs, "tau" );
	 if( !isnan( x ) ) tau=(int)x;
  }
  if( optarglist_has_key( optargs, "FAN" ) ){
	 x = optarglist_scalar_by_key( optargs, "FAN" );
	 if( !isnan( x ) ) FAN=(int)x;
  }
  if( optarglist_has_key( optargs, "noiseamp" ) ){
	 x = optarglist_scalar_by_key( optargs, "noiseamp" );
	 if( !isnan( x ) ) noiseamp=x;
  }

  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( n1, n2 );
  }

  dprintf("Using Parameters: m=%i, tau=%i, FAN=%i, noiseamp=%f\n",
			 m, tau, FAN, noiseamp );

  srand((long)time(NULL));
  p1 = phspace_init( m, tau, s1, n1 );
  p2 = phspace_init( m, tau, s2, n2 );

  R = recplot_init( n1, n2, FAN, RPLOT_FAN );
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

  free( p1 );
  free( p2 );
  recplot_free( R );

  return d;
}

/** Compute trial-to-trial distance matrix for all trials in eeg-struct.
	 Average over all channels in the EEG-struct.

	 \param eeg the EEG-data
	 \param f the function used to compare two ERPs
	 \param d user allocated memory, or NULL -> own memory is alloc'ed
	 \param optargs may contain:
	 - "progress=void*", a progress-bar function; default=NULL;
	 - <b> userdata is passed to the vectordist_*() function </b>
 */
double** eeg_distmatrix( EEG *eeg, VectorDistanceFunction f, double **d, OptArgList *optargs ){
  int i,j;
  int c;
  void *tmp;
  ProgressBarFunction progress=NULL;

  if( d==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 d = matrix_init( eeg->ntrials, eeg->ntrials );
  }

  if( optarglist_has_key( optargs, "progress" ) ){
	 tmp = optarglist_ptr_by_key( optargs, "progress" );
	 if( tmp ) progress = (ProgressBarFunction) tmp;
  }
  dprintf("progress=%p\n", progress );

  if( progress ){
	 progress( PROGRESSBAR_INIT, eeg->ntrials );
  }

  for( i=0; i<eeg->ntrials; i++ ){
	 d[i][i] = 0.0;	 
	 if( progress ){
		progress( PROGRESSBAR_CONTINUE_LONG, i );
	 }
	 for( j=i+1; j<eeg->ntrials; j++ ){
		d[i][j] = 0.0;
		for( c=0; c<eeg->nbchan; c++ ){
		  if( progress ){
			 progress( PROGRESSBAR_CONTINUE_SHORT, 0 );
		  }
		  d[i][j] += f( eeg->data[c][i], 
							 eeg->data[c][j],
							 eeg->n, optargs );
		  if( isnan( d[i][j] ) ){
			 errprintf("D[%i][%i] is nan\n", i, j);
		  }
		}
		d[i][j] /= (double)eeg->nbchan;
		d[j][i] = d[i][j];
	 }
  }
  if( progress ){
	 progress( PROGRESSBAR_FINISH, 0 );
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
  double det=0.0;
  double d=0.0;
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
