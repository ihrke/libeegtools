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
 */
double** vectordist_distmatrix( VectorDistanceFunction f, 
										  const double **X, int n, int p, double **D, 
										  void *userdata ){
  int i,j;
  
  if( D==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 D = matrix_init( n, n );
  }
  
  for( i=0; i<n; i++ ){
	 for( j=i+1; j<n; j++ ){
		D[i][j] = f( (double*)X[i], (double*)X[j], p, userdata );
		D[j][i] = D[i][j];
		if( isnan( D[j][i] ) ){
		  errprintf("D[%i][%i] is nan\n", j, i);
		}
	 }
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

/** compute the cumulated sum along the regularized warping path.
	 \param userdata (int)userdata[0] is nmarkers; 
	       ((int)userdata)[1] and following is an 2 x nmarkers matrix 
			 given the corresponding markers in the signals.
 */
double   vectordist_regularized_dtw( double *x1, double *x2, int p, void *userdata ){
  int nmarkers;
  int **markers;
  double **G, **D;
  double maxsigma;
  double Djk;

  nmarkers = *((int*)userdata);
  userdata = ((int*)userdata)+1;
  maxsigma = *((double*)userdata);
  userdata = ((double*)userdata)+1;
  markers = (int**)userdata;
  
  G = regularization_gaussian_line( markers[0], markers[1], nmarkers, p, maxsigma, NULL );
  /** TODO **/
  return Djk;
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

  theta1 = *((double*)userdata);
  userdata = ((double*)userdata)+1;
  theta2 = *((double*)userdata);

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
  settings = (SettingsHierarchicalDTW*)params;
  double theta[2];
  theta[0]=settings->theta1;
  theta[1]=settings->theta2;

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

/** compute trial-to-trial distance matrix for all trials in eeg-struct.
	 \todo implement that the called distance function incorporates time-markers or whatever is needed
	 via the userdata
	 \param eeg the EEG-data
	 \param f the function used to compare two ERPs
	 \param channel which electrode-channel
	 \param d user allocated memory, or NULL -> own memory is alloc'ed
 */
double** eegtrials_distmatrix_channel( EEGdata_trials *eeg, VectorDistanceFunction f, int channel, double **d ){
  int i,j;

  if( d==ALLOC_IN_FCT ){
	 warnprintf(" allocating matrix in fct\n");
	 d = matrix_init( eeg->ntrials, eeg->ntrials );
  }
  
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=i+1; j<eeg->ntrials; j++ ){
		d[i][j] = f( (eeg->data[i]->d[channel]), 
						 (eeg->data[j]->d[channel]), eeg->nsamples, NULL );
		d[j][i] = d[i][j];
		if( isnan( d[j][i] ) ){
		  errprintf("D[%i][%i] is nan\n", j, i);
		}
	 }
  }
  
  return d;
}

double pointdist_euclidean(double x, double y){
  return fabs(x-y);
}
