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
	 \param params is a SettingsPADTW struct that should contain valid values for theta1,theta2 weights for distance
	 \return d or NULL (error)
 */
double** eeg_distmatrix_euclidean_derivative_channel( const EEGdata *s1, const EEGdata *s2, 
																		int channel, double **d, 
																		void *params ){ 
  /* get params from void pointer */
  SettingsPADTW *settings;
  settings = (SettingsPADTW*)params;
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
	 \param s1,s2 
	 \param d matrix or NULL (alloc'd in function)
	 \param params is a settingsPADTW struct that should contain valid:
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
double** signaldist_stft( const double *s1, int n1, const double *s2, int n2, double **d, void *userdata ){
  int i, j;
  Spectrogram *sp1, *sp2;
  double *window;
  double *d1, *d2,
	 mean1, mean2;
  
  double sample_frequency;
  double*(*winfct)(double*,int);
  int winlength;
  int N_freq;
  int N_time;
  double corner_freqs[2];
  /* get params from void pointer */
  SettingsPADTW *settings;
  settings = (SettingsPADTW*)userdata;
  sample_frequency=settings->sampling_rate;
  winfct=settings->winfct;
  winlength=settings->winlength;
  N_freq=settings->N_freq;
  N_time=settings->N_time;
  corner_freqs[0]=settings->corner_freqs[0];
  corner_freqs[1]=settings->corner_freqs[1];

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
	 \param params is a settingsPADTW struct that should contain valid:
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

  d = signaldist_stft( s1->d[channel], s1->n, s2->d[channel], s2->n, d, params );
  return d;
}


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
  dprintf("path=%p,path->J=%i, path->K=%i\n",path, path->J, path->K);
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
  path->upath[idx] = j;
  path->spath[idx] = k;
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
	 /* 	path->upath[idx] = j+1; */
	 /* 	path->spath[idx] = k; */
	 /* 	idx--; */
	 /* 	path->upath[idx] = j; */
	 /* 	path->spath[idx] = k; */
	 /* } else { */
		path->upath[idx] = j;
		path->spath[idx] = k;
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


/* ---------------------------------------------------------------------------- 
   -- Timewarping                                                            -- 
   ---------------------------------------------------------------------------- */


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
  dprintf("path=%p,path->J=%i, path->K=%i\n",path, path->J, path->K);
  reset_warppath(path, J, K);


  /* Backtracking */
  j=J-1; k=K-1;

  idx = 1;
  path->upath[0] = j;
  path->spath[0] = k;
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


	 path->upath[idx] = j;
	 path->spath[idx] = k;
	 idx++;
	 /*dprintf("(j, k), idx = (%i,%i), %i\n", j, k, idx);*/
  }
  dprintf("leaving\n");
  return path;
}

/** same algorithm as get_warppath but does not do the backtracking
	 but only returns the distance according to DTW.
	 \param s is warped to match length of signal u (J)
	 \param theta are weights for the gradient part of the distance measure
	 (see (1))
	 \return Djk
 
*/
double DTW_get_warpdistance(const double *u, int J, const double *s, int K,
								double theta1, double theta2){
  double **d, Djk;
  int j;
 
  d = DTW_build_cumdistmatrix(u, J, s, K, theta1, theta2, NULL);
  Djk = d[J-1][K-1];

  for(j=0; j<J; j++) free(d[j]);
  free(d);
  return Djk;
}

/** same as get_warppath, but the Djk is set and 
	 the warpPath struct ist used as return value;
	 the path is allocated by this function.
 */
WarpPath* DTW_get_warppath2(const double *u, int J, 
									 const double *s, int K,	
									 double theta1, double theta2, double *Djk){
	WarpPath *path;
	double **d;
	int i;
  
	d = DTW_build_cumdistmatrix( u,J,s,K,theta1,theta2,NULL );
	path=DTW_path_from_cumdistmatrix( (const double**)d, J,K, NULL );

	for( i=0; i<J; i++) free(d[i]);
	free(d);

	return path;
}


/** Add magnitude of signals according to warppath P.
	 Data is resampled to have length (J+K)/2+1 using cspline interpolation.
 * \param path - contains warppath
 * \param weights - weights in average, for using it with hierarchical averaging
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* DTW_add_signals_by_path(const double *s1, int n1, 
										  const double *s2, int n2, 
										  const WarpPath *P, double *avg, 
										  const double weights[2]){
  double *x, *y, *xp;
  int idx, newn, i;
  
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

  x[0] = 0;
  y[0] = weights[0]*(s1[0]) + weights[1]*(s2[0]); /* first entry, because it is skipped in path */
  idx = 1;
  for(i=0; i<n1+n2; i++){
	 if( P->upath[i]==0 && P->spath[i]==0 ){
		continue;
	 }
	 x[idx] = (double)((P->upath[i])+(P->spath[i]))/2.0; /* average latencies */
	 y[idx] = weights[0]*(s1[P->upath[i]]) + weights[1]*(s2[P->spath[i]]); /* average magnitude */
	 dprintf( "x[%i]=(%f,%f), (%i, %i)\n", idx, x[idx], y[idx], (P->upath[i]), (P->spath[i]) );
	 idx++;
  }

  dprintf("idx=%i\n", idx );
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, idx);
  gsl_spline_init (spline, x, y, idx);

  for( i=0; i<newn; i++ ){
	 dprintf("xp[%i]=%f\n", i, xp[i] ); 
	 avg[i] = gsl_spline_eval( spline, xp[i], acc );
	 dprintf(" avg[%i] = (%f, %f)\n", i,xp[i], avg[i] );
  }
  /* double *tmp; */
  /* 	int i, idx; */

  /* 	tmp = (double*)calloc(n1+n2, sizeof(double)); */
  /* 	if(avg==NULL){ */
  /* 	  dprintf("Allocating own memory\n"); */
  /* 	  avg = (double*)calloc((n1+n2)/2+1, sizeof(double)); */
  /* 	} */
	
  /* 	idx = 0; */
  /* 	for(i=0; i<n1+n2; i++){ */
  /* 	  if( P->upath[i]==0 && P->spath[i]==0 ){ */
  /* 		 continue; */
  /* 	  } */
  /* 	  tmp[idx++] = weights[0]*(s1[P->upath[i]]) + weights[1]*(s2[P->spath[i]]); */
  /* 	} */
  /* 	dprintf("J=%i,K=%i, idx=%i\n", n1, n2, idx);	 */
  /* 	avg = resample_gsl( tmp, idx, (n1+n2)/2+1, avg, gsl_interp_cspline ); */

  /* 	free(tmp); */
  /* 	return avg; */ 

  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free( x );
  free( y );
  free( xp );

  return avg;
}

/** Warpaverage two signals together. Use method described in Picton.
	 Data is resampled to have length (J+K)/2+1 using cspline interpolation.
 * \param u,J - sig1
 * \param s,K - sig2
 * \param path - contains warppath
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* ADTW_from_path(const double *u, int J, 
							  const double *s, int K, 
							  const WarpPath *P, double *avg){
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
	dprintf("J=%i,K=%i, idx=%i\n", J, K, idx);	
	avg = resample_gsl( tmp, idx, (J+K)/2+1, avg, gsl_interp_cspline );

	free(tmp);
	return avg;
}

/** Warpaverage two signals (ADTW). Use method described in Picton.
 * \param u,J - sig1
 * \param s,K - sig2
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* ADTW(const double *s1, int n1, const double *s2, int n2, double *avg){
  WarpPath *p;
  double Djk;

  p = DTW_get_warppath2(s1, n1, s2, n2, 1.0, 1.0, &Djk);
  avg = ADTW_from_path(s1, n1, s2, n2, p, avg);

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
EEGdata* eegtrials_PADTW( EEGdata_trials *eeg_in, const double **distmatrix, int N, 
								  EEGdata *out, SettingsPADTW settings ){
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


  if( N>eeg_in->ntrials ){
	 errprintf(" distance matrix contains more trials than eegdata. this is fatal\n");
	 return NULL;
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
  while( trials_left >= 2 ){
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
	 G = settings.regularize( s1, s2, settings.sigma, G); /* we need G only 
																				once for all channels */

	 int nchan,chan;
	 /* loop this for all channels */
	 if(settings.channels==NULL){
		channels = linspace(0,num_chan-1);
		nchan = num_chan;
	 } else {
		channels = settings.channels;
		nchan = settings.num_channels;
	 }
	 for( c=0; c<nchan; c++){
		chan = channels[c];
		oprintf("Trials (%i,%i): Channel=%i\n", idx1, idx2, chan); 
		dprintf(" compute d\n");
		d = settings.trialdistance( s1, s2, chan, d, (void*)(&settings) );
		dprintf(" ...done\n");  
		dprintf(" Regularize d\n");
		maxdist = matrix_max( (const double**)d, n, n, NULL, NULL );
		scalar_minus_matrix( maxdist, d, n, n );
		matrix_dottimes_matrix( d, (const double**)G, n, n );
		scalar_minus_matrix( maxdist, d, n, n ); 
		dprintf(" ...done\n");  
		dprintf(" Compute path\n");  
		P = DTW_path_from_square_distmatrix( (const double**)d, n, P );
		dprintf(" ...done\n");  
		dprintf(" Warpavg\n");  
		new = eeg_DTW_add_signals_by_path( s1, s2, new, chan, P, weights );
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
  

  /* cleaning up */
  dprintf("Freeing Memory\n");
  matrix_free( d, n );
  matrix_free( G, n );
  dgram_free( T );
  free( indices );
  free_warppath( P );

  return out;
}


double DTW_distance_between_paths(const WarpPath *P1, const WarpPath *P2){
  int i, J, K;
  double dist=0.0;


  if(P1->J != P2->J || P1->K != P2->K){
	 errprintf( "P1 and P2 not comparable (P1->J, P2->J, P1->K, P2->K)=(%i,%i,%i,%i)\n", P1->J, P2->J, P1->K, P2->K);
  }
  J = P1->J;
  K = P1->K;

  for( i=0; i<(J + K); i++ ){
	 dist += ABS( P1->upath[i] - P2->upath[i] ) + ABS( P1->spath[i] - P2->spath[i] );
  }
  dist = (double)dist/(double)(J+K);

  return dist;
}
EEGdata* eeg_DTW_add_signals_by_path(const EEGdata *s1, const EEGdata *s2, 
												 EEGdata *target, int channel, 
												 const WarpPath *P, const double weights[2]){
  int i;

  if(target==NULL){
	 target = init_eegdata(s1->nbchan, (s1->n+s2->n)/2, 0);
	 warnprintf("ALLOC: allocated memory in function!\n");
  }
  DTW_add_signals_by_path( s1->d[channel], s1->n,
									s2->d[channel], s2->n,
									P, target->d[channel], weights );
  
  /* set time-markers */
  for( i=0; i<s1->nmarkers; i++ ){
	 target->markers[i] = (s1->markers[i]+s2->markers[i])/2;
	 dprintf( "New marker for target[%i]=%i\n", i, (int)target->markers[i] );
  }

  return target;
}


EEGdata* eeg_ADTW_from_path(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, const WarpPath *P){
  int i;

  if(target==NULL){
	 target = init_eegdata(s1->nbchan, (s1->n+s2->n)/2, 0);
	 warnprintf("ALLOC: allocated memory in function!\n");
  }
  ADTW_from_path( s1->d[channel], s1->n,
						s2->d[channel], s2->n,
						P, target->d[channel] );
  
  /* set time-markers */
  for( i=0; i<s1->nmarkers; i++ ){
	 target->markers[i] = (s1->markers[i]+s2->markers[i])/2;
	 dprintf( "New marker for target[%i]=%i\n", i, (int)target->markers[i] );
  }

  return target;
}


/** compute the  ADTW for one channel in s1,s2 and put it into target. Ignore time-markers.
	 \param s1,s2 sources
	 \param target output
	 \param channel index to use
	 \param theta restriction parameter (see \ref timewarping)
 */
void eeg_ADTW_channel(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta ){
  double **dist;
  WarpPath *P;

  dist = matrix_init( s1->n, s2->n );

  dist = DTW_build_restricted_cumdistmatrix( s1->d[channel], s1->n,  s2->d[channel], s2->n, theta, dist );
  P = DTW_path_from_cumdistmatrix( (const double**) dist, s1->n, s2->n, NULL );
  target->d[channel] = ADTW_from_path( s1->d[channel], s1->n, s2->d[channel], s2->n, P, target->d[channel] );
  
  matrix_free( dist, s1->n );
  free_warppath( P );
}

/** compute the ADTW for all channels in s1,s2 and put it into target. Ignore time-markers.
	 \param s1,s2 sources
	 \param target output
	 \param theta restriction parameter (see \ref timewarping)
 */
void eeg_ADTW(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, double theta ){
  int chan, nchan;

  if( eegdata_cmp_settings( s1, s2 ) ){
	 errprintf("cannot ADTW the to EEGdata structs, because they are too different\n");
	 return;
  }
  nchan = s1->nbchan;

  for( chan=0; chan<nchan; chan++ ){
	 dprintf("eeg_ADTW_channel for Channel: %i\n", chan);
	 eeg_ADTW_channel( s1, s2, target, chan, theta );
  }
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
void _settings_guess_from_eegtrials( const EEGdata_trials *eeg, SettingsPADTW *settings ){
  int n = eeg->nsamples;
  settings->winlength = MAX( SQR( sqrt(next_pow2( n ))-3 ), 4 );
  settings->N_time=n;
  settings->N_freq=(settings->winlength)*2;
  settings->sampling_rate=eeg->sampling_rate;
  settings->corner_freqs[0]=0.0;
  settings->corner_freqs[1]=settings->sampling_rate/2.0;
}
/** \endcond */

SettingsPADTW init_PADTW( const EEGdata_trials *eeg ){	 
  SettingsPADTW settings;
  settings.regularize=eeg_regularization_gaussian_line;
  settings.linkage=dgram_dist_completelinkage;
  settings.dont_touch_eeg=1;
  settings.trialdistance=eeg_distmatrix_stft_channel;
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
  _settings_guess_from_eegtrials( eeg, &settings );
  return settings;
}
