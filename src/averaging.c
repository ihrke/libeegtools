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

#include "averaging.h"

/* ---------------------------------------------------------------------------- 
   -- Timewarping                                                            -- 
   ---------------------------------------------------------------------------- */
/** construct the warppath for  u and s.
 * 
 *  - signal s is warped to match length of signal u (J)
 *  - theta are weights for the gradient part of the distance measure
 *    (see (1))
 *  - path must be K long and will contain the indices which warp s to u
 */
double DTW_get_warppath(const double *u, int J, const double *s, int K,
						  double theta1, double theta2, int* path){
  double **d, Djk;
  int j,k;
  double avgu, avgs, rmsu=0.0, rmss=0.0;
  double snorm, snormp; /* variable for a point of the normalized
									signal and the preceding point */
  double unorm, unormp;

  
  d = DTW_build_cumdistmatrix(u,  J, s, K, theta1, theta2, NULL);
  Djk = d[J-1][K-1];
  
  /* Backtracking */
  j=J-1; k=K-1;
  path[k]=j;
  while(j>0 || k>0){
    if(k==0){
      path[k]=0; break;
    } else if(j==0){
      while(k>=0){
		  path[k]=0; 
		  k--;
      }
      break;
    }
    if(d[j][k-1]<d[j-1][k-1] && d[j][k-1]<d[j-1][k]){
      k--;
      path[k]=j;
    } 
    else {
      if(d[j-1][k]<d[j-1][k-1]){
		  j--;
      } 
      else {
		  j--; k--;
		  path[k]=j;
      }
    }
  }
  
  for(j=0; j<J; j++) free(d[j]);
  free(d);
  //  printf("get_warppath: Djk=%f\n", Djk);
  return Djk;
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
  int i,j,k;
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
	 dprintf("Restriction R=%f too large, using 1.0\n");
	 R = 1.0;
  } else if( R<0 ){
	 dprintf("Restriction R=%f < 0, aborting\n");
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
  for( j=0; j<MAX( J, K ); j++ ){ // J>K
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
/** build the  pointwise (non-cumulated) dissimilarity matrix d
	 \param u,J 1st signal
	 \param s,K 2nd signal
	 \param theta1/theta2 - weights for metric
	 \param d -- if NULL, function allocates the memory
	 \return pointer to JxK matrix d
 */
double** DTW_build_distmatrix(const double *u, int J, const double *s, int K, 
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
  
  /* computing d_jk */
  for(j=0; j<J; j++){
    for(k=0; k<K; k++){
      snorm = (s[k]-avgs)/rmss;
      unorm = (u[j]-avgu)/rmsu;
      (k==0) ? (snormp = 0) : (snormp = ((s[k-1]-avgs)/rmss));
      (j==0) ? (unormp = 0) : (unormp = ((u[j-1]-avgu)/rmsu));
      d[j][k] = theta1 * fabs(unorm - snorm) + theta2 * fabs( (unorm-unormp) - (snorm-snormp) ); /* (1) */
    }
  }

  return d;
}

/** cumulate the pointwise distance-matrix d by choosing
	 \f[
	 \mathbf{D}_{jk} = \mathbf{d}_{jk}+\min{\{\mathbf{D}_{j,k-1}, \mathbf{D}_{j-1,k}, \mathbf{D}_{j-1, k-1}\}}
	 \f]
*/
void DTW_cumulate_distmatrix(double **d, int J, int K){
  int j,k;

  /* computing D_jk */
  for(j=0; j<J; j++){
    for(k=0; k<K; k++){
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

/** construct the warppath from distance matrix. Use DTW_build_cumdistmatrix() for that.
	 \param d - cumulated distmatrix
	 \param J,K - dims for d
	 \param path - pointer to WarpPath-struct of length (J+K) or NULL (allocated in function)
	 \return WarpPath struct
*/
WarpPath* DTW_path_from_cumdistmatrix(const double **d, int J, int K, WarpPath *path){
  int i, j, k;
  int idx;
  double left, down, downleft;

  if( path==NULL ){
	 path = init_warppath(J, K);
  }

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
  }

  return path;
}

/** insert markers into the distance matrix produced by \ref DTW_build_distmatrix()
	 \param d,J,K - JxK matrix of pointwise distances -- is overwritten
	 \param markers - (2xnmarkers) array of indices that are to be set to one;
                     markers[1][i] and markers[2][i] are corresponding points in the
							two eegdata-sequences.
 */
void DTW_markers_to_distmatrix(double **d, int J, int K, const unsigned long **markers, unsigned nmarkers){
  int i, m1_idx, m2_idx;

  for( i=0; i<nmarkers; i++ ){
	 m1_idx = markers[1][i];
	 m2_idx = markers[2][i];
	 d[m1_idx][m2_idx]=0;
  }

}
/** get warppath applying marker-based DTW
	 \param u,J 1st signal
	 \param s,K 2nd signal
	 \param theta1/theta2 - weights for metric
	 \param markers - (2xnmarkers) array of indices
	 \return pointer to WarpPath struct
 */
WarpPath* DTW_warppath_with_markers(const double *u, int J, const double *s, int K,
												double theta1, double theta2, 
												const unsigned long **markers, unsigned nmarkers){ 
  int i;
  double **d;
  WarpPath *P;

  d = DTW_build_distmatrix( u,J,s,K,theta1,theta2, NULL );
  DTW_markers_to_distmatrix( d, J, K, markers, nmarkers );
  DTW_cumulate_distmatrix( d, J, K );
  P = DTW_path_from_cumdistmatrix( (const double**) d, J, K, NULL );
  
  for( i=0; i<J; i++ )
	 free(d[i]);
  free(d);

  return P;
}
/** get warppath applying marker-based DTW
	 \param u,J 1st signal
	 \param s,K 2nd signal
	 \param theta1/theta2 - weights for metric
	 \param markers - (2xnmarkers) array of indices
	 \return Djk
 */
double    DTW_get_warpdistance_markers(const double *u, int J, const double *s, int K,
													double theta1, double theta2, const unsigned long **markers, unsigned nmarkers){
  int i;
  double **d;
  double Djk;

  d = DTW_build_distmatrix( u,J,s,K,theta1,theta2, NULL );
  DTW_markers_to_distmatrix( d, J, K, markers, nmarkers );
  DTW_cumulate_distmatrix( d, J, K );
  Djk = d[J-1][K-1];
  
  for( i=0; i<J; i++ )
	 free(d[i]);
  free(d);

  return Djk;
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
WarpPath* DTW_get_warppath2(const double *u, int J, const double *s, int K,	double theta1, double theta2, double *Djk){
	WarpPath *path;
	double **d;
	int i,j,k;
  
	d = DTW_build_cumdistmatrix( u,J,s,K,theta1,theta2,NULL );
	path=DTW_path_from_cumdistmatrix( (const double**)d, J,K, NULL );

	for( i=0; i<J; i++) free(d[i]);
	free(d);

	return path;
}

/** Warpaverage two signals together. Use method described in Picton.
 * \param u,J - sig1
 * \param s,K - sig2
 * \param path - contains warppath
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* ADTW_from_path(const double *u, int J, const double *s, int K, const WarpPath *P, double *avg){
	double *tmp;
	int i;
	tmp = (double*)calloc(K+J, sizeof(double));
	if(avg==NULL){
	  dprintf("Allocating own memory\n");
	  avg = (double*)calloc((K+J)/2+1, sizeof(double));
	}
	
	for(i=0; i<J+K; i++){
		tmp[i] = (u[P->upath[i]] + s[P->spath[i]])/2.0;
	}
	
	for(i=0; i<J+K; i+=2){
		if(i+1<J+K)
			avg[i/2]=(tmp[i]+tmp[i+1])/2.0;
		else
			avg[i/2]=tmp[i];
	}
	
	free(tmp);
	return avg;
}

/** Calculate a simple, pointwise average of two vectors. 
	 \f[ \hat{s}(t) = \frac{(s_1(t)+s_2(t))}{2} \f]
	 \param avg - user-allocated memory of length n; if NULL, the function 
	              allocates the memory
 */
double* simple_average_2v(const double *s1, const double *s2, int n, double *avg){
  int i;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++)
    avg[i] = (s1[i]+s2[i])/2.0;
  return avg;
}

/** Calculate a simple, pointwise average of n vectors
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N}s_i(t)
	 \f]
 * \param avg - user-allocated memory of length n; if NULL, the function 
 *              allocates the memory
 */
double* simple_average_nv(const double **s, int N, int n, double *avg){
  int i, j;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++) {
	 avg[i] = 0;
	 for( j=0; j<N; j++ ){
		avg[i] += s[i][j];
	 }
	 avg[i] /= (double) N;
  }
  return avg;
}

/** Calculate alternate average of n vectors
	 \f[
	 \hat{s}(t) = \frac{1}{N}\sum_{i=1}^{N} (-1)^{i-1}s_i(t)
	 \f]
 * \param avg - user-allocated memory of length n; if NULL, the function 
 *              allocates the memory
 */
double* alternate_average_nv(const double **s, int N, int n, double *avg){
  int i, j;
  if(!avg)
    avg = (double*)malloc(n*sizeof(double));
  for(i=0; i<n; i++) {
	 avg[i] = 0;
	 for( j=0; j<N; j++ ){
		avg[i] += pow(-1, j-1)*s[i][j];
	 }
	 avg[i] /= (double) N;
  }
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

/** Warpaverage two signals (ADTW). Do the warping for a time signal
 *  where stimulus onset (zero) and reaction time (sR) are given.
 *   0,...,zero -- simple average
 *   zero,...,sR-- warpavg
 *   sR,...,end -- warpavg.
 * \param s1 - sig1
 * \param s2 - sig2
 * \param zero  - stimulus onset in both signals
 * \param n     - num sampling points of both signals
 * \param sR1   - reaction time for signal 1 (sampling units)
 * \param sR2   - reaction time for signal 2
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* ADTW_signal(const double *s1, int sR1, 
						  const double *s2, int sR2, 
						  int zero, int n, double *avg){
  WarpPath *p;
  double Djk;
  dprintf("ADTW_signal: zero=%i, sR1=%i, sR2=%i\n", zero, sR1, sR2);
  if(!avg) avg = (double*)malloc(n*sizeof(double));
  avg = simple_average_2v( s1, s2, zero, avg );
  /*&(avg[zero]) = */
  ADTW(&(s1[zero]), sR1-zero+1, &(s2[zero]), sR2-zero+1, &(avg[zero]));
  /*  &(avg[((sR1-zero)+(sR2-zero))/2+1])=*/
  ADTW(&(s1[sR1]), n-sR1, &(s2[sR2]), n-sR2, &(avg[(sR1+sR2)/2]));

  return avg;
}


/** Warpaverage N signals (PADTW). Use method described in Ihrke Bachelor.
 * \param s,N,n - N segments of EEG-data, each n sampling points
 * \param zero  - the zero marker (stimulus onset) in sampling units
 * \param sR    - N reaction times in sampling units for each trial
 * \param wa    - pointer to caller-allocated memory (length n)
 *                if NULL is given, memory is allocated by the function.
 * \note NOT FUNCTIONAL YET!
 */  
double* PADTW(const double **s, int N, int n, int zero, int *sR, double *wa){
  WarpPath **Ps;
  double **sc;
  double **Djk;
  double **tmp;
  double *tmp2;
  int i,j,k, l;
  int msize;
  int curN;
  

  //  dprintf("PATDW: N=%i,n=%i,zero=%i,sR[0]=%i,sR[N-1]=%i,s[0][0]=%f,s[0][1]=%f,s[1][0]=%f,s[1][1]=%f\n",
  //	 N,n,zero,sR[0],sR[N-1],s[0][0],s[0][1],s[1][0],s[1][1]);
  /* initialization */
  sc = copy_double_ptrptr(s, N, n);
  Ps = (WarpPath**)malloc(N*sizeof(WarpPath*));
  tmp= (double**)malloc(N*sizeof(double*));
  Djk = (double**)malloc(N*sizeof(double*));
  for(i=0; i<N; i++)
    Djk[i] = (double*)calloc(N,sizeof(double));

  /* main loop -- iterate until top of pyramid */
  curN=N;
  while(curN>1){
	 dprintf("curN=%i\n", curN);
    /* calculate the distance matrix and pathes */
	 for(i=0; i<curN;i++)
		for(j=0; j<curN;j++)
		  Djk[i][j]=DBL_MAX;

    for(i=0; i<curN; i++){
      for(j=i+1; j<curN; j++){
		  dprintf("(%i,%i)\n", i,j);
		  DTW_get_warppath2(&(sc[i][zero]), sR[i],
										 &(sc[j][zero]), sR[j], 1.0,1.0,
										 &(Djk[i][j]));
		  //		  Djk[j][i]=Djk[i][j];
      }
    }

    msize = curN;
    l = 0;
    while(l<curN/2){
		matrix_print(Djk, curN, curN);
      matrix_min( (const double**)Djk, curN, curN, &i, &j);
		dprintf("argmin(Djk)=(%i,%i)\n",i,j);
      tmp[l]=(double*)malloc(n*sizeof(double));
      tmp[l] = ADTW_signal(sc[i], sR[i], sc[j], sR[j], zero, n, tmp[l]);
		for(k=0; k<curN; k++){ // delete row and col
		  Djk[i][k]=DBL_MAX; 
		  Djk[j][k]=DBL_MAX;		  
		}
      l++;
    }
	 for(i=0; i<l; i++){
		tmp2 = sc[i];
		sc[i]=tmp[i];
		free(tmp2);
	 }
	 curN = curN/2;
  }
  
  for(i=0; i<n; i++){
	 wa[i]=sc[0][i];
  }

  for(i=0; i<N; i++){
    free(sc[i]);
  }
  free(sc);
  free(Ps);
  return wa;
}

double DTW_distance_between_paths(const WarpPath *P1, const WarpPath *P2){
  int i, J, K;
  double dist=0.0;


  if(P1->J != P2->J || P1->K != P2->K){
	 errprintf( "P1 and P2 not comparable (P1->J, P2->J, P1->K, P2->K)\n", P1->J, P2->J, P1->K, P2->K);
  }
  J = P1->J;
  K = P1->K;

  for( i=0; i<(J + K); i++ ){
	 dist += ABS( P1->upath[i] - P2->upath[i] ) + ABS( P1->spath[i] - P2->spath[i] );
  }
  dist = (double)dist/(double)(J+K);

  return dist;
}


/** compute multiple warp-Pathes by computing restricted warppathes 
	 with restriction parameter \f$\theta\f$ for each of the segments 
	 from marker 0...i, i...i+1, ..., N-1...N (see \ref timewarping).
	 \param s1,s2 data
	 \param channel channel to use in the dataset
	 \param theta restriction parameter (Chiba-Band)
	 \param P pointer to Warppath-structs (nmarker many); if NULL, own memory is allocated
 */
WarpPath* eeg_DTW_get_paths_by_markers( const EEGdata *s1, const EEGdata *s2, int channel, double theta, WarpPath *P ){
  int i,j, N,n;
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
  mark1 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark2 = (unsigned long*) malloc( (N+2)*sizeof( unsigned long ) );
  mark1[0]   = 0;   mark2[0]   = 0;
  mark1[N+1] = n-1; mark2[N+1] = n-1;
  memcpy( mark1, s1->markers, N*sizeof(unsigned long) );
  memcpy( mark2, s2->markers, N*sizeof(unsigned long) );
  dist = matrix_init( n, n );

  N += 2;
  for( i=1; i<N; i++ ){
	 J = mark1[i] - mark1[i-1];
	 K = mark2[i] - mark2[i-1];

	 Pptr = &(P[i]);
	 /* prepare warppath */
	 if( P==NULL ){
		Pptr = init_warppath( J, K );
	 } else {
		reset_warppath( Pptr, J, K );
	 }

	 dist = DTW_build_restricted_cumdistmatrix( s1->d[channel], J, 
															  s2->d[channel], K, 
															  theta, dist );  
	 Pptr = DTW_path_from_cumdistmatrix( (const double**) dist, J, K, Pptr);
  }

  free(mark1); free(mark2);
  matrix_free(dist, n);

  return Pptr;
}

EEGdata* eeg_ADTW_from_path(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel, const WarpPath *P){
  if(target==NULL){
	 target = init_eegdata(s1->nbchan, (s1->n+s2->n)/2, 0);
	 dprintf("ALLOC: allocated memory in function!\n");
  }
  ADTW_from_path( s1->d[channel], s1->n,
						s2->d[channel], s2->n,
						P, target->d[channel] );

  return target;
}


void eeg_ADTW_markers_channel(const EEGdata *s1, const EEGdata *s2, EEGdata *target, int channel){
  WarpPath *P;
  unsigned long **markers;
  int i;
  
  markers = (unsigned long**) malloc( 2*sizeof(unsigned long*) );
  for( i=0; i<2; i++ ) 
	 markers[i] = (unsigned long*) malloc( s1->nmarkers*sizeof(unsigned long) );
  P = DTW_warppath_with_markers( s1->d[channel], s1->n, s2->d[channel], s2->n, 1, 1, markers, s1->nmarkers );
  ADTW_from_path( s1->d[channel], s1->n, s2->d[channel], s2->n, P, target->d[channel]);
}
