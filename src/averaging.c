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
double get_warppath(const double *u, int J, const double *s, int K,
	      double theta1, double theta2, int* path){
  double **d, Djk;
  int j,k;
  double avgu, avgs, rmsu=0.0, rmss=0.0;
  double snorm, snormp; /* variable for a point of the normalized
			   signal and the preceding point */
  double unorm, unormp;

  /* d is JxK */
  d = (double**)malloc(J*sizeof(double*));
  for(j=0; j<J; j++)
    d[j]=(double*)malloc(K*sizeof(double));

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
  return Djk;
}
/** same as get_warppath, but the warpPath struct ist used as return value;
 * the path is allocated by this function.
 */
WarpPath* get_warppath2(const double *u, int J, const double *s, int K, double theta1, double theta2){
	double Djk;
	WarpPath *path = get_warppath3(u, J, s, K, theta1, theta2, &Djk);
	return path;
}

/** same as get_warppath2, but the Djk is set
 */
WarpPath* get_warppath3(const double *u, int J, const double *s, int K,	double theta1, double theta2, double *Djk){
	WarpPath *path;
	int *p;
	int i,j,k;
  
	p = (int*)malloc(K*sizeof(int));
	path = (WarpPath*)malloc(sizeof(WarpPath));
	path->J=J;
	path->K=K;
	path->upath = (int*)calloc(J+K, sizeof(int));
	path->spath = (int*)calloc(J+K, sizeof(int));

	(*Djk) = get_warppath( u, J, s, K, theta1, theta2, p );

	/* handle base case */
	path->upath[0]=0;
	path->spath[0]=0;
	if(p[0]>0){
		for(i=0; i<p[0]; i++){
			path->spath[i]=0;
			path->upath[i]=i;
		}		  
	}
	k = p[0];
	for(i=1; i<K; i++){
		/* we are at index k in both arrays */
		if(p[i]==p[i-1]){ // we go in x-direction
			//dprintf("identical\n");
			path->spath[k]=i;
			path->upath[k]=p[i];
			k++;
		} else if(p[i]>p[i-1]){ // diagonal or step
			for(j=p[i-1]; j<=p[i]; j++){			  
				path->spath[k]=i;
				path->upath[k]=j;
			//	dprintf("p(%i)=%i; sp[%i]=%i, up[%i]=%i\n",i,p[i], k,path->spath[k], k, path->upath[k]);
				k++; 
			}
		}
	}
  
	free(p);
	return path;
}

int invert_path(int *path, int K){
  /* path must contain enough memory
   *  - returns length of new array (J)
   */
  int J, i;
  int *tmp;

  J = path[K-1];
  tmp = (int*)malloc(J*sizeof(int));
  for(i=0; i<J; i++){
    /*    tmp[i*/
  }
  
  path = memcpy(path, tmp, J*sizeof(int));
  free(tmp);
  return 0;
}

/** apply warp-path given as integer array to s.
 */
void warp_to_path(double *s, int K, int *path, int J){
  double *snew;
  int j, k, l;

  snew = (double*)calloc(J,sizeof(double));

  j=0; l=1;
  snew[0]=s[0];
  for(k=1; k<K; k++){   
    if(j==path[k]){
      snew[j]+=s[k]; l++;
    }
    else if(j<path[k]){
      snew[j]=snew[j]/((double)l);
      l=1; 
      for(j=j+1; j<=path[k]; j++){ 
	snew[j]+=s[k];
      }
      j=path[k];
    }
  }

  s = memcpy(s, snew, J*sizeof(double));
  free(snew);
}

/** Warpaverage two signals together. Use method described in Picton.
 * \param u,J - sig1
 * \param s,K - sig2
 * \param path - contains warppath
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* avgwarp_from_path(const double *u, int J, const double *s, int K, const WarpPath *P, double *avg){
	double *tmp;
	int i;
	tmp = (double*)calloc(K+J, sizeof(double));
	if(avg==NULL){
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

/** Warpaverage two signals (ADTW). Use method described in Picton.
 * \param u,J - sig1
 * \param s,K - sig2
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* ADTW(const double *s1, int n1, const double *s2, int n2, double *avg){
  WarpPath *p;
  double Djk;

  p = get_warppath3(s1, n1, s2, n2, 1.0, 1.0, &Djk);
  avg = avgwarp_from_path(s1, n1, s2, n2, p, avg);

  return avg;
}

/** Get the warpaverage of two signals computed as in Picton 1995.
 * \param u,s - two signal, both of length n
 * it will be warped between zero->sR and sR->n
 * \param wavg - caller allocated memory; if NULL, own memory is allocated
 * \return a (K+J)/2+1 long warpavg
 */
double* warpavg_two(const double *u, const double *s, int n, int zero, int sRu, int sRs, double *wavg){
	WarpPath *P;
	int i, J, K, len;
	const double *uptr, *sptr;
	
	if(wavg==NULL){
		wavg = (double*)malloc(n * sizeof(double));
	}
	
	
	for(i=0; i<zero; i++) wavg[i] = (u[i]+s[i])/2.0;
	/* first warp, between zero and sR */
	uptr = &(u[zero]);
	sptr = &(s[zero]);
	J = sRu-zero;
	K = sRs-zero;
	dprintf("  J=%i, K=%i, zero=%i, sRu=%i, sRs=%i\n", J, K, zero, sRu, sRs);
	P = get_warppath2(uptr, J, sptr, K, 1.0, 1.0);
	avgwarp_from_path(uptr, J, sptr, K, P, &(wavg[zero]));
	free_warppath(P);
	len = zero+(K+J)/2-1;
	
	/*for(i=zero+len; i<n; i++) wavg[i]=0;*/
	/* second warp, between sR and n */
	uptr = &(u[sRu]);
	sptr = &(s[sRs]);
	J = n-sRu;
	K = n-sRs;
	P = get_warppath2(uptr, J, sptr, K, 1.0, 1.0);
	avgwarp_from_path(uptr, J, sptr, K, P, &(wavg[len]));
	wavg[n-1]=wavg[n-2];
	free_warppath(P);	
	
	return wavg;
}

void idx_to_indices( int idx, int N, int *i, int *j){
	int t;
	t=0;
	for(*i=0; *i<N; (*i)++)
		for(*j=(*i)+1; *j<N; (*j)++)
			if(t==idx) return;
}

/** Warpaveraging N trials.
 * \todo Extend this to arbitrary many markers! TODO
 * \param ui Nxn, trials that are to be combined
 * \param markers Nx2 matrix, containing corresponding zero-markers and response-markers in the u_i and u_j
 * \param wa allocated by caller; if NULL, own memory is allocated.
 */
double* warpaverage(const double **ui, int N, int n, const int **markers, double *wa){
	WarpPath **Ps;
	double *Djk, *tmp;
	int *Djk_idx;
	int i, j, pidx;
	int zero, zeroj, J, K;
	int ii, ji, curN;
	int *seeni, *seenj;
	double **uis;
	double *u, *s;
	int *sRui;
	
	if(wa==NULL){
		wa = (double*)malloc(n*sizeof(double));
	}
/*	Ps = (WarpPath**)malloc((N*(N-1)/2) * sizeof(WarpPath*));
	Djk= (double*)malloc((N*(N-1)/2) * sizeof(double));
	Djk_idx= (int*)malloc((N*(N-1)/2) * sizeof(int));*/
	tmp   = (double*)malloc(n*sizeof(double));
/*	seeni = (int*)malloc(N*sizeof(int));
	seenj = (int*)malloc(N*sizeof(int));*/
	sRui  = (int*)calloc(N,sizeof(int));
	
	/* these are for holding the warped-together ui's */
	uis = (double**)malloc(N*sizeof(double*));
	for(i=0; i<N; i++){
		 uis[i]=(double*)malloc(n*sizeof(double));
		 for(j=0; j<n; j++){
			 uis[i][j] = ui[i][j];
		 }
	}
	
	zero = markers[0][0];
	for(i=0; i<N; i++){
		sRui[i]=markers[i][1];
	}
	
	curN = N;
	while(curN>1){ /* more than one trial left */
		for(i=1; i<curN; i+=2){
			u    = uis[i-1];
			s    = uis[i];
			dprintf("TW(%i, %i), i/2=%i, zero=%i, sRu=%i, sRs=%i, curN=%i\n", i-1, i, i/2, zero, sRui[i-1], sRui[i], curN);
			tmp = warpavg_two(u, s, n, zero, sRui[i-1], sRui[i], tmp);
			
			for(j=0; j<n; j++) uis[i/2][j]=tmp[j];
			sRui[i/2]=(sRui[i-1]+sRui[i])/2;
		}
		curN = curN/2;
		
		/*	
		pidx = 0;
		for(i=0; i<N; i++){
			for(j=i+1; j<N; j++){
				zeroi = markers[i][0];
				J = markers[i][1]-zeroi;
				zeroj = markers[j][0];
				K = markers[j][1]-zeroj;
				Ps[pidx] = get_warppath3( &(ui[i][zeroi]), J,
												&(ui[j][zeroj]), K, 1.0, 1.0, &(Djk[pidx]) );
			}
			seeni[i]=0;
			seenj[i]=0;
		}
		gsl_sort_index (Djk_idx, Djk, 1, (N*(N-1))/2);
		for(i=0; i<(N*(N-1))/2; i++){
			pidx = Djk_idx[i];
			idx_to_indices(pidx, N, &ii, &ji);
			if(seeni[ii]++ | seenj[ji]++) continue;
			zeroi = markers[ii][0];
			zeroj = markers[ji][0];
			J = markers[ii][1]-zeroi;
			K = markers[ji][1]-zeroj;
//			avgwarp_from_path(&(ui[ii][zeroi]), J, 
	//								  &(ui[ji][zeroj]), K, 1.0, 1.0, NULL, NULL);
	}*/
	}
	for(i=0; i<n; i++)
		wa[i] = uis[0][i];
	
	dprintf("1\n");
	/*free(Ps); free(Djk);*/
	dprintf("1\n");
	/*free(seeni); free(seenj);*/
	dprintf("1\n");
	free(tmp); 
	dprintf("1\n");
	for(i=0; i<N; i++) free(uis[i]);
	dprintf("1\n");
	free(uis);
	dprintf("sRui[%i]=%i\n", 4, sRui[4]);
	free(sRui);
	dprintf("1\n");
	return wa;
}

/** - signal s is warped to match length of signal u (J)
 *  - assumes that enough memory is allocated (if J>K)
 *  - theta are weights for the gradient part of the distance measure
 *    (see (1))
 * Return-value:
 *  - is the final Element of the cumulated path D_jk which serves as
 *    a similarity measure between the two curves
 */
double timewarp(const double *u, int J, double *s, int K, 
	      double theta1, double theta2){
  int *path;
  double *snew;
  int j, k, l;
  double Djk;

  path = (int*)   malloc(K*sizeof(int));
  snew = (double*)calloc(J,sizeof(double));

  Djk = get_warppath(u, J, s, K, theta1, theta2, path);
  warp_to_path( s, K, path, J );
// /*  
//   j=0; l=1;
//   snew[0]=s[0];
//   for(k=1; k<K; k++){   
// /*     printf("k=%i p(k)=%i, j=%i\n", k, path[k],j); */
//     if(j==path[k]){
//       snew[j]+=s[k]; l++;
//  /*      printf("(1) snew[%i]=%f, l=%i, j=%i, k=%i, p(k)=%i\n", j, snew[j], l, j, k, path[k]); */
//     }
//     else if(j<path[k]){
//       snew[j]=snew[j]/((double)l);
// /*       printf("(2) snew[%i]=%f, l=%i, j=%i, k=%i, p(k)=%i\n", j, snew[j], l, j, k, path[k]); */
//       l=1; 
//       for(j=j+1; j<=path[k]; j++){ 
// 			snew[j]+=s[k];
// /* 	printf("(3) snew[%i]=%f, l=%i, j=%i, k=%i, p(k)=%i\n", j, snew[j], l, j, k, path[k]); */
//       }
//       j=path[k];
//     }
//   }*/

//   s = memcpy(s, snew, J*sizeof(double));
//   free(snew);
  free(path);
  dprintf("db:Djk=%f\n", Djk);
  return Djk;
}

/** Estimate the hatu by computing the iterative WA after cleaning the si.
 * Function is a 'model' that can e.g. submitted to loocv.
 * \param d ModelData containing everything necessary
 * \param hatu allocated memory for the estimate hatu (returned as pointer)
 */
double* iterative_warpavg(const ModelData *d, double *hatu){
	int i, j;
	double *sRi;
	double *savg;
	double hatR, shatR;
	double zero;
	double **ui;
	
	/* deep copy, because we have to clean them */
	ui   = (double**)malloc(d->N*sizeof(double*));
	for(i=0; i<d->N; i++) {
		ui[i]=(double*)malloc(d->n*sizeof(double));
		for(j=0; j<d->n; j++){
			ui[i][j]=d->si[i][j];
		}
	}
	sRi  = (double*)malloc(d->N*sizeof(double));
	savg = (double*)malloc(d->n*sizeof(double));
	for(i=0; i<d->N; i++)
		sRi[i] = closest_index(d->times, d->n, d->Ri[i]);
	
	/* simple avg */
	for(i=0; i<d->n; i++){
		for(j=0; j<d->N; j++){
			savg[i] += d->si[j][i];
		}
		savg[i] /= (double)d->N;
	}

	/* start with simple avg as estimation for hatu */
	hatu = memcpy(hatu, savg, d->n*sizeof(double));
	hatR = gsl_stats_mean(d->Ri, 1, d->N);
	zero  = closest_index(d->times, d->n, 0);
	shatR = closest_index(d->times, d->n, hatR);
	
	dprintf("+-n=%i\n", d->n);
	extend_and_denoise(hatu, d->n, d->den_params->L, d->den_params->cleanfct,
							  d->den_params->eta, d->den_params->sigextfct);
	
	/* ------- Denoise everything ---------- */
	for(i=0; i<d->N; i++){
		extend_and_denoise(ui[i], d->n, d->den_params->L, d->den_params->cleanfct,
								 d->den_params->eta, d->den_params->sigextfct);		
	}
	
	/**\todo continue here and think of some smart way to get a nice estimate for hatu
	 */

	

	free(sRi); free(savg);
	for(i=0; i<d->N; i++)
		free(ui[i]);
	free(ui);
	return hatu;
}
/** Generate a NxN matrix of differences between single-trial data in m.
 * \param m the data
 * \param dm the difference-matrix allocated by the caller
 * \return a pointer to dm
 */
double** diffmatrix(ModelData *m, double **dm){
	int i, j;
	double **ui;
	int *dummy;
	
	print_modeldata(stderr, m);
	
	dummy = (int*)malloc(m->n*sizeof(int));
	/* deep copy, because we have to clean them */
	ui   = (double**)malloc(m->N*sizeof(double*));
	for(i=0; i<m->N; i++) {
		ui[i]=(double*)malloc(m->n*sizeof(double));
		for(j=0; j<m->n; j++){
			ui[i][j]=m->si[i][j];
		}
	}
	/* ------- Denoise everything ---------- */
//	for(i=0; i<m->N; i++){
	//	extend_and_denoise(ui[i], m->n, m->den_params->L, m->den_params->cleanfct,
		//						 m->den_params->eta, m->den_params->sigextfct);		
	//}
	
	for(i=0; i<m->N; i++){
		for(j=i; j<m->N; j++){
			dm[i][j]=get_warppath( ui[i], m->n, ui[j], m->n, 
									 m->tw_params->theta1, m->tw_params->theta2, dummy);
			dm[j][i]= dm[i][j];
		}
	}
	for(i=0; i<m->N; i++)
		free(ui[i]);
	free(ui);
	free(dummy);
	return dm;
}

