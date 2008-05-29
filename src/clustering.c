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

#include "clustering.h"


/** standardize matrix d (NxN) by computing
	 \f$
	 \hat{d}_{jk} = \frac{1}{\max d_{jk}} d_{jk}
	 \f$
*/
void diffmatrix_standardize(double **d, int N){
  double maxv;
  int i,j;

  maxv = -DBL_MAX;
  for( i=0; i<N; i++ ){
	 for( j=0; j<N; j++ ){
		if(d[i][j]>maxv)
		  maxv = d[i][j];
	 }
  }
  for( i=0; i<N; i++ ){
	 for( j=0; j<N; j++ ){
		d[i][j] /= maxv;
	 }
  }
    
}

/** return the distance of two ERPs from timewarping the complete ERP 
	 (no time-marker-based warping).
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_tw_complete(EEGdata *s1, EEGdata *s2, int channel){
  int *P;
  double Djk;
  P = (int*) malloc( s1->n*sizeof( int ) );
  Djk = get_warppath( s1->d[channel], s1->n, s2->d[channel], s2->n, 1.0, 1.0, P );
  free(P);

  return Djk;
}

/** return the distance of two ERPs from timewarping the  ERP 
	 based on the time-markers
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_tw_markers(EEGdata *s1, EEGdata *s2, int channel){
  int *P;
  double Djk;
  P = (int*) malloc( s1->n*sizeof( int ) );
  Djk = get_warppath( s1->d[channel], s1->n, s2->d[channel], s2->n, 1.0, 1.0, P );
  free(P);

  return Djk;
}
/** return the distance of two ERPs from euclidean distances.
	 \f$
	 \Delta(s_i, s_j) = ||s_i-s_j||^2 = \sum_t d_t(s_i(t), s_j(t))
	 \f$
	 \param s1,s2 ERP-signals
	 \param channel channel in the EEGdata-struct to be used
	 \return a number given the distance of the two ERPs
 */
double clusterdist_euclidean_pointwise(EEGdata  *s1, EEGdata *s2, int channel){
  int i;
  double dist;

  massert( s1->n!=s2->n, "n1 != n2, aborting\n" );
  dist = 0.0;
  for(i=0; i<s1->n; i++){
	 dist += SQR( s1->d[channel][i]-s2->d[channel][i] );
  }
  return dist;
}

/** Return a matrix of the differences between trials in 
	 the EEGdata_trials struct (NxN) for a given channel.
	 \param eeg 
	 \param dist - the distance function
	 \param channel - which electrode-channel? (index)
*/
double** eegtrials_diffmatrix_channel(EEGdata_trials *eeg, 
												  double(*dist)(EEGdata*,EEGdata*,int), 
												  int channel){
  double **dm;
  int i,j;

  dm = (double**) malloc( eeg->ntrials*sizeof( double* ) );
  for( i=0; i<eeg->ntrials; i++ ){
	 dm[i] = (double*) malloc( eeg->ntrials*sizeof( double ) );
	 dm[i][i] = 0;
  }
  
  for( i=0; i<eeg->ntrials; i++ ){
	 for( j=i+1; j<eeg->ntrials; j++ ){
		dprintf("getting distance %i x %i\n", i,j);
		dm[i][j] = dist( eeg->data[i], eeg->data[j], channel );
		dm[j][i] = dm[i][j];
	 }
  }
  return dm;
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

/** do K-Medoids clustering on the distance-matrix.
	 \param dist
	 \param K 
	 \param clusters - function-allocated memory. contains
	                   K cluster-structs.
*/
Clusters* kmedoids(const double **dist, int N, int K){
  Clusters *C, *Cnew;
  int *medoids;
  int i, j, k, r, minidx, num_not_changed;
  double sum, minsum, mindist, curdist;
  int count=0; 

  /* memory alloc */
  medoids = (int*) malloc( K*sizeof(int) );
  C    = init_cluster(K, N);
  Cnew = init_cluster(K, N);

  /* initialization step, random init */
  srandom( (unsigned int)time( (time_t *)NULL ) );
  /*random() / (RAND_MAX / N + 1); -- value from 0,...,N-1 */
  for( i=0; i<N; i++ ){
	 r = (random() / (RAND_MAX / K+1));
	 C->clust[r][C->n[r]] = i;
	 (C->n[r])++;
  }

/*   dprintf("initialized cluster\n");  */
/*   print_cluster(C); */

  num_not_changed = 0;
  while(num_not_changed<1){ /* iterate until convergence */
	 for( i=0; i<K; i++ ){ /* minimize distance to one of trials */
		minsum = DBL_MAX;
		minidx = 0;
		for( j=0; j<C->n[i]; j++ ){
		  sum = 0;
		  for( k=0; k<C->n[i]; k++){
			 sum += dist[C->clust[i][j]][C->clust[i][k]];
		  }
		  
		  if( sum<minsum ){
			 /* fprintf(stderr, "sum=%f, minsum=%f\n", sum, minsum);*/
			 minsum = sum;
			 minidx = C->clust[i][j];
		  }
		}
		medoids[i] = minidx;
	 }

	 /* reset cluster */
	 for( i=0; i<Cnew->K; i++ ){
		Cnew->n[i] = 0;
	 }
	 
	 for( i=0; i<N; i++ ){ /* reassign clusters */
		mindist = DBL_MAX;
		minidx = 0;
		for( j=0; j<K; j++ ){ /* look for closest center */
		  curdist = dist[medoids[j]][i];
		  if( curdist<mindist ){
			 mindist=curdist;
			 minidx=j;
		  }
		}
		/*dprintf("Assigning trial %i to cluster %i\n", i, minidx);*/
		Cnew->clust[minidx][Cnew->n[minidx]] = i;
		Cnew->n[minidx] += 1;
	 }
/* 	 if(count>=100000){ */
/* 		fprintf(stderr, "count=%i\n", count); */
/* 		fprintf(stderr, "medoids = (%i, %i)\n", medoids[0], medoids[1]); */
/* 		fprintf(stderr, "C=\n"); */
/* 		print_cluster(C); */
/* 		fprintf(stderr, "Cnew=\n"); */
/* 		print_cluster(Cnew); */
/* 		count = 0; */
/* 	 } else { */
/* 		count++; */
/* 	 } */
	 if(!compare_clusters(C, Cnew))
		num_not_changed++;
	 else
		copy_cluster(C, Cnew);
  }

/*   print_cluster(C); */
/*   fprintf(stderr, "medoids = (%i, %i)\n", medoids[0], medoids[1]); */
  free(medoids);
  free_cluster(Cnew);
  return C;
}

/** Compare whether two cluster structs are "identical", meaning
	 that each cluster contains the same arguments (sequence
	 does not matter). Superficial tests are first carried out, 
	 then intense comparison.
*/
int compare_clusters(const Clusters *c1, const Clusters *c2){
  int i, j, k, l, found, el;
  int c2i, prev_c2i;

  if( c1->K != c2->K ) return 1;
  /* clusters with equal number of trials? */
  for( i=0; i<c1->K; i++ ){
	 found = 0;
	 for( j=0; j<c2->K; j++ ){
		if(c1->n[i] == c2->n[j]){
		  found = 1;
		  break;
		}
	 }
	 if(!found)		
		return 1;
  }

  /* deep comparison */
  for( i=0; i<c1->K; i++ ){ /* each cluster from c1 */ 
	 dprintf("i = %i\n", i);
	 prev_c2i = c1->K+1;
	 c2i = 0;
	 for( j=0; j<c1->n[i]; j++ ){ /* each element from c1 */
		el = c1->clust[i][j];
		/*		dprintf("j = %i\n", j);*/
		for( k=0; k<c2->K; k++ ){
		  for( l=0; l<c2->n[k]; l++ ){
			 if( el == c2->clust[k][l] ){
				c2i = k;
				if(prev_c2i==c1->K+1)
				  prev_c2i = k; /* first run */
			 } 
		  }
		}
		if( c2->n[c2i] != c1->n[i] ){
		  /*		  dprintf("ab: c2->n[%i]=%i, c1->n[%i]=%i\n", c2i, c2->n[c2i], i, c1->n[i]);*/
		  return 1;
		}
		if(prev_c2i!=c2i){
		  /*		  dprintf("ab: prev_c2i=%i, c2i=%i\n", prev_c2i, c2i); */
		  return 1;
		}
	 }
  }
  return 0;
}

/** assumes that dest is already allocated */
void copy_cluster(Clusters *dest, const Clusters *src){
  int i, j;

  dest->K = src->K;
  for( i=0; i<src->K; i++ ){
	 dest->n[i] = src->n[i];
	 for( j=0; j<src->n[i]; j++ ){
		dest->clust[i][j] = src->clust[i][j];
	 }
  }
  
}

void free_cluster(Clusters *c){
  int i;
  free( c->n );
  for( i=0; i<c->K; i++)
	 free(c->clust[i]);
  free(c->clust);
  free(c);
}

Clusters* init_cluster(int K, int maxN){
  Clusters *c;
  int i;

  c = (Clusters*) malloc( sizeof(Clusters) );
  c->K=K;
  c->n = (int*) malloc( K*sizeof(int) );
  c->clust = (int**) malloc( K*sizeof(int*) );
  for( i=0; i<K; i++ ){
	 c->clust[i] = (int*) malloc( maxN*sizeof(int) );
	 c->n[i] = 0;
  }

  return c;
}

void print_cluster(Clusters *c){
  FILE *o;
  int i, j;
  o = stderr;
  fprintf(o, "Cluster with %i partitions\n", c->K);
  for( i=0; i<c->K; i++ ){
	 fprintf(o, " C[%i]=\{ ", i+1);
	 for( j=0; j<c->n[i]; j++ ){
		fprintf( o, "%i ", c->clust[i][j] );
	 }
	 fprintf( o, "}\n" );
  }
  fprintf( o, "\n" );
}


int      gap_get_K(const double *gapstat, int k);
double*  gap_get_gapstat(const double *Wk, const double **Wkref, int B, int k, double *gapstat);
double** gap_get_reference_distribution(const double **d, int N, int n, double **ref);

/** compute \f$ W_k = \sum_{r=1}^{k}\frac{1}{2n_r}D_r\f$
	 where
	 \f$
	 D_r = \sum_{i,i'\in C_r}d_{ii'}
	 \f$
	 for a given cluster assignment.
	 \param d - difference matrix
	 \param N - number of entries in d (NxN)
	 \param c - cluster-assigment (Clusters - struct)
 */
double   gap_get_within_scatter(const double **d, int N, const Clusters *c){
  double *sums;
  double W;
  int i, j, k;
  print_cluster(c);
  sums = (double*) malloc( c->K*sizeof(double) );
  W = 0;
  for( i=0; i<c->K; i++ ){
	 sums[i] = 0;
	 for( j=0; j<(c->n[i])-1; j++ ){
		for( k=j+1; k<c->n[i]; k++){	 
	/* 	  dprintf("sums[i]=%f\n", sums[i]); */
/* 		  dprintf("i,j,k=(%i,%i,%i), c->clust[i][j]=%i, c->clust[i][k]=%i, d[]=%f\n", */
/* 					 i,j,k, c->clust[i][j], c->clust[i][k], d[c->clust[i][j]][c->clust[i][k]]); */
		  sums[i] += d[c->clust[i][j]][c->clust[i][k]];
		}
	 }
	 sums[i] /= (double)2*c->n[i];
	 dprintf("sums[i]=%f\n", sums[i]);
	 W += sums[i];
  }

  free(sums);

  return W;
}
double*  gap_get_within_scatter_distribution(const double **d, int N, int k, const Clusters **c, double *Wk);
