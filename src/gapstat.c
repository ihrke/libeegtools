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
/** \experimental_do_not_document */

#include "gapstat.h"
#include "linalg.h"
#include "clustering.h"
#include <time.h>
#include "eeg.h"


#ifdef LIBEEGTOOLS_FIXME

/** initialize GapStatistic struct.
	 \param g - a pointer to allocated mem, or NULL -> own memory is allocated
	 \param K - max num clusters
	 \param B - num rep ref dist
 */
GapStatistic* gapstat_init( GapStatistic *g, int K, int B ){
  int i;

  if( g==ALLOC_IN_FCT ){
	 g = (GapStatistic*)malloc( sizeof(GapStatistic) );
  }
  g->K = K;
  g->B = B;
  g->khat = -1;
  dprintf("allocating memory\n");
  g->gapdistr = (double*)malloc( K*sizeof(double) );
  g->sk = (double*)malloc( K*sizeof(double) );
  g->Wk= (double*)malloc( K*sizeof(double) );
  g->Wkref= (double**) malloc( B*sizeof(double*) );
  for( i=0; i<B; i++ )
	 g->Wkref[i] = (double*) malloc( K*sizeof(double) );
  g->progress=NULL;

  return g;
}

void gapstat_free( GapStatistic *g ){
  int i;
  free( g->Wk );
  for( i=0; i<g->B; i++ )
	 free( g->Wkref[i] );
  free( g->Wkref );
  free( g->gapdistr );
  free( g->sk );
  free( g );
}

/** calculate the gap-statistic. the struct contains all important information.
	 \todo fix this, there appears to be an NaN error
	 \param gap
	 \param X is nxp where n is number of observations and p num features (e.g. 
	          trials x voltage
	 \param D is the distance matrix for X (nxn)
 */
void gapstat_calculate( GapStatistic *gap, double **X, int n, int p, 
								VectorDistanceFunction distfunction, const double** D ){
  int i, k, b;
  Clusters *C;
  double **Xr; /* reference distribution */
  double **Dr; /* distance in ref dist */
  double l_bar; /* (1/B sum_b( log(Wkref[b]) )) */

  /* init */
  Xr = (double**) malloc( n*sizeof( double* ) );
  Dr = (double**) malloc( n*sizeof( double* ) );
  for( i=0; i<n; i++ ){
	 Xr[i] = (double*) malloc( p*sizeof( double ) );
	 Dr[i] = (double*) malloc( n*sizeof( double ) );
  }

  if( gap->progress ){
	 gap->progress( PROGRESSBAR_INIT, (gap->B+1)*gap->K );
  }

  /* calculation */
  for( k=1; k<=gap->K; k++ ){ /* data within-scatter */
	 if( gap->progress )
		gap->progress( PROGRESSBAR_CONTINUE_LONG, k );
	 C = kmedoids_repeat( D, n, k, 50 );
	 //	 cluster_print( stderr, C );
	 gap->Wk[k-1] = get_within_scatter( D, n, C );
	 cluster_free( C );
  }

  for( b=0; b<gap->B; b++ ){ /* monte Carlo for ref-data*/
	 for( k=1; k<=gap->K; k++ ){
		if( gap->progress ){
		  gap->progress( PROGRESSBAR_CONTINUE_LONG, (b+1)*gap->K+k );
		}
		
		Xr = gap_get_reference_distribution_simple( (const double **) X, n, p, Xr );
		Dr = vectordist_distmatrix( distfunction, (const double**)Xr, n, p, Dr, NULL, NULL );
		C = kmedoids_repeat( (const double**)Dr, n, k, 50 );
		//		cluster_print(stderr, C );
		gap->Wkref[b][k-1] = get_within_scatter( (const double**)Dr, n, C );
		cluster_free( C );
	 }
  }
  
  for( k=1; k<=gap->K; k++ ){ /* gap distr */
	 gap->gapdistr[k-1]=0;
	 for( b=0; b<gap->B; b++ ){
		gap->gapdistr[k-1] += log( gap->Wkref[b][k-1] );
	 }
	 l_bar = ((gap->gapdistr[k-1])/(double)gap->B);
	 gap->gapdistr[k-1] = l_bar - log( gap->Wk[k-1] );

	 /* need std and sk */
	 gap->sk[k-1] = 0;
	 for( b=0; b<gap->B; b++ ){
		gap->sk[k-1] += SQR( log( gap->Wkref[b][k-1] ) - l_bar );
	 }
	 gap->sk[k-1] = sqrt( gap->sk[k-1]/(double)gap->B );
	 gap->sk[k-1] = gap->sk[k-1] * sqrt( 1.0 + 1.0/(double)gap->B );
  }
  
  /* find best k */
  for( k=0; k<=gap->K-1; k++ ){
	 if( gap->gapdistr[k] >= ((gap->gapdistr[k+1])-(gap->sk[k+1])) ){
		gap->khat = k+1;
		break;
	 }
  }
  if( gap->progress )
	 gap->progress( PROGRESSBAR_FINISH, 0 );

  /* free */
  for( i=0; i<n; i++ ){
	 free(Xr[i]);
	 free(Dr[i]);
  }
  free( Xr );
  free( Dr );
}

void gapstat_print( FILE *out, GapStatistic *g ){
  int i, j;
  double mean; 

  fprintf(stderr, "GapStatistic:\n"
			 " K=%i\n"
			 " B=%i\n"
			 " khat=%i\n", g->K, g->B, g->khat );
  fprintf(stderr, " Gap       = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.2f) ", i+1, g->gapdistr[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " sk        = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.4f) ", i+1, g->sk[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " Wk        = [ ");
  for(i=0; i<g->K; i++ )
	 fprintf(stderr, "(%i, %.4f) ", i+1, g->Wk[i]);
  fprintf(stderr, "]\n");

  fprintf(stderr, " <Wkref>_B = [ ");
  for(i=0; i<g->K; i++ ){
	 mean = 0;
	 for(j=0; j<g->B; j++){
		mean += g->Wkref[j][i];
	 }
	 mean /= (double)g->B;
	 fprintf(stderr, "(%i, %.4f) ", i+1, mean);
  }
  fprintf(stderr, "]\n");

}

/** Compute a reference distribution for data X.
	 \param X original data
	 \param n,p dimensions of X
	 \param Xr n x p matrix with reference data
 */
double** gap_get_reference_distribution_simple( const double **X, int n, int p, double **Xr ){
  int i, j;
  double minf, maxf;

  for( i=0; i<p; i++ ){ /* features */
	 minf=DBL_MAX;
	 maxf=DBL_MIN;
	 for( j=0; j<n; j++ ){ 		  /* find min/max for feature i */
		if( X[j][i]<minf )
		  minf = X[j][i];
		if( X[j][i]>maxf )
		  maxf = X[j][i];
	 }
	 
	 for( j=0; j<n; j++ ){ /* uniform random values in [minf,maxf] */
		Xr[j][i] = (drand48()*maxf)+minf;
		if( isnan(Xr[j][i]) ){
		  errprintf(" Xr[%i][%i] is nan\n", j, i);
		}
	 }
  }

  return Xr;
}

/** Compute a reference distribution for data X. 
	 CAUTION, not tested.
	 \param X original data
	 \param n,p dimensions of X
	 \param Xr n x p matrix with reference data
	 \todo does not yet work for n<p, since libgsl does not support this for SVD
 */
double** gap_get_reference_distribution_svd( const double **X, int n, int p, double **Xr ){
  int i, j;
  double minf, maxf;
  gsl_matrix *gX, *gXd, *V;
  gsl_vector *work, *S;

  /* convert to gsl */
  gX = gsl_matrix_alloc( n, p );
  gXd= gsl_matrix_alloc( n, p );
  V  = gsl_matrix_alloc( n, p );
  S  = gsl_vector_alloc( p );
  work=gsl_vector_alloc( p );

  for( i=0; i<n; i++ ){
	 for( j=0; j<p; j++ ){
		gsl_matrix_set( gX, i, j, X[i][j] );
	 }
  }
  gsl_matrix_memcpy( gXd, gX );

  gsl_linalg_SV_decomp( gXd, V, S, work);

  gsl_blas_dgemm(CblasNoTrans, /* matrix mult */
					  CblasNoTrans,
					  1.0, gX, V, 0.0, gXd);

  for( i=0; i<p; i++ ){ /* features */
	 minf=DBL_MAX;
	 maxf=DBL_MIN;
	 for( j=0; j<n; j++ ){ 		  /* find min/max for feature i */
		if( gsl_matrix_get( gXd, j, i )<minf )
		  minf = gsl_matrix_get( gXd, j, i );
		if( gsl_matrix_get( gXd, j, i )>maxf )
		  maxf = gsl_matrix_get( gXd, j, i );
	 }
	 
	 for( j=0; j<n; j++ ){ /* uniform random values in [minf,maxf] */
		gsl_matrix_set( gX, j, i, ((drand48()*maxf)+minf) );
	 }
  }

  gsl_matrix_transpose_memcpy( gXd, V );
  gsl_blas_dgemm(CblasNoTrans,
					  CblasNoTrans,
					  1.0, gX, gXd, 0.0, V);

  for( i=0; i<n; i++ ){
	 for( j=0; j<p; j++ ){
		Xr[j][i] = gsl_matrix_get( V, i, j );
	 }
  }

  gsl_matrix_free( gX );
  gsl_matrix_free( gXd );
  gsl_matrix_free( V );
  gsl_vector_free( work );
  gsl_vector_free( S );

  return Xr;
}

#endif
