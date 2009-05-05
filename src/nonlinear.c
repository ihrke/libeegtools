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

#include "nonlinear.h"

/**\cond PRIVATE */
double make_cond_entropy(long t, long *h1, long *h11, long **h2, long *array, int partitions, int length)
{
  long i,j,hi,hii,count=0;
  double hpi,hpj,pij,cond_ent=0.0,norm;

  for (i=0;i<partitions;i++) {
    h1[i]=h11[i]=0;
    for (j=0;j<partitions;j++){
      h2[i][j]=0;
	 }
  }

  for (i=0;i<length;i++){
    if (i >= t) {
      hii=array[i];
      hi=array[i-t];
      h1[hi]++;
      h11[hii]++;
      h2[hi][hii]++;
      count++;
    }
  }

  norm=1.0/(double)count;
  cond_ent=0.0;

  for (i=0;i<partitions;i++) {
    hpi=(double)(h1[i])*norm;
    if (hpi > 0.0) {
      for (j=0;j<partitions;j++) {
		  hpj=(double)(h11[j])*norm;
		  if (hpj > 0.0) {
			 pij=(double)h2[i][j]*norm;
			 if (pij > 0.0){
				cond_ent += pij*log(pij/hpj/hpi);
			 }
		  }
      }
    }
  }

  return cond_ent;
}
/** \endcond */


/** estimate time-lag based on the first local minimum of the 
	 mutual information.
	 
	 This function is "inspired" (stolen and adapted) from TISEAN 3.1,
	 http://www.mpipks-dresden.mpg.de/~tisean

	 Ref: 
	 R. Hegger, H. Kantz, and T. Schreiber, Practical
	 implementation of nonlinear time series methods: The TISEAN
	 package, CHAOS 9, 413 (1999)

	 \param p the phase-space rep of the signal
	 \param partitions number of partitions to use (if <0, 16 is used as default)
	 \param corrlength maximum corrlength (if <0, 20 is used as default)
	 \param mutual if != NULL, the array is filled with the mutual information for each time-lag (corrlength-long)
	 \return the minimum of mutual
*/
int         phspace_estimate_timelag_mutual( PhaseSpace *p, long partitions, long corrlength, double *mutual ){
  long *array,*h1,*h11,**h2;
  long tau,i,length;
  double *series,min,interval,shannon;
  int free_mutual=0;

  if(partitions<=0)
	 partitions=16;
  if(corrlength<=0)
	 corrlength=20;
  
  series = p->x;
  length = p->xn;

  h1 =(long *)malloc(sizeof(long) *partitions);
  h11=(long *)malloc(sizeof(long) *partitions);
  h2 =(long**)malloc(sizeof(long*)*partitions);
  for (i=0;i<partitions;i++) 
    h2[i]=(long *)malloc(sizeof(long)*partitions);
  array=(long *)malloc(sizeof(long)*length);

  /* rescaling data */
  min=interval=series[0];
  for (i=1;i<length;i++) {
    if (series[i] < min) 
		min=series[i];
    if (series[i] > interval) 
		interval=series[i];
  }
  if (interval != 0.0) {
    for (i=0;i<length;i++){
      series[i]=(series[i]- min)/ (double)interval;
	 }
  }


  for (i=0;i<length;i++){
    if (series[i] < 1.0){
      array[i]=(long)(series[i]*(double)partitions);
	 } else {
      array[i]=partitions-1;
	 }
  }
  
  shannon=make_cond_entropy(0, h1, h11, h2, array, partitions, length );
  if (corrlength >= length){
    corrlength=length-1;
  }
  
  if( !mutual ){
	 free_mutual=1;
	 mutual = (double*) malloc( (corrlength+1)*sizeof(double) );
  }

  mutual[0] = shannon;
  for (tau=1;tau<=corrlength;tau++) {
	 mutual[tau] = make_cond_entropy(tau, h1, h11, h2, array, partitions, length );
  }

  for( i=0; i<corrlength; i++ ) {
	 dprintf("mutual[%li,%i] = %f,%f\n", i,i+1, mutual[i], mutual[i+1]); 
	 if( mutual[i+1]>mutual[i] ){
		tau = i;
		break;
	 }
  }

  //  vector_min( mutual, corrlength+1, &tau );

  free( h1 );
  free( h11 );
  for( i=0; i<partitions; i++ )
	 free( h2[i] );
  free( h2 );
  free( array );
  if( free_mutual )
	 free( mutual );

  return tau;
}


/** estimate time-lag based on the first zero of the autorcorrelation function.
	 
*/
int         phspace_estimate_timelag_autocorr( PhaseSpace *p ){
  int tau=-1;
  int i,j;
  double Rxx;
  for( i=0; i<p->xn; i++ ){
	 Rxx = 0.0;
	 for( j=0; j<p->xn; j++ ){
		Rxx += p->x[j]*p->x[ABS(j-i)];
	 }
	 dprintf("Rxx(%i)=%f\n", i, Rxx);
	 if( Rxx<=0 ){
		tau=i;
		break;
	 }
  }
  if( tau<0 ){
	 warnprintf("Sorry, have not been able to compute tau properly."
					"There was no zero crossing in the autocorrelation fucntion\n");
  }
  return tau;
}

/** Estimate the embedding dimension for the phase-space from the
	 number of false-nearest-neighbours as proposed in

	 KENNEL et al. DETERMINING EMBEDDING DIMENSION FOR PHASE-SPACE
	 RECONSTRUCTION USING A GEOMETRICAL CONSTRUCTION. Physical Review A
	 (1992) vol. 45 (6) pp. 3403-3411

	 \param 
	 \return 
 */
int         phspace_estimate_dimension_fnn( PhaseSpace *p, double Rtol, double Atol, int num_dim  ){
  
}

/** Estimate number of false-nearest-neighbours as proposed in

	 KENNEL et al. DETERMINING EMBEDDING DIMENSION FOR PHASE-SPACE
	 RECONSTRUCTION USING A GEOMETRICAL CONSTRUCTION. Physical Review A
	 (1992) vol. 45 (6) pp. 3403-3411


	 by applying two criteria.
	 
	 (1) if 
	 \f[
	 \frac{ |x_{i+m\tau} - x^{r}_{i+m\tau} |} 
	 {} > R_{tol}
	 \f]
	 then state i is

 */
double  phspace_fnn_ratio( PhaseSpace *p, double Rtol, double Atol ){
  int i, j;
  int i_nn;							  /* index nearest_neighbour */
  double d_nn;						  /* dist nn */
  int nfnn;							  /* number of false nearest neighbours */
  double ratio_fnn;
  double **d;
  double **X;
  PhaseSpace pnext;
  double Ra;
  double crit1,crit2;


  if( p->m <= 0){
	 errprintf("p->m=%i is <=0", p->m);
	 return -1;
  }

  Ra = phspace_attractor_size( p );
  dprintf( "Ra=%f\n", Ra);
  /* create distance matrix between points in phase-space */
  X = matrix_init( p->xn, p->m );
  for( i=0; i<p->xn; i++ ){
	 phspace_index_i( p, i, X[i] );
  }
  d = vectordist_distmatrix( vectordist_euclidean, X, p->xn, p->m, ALLOC_IN_FCT, NULL, NULL );
  dprintf( "dmin,dmax=%f,%f\n", matrix_min(d,p->xn,p->xn,NULL,NULL), matrix_max(d,p->xn,p->xn,NULL,NULL) );
  /* for each column, find minimum entry in d -> nearest neighbour;
	  for this neighbour, determine if it is a false one, or not 
  */
  nfnn = 0;
  /* go up one dimension */
  pnext.m = p->m+1;
  pnext.tau=p->tau;
  pnext.x = p->x;
  pnext.xn = p->xn;

  for( i=0; i<p->xn; i++ ){
	 d_nn = DBL_MAX;
	 for( j=0; j<p->xn; j++ ){
		if( i==j ) continue;
		if( d[i][j]<d_nn ){
		  d_nn = d[i][j];
		  i_nn = j;
		}
	 }
	 dprintf( "d_nn = %f, i_nn = %i\n", d_nn, i_nn );
	 
	 /* crit (1) */
	 crit1 = (ABS( phspace_index_ij(&pnext, i, p->m) - phspace_index_ij(&pnext, i_nn, p->m) )/d_nn);
	 dprintf( "crit1=%f, Rtol=%f\n", crit1, Rtol );
	 if( crit1 > Rtol ){
		dprintf(" Crit 1 for %i failed\n", i);
		nfnn++;
		continue;
	 }

	 /* crit (2) */
	 dprintf( "crit2=%f, Rtol=%f\n", crit2, Rtol );
	 crit2 = sqrt( SQR( d_nn ) + 
						SQR( phspace_index_ij(&pnext, i, p->m) - 
							  phspace_index_ij(&pnext, i_nn, p->m) ) )/Ra;
	 if( crit2 > Atol ){
		dprintf(" Crit 2 for %i failed\n", i);
		nfnn++;
		continue;
	 }
  }
  ratio_fnn = nfnn/(double)p->xn*100.0;
  dprintf(" nfnn=%i, ratio=%f\n", nfnn, ratio_fnn );

  matrix_free( X, p->xn );
  matrix_free( d, p->xn );

  return ratio_fnn;
}

/** estimate size of the attractor from a scalar time-series using 
	 \f[
	 R = \sqrt{ \frac{1}{N}\sum_{i=1}^{N} (x_i - \langle x\rangle_i)^2 }
	 \f]
	 where
	 \f[
	 \langle x\rangle_i = \frac{1}{N}\sum_{i=1}^{N} x_i
	 \f]
 */
double      phspace_attractor_size( PhaseSpace *p ){
  double R=0.0;
  double mean=0.0;
  int i;
  
  for( i=0; i<p->xn; i++ ){
	 mean += p->x[i];
  }
  mean /= (double)p->xn;

  for( i=0; i<p->xn; i++ ){
	 R += SQR ( p->x[i] - mean );
  }
  R /= (double)p->xn;

  return sqrt( R );
}

/** Get an element from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param i,j as in the formula above
 */
double phspace_index_ij( PhaseSpace *p, int i, int j ){
  int    idx;
  
  if( i<0 || i>p->xn || j<0 || j>=p->m ){
	 errprintf("i=%i, p->xn=%i\n, j=%i, p->m=%i\n", i, p->xn, j, p->m);
	 return 0;
  }
  idx = i - (((p->m) - (j+1)) * (p->tau));
  if( idx<0 ){
	 idx = (p->xn)+idx;
  }
  /* dprintf("idx=%i\n", idx); */
  return(p->x[idx]);
}

/** Get an element from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param i as in the formula above
	 \param x output (m long vector containing the final vector)
 */
void phspace_index_i( PhaseSpace *p, int i, double *x){
  int    idx, j;

  if( i<0 || i>p->xn ){
	 errprintf("i=%i, p->xn=%i\n", i, p->xn);
	 return;
  }
  for( j=0; j<p->m; j++ ){
	 idx = i - (((p->m) - (j+1)) * (p->tau));
	 if( idx<0 ){
		idx = (p->xn)+idx;
	 }
	 x[j] = p->x[idx];
  }
}
/** Get a dimension from phase-space reconstruction:
	 \f[
	 \vec{x}_i = \sum_{j=1}^{m}s_{i+(j-1)\tau}\vec{e}_j
	 \f]
	 \param p is the phase-space rep. of the data
	 \param j as in the formula above
	 \param x output (n long vector containing the final vector)
 */
void phspace_index_j( PhaseSpace *p, int j, double *x){
  int    idx, i;

  if( j<0 || j>=p->m ){
	 errprintf("j=%i, p->m=%i\n", j, p->m);
	 return;
  }
  for( i=0; i<p->xn; i++ ){
	 idx = i - (((p->m) - (j+1)) * (p->tau));
	 if( idx<0 ){
		idx = (p->xn)+idx;
	 }
	 x[i] = p->x[idx];
  }
}

PhaseSpace* phspace_init( int m, int tau, double *x, int n ){
  PhaseSpace *p;
  p = (PhaseSpace*)malloc(sizeof(PhaseSpace));
  p->m = m;
  p->tau = tau;
  p->x = x;
  p->xn= n;
  return p;
}

void phspace_free( PhaseSpace *p ){
  free(p);
}

void phspace_print( FILE *out, PhaseSpace *p){
  fprintf( out, "PhaseSpace '%p':\n"
			  " m   = %i\n"
			  " tau = %i\n"
			  " x   = %p\n"
			  " xn  = %i\n", p, p->m, p->tau, p->x, p->xn );
}
