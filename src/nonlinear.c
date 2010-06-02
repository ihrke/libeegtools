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
double make_cond_entropy(long t, long *h1, long *h11, long **h2, long *array, 
								 int partitions, int length)
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

/** \brief simple nonlinear predicition.

	 Using the Simple-Prediction algorithm from 

	 Kantz & Schreiber 1997 Nonlinear Time-Series Analysis. Cambridge University Press

	 this function calculates the current+npredict's sample of the time series from which
	 sample has been taken based on the time-series in p. It's basically an average
	 over all points that are npredict time-steps after those points in p that are closer to 
	 sample than epsilon.

	 It's also described here:

	 R. Hegger, H. Kantz, and T. Schreiber, Practical
	 implementation of nonlinear time series methods: The TISEAN
	 package, CHAOS 9, 413 (1999)

	 \param p prediction time-series
	 \param sample predict time-points after this sample (must be of same dimension as p->m)
	 \param npredict predict n time-steps
	 \param epsilon epsilon-ball placed around all points in p to look for sample
	 \return the predicted sample in the time-series
 */
double      tdelay_predict_simple( TimeDelayReconstruction *p, double *sample, int npredict, double epsilon ){
  double predict=0;
  int nfound;
  int i;
  double *s; /* compare to sample */
  double dist; 

  s = dblp_init( NULL, p->m, 0.0 );
  for( nfound=0; nfound==0; epsilon*=1.2 ){ /* increase epsilon if no matching point is found */
	 dprintf("epsilon=%f\n", epsilon );
	 for( i=0; i<p->xn-npredict; i++ ){
		tdelay_index_i( p, i, s );
		dist = vectordist_euclidean( s, sample, p->m, NULL );
		if( dist<epsilon ){
		  nfound++;
		  predict+=p->x[i+npredict];
		}
	 }
  }
  dprintf("nfound after alg=%i\n", nfound);
  free( s );

  return predict/(double)nfound;
}

/** Calculates the root mean square prediction error averaged over all
	 samples in y. I.e. the function predicts the npredict'th following sample
	 for each point in y and sums up the RMSE for all these predictions. The result
	 is divided by yn.
	 \f[
	 result = \frac{1}{yn}\sum_{i=1}^{i=yn} RMSE( predict(y,npredict), y[i+npredict] )
	 \f] 
	 with 
	 \f$ 
    RMSE(r,d) = \sqrt{ \frac{1}{\#r}  \sum{ (r-d)^2)}} 
    \f$ 
	 For prediction, the function uses tdelay_predict_simple().

	 \param reference is the reference signal used for prediction 
	 \param y signal to predict for
	 \param yn length of y
	 \param npredict prediction for npredict timesteps ahead
	 \param epsilon defines the initial neighbourhood
	 \return root mean square prediction error cumulated (see above)
 */
double tdelay_simple_nonlinear_prediction_error( TimeDelayReconstruction *reference, double *y, int yn, 
																  int npredict, double epsilon ){
  int i;
  double crmse=0.0;
  TimeDelayReconstruction *Y;
  double *ysample;
  double pred; 
  Y = tdelay_init( reference->m, reference->tau, y, yn );
  ysample = dblp_init( NULL, reference->m, 0.0 );

  for( i=0; i<yn-npredict; i++ ){
	 tdelay_index_i( Y, i, ysample );
	 pred = tdelay_predict_simple( reference, ysample, npredict, epsilon );
	 crmse += rmse( &pred, &(y[i+npredict]), 1 );
  }

  free( ysample );
  tdelay_free( Y );
  return crmse/(double)(yn-npredict);
}

/** Calculate a trial x trial matrix M_ij that contains the prediction error when predicting
	 time-series points in trial j using trial i as the reference.

	 The function uses tdelay_simple_nonlinear_prediction_error() for prediction.

	 \param eeg
	 \param embedding_dim embedding dimensionm for phase space reconstruction
	 \param time_lag time-lag in sampling units used for phase space reconstruction
	 \param npredict prediction for npredict timesteps ahead
	 \param epsilon defines the initial neighbourhood (e.g. 1/4 of the variance in the data points)
	 \param allocated memory of eeg->ntrials x eeg->ntrials or ALLOC_IN_FCT
	 \param optargs may contain:
	 - "channel=int" calculate the trial x trial matrix for this channel, default=0
	 - "progress=void*" progressbar callback-function; default=NULL
	 \return the trial x trial matrix contain the prediction errors; NULL if an error occured
 */
double** eeg_nonlinear_prediction_error( const EEG *eeg, int embedding_dim, int time_lag,
													  int npredict, double epsilon, 
													  double** output, OptArgList *optargs ){
  int channel=0;
  double x;
  int i, j;
  void *ptr;
  ProgressBarFunction progress=NULL;

  /* optarg parsing */
  if( optarglist_has_key( optargs, "channel" ) ){
	 x = optarglist_scalar_by_key( optargs, "channel" );
	 if( !isnan( x ) ) channel=(int)x;
  }
  if( eeg->nbchan<=channel ){
	 errprintf( "Channel '%i' is not in EEG-set which holds only '%i' channels\n", 
					channel, eeg->nbchan );
	 return NULL;
  }
  if( optarglist_has_key( optargs, "progress" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "progress" );
	 if( ptr ) progress = (ProgressBarFunction)ptr;
  }

  /* output allocation */
  if( !output ){
	 output = dblpp_init( eeg->ntrials, eeg->ntrials );
  }

  /* calculation */
  TimeDelayReconstruction *reference;
  reference = tdelay_init( embedding_dim, time_lag, NULL, eeg->n );  

  if( progress ){
	 progress( PROGRESSBAR_INIT, eeg->ntrials );
  }
  for( i=0; i<eeg->ntrials; i++ ){
	 /* fill phase-space with data */
	 reference->x = eeg->data[channel][i];	
	 if( progress ) progress( PROGRESSBAR_CONTINUE_LONG, i );
	 for( j=0; j<eeg->ntrials; j++ ){ /* not symmetric */
		output[i][j] = tdelay_simple_nonlinear_prediction_error( reference, eeg->data[channel][j], eeg->n,
																					 npredict, epsilon );
		if( progress )	progress( PROGRESSBAR_CONTINUE_SHORT, j );
	 }
  }

  if( progress ) progress( PROGRESSBAR_FINISH, 0 );
  return output;
}


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
int         tdelay_estimate_timelag_mutual( TimeDelayReconstruction *p, long partitions, long corrlength, double *mutual ){
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
	 dprintf("mutual[%li,%li] = %f,%f\n", i,i+1, mutual[i], mutual[i+1]); 
	 if( mutual[i+1]>mutual[i] ){
		tau = i;
		break;
	 }
  }

  //  dblp_min( mutual, corrlength+1, &tau );

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
int         tdelay_estimate_timelag_autocorr( TimeDelayReconstruction *p ){
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
	 \todo implement this (it's easy!)
	 \return 
 */
int         tdelay_estimate_dimension_fnn( TimeDelayReconstruction *p, double Rtol, double Atol, int num_dim  ){
  errprintf("Not implemented ?!?!?\n");
  return -1;
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
double  tdelay_fnn_ratio( TimeDelayReconstruction *p, double Rtol, double Atol ){
  int i, j;
  int i_nn=0;							  /* index nearest_neighbour */
  double d_nn;						  /* dist nn */
  int nfnn;							  /* number of false nearest neighbours */
  double ratio_fnn;
  double **d;
  double **X;
  TimeDelayReconstruction pnext;
  double Ra;
  double crit1,crit2;


  if( p->m <= 0){
	 errprintf("p->m=%i is <=0", p->m);
	 return -1;
  }

  Ra = tdelay_attractor_size( p );
  dprintf( "Ra=%f\n", Ra);
  /* create distance matrix between points in phase-space */
  X = dblpp_init( p->xn, p->m );
  for( i=0; i<p->xn; i++ ){
	 tdelay_index_i( p, i, X[i] );
  }
  d = vectordist_distmatrix( vectordist_euclidean, (const double**)X, p->xn, p->m, ALLOC_IN_FCT, NULL, NULL );
  dprintf( "dmin,dmax=%f,%f\n", 
			  dblpp_min((const double**)d,p->xn,p->xn,NULL,NULL), 
			  dblpp_max((const double**)d,p->xn,p->xn,NULL,NULL) );
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
	 crit1 = (ABS( tdelay_index_ij(&pnext, i, p->m) - tdelay_index_ij(&pnext, i_nn, p->m) )/d_nn);
	 dprintf( "crit1=%f, Rtol=%f\n", crit1, Rtol );
	 if( crit1 > Rtol ){
		dprintf(" Crit 1 for %i failed\n", i);
		nfnn++;
		continue;
	 }

	 /* crit (2) */
	 crit2 = sqrt( SQR( d_nn ) + 
						SQR( tdelay_index_ij(&pnext, i, p->m) - 
							  tdelay_index_ij(&pnext, i_nn, p->m) ) )/Ra;
	 dprintf( "crit2=%f, Rtol=%f\n", crit2, Rtol );

	 if( crit2 > Atol ){
		dprintf(" Crit 2 for %i failed\n", i);
		nfnn++;
		continue;
	 }
  }
  ratio_fnn = nfnn/(double)p->xn*100.0;
  dprintf(" nfnn=%i, ratio=%f\n", nfnn, ratio_fnn );

  dblpp_free( X, p->xn );
  dblpp_free( d, p->xn );

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
double      tdelay_attractor_size( TimeDelayReconstruction *p ){
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
double tdelay_index_ij( TimeDelayReconstruction *p, int i, int j ){
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
void tdelay_index_i( TimeDelayReconstruction *p, int i, double *x){
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
void tdelay_index_j( TimeDelayReconstruction *p, int j, double *x){
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

TimeDelayReconstruction* tdelay_init( int m, int tau, double *x, int n ){
  TimeDelayReconstruction *p;
  p = (TimeDelayReconstruction*)malloc(sizeof(TimeDelayReconstruction));
  p->m = m;
  p->tau = tau;
  p->x = x;
  p->xn= n;
  return p;
}

void tdelay_free( TimeDelayReconstruction *p ){
  free(p);
}

void tdelay_print( FILE *out, TimeDelayReconstruction *p){
  fprintf( out, "TimeDelayReconstruction '%p':\n"
			  " m   = %i\n"
			  " tau = %i\n"
			  " x   = %p\n"
			  " xn  = %i\n", p, p->m, p->tau, p->x, p->xn );
}
