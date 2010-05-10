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
#include "optarg.h"
#include "array.h"
#include "linalg.h"
#include "eeg.h"

/** \cond PRIVATE

	 apply slope constraint by using cumulation along different
	 possible pathes. See Skaoe & Chiba, 1978.
*/

double _slope_constraint( SlopeConstraint constraint, Array *d, Array *D, int i, int j ){
  double r;

  if( (i==0 || j==0) ){
	 errprintf("i and j must be larger than zero '%i'\n", (int)constraint);
	 return NAN;
  }

  r = mat_IDX( D, i-1, j-1 ) + 2*mat_IDX( d, i, j );
  switch( constraint ){
  case SLOPE_CONSTRAINT_NONE:
	 r = MIN( r, mat_IDX( D, i, j-1) + mat_IDX( d, i, j ) );
	 r = MIN( r, mat_IDX( D, i-1, j) + mat_IDX( d, i, j ) );
	 break;
  case SLOPE_CONSTRAINT_LAX:
	 if( j>=3 )
		r = MIN( r, mat_IDX( D, i-1, j-3 ) + 2*mat_IDX( d, i, j-2 ) + mat_IDX( d, i, j-1 ) + mat_IDX( d, i, j ) );
	 if( j>=2 ){
		r = MIN( r, mat_IDX( D, i-1, j-2 ) + 2*mat_IDX( d, i, j-1 ) + mat_IDX( d, i, j ) );
	 }
	 if( i>=2 ){
		r = MIN( r, mat_IDX( D, i-2, j-1 ) + 2*mat_IDX( d, i-1, j ) + mat_IDX( d, i, j ) );
	 }
	 if( i>=3 ){
		r = MIN( r, mat_IDX( D, i-3, j-1 ) + 2*mat_IDX( d, i-2, j ) + mat_IDX( d, i-1, j ) + mat_IDX( d, i, j ) );
	 }
	 break;
  case SLOPE_CONSTRAINT_MEDIUM:
	 if( j>=2 ){
		r = MIN( r, mat_IDX( D, i-1, j-2 ) + 2*mat_IDX( d, i, j-1 ) + mat_IDX( d, i, j ) );
	 }
	 if( i>=2 ){
		r = MIN( r, mat_IDX( D, i-2, j-1 ) + 2*mat_IDX( d, i-1, j ) + mat_IDX( d, i, j ) );
	 }
	 break;
  case SLOPE_CONSTRAINT_SEVERE:
	 if( i>=2 && j>=3 ){
		r = MIN( r, mat_IDX( D, i-2, j-3 ) + 2*mat_IDX( d, i-1, j-2 ) + 2*mat_IDX( d, i, j-1 ) + mat_IDX( d, i, j ) );
	 }
	 if( j>=2 && i>=3 ){
		r = MIN( r, mat_IDX( D, i-3, j-2 ) + 2*mat_IDX( d, i-2, j-1 ) + 2*mat_IDX( d, i-1, j ) + mat_IDX( d, i, j ) );
	 }
	 break;
  default:
	 errprintf("Slope constraint not supported '%i'\n", (int)constraint);
	 r = NAN;
	 break;
  }
  return r;
}
/**\endcond*/

/** \brief cumulate a distance matrix d for Dynamic Time-Warping.

	 \f[
	 D_{jk} = \mathbf{d}_{jk}+\min{\{D_{j,k-1}, D_{j-1,k}, D_{j-1, k-1}\}}
	 \f]
	 The formula is modified, depending on the slope constraint.

	 \param mat distance matrix
	 \param alloc if TRUE, output matrix is freshly allocated; else mat is 
	        overwritten
	 \param optargs may contain:						 
	 - "slope_constraint=int" slope constraint, one of SLOPE_CONSTRAINT_*; 
	     default=SLOPE_CONSTRAINT_NONE
 */
Array*      matrix_dtw_cumulate ( Array *mat, bool alloc, OptArgList *optargs ){
  /* j = 0,...,M-1
	  k = 0,...,N-1
  */
  int j, k, M, N;
  double x;
  Array *D;
  SlopeConstraint constraint = SLOPE_CONSTRAINT_NONE;

  bool ismatrix;
  matrix_CHECK( ismatrix, mat );
  if( !ismatrix ) return NULL;

  if( optarglist_has_key( optargs, "slope_constraint" ) ){
	 x = optarglist_scalar_by_key( optargs, "slope_constraint" );
	 if( !isnan( x ) ) constraint=(SlopeConstraint)x;
  }
  dprintf("using slope constraint '%i'\n", constraint );

  D = array_copy( mat, TRUE );

  M = mat->size[0];
  N = mat->size[1];
   /* computing D_jk */
  for(j=0; j<M; j++){
    for(k=0; k<N; k++){
      if(k==0 && j==0) ;
      else if(k==0 && j>0){
		  mat_IDX( D, j,k) = mat_IDX( D, j-1, k ) + mat_IDX( mat, j, k );
		} else if(k>0 && j==0){
		  mat_IDX( D, j, k) = mat_IDX( D, j, k-1 ) + mat_IDX( mat, j, k );
		} else { /* j,k > 0 */
		  mat_IDX( D, j, k) = _slope_constraint( constraint, mat, D, j, k ); 
		}
    }
  }
  if( !alloc ){
	 memcpy( mat->data, D->data, D->nbytes );
	 array_free( D );
	 D = mat;
  }
	 
  return D;
}

/** \brief calculate the warping path.
	 
	 \param d is the cumulated distances matrix (usually output from matrix_dtw_cumulate())
	 \return the warp-Path (2D uint array, 2 x N)
 */
Array* matrix_dtw_backtrack ( const Array *d ){ 
  int i,j,k, M, N; /* j = 0,...,M-1
						  k = 0,...,N-1  */
  int idx;
  double left, down, downleft;

  bool isthing;
  matrix_CHECK( isthing, d );
  if( !isthing ) return NULL;

  M = d->size[0];
  N = d->size[1];

  Array *path,
	 *path2 = array_new2( UINT, 2, 2, M+N ); /* dummy path */

  /*------------ computation -------------------------------*/
  /* Backtracking */
  j=M-1; k=N-1;

  idx = 1;
  array_INDEX2( path2, uint, 0, 0 ) = j;
  array_INDEX2( path2, uint, 1, 0 ) = k;
  while( j>0 || k>0 ){ 
	 if( k==0 ){ /* base cases */
  		j--;
  	 } else if( j==0 ){
  		k--;
  	 } else { /* min( d[j-1][k], d[j-1][k-1], d[j][k-1] ) */
  		left     = mat_IDX( d, j-1, k );
  		down     = mat_IDX( d, j  , k-1 );
  		downleft = mat_IDX( d, j-1, k-1 );

  		(isnan( mat_IDX( d, j-1,k  ) ) )?(left=DBL_MAX):(left=mat_IDX( d, j-1, k ) );
  		(isnan( mat_IDX( d, j  ,k-1) ) )?(down=DBL_MAX):(down=mat_IDX( d, j  ,k-1) );
  		(isnan( mat_IDX( d, j-1,k-1) ) )?(downleft=DBL_MAX):(downleft=mat_IDX(d,j-1,k-1));

  		/* this order of comparisons was chosen to ensure, that diagonal
  			steps are preferred in case of equal distances
  		*/
  		if( downleft<=left ){	  /* downleft or down */
  		  if( downleft<=down ){
  			 k--; j--;			 	  /* downleft */
  		  } else {
  			 k--;						  /* down */
  		  }
  		} else {						  /* left or down */
  		  if( down<=left ){
  			 k--;						  /* down */
  		  } else {
  			 j--; 					  /* left */
  		  }
  		}


  	 }	/* if */

  	 array_INDEX2( path2, uint, 0, idx ) = j;
  	 array_INDEX2( path2, uint, 1, idx ) = k;
  	 idx++;
  }
  
  /* now the warppath has wrong order, reverse */ 
  path = array_new2( UINT, 2, 2, idx );

  for( i=0; i<idx; i++ ){
	 array_INDEX2( path, uint, 0, idx-i-1 )=array_INDEX2( path2, uint, 0, i );
	 array_INDEX2( path, uint, 1, idx-i-1 )=array_INDEX2( path2, uint, 1, i );
  }
  /*------------ /computation -------------------------------*/

  array_free( path2);
  return path;
}

/** Adjust time-markers according to
	 \f[
	 \hat{\tau} = \left(\frac{\tau_1}{\omega_1} + \frac{\tau_2}{\omega_2}\right)
	 \f]
	 where \f$ \omega_1,\omega_2\f$ are the weights.
	 \param m1,m2,nmarkers are the markers to be adjusted
	 \param outmarkers is memory or NULL to store the new adjusted markers,
	 \param weights are the weights according to the above equation
	 \return outmarkers
 */
unsigned int* warp_adjust_time_markers(const unsigned int *m1, const unsigned int *m2, 
											  int nmarkers, unsigned int *outmarkers,
											  const double weights[2] ){
  int i;
  if(!outmarkers){
	 outmarkers = (unsigned int*)malloc(nmarkers*sizeof(unsigned int) );
  }
  for( i=0; i<nmarkers; i++ ){
	 outmarkers[i] = m1[i]*weights[0] + m2[i]*weights[1];
  }
  
  return outmarkers;
}

/** Add magnitude of signals according to warppath P.
	 Data is resampled to have length (J+K)/2+1 using cspline interpolation.
 * \param path - contains warppath
 * \param weights - weights in average, for using it with hierarchical averaging
 * \param avg - pointer to caller-allocated memory (length (K+J)/2+1); 
 *              if NULL is given, memory is allocated by the function.
 */
double* warp_add_signals_by_path(const double *s1, int n1, 
											const double *s2, int n2, 
											const WarpPath *P, double *avg, 
											const double weights[2]){
  double *x, *y, *xp;
  int newn, i;
  
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


  for(i=0; i<P->n; i++){
	 x[i] = (double)((P->t1[i])*weights[0]+(P->t2[i])*weights[1]); /* average latencies */
	 y[i] = weights[0]*(s1[P->t1[i]]) + weights[1]*(s2[P->t2[i]]); /* average magnitude */
	 //	 dprintf( "x[%i]=(%f,%f), (%i, %i)\n", i, x[i], y[i], (P->t1[i]), (P->t2[i]) );
	 if( i>0 && x[i]<=x[i-1] ){
		warnprintf("x not monotonic at P=%p x[%i]=(%f,%f), (%i, %i)\n", P, i, x[i], y[i], (P->t1[i]), (P->t2[i]) );
	 }
  }

  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, P->n);
  gsl_spline_init (spline, x, y, P->n);

  for( i=0; i<newn; i++ ){
	 //	 dprintf("xp[%i]=%f\n", i, xp[i] ); 
	 avg[i] = gsl_spline_eval( spline, xp[i], acc );
	 //	 dprintf(" avg[%i] = (%f, %f)\n", i,xp[i], avg[i] );
  }

  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free( x );
  free( y );
  free( xp );

  return avg;
}


/** Hierarchically warp/average the trials within the EEG-struct.

	 \param eeg_in input data
	 \param distmatrix Distance matrix between the trials in eeg_in (must be eeg_in->ntrials x eeg_in->ntrials matrix)
	 \param out memory for the final average or NULL (own memory is allocated)
	 \param optargs may contain:
	 - <b>Directly used by this function:</b>
	 - <tt>regularize=void*</tt> regularization function, default=\c NULL
	 - <tt>linkage=void*</tt> rule for building the dendrogram, default=\c dgram_dist_completelinkage
	 - <tt>pointdistance=void*</tt> pointwise distance function for two time-series, default=\c signaldist_euclidean
	 - <tt>pass_trial_markers_in_optargs=int</tt> if true, replace (or create) the "markers=int**" and "nmarkers=int" fields in optargs to match the currently used trials (in each iteration, two trials are combined and their markers are sometimes required by the distance functions);  the function uses an own copy of the optargs handed in to the function, such that the original optargs is not modified!, default=\c TRUE
	 - <tt>progress=void*</tt> progress-bar function, called every now and then.
	 The function is called with PROGRESSBAR_INIT and 
	 the max number of calls once. Then it is called
	 for a major step with PROGRESSBAR_CONTINUE_LONG
	 and the num of the step as second arg. For small
	 steps it is called with PROGRESSBAR_CONTINUE_SHORT.
	 If NULL, nothing is done., default=\c NULL
	 - <tt>use_eeg_in=int</tt> if provided, use the memory in eeg_in. default is to allocate new memory
	 - <b> Depending on which metric/regularization etc you chose, you might want to provide further optional arguments </b> 
	 - The optargs list is passed to the regularization, linkage, dtw_cumulate_matrix and metric functions
	 \return pointer to final average (EEG-struct with one trial which holds the average)
 */  
EEG* eeg_dtw_hierarchical( EEG *eeg_in, const double **distmatrix,
									EEG *out, OptArgList *optargs ){
#if 0
  Dendrogram *T, *Tsub; 
  EEG *eeg;
  WarpPath *P;
  int nmarkers, n, N;
  int i, idx1, idx2, c,
	 trials_left;
  double *s1, *s2, *s1s2;
  unsigned int    *s1s2_markers;
  unsigned int    **markers=NULL;
  double **G, **d;
  double weights[2]={1.0,1.0}; 			  /* this is for recursive averaging */
  int    *indices;				  /* n-array containing number averagings for each trial */
  void *ptr;
  double x;

  /* set defaults */
  RegularizationFunction regularize=NULL;
  double sigma_max=0.2;
  LinkageFunction linkage=dgram_dist_completelinkage;
  PointwiseDistanceFunction pointdistance=signaldist_euclidean;
  ProgressBarFunction progress=NULL;
  int use_eeg_in=FALSE;
  int pass_trial_markers_in_optargs = TRUE;
  OptArgList *new_optargs, *tmp_arglist;

  /* override params */
  if( optarglist_has_key( optargs, "regularize" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "regularize" );
	 if( ptr ) regularize = (RegularizationFunction)ptr;
  }
  if( optarglist_has_key( optargs, "linkage" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "linkage" );
	 if( ptr ) linkage = (LinkageFunction)ptr;
  }
  if( optarglist_has_key( optargs, "pointdistance" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "pointdistance" );
	 if( ptr ) pointdistance = (PointwiseDistanceFunction)ptr;
  }
  if( optarglist_has_key( optargs, "progress" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "progress" );
	 if( ptr ) progress = (ProgressBarFunction)ptr;
  }
  if( optarglist_has_key( optargs, "use_eeg_in" ) ){
	 x = optarglist_scalar_by_key( optargs, "use_eeg_in" );
	 if( !isnan( x ) ) use_eeg_in=(int)x;
  }
  if( optarglist_has_key( optargs, "sigma_max" ) ){
	 x = optarglist_scalar_by_key( optargs, "sigma_max" );
	 if( !isnan( x ) ) sigma_max=x;
  }
  if( optarglist_has_key( optargs, "pass_trial_markers_in_optargs" ) ){
	 x = optarglist_scalar_by_key( optargs, "pass_trial_markers_in_optargs" );
	 if( !isnan( x ) ) pass_trial_markers_in_optargs=(int)x;
  }

  /* what about the EEG memory ?*/
  if( !use_eeg_in ){
	 dprintf("Cloning EEGData\n");
	 eeg = eeg_clone( eeg_in, EEG_CLONE_ALL );
	 dprintf("...done\n");
  } else {
	 warnprintf("using eeg_in and deleting it because you want it that way!\n");
	 eeg = eeg_in; 
  }

  /* shortcuts */
  N = eeg->ntrials;
  n = eeg->n;
  nmarkers = eeg->nmarkers[0]+2; /* trivial markers */

  /* what about optional arguments? */
  OptArg *tmp_arg=NULL;
  if( pass_trial_markers_in_optargs ){ /* create own optarg memory*/
	 markers = (unsigned int**) malloc( 2*sizeof(unsigned int*));
	 markers[0] = (unsigned int*) malloc( nmarkers*sizeof(unsigned int) ); 
	 markers[1] = (unsigned int*) malloc( nmarkers*sizeof(unsigned int) );
	 markers[0][0] = 0; /* trivial markers */
	 markers[0][nmarkers-1] = n-1;
	 markers[1][0] = 0;
	 markers[1][nmarkers-1] = n-1;

	 new_optargs = (OptArgList*) malloc( optargs->nargs*sizeof(OptArg) );
	 memcpy( new_optargs, optargs,  optargs->nargs*sizeof(OptArg) );
	 if( !optarglist_has_key( new_optargs, "nmarkers" ) ){
		tmp_arg = optarg_scalar( "nmarkers", (double) nmarkers );
	 }
	 tmp_arglist = optarglist_append_arg( new_optargs, tmp_arg );
	 optarglist_free( new_optargs );
	 new_optargs = tmp_arglist;
	 free( tmp_arg );
	 if( !optarglist_has_key( new_optargs, "markers" ) ){
		tmp_arg = optarg_ptr( "markers", markers );
	 }
	 tmp_arglist = optarglist_append_arg( new_optargs, tmp_arg );
	 optarglist_free( new_optargs );
	 new_optargs = tmp_arglist;
	 free( tmp_arg );
  } else {
	 new_optargs = optargs;
  }

#ifdef DEBUG
  optarglist_print( new_optargs, stderr );
#endif

  /* build hierarchical dendrogram */
  T = agglomerative_clustering( (const double**)distmatrix, eeg->ntrials, linkage );


  G = dblpp_init( n, n );
  d = dblpp_init( n, n );  
  indices = (int*)malloc( N*sizeof(int) );
  for( i=0; i<N; i++ ){
	 indices[i] = 1;
  }
  if( !out ){
	 out = eeg_init( eeg->nbchan, 1, n );
	 out = eeg_init_markers( nmarkers-2, out );
	 out->times = (double*)malloc( n*sizeof(double) );
	 memcpy( out->times, eeg->times, n*sizeof(double) );
  } else {
	 warnprintf( "overwriting your input eeg!\n");
  }
  
  P = init_warppath( ALLOC_IN_FCT, n, n );

  idx1=0; 
  idx2=0;
  /* now walk the tree to find pairs of trials to match */
  trials_left = N;
  if( progress ){
	 progress( PROGRESSBAR_INIT, N );
  }
  while( trials_left >= 2 ){
	 if( progress ){
		progress( PROGRESSBAR_CONTINUE_LONG, N-trials_left );
	 }
	 
	 Tsub = dgram_get_deepest( T );
	 idx1 = Tsub->left->val;
	 idx2 = Tsub->right->val;
	 dprintf("trials_left = %i, Tsub=%p\n", trials_left, Tsub);
	 dprintf("Tsub=%p, val=%i, lval=%i, rval=%i, d[%i,%i]=%f\n", Tsub, Tsub->val, Tsub->left->val, 
				Tsub->right->val,  idx1,  idx2, d[idx1][idx2]);

	 /* fill new markers to optargs */
	 if( pass_trial_markers_in_optargs ){
		for( i=1; i<=nmarkers-2; i++ ){
		  markers[0][i] = eeg->markers[idx1][i-1];
		  markers[1][i] = eeg->markers[idx2][i-1];
		}
	 }

	 weights[0] = indices[idx1]/(double)(indices[idx1]+indices[idx2]);
	 weights[1] = indices[idx2]/(double)(indices[idx1]+indices[idx2]);
	 dprintf("indices[%i,%i]=(%i,%i)\n", idx1, idx2, indices[idx1], indices[idx2]);
	 dprintf("weights=(%f,%f)\n", weights[0], weights[1] );
	 	 

	 dprintf(" compute regularization matrix G\n");
	 if( regularize ){
		G = regularize( G, n, new_optargs ); /* we need G only 
														once for all channels */
		if(!G){
		  errprintf("Regularization did not work for some reason\n");
		  return NULL;
		}
	 }

  
	 /* loop this for all channels */
	 for( c=0; c<eeg->nbchan; c++){	
		/* prepare average */
		s1 = eeg->data[c][idx1];
		s2 = eeg->data[c][idx2];
		if( progress ){
		  //		  oprintf("Trials (%i,%i): Channel=%i\n", idx1, idx2, chan); 
		  progress( PROGRESSBAR_CONTINUE_SHORT, 0 );
		}
		
		dprintf(" compute d for trial %i,%i\n", idx1, idx2);
		if(!(d = pointdistance( s1, n, s2, n, d, new_optargs ) ) ){
		  errprintf("pointdistance faulty\n");
		  return NULL;
		}
		dprintf(" ...done\n");  
		if( regularize ){
		  dprintf(" Regularize d\n");
		  //		  dtw_regularize_matrix( d, (const double**)G, n, n );
		  dprintf(" ...done\n");  
		}
		dprintf(" Compute path\n");  
		dtw_cumulate_matrix( d, n, n, new_optargs );
		P = dtw_backtrack( (const double**) d, n, n, P );
		dprintf(" ...done\n");  
		dprintf(" Warpavg\n");  
		s1s2 = warp_add_signals_by_path( s1, n, s2, n, P, ALLOC_IN_FCT, weights );
		s1s2_markers = warp_adjust_time_markers( eeg->markers[idx1], 
															  eeg->markers[idx2], 
															  nmarkers-2, ALLOC_IN_FCT,
															  (const double*)weights );
		dprintf(" ...done\n");  

		/* replace trials with ADTW */
		free( s1 );
		//		free( s2 );
		eeg->data[c][idx1] = s1s2; /* ADTW goes to idx1 */
		//		eeg->data[c][idx2] = NULL; /* do not touch this again */
		free( eeg->markers[idx1] );
		//		free( eeg->markers[idx2] );
		eeg->markers[idx1] = s1s2_markers;
		//		eeg->markers[idx2] = NULL;
	 }

	 indices[idx1] += indices[idx2];
	 indices[idx2]=-1; 			  /* never to be used again */
	 dprintf("new indices[%i,%i]=(%i,%i)\n", idx1, idx2, 
				indices[idx1], indices[idx2]);

	 /* replace node by leaf representing ADTW(idx1, idx2) */
	 Tsub->val = idx1;
	 Tsub->left = NULL;
	 Tsub->right = NULL;

	 trials_left--;
  }
  for( i=0; i<N; i++ ){
	 dprintf("indices[%i]=%i\n", i, indices[i]);
  }

  /* fill output data */
  for( c=0; c<eeg->nbchan; c++ ){
	 memcpy( out->data[c][0], eeg->data[c][idx1], n*sizeof(double) );
	 memcpy( out->markers[0], eeg->markers[idx1], (nmarkers-2)*sizeof(int) );
	 //	 free( eeg->data[c][idx1] );
	 //	 free( eeg->markers[idx1] );
  }
  eeg_free( eeg );


  if( progress ){
	 progress( PROGRESSBAR_FINISH, 0 );
  }

  /* cleaning up */
  dprintf("Freeing Memory\n");
  if( pass_trial_markers_in_optargs ){
	 optarglist_free( new_optargs );
  }
  dblpp_free( d, n );
  dblpp_free( G, n );
  dgram_free( T );
  free( indices );
  free_warppath( P );

  return out;
#endif
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
	 \param stmulus_marker gives the index indicating which of the markers within eeg_in  
	                       is the stimulus-onset
	 \param response_marker gives the index indicating which of the markers within eeg_in  
	                       is the response-onset						  
	 \param k parameter for gibbon's method
	 \return pointer to  newly allocated memory
*/ 
EEG* eeg_gibbons( EEG *eeg, int stimulus_marker, int response_marker, double k ){
  double *times, *new_times;
  double meanRT=0, curRT;
  double step;
  int n, N, i, j, l;
  int numchan, c;
  int stim_onset=0, resp_onset=0;
  EEG *eeg_out;

  for( i=0; i<eeg->ntrials; i++ ){
	 if( stimulus_marker>=eeg->nmarkers[i] ||
		  response_marker>=eeg->nmarkers[i] ||
		  stimulus_marker>=response_marker ){
		errprintf( "stimulus or response marker not correct (%i, %i) in trial %i, abort fct\n", 
					  stimulus_marker, response_marker, i );
		return NULL;
	 }
  }

  eeg_out = eeg_init( eeg->nbchan, 1, eeg->n );
  eeg_out->sampling_rate=eeg->sampling_rate; 
  if( eeg->times ){
	 eeg_out->times = (double*) malloc( eeg->n*sizeof(double) );
	 memcpy( eeg_out->times, eeg->times, eeg->n*sizeof(double) );
  }
  if( eeg->chaninfo ){
	 eeg_out->chaninfo = (ChannelInfo*) malloc( eeg->nbchan*sizeof(ChannelInfo) );
	 memcpy( eeg_out->chaninfo, eeg->chaninfo,  eeg->nbchan*sizeof(ChannelInfo) );
  }
  eeg_append_comment( eeg_out, "output from eeg_gibbons()\n");
  eeg_init_markers( 2, eeg_out );

  /* for convenience */
  times = eeg->times;  
  N = eeg->ntrials;
  n = eeg->n;
  numchan = eeg->nbchan;

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
  for( i=0; i<N; i++ ){
	 memcpy( new_times, times, n*sizeof(double) ); /* copy times */
	 stim_onset = eeg->markers[i][stimulus_marker];
	 resp_onset = eeg->markers[i][response_marker];
	 dprintf("markers = (%i,%i)\n", stim_onset, resp_onset);
	 curRT = times[resp_onset]-times[stim_onset];

	 for( j=stim_onset; j<resp_onset; j++ ){
		new_times[j] = times[j] + pow( times[j], k )
		  /pow( curRT, k )*( meanRT-curRT ); /* gibbons eq. */
	 }
	 /* warp the rest of the segments linearly */
	 step =(times[n-1]-curRT)/(double)(n-1-resp_onset);
	 l = 0;
	 for( j=resp_onset; j<n; j++ ){
		new_times[j]=meanRT + (l++)*step;
	 }

	 for( c=0; c<numchan; c++ ){ /* loop channels */
		gsl_spline_init (spline, new_times, eeg->data[c][i], n);
		for( j=0; j<n; j++ ){	  /* and samples */
		  eeg_out->data[c][0][j] += gsl_spline_eval ( spline, times[j], acc );
		}
	 }
  }
  eeg_out->markers[0][0] = stim_onset;
  eeg_out->markers[0][1] = meanRT;
  eeg_out->marker_labels[0][0] = create_string( "stimulus" );
  eeg_out->marker_labels[0][1] = create_string( "response" );

  for( c=0; c<numchan; c++ ){
	 for( j=0; j<n; j++ ){
		eeg_out->data[c][0][j] /= (double)N; /* average */
	 }
  }
  
  /* free */
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free( new_times );

  return eeg_out;
}

/*****************************************************************************
GOING TO BE OBSOLETE 
******************************************************************************/

WarpPath* init_warppath(WarpPath *path, int n1, int n2){
  if(path==NULL){
	 path = (WarpPath*)malloc(sizeof(WarpPath));
	 path->t1=NULL;
	 path->t2=NULL;
  }
  path->n = (n1+n2);
  path->n1= n1;
  path->n2= n2;

  if( path->t1 ){
	 free( path->t1 );
  } 
  if( path->t2 ){
	 free( path->t2 );
  }
  path->t1 = (int*)calloc(n1+n2, sizeof(int));
  path->t2 = (int*)calloc(n1+n2, sizeof(int));

  return path;
}

void      reset_warppath(WarpPath *P, int n1, int n2){
  memset( P->t1, 0, P->n*sizeof(int) );
  memset( P->t2, 0, P->n*sizeof(int) );
  P->n1 = n1; P->n2 = n2;
  P->n = 0;
  memset( P->t1, 0, n1*sizeof(int) );
  memset( P->t2, 0, n2*sizeof(int) );
}

void free_warppath(WarpPath *p){
	free(p->t1);
	free(p->t2);
	free(p);
}

void print_warppath( FILE *out, WarpPath *P ){
  int howmany, i;
  howmany=3;

  fprintf( out, "WarpPath '%p':\n"
			  " n1 = %i\n"
			  " n2 = %i\n"
			  " n  = %i\n"
			  " P  = ", P, P->n1, P->n2, P->n );
  if( P->n>=howmany ){
	 for( i=0; i<howmany; i++ ){
		fprintf( out, " (%i -> %i),", P->t1[i], P->t2[i] );
	 }
  } else if( P->n>0 ){
	 for( i=0; i<P->n; i++ ){
		fprintf( out, " (%i -> %i),", P->t1[i], P->t2[i] );
	 }
  } else {
	 fprintf( out, "<empty>");
  }
  fprintf( out, "\n" );
}

/** cumulate a distance matrix d to give
	 \f[
	 D_{jk} = \mathbf{d}_{jk}+\min{\{D_{j,k-1}, D_{j-1,k}, D_{j-1, k-1}\}}
	 \f]
	 \param d input/output matrix
	 \param M,N dimensions of d
	 \param optargs may contain:						 
	 - "slope_constraint=int" slope constraint, one of SLOPE_CONSTRAINT_*; default=SLOPE_CONSTRAINT_NONE
 */
void      dtw_cumulate_matrix  ( double **d, int M, int N, OptArgList *optargs ){
  /* j = 0,...,M-1
	  k = 0,...,N-1
  */
  int j, k;
  double x;
  double **D;
  SlopeConstraint constraint = SLOPE_CONSTRAINT_NONE;
  
  if( optarglist_has_key( optargs, "slope_constraint" ) ){
	 x = optarglist_scalar_by_key( optargs, "slope_constraint" );
	 if( !isnan( x ) ) constraint=(SlopeConstraint)x;
  }
  dprintf("using slope constraint '%i'\n", constraint );

  D = dblpp_init( M, N );
  dblpp_copy( (const double**)d, D, M, N );

   /* computing D_jk */
  for(j=0; j<M; j++){
    for(k=0; k<N; k++){
      if(k==0 && j==0) ;
      else if(k==0 && j>0){
		  D[j][k] = D[j-1][k] + d[j][k];
		} else if(k>0 && j==0){
		  D[j][k] = D[j][k-1] + d[j][k];
		} else { /* j,k > 0 */
		  D[j][k] = _slope_constraint( constraint, d, D, j, k ); // MIN(MIN(d[j][k-1], d[j-1][k]), d[j-1][k-1]);
		}
    }
  }
  
  dblpp_copy( (const double**)D, d, M, N );
  dblpp_free( D, M );
}
/**

 */
WarpPath* dtw_backtrack        ( const double **d, int M, int N, WarpPath *P ){ 
  int j,k; /* j = 0,...,M-1
				  k = 0,...,N-1  */
  int idx;
  double left, down, downleft;

  if( P==NULL ){
	 P = init_warppath( NULL, M, N );
  }

  /* Backtracking */
  j=M-1; k=N-1;

  idx = 1;
  P->t1[0] = j;
  P->t2[0] = k;
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


		/* this order of comparisons was chosen to ensure, that diagonal
			steps are preferred in case of equal distances
		*/
		if( downleft<=left ){	  /* downleft or down */
		  if( downleft<=down ){  
			 k--; j--;			 	  /* downleft */
		  } else {
			 k--;						  /* down */
		  } 
		} else {						  /* left or down */
		  if( down<=left ){
			 k--;						  /* down */
		  } else {
			 j--; 					  /* left */
		  }
		}


	 }	/* if */

	 P->t1[idx] = j;
	 P->t2[idx] = k;
	 //	 dprintf("(j, k), idx = (%i,%i), %i\n", j, k, idx);
	 if(idx>0 && P->t1[idx]>=P->t1[idx-1] && P->t2[idx]>=P->t2[idx-1]  ){
		warnprintf("P=%p: idx=%i\n", P, idx);
	 }
	 idx++;
  }
  P->n=idx-1;
  //  print_warppath( stderr, P );

  /* now the warppath has wrong order, reverse */
  for( j=0; j<(P->n)/2; j++ ){
	 swap2i( &(P->t1[j]), &(P->t1[P->n-1-j]) );
	 swap2i( &(P->t2[j]), &(P->t2[P->n-1-j]) );
  }
  //  print_warppath( stderr, P );

  return P;
}
