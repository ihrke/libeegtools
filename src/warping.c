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
#include <math.h>

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

	 \ingroup grpdtw
	 \f[
	 D_{jk} = d_{jk}+\min{\{D_{j,k-1}, D_{j-1,k}, D_{j-1, k-1}\}}
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
	 
	 \ingroup grpdtw
	 \param d (M x N) is the cumulated distances matrix (usually output from matrix_dtw_cumulate())
	 \return the warp-Path (2D uint array, 2 x N); row 1 has the column-indices, row 2 has the row indices
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

/** \brief Add signals according to a warppath and downsample to specified time-values.

  \warning this works currently only for 1D signals!

	 \ingroup grpdtw

	 If called from hierarchical averaging routines, you might
	 want to pass a "weights" field in the optional arguments.

	 Both arrays must be of equal length N x p.

	 \warning p=1, currently!

	 \param s1 N x p (DOUBLE) array; first signal
	 \param s2 N x p (DOUBLE) array; second signal
	 \param path - contains warppath (2xN INT);
	 \param times is the common time-scale of s1 and s2
	 \param opts may contain:
	 - "weights=double*" - weights in average, for using it with hierarchical averaging;
					  should be 2 double values with w[0]+w[1]=1.0
	 \return the warped average of the signals; Nxp samples
 */
Array*  dtw_add_signals_downsample( const Array *s1, const Array *s2, const Array *times, const Array *path, OptArgList *opts ){
	int ispath;
	warppath_CHECK( ispath, path );
	if( !ispath ) return NULL;

	bool isvector;
	vector_CHECK( isvector, times );
	if( !isvector ){ errprintf("times not a vector\n");	return NULL; }
	vector_CHECK( isvector, s1 );
	if( !isvector ){ errprintf("s1 not a vector\n");	return NULL; }
	vector_CHECK( isvector, s2 );
	if( !isvector ){ errprintf("s2 not a vector\n");	return NULL; }


	int N = s1->size[0];
	int Nn = path->size[1];
	dprintf("Signals have %i samples\n", N);

	double defaultw[2]={0.5,0.5};
	void *tmp;
	double *w=defaultw;
	optarg_PARSE_PTR( opts, "weights", w, double*, tmp );

	dprintf("weights=(%f,%f)\n", w[0], w[1]);
	/* --------------------- computation ------------------------*/
	Array *wa = array_new2( DOUBLE, 1, Nn );
	Array *ntimes=array_new2( DOUBLE, 1, Nn );
	int i;
	for( i=0; i<Nn; i++ ){
		vec_IDX( ntimes, i)=w[0]*vec_IDX( times, array_INDEX2(path,uint,0,i) )+
							w[1]*vec_IDX( times, array_INDEX2(path,uint,1,i) );
		vec_IDX( wa, i) = w[0]*vec_IDX( s1, array_INDEX2(path,uint,0,i) ) +
						  w[1]*vec_IDX( s2, array_INDEX2(path,uint,1,i) );
	}
	Array *wavg=vector_interp1( ntimes, wa, times, NULL );
	/* --------------------- /computation ------------------------*/

	array_free( ntimes );
	array_free( wa );

	return wavg;
}

/** \brief Add signals according to a warppath.

	 \ingroup grpdtw

	 If called from hierarchical averaging routines, you might
	 want to pass a "weights" field in the optional arguments.

	 \param s1 N1 x p (DOUBLE) array; first signal
	 \param s2 N2 x p (DOUBLE) array; second signal
	 \param path - contains warppath (2xN INT);
	 \param opts may contain:
	 - "weights=double*" - weights in average, for using it with hierarchical averaging;
	                  should be 2 double values with w[0]+w[1]=1.0
	 \return the warped average of the signals; N1+N2 samples
 */
Array*  dtw_add_signals( const Array *s1, const Array *s2, const Array *path, OptArgList *opts ){
  int ispath;
  warppath_CHECK( ispath, path );
  if( !ispath ) return NULL;
  
  if( s1->dtype!=DOUBLE || s2->dtype!=DOUBLE ){
	 errprintf("Input signals must be double-arrays\n");
	 return NULL;
  } 
  if( s1->ndim>2 || s2->ndim>2 || s1->ndim != s2->ndim ){
	 warnprintf("Input signals should be 1D or 2D... continue at your own risk (%i,%i)\n",
					s1->ndim, s2->ndim );
  }
  /* thin matrix-wrapper (in case of 1D data */
  Array *s1m, *s2m;
  s1m = array_fromptr2( DOUBLE, 2, s1->data, s1->size[0], (s1->ndim>1)?(s1->size[1]):1 );
  s2m = array_fromptr2( DOUBLE, 2, s2->data, s2->size[0], (s2->ndim>1)?(s2->size[1]):1 );

  if( s1m->size[1] != s2m->size[1] ){
	 errprintf("Signals must have the same dimension 2\n");
	 array_free( s1m );
	 array_free( s2m );
	 return NULL; 
  }
  int N1 = s1m->size[0];
  int N2 = s2m->size[0];
  int Nn = path->size[1];
  int p = s1m->size[1];

  dprintf("Signals have %i and %i samples and are of dimension p=%i\n", N1, N2, p);

  double defaultw[2]={0.5,0.5};
  void *tmp;
  double *w=defaultw;
  optarg_PARSE_PTR( opts, "weights", w, double*, tmp );
  

  /* --------------------- computation ------------------------*/
  Array *wa = array_new2( DOUBLE, 2, Nn, p );
  int i,j;
  for( i=0; i<Nn; i++ ){
	 for( j=0; j<p; j++ ){
		mat_IDX( wa, i, j )=w[0]*mat_IDX( s1m, array_INDEX2(path,uint,0,i),j )
		  + w[1]*mat_IDX( s2m, array_INDEX2(path,uint,1,i),j );
	 }
  }
  /* --------------------- /computation ------------------------*/

  array_free( s1m );
  array_free( s2m );
  return wa;
}

/** \brief Warp-average according to Gibbons+Stahl 2007.
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
	 \param data C x N x n data input (DOUBLE, channels x trials x samples)
	 \param times n-vector giving the sampling times of the data
	 \param markers N x M time-markers (INT) in sampling-point units
	 \param stmulus_marker gives the index indicating which of the markers within markers
						   is the stimulus-onset; 0<m<M
	 \param response_marker gives the index indicating which of the markers within markers
						   is the response-onset
	 \param k parameter for gibbon's method
	 \param alloc if TRUE, new memory is allocated, else it is overwritten
	 \return pointer to  newly allocated memory
*/
Array* array_warp_gibbons( Array *data, Array *markers, int stimulus_marker, int response_marker, double k, bool alloc ){
	if( data->dtype!=DOUBLE || data->ndim!=3 ){
		errprintf("Data must be 3D (C x N x n) and DOUBLE precision\n");
		return NULL;
	}
	int C,N,n;
	C=data->size[0]; N=data->size[1]; n=data->size[2];
	if( markers->dtype!=INT || markers->ndim!=2 || markers->size[0]!=N ){
		errprintf("Markers must be 2D (N x M) and INT\n");
		return NULL;
	}
	int M=markers->size[0];
	if( stimulus_marker<0 || stimulus_marker>=M || response_marker<0 || response_marker>=M ){
		errprintf("Stimulus and Response-Markers must be between 0 and %i", M);
		return NULL;
	}
	if( k<0 || k>4 ){
		errprintf("k must be between 0 and 4\n");
		return 0;
	}

	Array *out;
	if( alloc ){
		out = array_copy( data, TRUE );
	} else {
		out=data;
	}

	/* gsl interpolation allocation */
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n); /* gibbons uses linear */

	/* computing mean reaction time */
	double meanRT=0.0;
	int c,i,j;
	for( i=0; i<N; i++ ){
		meanRT+=array_INDEX2( markers, int, i, response_marker )-array_INDEX2( markers, int, i, stimulus_marker );
	}
	meanRT /= (double)N;

	/* compute warping for each trial */
	Array *ntimes = array_new2( DOUBLE, 1, n );
	int stim_onset;
	int resp_onset;
	double curRT;
	double step;
	int l;
	for( i=0; i<N; i++ ){
		for( j=0; j<n; j++ ){
			vec_IDX( ntimes, j )=(double)j;
		}
	   stim_onset = array_INDEX2( markers, int, i, stimulus_marker );
		resp_onset = array_INDEX2( markers, int, i, response_marker );
	   dprintf("markers = (%i,%i)\n", stim_onset, resp_onset);
	   curRT = resp_onset-stim_onset;

	   for( j=stim_onset; j<resp_onset; j++ ){
		  vec_IDX( ntimes, j) = (double)j + pow( (double)j, k )
			/pow( curRT, k )*( meanRT-curRT ); /* gibbons eq. */
	   }
	   /* warp the rest of the segments linearly */
	   step =((double)(n-1)-curRT)/(double)(n-1-resp_onset);
	   l = 0;
	   for( j=resp_onset; j<n; j++ ){
		  vec_IDX( ntimes,j)=meanRT + (l++)*step;
	   }

	   for( c=0; c<C; c++ ){ /* loop channels */
		  gsl_spline_init (spline, ntimes->data, array_INDEXMEM3( data, c, i, 0), n);
		  for( j=0; j<n; j++ ){	  /* and samples */
			array_INDEX3( data, double, c, i, j ) = gsl_spline_eval ( spline, (double)j, acc );
		  }
	   }
	}

	/* free */
	gsl_spline_free (spline);
	gsl_interp_accel_free (acc);
	array_free( ntimes );

	return out;
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
#ifdef FIXEEG
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
#endif
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
	 D_{jk} = d_{jk}+\min{\{D_{j,k-1}, D_{j-1,k}, D_{j-1, k-1}\}}
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


#ifdef EXPERIMENTAL

/**
	return the sum of the individual pairwise
	distances (euclidean) divided by N(n-1)/2
 */
double multidist( Array *dat, uint *idx ){  
  bool ismat;
  matrix_CHECK( ismat, dat );
  if( !ismat ) return NAN;
  
  double dist=0.0;
  int N=dat->size[0];
  int i,j;

  for( i=0; i<N-1; i++ ){
	 for( j=i+1; j<N; j++ ){
		//		dprintf("(i,j)=%i,%i, idx=%i,%i\n",i,j, idx[i], idx[j]);
		dist += 0.9* ABS( mat_IDX( dat, i, idx[i] )-
								mat_IDX( dat, j, idx[j] ) ) +
		  0.1* ABS( idx[i]-idx[j] );
	 }
  }
  dist /= (double)((N*(N-1))/2);

  //  dprintf("d(i,j)=%f\n", dist );
  return dist;
}

/** \brief Calculate the N-dimensional warping path (brute force).


	 \warning this function uses A LOT of memory and CPU. Memory
	 and CPU usage scales O(n^N) where N is the number of trials
	 and n the number of sampling points. In practice, that means that 
	 it is usually only applicable for 3-4 trials with limited sampling
	 rate (depending on the kind of computer you use).


	 \note works for now only for 1D signals of equal length

	 Example: 3 trials with the 3D warping path

	 \image html multiwarp.jpg

	 \param in 
	 \param optargs may contain:
	 - <tt>progress=void*</tt> progress-bar function, called every now and then.
	 \return the N-dimensional warp-path
 */
Array* multiwarp( Array *in, OptArgList *optargs ){
  bool ismat;
  matrix_CHECK( ismat, in );
  if( !ismat ) return NULL;
  void *ptr;
  ProgressBarFunction progress=NULL;  
  optarg_PARSE_PTR( optargs, "progress", progress, ProgressBarFunction, ptr );
  if( progress )
	 dprintf("found progressbar\n");

  int i, j;  
  int N=in->size[0]; /* num signals */
  int n=in->size[1]; /* num samples per signal */

  dprintf("N=%i,n=%i\n", N, n );

  /* allocate N-dimensional distance matrix */
  uint *dims;
  MALLOC( dims, N, uint );
  for( i=0; i<N; i++ )
	 dims[i]=n;
  Array *distmat=array_new( DOUBLE, N, dims );
  free( dims );

  /* calculate distmat and cumulate at the same time */
  uint *idx, *cidx, *tidx;
  MALLOC( idx, N, uint );
  MALLOC( cidx, N, uint );
  MALLOC( tidx, N, uint );

  ulonglong lidx=0;
  Array *minvec=array_new2( DOUBLE, 1, (1<<N)-1 );

  if( progress )  progress( PROGRESSBAR_INIT, array_NUMEL( distmat ) );
  for( lidx=0; lidx<array_NUMEL( distmat ); lidx++ ){	 
	 /*--- distance calculation */
	 array_calc_rowindex( lidx, distmat->size, distmat->ndim, idx );
	 (*((double*)array_index( distmat, idx )))=multidist( in, idx );
	 //	 dprintf("lidx=%li, %i, %i, (1<<N)-1=%i\n", lidx, idx[0], idx[1], (1<<N)-1 );	 

	 /*--- cumulate the matrix */
	 /* get all elements "surrounding" idx */
	 for( i=0; i<( 1<<N )-1; i++ ){
		for( j=0; j<N; j++ ){ /* calculate index */
		  cidx[j] = idx[j];
		  if( CHECK_BIT( i+1, j ) ){
			 cidx[j] = (cidx[j]==0)?(0):(cidx[j]-1);
		  }
		}
		vec_IDX( minvec, i ) = (*(double*)array_index( distmat, cidx ));
	 }

	 /* cumulating */
	 (*((double*)array_index( distmat, idx ))) += (*(double*)array_min( minvec ));

	 if( progress && (lidx % n==0) & lidx/n % n==0 ) 
		progress( PROGRESSBAR_CONTINUE_LONG, lidx );
  }
  if( progress ) progress( PROGRESSBAR_FINISH, 0 );

  /*--- backtracking-phase */
  Array *p1=array_new2( UINT, 2, N, N*n ); /* dummy of maximum length */

  dprintf("Backtrack\n");
  lidx=0;
  for( i=0; i<N; i++ ){
	 idx[i]=n-1; /* start in the upper right */
	 (*(uint*)array_index2( p1, i, lidx ))=idx[i];
  }
  lidx=1;
  bool continue_flag, dontcheck_flag;
  double dmin;
  while(1){
	 continue_flag=FALSE;
	 dontcheck_flag=FALSE;

	 /* get all elements "surrounding" idx and put minimal index in idx */
	 dmin=DBL_MAX;
	 for( i=0; i<( 1<<N )-1; i++ ){
		for( j=0; j<N; j++ ){ /* calculate index */
		  cidx[j] = idx[j];
		  if( CHECK_BIT( i+1, j ) ){
			 if( cidx[j]==0 ) 
				dontcheck_flag=TRUE;
			 else 
				cidx[j] = cidx[j]-1;
		  }
		}
		if( !dontcheck_flag ){
		  //		  dprintf("Check field: (%i,%i,%i), dmin=%f\n", cidx[0], cidx[1], cidx[2], (*(double*)array_index( distmat, cidx )) );

		  if(  (*(double*)array_index( distmat, cidx )) < dmin ){
			 dmin = (*(double*)array_index( distmat, cidx ));
			 for( j=0; j<N; j++ )
				tidx[j]=cidx[j];
		  }
		} else {
		  dontcheck_flag=FALSE;
		}
	 }	 
	 for( j=0; j<N; j++ )
		idx[j]=tidx[j];

	 //	 dprintf("Select field: (%i,%i,%i), dmin=%f\n", idx[0], idx[1], idx[2], dmin );

	 /* advance path */
	 for( j=0; j<N; j++ ){
		(*(uint*)array_index2( p1, j, lidx ))=idx[j];
	 }
	 lidx++;

	 /* done? */
	 for( i=0; i<N; i++ ){
		if( idx[i]>0 ){
		  continue_flag=TRUE;
		  break;
		}
	 }

	 if( !continue_flag )
		break;
  }

  /* reverse warppath and trim to correct size */
  Array *p2=array_new2( UINT, 2, N, lidx );
  for( i=0; i<lidx; i++ ){
	 for( j=0; j<N; j++ ){
		array_INDEX2( p2, uint, j, lidx-i-1 ) = array_INDEX2( p1, uint, j, i );
	 }
  }

  /* cleaning up */
  array_free( minvec );
  free( cidx );
  free( idx );
  free( tidx );
  array_free( distmat );
  array_free( p1 );

  return p2;
}

Array* multiwarp_add( Array *in, Array *times, Array *path ){
  bool ismat, isvec;
  matrix_CHECK( ismat, in );
  if( !ismat ) return NULL;
  vector_CHECK( isvec, times );
  if( !isvec ) return NULL;

  int i, j;  
  int N=in->size[0]; /* num signals */
  int n=in->size[1]; /* num samples per signal */
  int Nn=path->size[1];

  dprintf("N=%i,n=%i,Nn=%i\n", N, n, Nn );

  double *w;
  MALLOC( w, N, double );
  for( i=0; i<N; i++ ){
	 w[i] = 1.0/(double)N;
  }
  dprintf("w=%f\n", w[0] );

  /* return times array in first dim */
  Array *wa = array_new2( DOUBLE, 2, 2, Nn );
  for( i=0; i<Nn; i++ ){
	 for( j=0; j<N; j++ ){
		mat_IDX( wa, 0, i) += w[j]*vec_IDX( times, array_INDEX2(path,uint,j,i));
		mat_IDX( wa, 1, i) += w[j]*mat_IDX( in, j, array_INDEX2(path,uint,j,i));
	 }
  }
  
  free( w );
  
  return wa;
}


#endif
