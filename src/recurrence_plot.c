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

#include "recurrence_plot.h"

/** \brief Calculate a (cross-)recurrence plot. 

	 Given the input multivariate time-series \f$\vec{s}_1(t)\f$
	 and \f$\vec{s}_2(t)\f$, calculate
	 \f[
	 \vec{C}^{\vec{s}_1,\vec{s}_2}(t_1,t_2) =  \Theta(\varepsilon - 
	  ||\vec{s}_1(t_1)-\vec{s}_2(t_2)||)
	 \f]

	 If \f$\vec{s}_1(t)=\vec{s}_2(t)\f$, it is a usual recurrence-plot,
	 else it is a cross-recurrence plot (CRP).

	 \param s1,s2 (multivariate) signals used to construct the 
	       recurrence plot (1D or 2D arrays)
	 \param out the recurrence plot or NULL -> memory is alloc'ed; for convenience, 
	       the recurrence plot is a double-matrix
	 \param epsilon; either a number given a fixed criterion for the 
	        distance between points in phase space
			  or: a fixed amount of nearest neighbours, if flags 
			  contains RPLOT_FAN
	 \param optargs  may contain:
	 <ul>
	 <li> <tt>fan=int</tt> flag, whether to use a "fixed-amount of nearest neighbours"
	      instead of the epsilon-ball. In this case, the epsilon argument is ignored
	 <li> optargs for the distance-function (optargs is passed 'as is' to this function)
	 </ul>
	 \return out (an N x N double matrix)
 */
Array* recplot( const Array *s1, const Array *s2, Array *out, double epsilon, 
					 OptArgList *optargs ){
  int fan=-1; /* fixed amount of neighbours? */
  double x;
  int i,j;
  if( optarglist_has_key( optargs, "fan" ) ){
	 x = optarglist_scalar_by_key( optargs, "fan" );
	 if( !isnan( x ) )
		fan=(int)x;
  }
  if( !array_comparable( s1, s2 ) ){
	 return NULL;
  }
  if( s1->ndim>2 || s2->ndim>2 || s1->dtype!=DOUBLE || s2->dtype!=DOUBLE){
	 errprintf("Input arrays must be vectors or matrices\n");
	 return NULL;
  }
  
  /* thin matrix-wrapper (in case of 1D data */
  Array *s1m, *s2m;
  s1m = array_fromptr2( DOUBLE, 2, s1->data, s1->size[0], (s1->ndim>1)?(s1->size[1]):1 );
  s2m = array_fromptr2( DOUBLE, 2, s2->data, s2->size[0], (s2->ndim>1)?(s2->size[1]):1 );

  int N = s1m->size[0];
  int p = s1m->size[1];
  if( !out ){
	 out = array_new2( DOUBLE, 2, N, N );
  } else {
	 if( out->ndim!=2 || out->dtype!=DOUBLE || out->size[0]<N || out->size[1]<N ){
		errprintf("Output must be a N x N matrix (N=%i)\n", N);
		return out;
	 }
  }

  /* init FAN */
  double *faneps=NULL;
  if( fan>0 ){
	 faneps=recplot_calculate_epsilons( s1m, s2m, NULL, fan );
  }

  /*---------- computation ----------- */
  double eps;
  for( i=0; i<N; i++ ){
	 if( faneps ){	 /* fixed amount of neighbours? */
		eps = faneps[i];
	 } else {
		eps = epsilon;
	 }

	 for( j=0; j<N; j++ ){
		if( eps - vectordist_euclidean( (double*)array_INDEXMEM2( s1m, i, 0 ), 
												  (double*)array_INDEXMEM2( s2m, j, 0 ), 
												  p, optargs ) > 0 ){
		  mat_IDX(out, i, j) = 1.0;
		} else {
		  mat_IDX(out, i, j) = 0.0;
		}
	 }
  }
  /*---------- /computation ----------- */

  /* clean up */
  array_free( s1m ); 
  array_free( s2m ); 
  if( faneps ) free( faneps );
  return out;
}

/** \brief calculate epsilons for each i such that a fixed amount of neighbours
	 is included.
	 
	 Currently it takes O(n^2) but could probably be done faster.

	 \param s1 2D-array (matrix), N x p
	 \param s2 2D-array (matrix), N x p
	 \param eps N-doubles or NULL
	 \param fan number of neighbours
	 \return eps (N doubles)
 */
double* recplot_calculate_epsilons( Array *s1, Array *s2, double *eps, int fan ){
  int i,j;
  bool ismatrix;
  if( !array_comparable( s1, s2 ) ){
	 return NULL;
  }
  matrix_CHECK( ismatrix, s1 );
  if( !ismatrix ) return NULL;
  matrix_CHECK( ismatrix, s2 );
  if( !ismatrix ) return NULL;

  int N=s1->size[0];
  int p=s1->size[1];

  if( !eps ){
	 eps = (double*)malloc( N*sizeof(double) );
  }
  double *tmp = (double*)malloc( N*sizeof(double) );
  double *epst= (double*)malloc( fan*sizeof(double) );

  /*---------- computation ----------- */
  for( i=0; i<N; i++ ){
	 for( j=0; j<N; j++ ){
		tmp[j] = vectordist_euclidean( (double*)array_INDEXMEM2( s1, i, j ), 
												 (double*)array_INDEXMEM2( s1, i, j ), 
												 p, NULL );
	 }
	 /*int gsl_sort_smallest (double * dest, size_t k, const double *
		src, size_t stride, size_t n) 
		This function copies the k smallest elements of the array src,
		of size n and stride stride, in ascending numerical order into
		the array dest. The size k of the subset must be less than or
		equal to n. The data src is not modified by this operation. */
	 gsl_sort_smallest( epst, fan, tmp, 1, N );
	 eps[i] = epst[fan-1];
  } 
  /*---------- /computation ----------- */
  
  free( tmp );
  free( epst);
  return eps;
}

/** calculate epsilons for each i such that a fixed amount of neighbours
	 is included.
	 \ref obsolete
 */
void recplot_fan( RecurrencePlot *R, PhaseSpace *p1, PhaseSpace *p2 ){
  int i,j;
  double *d;
  double *x1, *x2;
  double *eps;

  dprintf("Calcluating FAN\n");
  x1 = (double*) malloc ( p1->m*sizeof(double) );
  x2 = (double*) malloc ( p2->m*sizeof(double) );
  
  eps = (double*) malloc( R->fan*sizeof(double) );
  d   = (double*) malloc( p2->xn*sizeof(double) );
  R->nepsilon = p1->xn;
  R->epsilon = (double*) malloc( R->nepsilon * sizeof(double) );

  for( i=0; i<p1->xn; i++ ){
	 for( j=0; j<p2->xn; j++ ){
		phspace_index_i( p1, i, x1 );
		phspace_index_i( p2, j, x2 );
		d[j] = vectordist_euclidean( x1, x2, p1->m, NULL );
	 }
	 gsl_sort_smallest( eps, R->fan, d, 1, p2->xn );
	 R->epsilon[i] = eps[(R->fan)-1];
  }
  free( eps );
  free( d );
}

/** prepare a recurrence plot. 
	 \param m,n dimensions of the plot
	 \param epsilon; either a number given a fixed criterion for the 
	        distance between points in phase space
			  or: a fixed amount of nearest neighbours, if flags 
			  contains RPLOT_FAN
    \param flags flags beginnging with RPLOT_	 
	 \ref obsolete
 */
RecurrencePlot* recplot_init( int m, int n, double epsilon, int flags ){
  RecurrencePlot *R;
  int i;

  R = (RecurrencePlot*)malloc( sizeof(RecurrencePlot) );
  R->m = m;
  R->n = n;
  R->nepsilon = 0;
  R->epsilon = NULL;
  R->fixed_epsilon = -1;
  R->fan = 0;

  if( flags & RPLOT_FAN ){
	 R->fan = (int) epsilon;
  } else {
	 R->fixed_epsilon = epsilon; 
  }

  R->R = (double**)malloc( m*sizeof(double*) );
  for( i=0; i<m; i++ ){
	 R->R[i] = (double*)malloc( n*sizeof(double) );
  }

  return R;
}

/**
	 \ref obsolete
*/
void recplot_print( FILE *out, RecurrencePlot *R ){
  int i;
  fprintf( out, "RecurrencePlot '%p':\n"
			  " m   = %i\n"
			  " n   = %i\n"
			  " fan = %i\n", R, R->m, R->n, R->fan );
  if( R->fan ){
	 fprintf( out, " e   = {");
	 for( i=0; i<R->nepsilon; i++ ){
		fprintf( out, "%.2f ", R->epsilon[i] );
	 }
	 fprintf( out, "}\n");
  } else {
	 fprintf( out, " e = %f\n", R->fixed_epsilon );
  }
}

/**
	 \ref obsolete
*/
void recplot_free( RecurrencePlot *R ){
  int i;
  for( i=0; i<R->m; i++ ){
	 free( R->R[i] );
  }
  free( R->R );
  free( R->epsilon );
  free( R );
}

/** calculate the recurrence plot from two signals in phase-space representation
	 p1,p2. If p1!=p2, it is the cross-recurrence plot.
	 \f[
	 R_{ij} = \Theta\left( \epsilon - ||\vec{x}_i - \vec{x}_j||^2  \right)
	 \f]	
	 Depending on the settings in R, epsilon is chosen for each i or a fixed
	 epsilon is used.
	 \ref obsolete

	 \param p1,p2
	 \param R is the recurrence plot. if NULL, own memory is allocated.
*/
void recplot_calculate( RecurrencePlot *R, PhaseSpace *p1, 
								PhaseSpace *p2 ){
  int i, j;
  double *x1, *x2;
  double epsilon;

  x1 = (double*) malloc ( p1->m*sizeof(double) );
  x2 = (double*) malloc ( p2->m*sizeof(double) );

  /* init FAN */
  if( R->fan ){
	 recplot_fan( R, p1, p2 );
  }

  for( i=0; i<R->m; i++ ){
	 /* fixed amount of neighbours? */
	 if( R->fan ){
		epsilon = R->epsilon[i];
	 } else {
		epsilon = R->fixed_epsilon;
	 }

	 /* calculation main loop */
	 for( j=0; j<R->n; j++ ){
		phspace_index_i( p1, i, x1 );
		phspace_index_i( p2, j, x2 );
		if( epsilon - vectordist_euclidean( x1, x2, p1->m, NULL ) > 0 ){
		  R->R[i][j] = 1.0;
		} else {
		  R->R[i][j] = 0.0;
		}
	 }
  }

  /* cleaning up */
  free( x1 );
  free( x2 );
}

/** Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot, following
	 the algorithm in Marwan et al. Cross recurrence plot based synchronization
	 of time series. Nonlinear Processes in Geophysics (2002) vol. 9 (3-4) pp. 325-331.

	 Algorithm:
	 \code
	 Input: Parameters dx, dy
	 1. find recurrence point next to origin (i,j)
	 2. square window of size w=2 starting from (i,j);
	    if another point in this window, goto 3.
		 else w++, goto 2.
	 3. if point is found in x-direction, increase w in x-dir
	    if point is found in y-direction, increase w in y-dir
		 until: either no point is found or window is of 
		        size (w+dx, w+dy)
		 after this step, window is of size (w+deltawx, w+deltawy)
	 4. set i^ = i + (w+deltawx)/2
	    set j^ = j + (w+deltawy)/2
	 5. (i,j) = (i^,j^), goto 2.
	 \endcode
 */
WarpPath* recplot_los_marwan( const RecurrencePlot *R, int dx, int dy ){
  WarpPath *P,*P2;
  int ii, ij,   /* current LOS-point */
	 xcont,ycont,/* continue in x/y direction flag */
	 deltawi, deltawj, /* addition in x/y direction to (wi,wj) */
	 wi, wj;     /* dimensions of the window */
  int i,j, pidx, breakflag;
  int no_recpoint_found, reached_border;

  P = init_warppath( NULL, R->m, R->n );
  
  breakflag=0;
  ii=0; ij=0;
  for( i=0; MIN( R->m, R->n ); i++ ){ /* find first point (1) */
	 for( j=0; j<=i; j++ ){
		if( R->R[i][j] ){
		  ii = i; ij = j;
		  breakflag=1;
		} else if( R->R[j][i] ){
		  ii = j; ij = i;
		  breakflag=1;
		}
		if(breakflag) 
		  break;
	 }
	 if(breakflag)
		break;
  }
  /* dprintf("(ii,ij) = (%i,%i)\n", ii,ij); */
  pidx = 0;
  P->t1[pidx]=ii;
  P->t2[pidx]=ij;
  pidx++;

  reached_border=0;
  while( !reached_border ){
	 wi = 2;
	 wj = 2;
	 no_recpoint_found=1;
	 xcont = 0;
	 ycont = 0;
	 while( no_recpoint_found ){ /* initial window sweep (2) */
		if( ii+wi>=R->m || ij+wj>=R->n ){
		  reached_border=1;
		  break;
		}

		for( i=ii; i<ii+wi; i++ ){
		  if( R->R[i][ij+wj] ){
			 no_recpoint_found=0;
			 xcont = 1;
			 break;
		  }
		} /* for i */

		for( j=ij; j<ij+wj; j++ ){
		  if( R->R[ii+wi][j] ){
			 no_recpoint_found=0;
			 ycont = 1;
			 break;
		  }
		} /* for j */

		wi++;	wj++;
	 } /* while(no_recpoint_found) */
	 if( reached_border )
		break;

	 /* dprintf("w=(%i,%i), cont=(%i,%i)\n", wi,wj,	xcont, ycont ); */
	 deltawi=0;
	 deltawj=0;
	 while( xcont || ycont ){	  /* (3) */
		if( ij+wj+deltawj>=(R->n-1) || ii+wi+deltawi>=(R->m-1) ){
		  reached_border=1;
		  break;
		}
		if( xcont ){ 				  /* check border in x */
		  deltawi++;
		  for( j=ij; j<ij+wj+deltawj; j++ ){
			 if( R->R[ii+wi+deltawi][j] ){
				xcont = 1;
			 } else {
				xcont = 0;
			 }
		  } /* for j */
		}

		if( ycont ){				  /* check border in y */
		  deltawj++;
		  for( i=ii; i<ii+wi+deltawi; i++ ){
			 if( R->R[i][ij+wj+deltawj] ){
				ycont=1;
			 } else {
				ycont=0;
			 }
		  } /* for i */
		}
		if( deltawi >= dx || deltawj >= dy )
		  break;
	 }	/* while xcont */
	 /* dprintf("deltaw=(%i,%i)\n", deltawi, deltawj ); */
	 if( reached_border )
		break;

	 ii = ii + (wi+deltawi)/2;
	 ij = ij + (wj+deltawj)/2;
	 P->t1[pidx]=ii;
	 P->t2[pidx]=ij;
	 pidx++;
  } /* while !reached_border */
  P->t1[pidx]=R->m;
  P->t2[pidx]=R->n;
  P->n = ++pidx;

  /* fill gaps with bresenham */
  P2 = init_warppath( NULL, R->m, R->n );
  pidx = 0;
  int nump;
  int *line;
  line = (int*)malloc( R->m*R->n*2*sizeof(int)); /* much too much */

  for( i=0; i<P->n-1; i++ ){
	 nump = bresenham_howmany_points( P->t1[i], P->t2[i], P->t1[i+1], P->t2[i+1] );
	 line = bresenham( P->t1[i], P->t2[i], P->t1[i+1], P->t2[i+1], line );
	 dprintf(" (%i,%i)->(%i,%i)\n", P->t1[i], P->t2[i], P->t1[i+1], P->t2[i+1] );
	 for( j=0; j<nump; j++ ){
		dprintf(" line=(%i,%i)\n", line[(2*j)+0], line[(2*j)+1]);
		P2->t1[pidx++]=line[(2*j)+0];
		P2->t2[pidx  ]=line[(2*j)+1];
	 }
  }
  P2->n = pidx;

  free_warppath( P );
  free( line );

  return P2;
}

/** Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot using 
	 a dynamic time-warping strategy.

	 Algorithm:
	 \code
	 Input: RecurrencePlot R (1 where recurrence, 0 otherwise)
	 1. calculate d = 1-R;
	 2. cumulate d, such that D[jk] = min{ D[j-1k], D[jk-1], D[j-1k-1] }
	 3. Backtrack
	 \endcode
	 (see \ref dtw for details)
 */
WarpPath* recplot_los_dtw( const RecurrencePlot *R ){
  double **d;
  WarpPath *P;


  d = dblpp_init( R->m, R->n );
  dblpp_copy( (const double**)R->R, d, R->m, R->n ); 

  /* flip binary matrix */
  scalar_minus_dblpp( 1.0, d, R->m, R->n );
  /* add a bias such that diagonal is preferred */
  
  
  dtw_cumulate_matrix( d, R->m, R->n, NULL );
  //dblpp_normalize_by_max( D, R->m, R->n );
  P = dtw_backtrack( (const double**) d, R->m, R->n, NULL );

  dblpp_free( d, R->m );

  return P;
}

/** Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot using 
	 a regularized dynamic time-warping strategy.

	 Algorithm:
	 \code
	 Input: RecurrencePlot R (1 where recurrence, 0 otherwise)
	 1. calculate d = 1-R;
	 2. get the regularization function G that is the distance transform of the
	    linear interpolation between corresponding points in time 
	 3. d = d*G;
	 4. cumulate d, such that D[jk] = min{ D[j-1k], D[jk-1], D[j-1k-1] }
	 5. Backtrack
	 \endcode
	 (see \ref dtw for details)

	 \param R the recurrence plot
	 \param markers 2xnmarkers array; markers[0] is markers of 1st signal,
	                markers[1] for second; if NULL, { {0,n-1}, {0,n-1} } is used
	 \param nmarkers number of markers
	 \return the warping function
 */
WarpPath* recplot_los_dtw_markers( const RecurrencePlot *R, int **markers, int nmarkers ){
  double **d, **G;
  WarpPath *P;
  int mflag=0;
  int n;

  n = R->m;
  d = dblpp_init( R->m, R->n );
  dblpp_copy( (const double**)R->R, d, R->m, R->n ); 

  /* flip binary matrix */
  scalar_minus_dblpp( 1.0, d, R->m, R->n );

  /* add a bias such that f is preferred */
  if( !markers ){
	 mflag = 1;
	 markers = (int**)malloc( 2*sizeof( int* ) );
	 markers[0] = (int*)malloc( 2*sizeof( int ) );
	 markers[1] = (int*)malloc( 2*sizeof( int ) );
	 markers[0][0] = 0; markers[1][0] = 0;
	 markers[0][1] = n-1; markers[1][1] = n-1;
  }
  OptArgList *opts = optarglist( "markers=int**,nmarkers=int",
											markers, nmarkers );
  G = regularization_linear_points( ALLOC_IN_FCT, R->m, opts );
  optarglist_free( opts );
  dblpp_normalize_by_max( G, R->m, R->n );

  /* apply */
  dblpp_add_dblpp( d, (const double**)G, R->m, R->n );
  
  /* dtw */
  dtw_cumulate_matrix( d, R->m, R->n, NULL );
  P = dtw_backtrack( (const double**) d, R->m, R->n, NULL );

  /* clean */
  dblpp_free( d, R->m );
  dblpp_free( G, R->m );
  if( mflag ){
	 free( markers[0] );
	 free( markers[1] );
	 free( markers );
  }

  return P;
}

/** Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot using 
	 an algorithm based on the distance transform of the plot.

	 Algorithm:
	 \code
	 Input: RecurrencePlot R (1 where recurrence, 0 otherwise)
	 1. calculate d = DT(R) \ref disttransform
	 2. cumulate d, such that D[jk] = min{ D[j-1k], D[jk-1], D[j-1k-1] }
	 3. Backtrack
	 \endcode
	 (see \ref dtw for details)

	 \param R the recurrence plot
	 \return the warping function
 */
WarpPath* recplot_los_disttransform( const RecurrencePlot *R ){
  WarpPath *P;
  int **I;
  double **d;
  int i,j;

  /* dummy, needed for distance transform */
  I = (int**) malloc( R->m*sizeof( int* ) );
  for( i=0; i<R->m; i++ ){
	 I[i] = (int*) malloc( R->n*sizeof( int ) );
	 for( j=0; j<R->n; j++ ){
		I[i][j] = 0;
		if( R->R[i][j]>0 ){
		  I[i][j] = 1;
		}
	 }
  }

  d = disttransform_deadreckoning( I, R->m, R->n, ALLOC_IN_FCT );
  dtw_cumulate_matrix( d, R->m, R->n, NULL );
  P = dtw_backtrack( (const double**) d, R->m, R->n, NULL );

  /* free */
  dblpp_free( d, R->m );
  for( i=0; i<R->m; i++ )
	 free( I[i] );
  free( I );

  return P;
}

/** Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot using 
	 a dynamic time-warping strategy.

	 Algorithm:
	 \code
	 Input: RecurrencePlot R (1 where recurrence, 0 otherwise)
	 1. calculate d = 2-R;
	 2. add a small amount of noise:
	     if d_jk = 1, d_jk = d_jk + epsilon
	 2. cumulate d, such that D[jk] = min{ D[j-1k], D[jk-1], D[j-1k-1] }
	 3. Backtrack
	 \endcode
	 (see \ref dtw for details)
 */
WarpPath* recplot_los_dtw_noise( const RecurrencePlot *R ){
  double **d;
  WarpPath *P;
  double noise;
  double noiseamp = 0.01;
  int i,j;

  srand((long)time(NULL));
  d = dblpp_init( R->m, R->n );
  dblpp_copy( (const double**)R->R, d, R->m, R->n ); 

  /* flip binary matrix */
  scalar_minus_dblpp( 1.0, d, R->m, R->n );
  
  /* add small amount of noise */
  for( i=0; i<R->m; i++ ){
	 for( j=0; j<R->n; j++ ){
		noise = noiseamp*(((double)rand()) / RAND_MAX);		
		//		if( d[i][j]<2 )
		  d[i][j] += noise;
	 }
  }
  
  OptArgList *opt = optarglist( "slope_constraint=int", SLOPE_CONSTRAINT_LAX );
  dtw_cumulate_matrix( d, R->m, R->n, opt );
  optarglist_free( opt );

  P = dtw_backtrack( (const double**) d, R->m, R->n, NULL );

  dblpp_free( d, R->m );

  return P;
}
