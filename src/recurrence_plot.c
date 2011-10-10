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
#include "imageproc.h"
#include "linalg.h"

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
	 \return out (an N (s1) x N (s2) double matrix);
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
	 dprintf("FAN=%i\n", fan);
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
  s1m = array_fromptr2( DOUBLE, 2, s1->data, s1->size[0],
								(s1->ndim>1)?(s1->size[1]):1 );
  s2m = array_fromptr2( DOUBLE, 2, s2->data, s2->size[0],
								(s2->ndim>1)?(s2->size[1]):1 );

  int N = s1m->size[0];
  int p = s1m->size[1];
  if( !out ){
	 out = array_new2( DOUBLE, 2, N, N );
  } else {
	 if( out->ndim!=2 || out->dtype!=DOUBLE ||
		  out->size[0]<N || out->size[1]<N ){
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
	 /*dprintf("Current Epsilon=%f\n", eps );*/

	 for( j=0; j<N; j++ ){
		if( vectordist_euclidean( (double*)array_INDEXMEM2( s1m, i, 0 ),
								  (double*)array_INDEXMEM2( s2m, j, 0 ),
								  p, optargs ) < eps ){
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
		tmp[j] = vectordist_euclidean( (double*)array_INDEXMEM2( s1, i, 0 ),
												 (double*)array_INDEXMEM2( s2, j, 0 ),
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


/**\cond PRIVATE
  */
void local_centre_of_mass( const Array *m, int ii, int ij, int isize,
									int jsize, int *sx, int *sy ){
	 int i, j;
	 double cx=0.0, cy=0.0;
	 double mass=0.0;
	 for( i=0; i<isize; i++ ){
		  for( j=0; j<jsize; j++ ){
				mass += mat_IDX( m, ii+i, ij+j );
				cx += i*mat_IDX( m, ii+i, ij+j );
				cy += j*mat_IDX( m, ii+i, ij+j );
		  }
	 }
	 // dprintf("cx,cy=%f,%f, mass=%f\n", cx, cy, mass);
	 *sx = (int)(cx/mass);
	 *sy = (int)(cy/mass);
}

/** \brief Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot, following
	 the algorithm in Marwan et al.

	 Reference:

	 Marwan et al. Cross recurrence plot based synchronization
	 of time series. Nonlinear Processes in Geophysics (2002)
	 vol. 9 (3-4) pp. 325-331.

	 \todo REPAIR THIS ALGORITHM!

	 Algorithm:
	 \code
	 Input: Recurrence Plot R, Parameters dx, dy

	 \endcode

	 \param R the recurrence plot
	 \return a warp-path (2 x M UINT array giving corresponding points of the LOS)
 */
Array* recplot_los_marwan( const Array *R, int dx, int dy ){
	 bool ismat;
	 matrix_CHECK( ismat, R );
	 if( !ismat ) return NULL;

	 Array *los=array_new2( UINT, 2, 2, R->size[0]+R->size[1] );

	 /* find first point on LOS (step [1]) */
	 int ii, ij, i, j;
	 int pidx=0; /* index in LOS */
	 for( ii=0; ii<R->size[0]; ii++ ){
		  if( mat_IDX(R, ii, 0)>0 )
				break;
	 }
	 for( ij=0; ij<R->size[1]; ij++ ){
		  if( mat_IDX(R, 0, ij)>0 )
				break;
	 }		  
	 if( ii<ij ){
		  for( i=0; i<=ii; i++ ){
				array_INDEX2( los, uint, 0, pidx )=0;
				array_INDEX2( los, uint, 1, pidx )=i;
				pidx++;
		  }
	 } else {
		  for( i=0; i<=ij; i++ ){
				array_INDEX2( los, uint, 0, pidx)=i;
				array_INDEX2( los, uint, 1, pidx)=0;
				pidx++;
		  }
	 }

	 int iip, ijp; /* previous points */
	 iip = array_INDEX2( los, uint, 0, pidx-1 );
	 ijp = array_INDEX2( los, uint, 1, pidx-1 );
	 ii=-1; ij=-1;

	 /* next recurrence point (step [2])*/
	 int w;
	 int sx, sy, kx, ky;
	 int cur[2],prev[2]; /* for bresenham */
	 Array *bres;

	 /* set the next recurrence point to (ii,ij) */
	 do {
		  w=2; ii=-1; ij=-1;
		  while( ii<0 || ij<0 ){
				for( i=w-1; i>=0; i-- ){
					 if( mat_IDX( R, iip+w-1, ijp+i )>0 ){
						  ii=iip+w-1;
						  ij=ijp+i;
						  break;
					 }
				}
				w++;
		  }
		  dprintf("(iip,ijp)->(ii,ij)=(%i,%i)->(%i,%i)\n", iip,ijp,ii, ij);

		  /* centre of mass from this recurrence point */
		  local_centre_of_mass( R, ii, ij, dx, dy, &sx, &sy );

		  ky=sy; kx=sx;
		  for( i=0; i<sy; i++ )
				if( mat_IDX( R, ii+i, ij )<=0 ) /* non-rec point */
					 ky=i;
		  for( j=0; j<sx; j++ )
				if( mat_IDX( R, ii, ij+j )<=0 ) /* non-rec point */
					 kx=j;

		  ii+=kx/2;
		  ij+=ky/2;

		  /* (iip,ijp) is the previous LOSpoint, (ii,ij) the current one */
		  /* linear interpolation between (iip,ijp)->(ii,ij) */
		  prev[0]=iip; prev[1]=ijp; cur[0]=ii; cur[1]=ij;
		  bres=bresenham_line( prev, cur );
		  for( i=0; i<bres->size[1]; i++ ){
				array_INDEX2( los, uint, 0, pidx ) =
						  array_INDEX2( bres, uint, 0, i );
				array_INDEX2( los, uint, 1, pidx ) =
						  array_INDEX2( bres, uint, 1, i );
				pidx++;
		  }

		  iip=ii; ijp=ij;
		  array_free( bres );
		  dprintf("pidx=%i, (%i,%i)\n", pidx, ii, ij);
	 }	 while( ii<R->size[0]-1 && ij<R->size[1]-1 );


	 return los;
}

 /** \brief Calculate the Line-Of-Synchrony of a Cross-Recurrence-Plot using
	 a dynamic time-warping strategy.

	 Algorithm:
	 \code
	 Input: RecurrencePlot R (1 where recurrence, 0 otherwise)
	 1. calculate d = 2-R;
	 2. add a small amount of noise:
		  d_jk = d_jk + epsilon
	 2. cumulate d, such that D[jk] = min{ D[j-1k], D[jk-1], D[j-1k-1] }
	 3. Backtrack
	 \endcode
	 (see \ref dtw for details)

	 \param R the recurrence plot
	 \return a warp-path (2 x M UINT array giving corresponding points of the LOS)
 */

Array* recplot_los_dtw_noise( const Array *R ){
  double noise;
  double noiseamp = 0.000000001;

  srand((long)time(NULL));

  bool ismat;
  matrix_CHECK( ismat, R );
  if( !ismat ) return NULL;

  Array *d = array_copy( R, TRUE );

  ulong i;
  for( i=0; i<array_NUMEL(d); i++ ){
		/* flip binary matrix and add noise*/
		noise = noiseamp*(((double)rand()) / RAND_MAX);
		array_INDEX1( d, double, i ) = (1.0-array_INDEX1( d, double, i ))+noise;
  }


  OptArgList *opt = optarglist( "slope_constraint=int", SLOPE_CONSTRAINT_LAX );
  matrix_dtw_cumulate( d, FALSE, opt );
  optarglist_free( opt );

  Array *P = matrix_dtw_backtrack( d );

  array_free( d );

  return P;
}
