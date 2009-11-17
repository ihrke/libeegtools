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

#include "regularization.h"
#include "optarg.h"

/** Calculate the regularization function that is the distance
	 transform of $f$, the piecwise linear interpolation between
	 event-markers (approximated with bresenham-alg). This is a linear
	 fall-off away from the line passing through all time-marker pairs.

	 \note if you don't pass any arguments, the regularization is done 
	 along the main diagonal of the DTW-matrix.

	 \param d matrix or NULL (alloc'd in function)
	 \param n number of points in the returned matrix
	 \param optargs may contain:
	 - <tt>markers=int**</tt>time-markers within the nsignal x nsignal matrix, default=<tt>( (0,n-1), (0,n-1) )</tt>
	 - <tt>nmarkers=int</tt>number of time-marker pairs; REQUIRED if markers is set, default=\c 2
	 \return d or NULL (error)
*/
double** regularization_linear_points( double **d, int n, OptArgList *optargs ){ 
  int i,j;
  int maxmem, npoints;
  double x;
  int *points;
  int **I;
  int **markers;
  int nmarkers;
  void *ptr;

  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( n, n );
  } 

  /* defaults */
  nmarkers = 2;
  markers = NULL;

  /* override */
  if( optarglist_has_key( optargs, "nmarkers" ) ){
	 x = optarglist_scalar_by_key( optargs, "nmarkers" );
	 if( !isnan( x ) ) nmarkers=(int)x;
  }
  if( optarglist_has_key( optargs, "markers" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "markers" );
	 if( ptr ) markers = (int**)ptr;
  }

  /* allocate own markers if not provided via optargs */
  if( !markers ){
	 nmarkers = 2;
	 markers = (int**) malloc( 2*sizeof( int* ) );
	 markers[0] = (int*) malloc( nmarkers*sizeof( int ) );
	 markers[1] = (int*) malloc( nmarkers*sizeof( int ) );
	 markers[0][0] = 0;
	 markers[1][0] = 0;
	 markers[0][1] = n-1;
	 markers[1][1] = n-1;
  }
	 
  /* memory for distance transform */
  I = (int**) malloc( n*sizeof( int* ) );
  for( i=0; i<n; i++ ){
	 I[i] = (int*) malloc( n*sizeof( int ) );
	 memset( I[i], 0, n*sizeof( int ) );
  }

  /* compute memory needed for bresenham */
  maxmem=0;
  for( i=0; i<nmarkers-1; i++ ){
	 npoints = bresenham_howmany_points( markers[0][i], markers[1][i], markers[0][i+1], markers[1][i+1] );
	 dprintf("npoints=%i\n", npoints );
	 maxmem  = MAX( npoints, maxmem );
  }
  dprintf("2*(maxmem+1)=%i\n", 2*(maxmem+1));
  points = (int*) malloc( (2*(maxmem+1))*sizeof(int) );
  
  /* compute bresenham for the line segments */
  for( i=0; i<nmarkers-1; i++ ){
	 npoints = bresenham_howmany_points( markers[0][i], markers[1][i], markers[0][i+1], markers[1][i+1] );
	 dprintf("Line from (%i,%i)->(%i,%i) with %i points\n", 
				markers[0][i], markers[1][i], markers[0][i+1], markers[1][i+1], npoints);
	 points  = bresenham( markers[0][i], markers[1][i], markers[0][i+1], markers[1][i+1], points );
	 for( j=0; j<2*npoints; j+=2 ){ /* draw the line */
		I[points[j]][points[j+1]] = 1;
	 }
  }
	 
  /* distance transform of line-segments */
  d = disttransform_deadreckoning( I, n, n, d );

  /* free */
  for( i=0; i<n; i++ ){
	 free( I[i] );
  }
  free(I);
  free( points );

  return d;
}

/** Calculate regularization function
	 \f[
	 G_f(x,y; \sigma) = \frac{1}{\sigma 2\pi} \exp{\left( -\frac{\min_{\xi}\sqrt{(\xi-x)^2+(y-f(\xi))^2}}{2\sigma^2} \right)}
	 \f]
	 where $f$ is piecwise linear (approximated with bresenham-alg) and the minimization
	 is approximated with distance-transform (deadreckoning).

	 \todo there are artifacts here, check for numerical errors

	 \param n number of points in the returned matrix
	 \param d matrix or NULL (alloc'd in function)
	 \param optargs may contain:
	 - <tt>markers=int**</tt> time-markers within the nsignal x nsignal matrix, default=<tt>( (0,n-1), (0,n-1) )</tt>
	 - <tt>nmarkers=int</tt> number of time-marker pairs; REQUIRED if markers is set, default=\c 2
	 - <tt>max_sigma=double</tt> std of a gaussian applied to the regularization matrix as returned from regularization_linear_points(); default=\c 0.2	 \return d or NULL (error)
 */
double** regularization_gaussian_line( double **d, int n, OptArgList *optargs ){
  int i,j,k;
 
  double max_sigma;
  double sigma;
  double maxdist, dist, closest_dist, normgauss;
  int flag;
  double pointdist=0.1;
  int nmarkers;
  int **markers;
  double x;
  void *ptr;

  /* get distance transform of the linear interpolation between markers */
  d = regularization_linear_points( d, n, optargs );

  /* defaults */
  nmarkers = 2;
  markers = NULL;

  /* override */
  if( optarglist_has_key( optargs, "nmarkers" ) ){
	 x = optarglist_scalar_by_key( optargs, "nmarkers" );
	 if( !isnan( x ) ) nmarkers=(int)x;
  }
  if( optarglist_has_key( optargs, "markers" ) ){
	 ptr = optarglist_ptr_by_key( optargs, "markers" );
	 if( ptr ) markers = (int**)ptr;
  }

  /* allocate own markers if not provided via optargs */
  if( !markers ){
	 nmarkers = 2;
	 markers = (int**) malloc( 2*sizeof( int* ) );
	 markers[0] = (int*) malloc( nmarkers*sizeof( int ) );
	 markers[1] = (int*) malloc( nmarkers*sizeof( int ) );
	 markers[0][0] = 0;
	 markers[1][0] = 0;
	 markers[0][1] = n-1;
	 markers[1][1] = n-1;
  }
	 
  /* apply gaussian with varying sigma */
  flag = 0;
  maxdist = sqrt(2.0)*(double)(n-1);
  for( i=0; i<n; i++ ){
  	 for( j=0; j<n; j++ ){
		closest_dist = DBL_MAX;
  		for( k=0; k<nmarkers; k++ ){
  		  dist = ( SQR( (double)(markers[0][k]-i) )+SQR( (double)(markers[1][k]-j) ) );
  		  dist = sqrt( dist - SQR( d[i][j] ) );
		  dist /= maxdist;
		  if(dist<closest_dist){
			 closest_dist=dist;
		  }
  		  if( dist<(pointdist) )
  			 flag=1;
  		}

		if( closest_dist==0 ){
		  closest_dist = 1e-10;
		}
		sigma = max_sigma*maxdist;
  		if( flag )
		  sigma = sigma*closest_dist/pointdist;

		normgauss = gaussian( 0.0, sigma, 0 );
		if(normgauss > 10000000) {
		  d[i][j]=1;
		} else {
		  d[i][j] = gaussian( d[i][j], sigma, 0 )/normgauss;
		}
		
  		flag = 0;
  	 }
  }

  return d;
}

