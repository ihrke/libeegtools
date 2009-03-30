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


/** Calculate regularization function
	 \f[
	 G_f(x,y; \sigma) = \frac{1}{\sigma 2\pi} \exp{\left( -\frac{\min_{\xi}\sqrt{(\xi-x)^2+(y-f(\xi))^2}}{2\sigma^2} \right)}
	 \f]
	 where $f$ is piecwise linear (approximated with bresenham-alg) and the minimization
	 is approximated with distance-transform (deadreckoning).
	 \param markers1,2 time-markers within the nsignal x nsignal matrix
	 \param nsignal number of points in the returned matrix
	 \param maxsigma
	 \param d matrix or NULL (alloc'd in function)
	 \return d or NULL (error)
 */
double** regularization_gaussian_line( const int *markers1, const int *markers2, int nmarkers, int nsignal, double maxsigma, double **d ){
  int i,j,k;
  int maxmem, npoints;
  int *points;
  int **I;
  double sigma;
  double maxdist, dist, closest_dist, normgauss;
  int flag;
  double pointdist=0.1;
  int n;

  /* convenience */
  n = nsignal;

  if( d==ALLOC_IN_FCT ){
	 d=matrix_init( nsignal, nsignal );
  } 

  I = (int**) malloc( n*sizeof( int* ) );
  for( i=0; i<n; i++ ){
	 I[i] = (int*) malloc( n*sizeof( int ) );
	 memset( I[i], 0, n*sizeof( int ) );
  }


 /* compute memory needed for bresenham */
  maxmem=0;
  for( i=0; i<nmarkers+1; i++ ){
	 npoints = bresenham_howmany_points( markers1[i], markers2[i], markers1[i+1], markers2[i+1] );
	 dprintf("npoints=%i\n", npoints );
	 maxmem  = MAX( npoints, maxmem );
  }
  dprintf("2*(maxmem+1)=%i\n", 2*(maxmem+1));
  points = (int*) malloc( (2*(maxmem+1))*sizeof(int) );

  /* compute bresenham for the line segments */
  for( i=0; i<nmarkers+1; i++ ){
	 npoints = bresenham_howmany_points( markers1[i], markers2[i], markers1[i+1], markers2[i+1] );
	 dprintf("Line from (%i,%i)->(%i,%i) with %i points\n", markers1[i], markers2[i], markers1[i+1], markers2[i+1], npoints);
	 points  = bresenham( markers1[i], markers2[i], markers1[i+1], markers2[i+1], points );
	 for( j=0; j<2*npoints; j+=2 ){ /* draw the line */
		I[points[j]][points[j+1]] = 1;
	 }
  }
	 
  /* distance transform of line-segments */
  d = disttransform_deadreckoning( I, n, n, d );
  
  /* apply gaussian with varying sigma */
  flag = 0;
  maxdist = sqrt(2.0)*(double)(n-1);
  for( i=0; i<n; i++ ){
  	 for( j=0; j<n; j++ ){
		closest_dist = DBL_MAX;
  		for( k=0; k<nmarkers+2; k++ ){
  		  dist = ( SQR( (double)(markers1[k]-i) )+SQR( (double)(markers2[k]-j) ) );
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
		sigma = maxsigma*maxdist;
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
	 \param s1,s2 
	 \param maxsigma
	 \param d matrix or NULL (alloc'd in function)
	 \return d or NULL (error)
 */
double** eeg_regularization_gaussian_line( const EEGdata *s1, const EEGdata *s2, 
														 double maxsigma,
														 double **d ){
  int i;
  int *m1, /* markers for signal 1 */
	 *m2;   /* markers for signal 2 */
  int nmarkers,n;

  if( eegdata_cmp_settings( s1, s2 ) ){
	 errprintf(" ERROR: two EEGdata sets are not similar enough\n");
	 return NULL;
  }
  nmarkers=s1->nmarkers;
  n = s1->n;


  /* for convenience, add (0,0) and (n,n) as markers */
  m1 = (int*) malloc( (nmarkers+2)*sizeof( int ) );
  m2 = (int*) malloc( (nmarkers+2)*sizeof( int ) );
  m1[0]=0; m2[0]=0;
  m1[nmarkers+1]=n-1; m2[nmarkers+1]=n-1;
  for( i=1; i<nmarkers+1; i++ ){
	 m1[i] = s1->markers[i-1];
	 m2[i] = s2->markers[i-1];
	 dprintf("m=%i,%i, nm=%i\n", m1[i], m2[i], nmarkers);
  }

  d = regularization_gaussian_line( m1, m2, nmarkers, n, maxsigma, d );
 
  free(m1); free(m2);

  return d;
}
