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
#include "linalg.h"
#include "imageproc.h"

/** \brief Calculate a regularization function that is the distance
	 transform of the piecwise linear interpolation between
	 points.
	 
	 (approximated with bresenham-alg). This is a linear
	 fall-off away from the line passing through all points.
	 This function is used e.g. to construct a regularization matrix
	 that punishes deviations from time-markers (stimulus onset,
	 response onset etc).

	 Example:
	 \image html regularization_line.jpg

	 \note if you don't pass any arguments, the regularization is done 
	 along the main diagonal
	 
	 \param points defining the piecewise linear function through the regularization matrix;
	        this is a 2 x M dimensional INT-array, all points must be within the 
			  dimensions; if NULL is passed, the function assumes (0,0),(dims[0]-1,dims[1]-1) 
			  as points, i.e. the regularization is done along the main diagonal
    \param dims the dimensions of the output matrix (rows x cols)
	 \param m the output matrix or NULL -> allocate in function 
	 \return the regularization matrix or NULL (error)
*/
Array* regularization_linear_points( const Array *points, uint dims[2], Array *m ){
  if( points->dtype!=INT ){
	 errprintf("Need INT-array for coordinates\n"); return NULL;
  }
  if( points->ndim!=2 ){
	 errprintf("Need 2D-array for coordinates\n"); return NULL;
  }
  if( points->size[1]<2 ){
	 errprintf("Need 2 coordinates per point\n"); return NULL;
  }
  if( m!=NULL ){
	 bool ismatrix;
	 matrix_CHECK( ismatrix, m );
	 if( !ismatrix) return NULL;
	 if( m->size[0]!=dims[0] || m->size[1]!=dims[1] ){
		errprintf( "output matrix must be of same dimension as dims[2]\n");
		return NULL;
	 }
  } else {
	 m=array_new2(DOUBLE, 2, dims[0], dims[1] );
  }

  /* bresenham bitmap-image */
  Array *mask=array_new2( INT, 2, dims[0], dims[1] );
  
  /* get line segments */
  Array *bres;
  if( points==NULL ){
	 Array *defpoints=array_new2( INT, 2, 2, 2 );
	 array_INDEX2( points,int, 0, 0)=0;
	 array_INDEX2( points,int, 1, 0)=0;
	 array_INDEX2( points,int, 0, 1)=dims[0]-1;
	 array_INDEX2( points,int, 1, 1)=dims[1]-1;
	 bres=bresenham_linesegments( defpoints );
	 array_free( defpoints );
  } else {
	 bres=bresenham_linesegments( points );
  }

  /* set to mask */
  int i;
  for( i=0; i<bres->size[1]; i++ ){
	 array_INDEX2(mask,int, array_INDEX2(bres,int,0,i),
					  array_INDEX2(bres,int,1,i))=1;
  }
  m=disttransform_deadreckoning( mask, m );

  array_free( mask );
  
  return m;
}

/** \brief Calculate a ''gaussian corridor''.
	 
	 \f[
	 G_f(x,y; \sigma) = \frac{1}{\sigma 2\pi} \exp{\left( -\frac{\min_{\xi}\sqrt{(\xi-x)^2+(y-f(\xi))^2}}{2\sigma^2} \right)}
	 \f]
	 where $f$ is piecwise linear (approximated with bresenham-alg) and the minimization
	 is approximated with distance-transform (deadreckoning).

	 \todo there are artifacts here, check for numerical errors

	 Example:
	 \image html regularization_gauss.jpg

	 \note if you don't pass any arguments, the regularization is done 
	 along the main diagonal
	 
	 \param points defining the piecewise linear function through the regularization matrix;
	        this is a 2 x M dimensional INT-array, all points must be within the 
			  dimensions; if NULL is passed, the function assumes (0,0),(dims[0]-1,dims[1]-1) 
			  as points, i.e. the regularization is done along the main diagonal
    \param dims the dimensions of the output matrix (rows x cols)
	 \param m the output matrix or NULL -> allocate in function 
	 \param max_sigma the regularization parameter
	 \return the regularization matrix or NULL (error)
*/
Array* regularization_gaussian_corridor( const Array *points, uint dims[2], Array *m, double max_sigma ){
  if( points->dtype!=INT ){
	 errprintf("Need INT-array for coordinates\n"); return NULL;
  }
  if( points->ndim!=2 ){
	 errprintf("Need 2D-array for coordinates\n"); return NULL;
  }
  if( points->size[1]<2 ){
	 errprintf("Need 2 coordinates per point\n"); return NULL;
  }
  if( m!=NULL ){
	 bool ismatrix;
	 matrix_CHECK( ismatrix, m );
	 if( !ismatrix) return NULL;
	 if( m->size[0]!=dims[0] || m->size[1]!=dims[1] ){
		errprintf( "output matrix must be of same dimension as dims[2]\n");
		return NULL;
	 }
  } else {
	 m=array_new2(DOUBLE, 2, dims[0], dims[1] );
  }
  
  /* get distance transform of the linear interpolation between markers */
  m = regularization_linear_points( points, dims, m );

  int i,j,k;
  double sigma;
  double maxdist, dist, closest_dist, normgauss;
  int flag;
  double pointdist=0.1;

  /* apply gaussian with varying sigma */
  flag = 0;
  maxdist = sqrt(2.0)*(double)(dims[0]-1);

  dprintf("maxdist=%f, max_sigma=%f\n", maxdist, max_sigma);
  for( i=0; i<dims[0]; i++ ){
  	 for( j=0; j<dims[1]; j++ ){
		closest_dist = DBL_MAX;
  		for( k=0; k<points->size[1]; k++ ){
  		  dist = ( SQR( (double)(array_INDEX2(points,int,0,k)-i) )
					  + SQR( (double)(array_INDEX2(points,int,1,k)-j) ) );
  		  dist = sqrt( dist - SQR( mat_IDX(m,i,j) ) );
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
		  mat_IDX(m,i,j)=1;
		} else {
		  mat_IDX(m,i,j)=gaussian( mat_IDX(m,i,j), sigma, 0 )/normgauss;
		}
		
  		flag = 0;
  	 }
  }

  return m;
}

